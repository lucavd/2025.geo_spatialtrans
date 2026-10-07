# Revisione avversariale S1.1 — 06_moore.py
# CP-1: contorno sui CENTRI dei pixel, implementazione propria.
#  (a) Moore-neighbour tracing (8 vicini, orario con y verso il basso, criterio d'arresto di Jacob:
#      si torna al pixel di partenza entrando dalla stessa direzione), area di Gauss del bordo esterno;
#  (b) marching squares (skimage.measure.find_contours a livello 0.5) con buchi sottratti;
#  (c) Moore esterno meno l'area (Moore) dei buchi di fondo, per separare i due effetti di segno opposto.
import sys, json, numpy as np, pandas as pd
from scipy import ndimage as ndi
from skimage import measure
EXP = "/mnt/micron/geo_spatialtrans/S1.1/review/export"
OUT = "/mnt/micron/geo_spatialtrans/S1.1/review/py_out"
S4 = np.array([[0,1,0],[1,1,1],[0,1,0]]); S8 = np.ones((3,3), int)
OFF = [(0,-1),(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1)]   # W NW N NE E SE S SW (orario, y giu')
IDX = {o: k for k, o in enumerate(OFF)}

def shoelace(P):
    r, c = P[:, 0], P[:, 1]
    return abs(np.dot(c, np.roll(r, -1)) - np.dot(np.roll(c, -1), r)) / 2

def moore_outer(B):
    B = np.pad(B, 1)
    rr, cc = np.nonzero(B)                       # ordine di riga: il primo e' il piu' in alto a sinistra
    if len(rr) <= 1: return 0.0
    s = (int(rr[0]), int(cc[0]))
    p = s; bdir = 0                              # il vicino ovest e' fondo
    path = [p]; q1 = None; closed = False
    for _ in range(8 * len(rr) + 16):
        found = False
        for k in range(1, 9):
            d = (bdir + k) % 8
            q = (p[0] + OFF[d][0], p[1] + OFF[d][1])
            if B[q]: found = True; break
        if not found: break
        if q1 is None: q1 = q
        elif p == s and q == q1: closed = True; break   # arresto: si ripete la prima mossa (s -> q1)
        bp = (p[0] + OFF[(d - 1) % 8][0], p[1] + OFF[(d - 1) % 8][1])
        bdir = IDX[(bp[0] - q[0], bp[1] - q[1])]
        p = q; path.append(p)
    assert closed, "Moore non chiuso"
    P = np.array(path, float)
    if len(P) > 1 and tuple(P[-1]) == s: P = P[:-1]
    return 0.0 if len(P) < 3 else shoelace(P)


def ms_area(B):
    Bp = np.pad(B, 1).astype(float)
    a = 0.0
    for C in measure.find_contours(Bp, 0.5):
        r, c = C[:, 0], C[:, 1]
        a += (np.dot(c, np.roll(r, -1)) - np.dot(np.roll(c, -1), r)) / 2   # con segno: buchi opposti
    return a

def comp_rows(id_):
    df = pd.read_parquet(f"{EXP}/{id_}.parquet")
    x = df.x.astype(int).values; y = df.y.astype(int).values
    lv, code = np.unique(df.cluster.values, return_inverse=True)
    G = np.zeros((y.max() - y.min() + 1, x.max() - x.min() + 1), np.int32); G[y - y.min(), x - x.min()] = code + 1
    rows = []
    for k in range(1, len(lv) + 1):
        L, n = ndi.label(G == k, structure=S4)
        for sl_i, sl in enumerate(ndi.find_objects(L)):
            B = L[sl] == sl_i + 1
            npx = int(B.sum())
            mo = moore_outer(B)
            # buchi: fondo (non-B) 8-connesso... il complemento di una regione 4-connessa e' 8-connesso
            Bp = np.pad(B, 1); H, nh = ndi.label(~Bp, structure=S8)
            outside = H[0, 0]
            hole_px = int(((H != outside) & (H > 0)).sum())
            msa = ms_area(B)
            rows.append((id_, lv[k - 1], npx, mo, hole_px, msa))
    return pd.DataFrame(rows, columns=["id", "cluster", "n_px", "area_moore_outer", "hole_px", "area_ms_signed"])

ids = sys.argv[1:]
allr = pd.concat([comp_rows(i) for i in ids], ignore_index=True)
allr.to_csv(f"{OUT}/moore_components_{ids[0]}.csv", index=False)
summ = allr.groupby("id").apply(lambda r: pd.Series({
    "n_comp": len(r), "sum_npx": r.n_px.sum(),
    "relerr_moore_outer": (r.area_moore_outer.sum() - r.n_px.sum()) / r.n_px.sum(),
    "relerr_moore_vs_npx_plus_holes": (r.area_moore_outer.sum() - (r.n_px + r.hole_px).sum()) / r.n_px.sum(),
    "relerr_marching_squares": (r.area_ms_signed.abs().sum() - r.n_px.sum()) / r.n_px.sum(),
    "frac_hole_px": r.hole_px.sum() / r.n_px.sum()}))
summ.to_csv(f"{OUT}/moore_summary_{ids[0]}.csv")
print(summ.to_string())

