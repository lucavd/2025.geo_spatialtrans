#!/usr/bin/env python
"""tools/S1.4c_raster.py — S1.4c: operatore di osservazione unico sulle maschere, per nuclei reali e simulati.

Convenzione (verificata sui reali in sviluppo, scarto 0 su A1 r1 cellpose e A4 f1 StarDist): il pixel (r, c) ha
coordinate x = c * upp, y = r * upp; area = n_pixel * upp^2. Un poligono simulato copre il pixel (r, c) se il punto
(r, c), in unita' di pixel, cade dentro (skimage.draw.polygon). Misure per etichetta:
  n_px, area, cx, cy (centroide dei pixel), sol (solidita' di skimage regionprops: area / area dell'inviluppo
  convesso dei pixel), n_bnd (pixel con un 4-vicino di etichetta diversa, fuori immagine = 0), n_c1, n_c2 (pixel di
  bordo con un'altra etichetta > 0 entro distanza di Chebyshev 1, 2 px), frac_out (D-1: frazione dei pixel il cui
  centroide piu' vicino non e' il proprio = Voronoi dei centroidi, non ritagliato).
Mutanti (variabile S14C_MUTANT): MU3 raster sfasato di mezzo pixel; MU4 delta = 0 (nessun contatto).
Uso (root del repo):
  .venv/bin/python tools/S1.4c_raster.py test
  .venv/bin/python tools/S1.4c_raster.py sim <in_dir> <out_dir> [ncores]
  .venv/bin/python tools/S1.4c_raster.py real <out_dir>
  .venv/bin/python tools/S1.4c_raster.py cp2 <out_dir> [ncores]
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
import shapely
from skimage.draw import polygon as sk_polygon
from skimage.measure import regionprops

MUT = os.environ.get("S14C_MUTANT", "none")
R3 = "/mnt/micron/geo_spatialtrans/R3"
MASKS = "/mnt/micron/geo_spatialtrans/R2/masks"
SHAPE_SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he",
             "A5": "cellpose_rgb", "A6": "stardist_he"}
BIG = np.iinfo(np.int64).max
CP2_MARGIN_UM = 10.0


def rasterize(geoms, ids, shape, upp):
    """Etichette per pixel; restituisce (lab, n_conflict = pixel coperti da piu' di un poligono)."""
    lab = np.zeros(shape, np.int64); cnt = np.zeros(shape, np.uint16)
    off = 0.5 if MUT == "MU3" else 0.0
    H, W = shape
    for g, i in zip(geoms, ids):
        if g is None or shapely.is_empty(g):
            continue
        for p in shapely.get_parts(g):
            if p.geom_type != "Polygon" or p.is_empty:
                continue
            ex = np.asarray(p.exterior.coords)
            rr = ex[:, 1] / upp + off; cc = ex[:, 0] / upp + off
            r0 = max(int(np.floor(rr.min())), 0); r1 = min(int(np.ceil(rr.max())) + 1, H)
            c0 = max(int(np.floor(cc.min())), 0); c1 = min(int(np.ceil(cc.max())) + 1, W)
            if r1 <= r0 or c1 <= c0:
                continue
            loc = np.zeros((r1 - r0, c1 - c0), bool)
            a, b = sk_polygon(rr - r0, cc - c0, loc.shape); loc[a, b] = True
            for h in p.interiors:
                hc = np.asarray(h.coords)
                a, b = sk_polygon(hc[:, 1] / upp + off - r0, hc[:, 0] / upp + off - c0, loc.shape); loc[a, b] = False
            sub = lab[r0:r1, c0:c1]; sub[loc] = i
            cnt[r0:r1, c0:c1] += loc
    return lab, int(np.sum(cnt > 1))


def boundary(lab):
    P = np.pad(lab, 1); core = P[1:-1, 1:-1]
    return (core > 0) & ((P[:-2, 1:-1] != core) | (P[2:, 1:-1] != core) | (P[1:-1, :-2] != core) | (P[1:-1, 2:] != core))


def other_within(lab, d):
    """Pixel con un'altra etichetta > 0 entro Chebyshev <= d."""
    dd = 0 if MUT == "MU4" else d
    if dd == 0:
        return np.zeros(lab.shape, bool)
    size = 2 * dd + 1
    S = int(lab.max()) + 1                                   # sentinella piccola: cval di ndimage passa per double
    mx = ndi.maximum_filter(lab, size=size, mode="constant", cval=0)
    mn = ndi.minimum_filter(np.where(lab > 0, lab, S), size=size, mode="constant", cval=S)
    return (lab > 0) & ((mx > lab) | (mn < lab))


def measure(lab, upp, deltas=(1, 2), solidity=True, frac_out=True):
    lab = lab.astype(np.int64)
    r, c = np.nonzero(lab); l = lab[r, c]
    ids = np.unique(l); M = int(lab.max()) + 1
    n = np.bincount(l, minlength=M)
    sx = np.bincount(l, c.astype(np.float64), M); sy = np.bincount(l, r.astype(np.float64), M)
    out = pd.DataFrame({"label": ids, "n_px": n[ids], "area": n[ids] * upp ** 2,
                        "cx": sx[ids] / n[ids] * upp, "cy": sy[ids] / n[ids] * upp})
    bnd = boundary(lab)
    out["n_bnd"] = np.bincount(lab[bnd], minlength=M)[ids]
    for d in deltas:
        oth = other_within(lab, d)
        out[f"n_c{d}"] = np.bincount(lab[bnd & oth], minlength=M)[ids]
    if solidity:
        sol = np.full(M, np.nan)
        for p in regionprops(lab):
            sol[p.label] = p.solidity
        out["sol"] = sol[ids]
    if frac_out:
        tree = cKDTree(np.column_stack([out["cx"].to_numpy(), out["cy"].to_numpy()]))
        _, near = tree.query(np.column_stack([c * upp, r * upp]), k=1, workers=4)
        pos = np.searchsorted(ids, l)
        bad = np.bincount(pos[near != pos], minlength=len(ids))
        out["frac_out"] = bad / out["n_px"].to_numpy()
    return out


def contact_pairs(lab, d=2):
    """Coppie non ordinate (i < j) di etichette > 0 con pixel entro Chebyshev <= d."""
    H, W = lab.shape; M = int(lab.max()) + 1; keys = []
    for dr in range(0, d + 1):
        for dc in range(-d, d + 1):
            if dr == 0 and dc <= 0:
                continue
            c0 = max(0, -dc); c1 = W - max(0, dc)
            a = lab[0:H - dr, c0:c1]; b = lab[dr:H, c0 + dc:c1 + dc]
            m = (a > 0) & (b > 0) & (a != b)
            lo = np.minimum(a[m], b[m]); hi = np.maximum(a[m], b[m])
            keys.append(np.unique(lo.astype(np.int64) * M + hi))
    k = np.unique(np.concatenate(keys)) if keys else np.zeros(0, np.int64)
    return pd.DataFrame({"label_i": k // M, "label_j": k % M})


def rois_table():
    return pd.read_csv("results/R2/R2_rois_checked.csv")


# ---------------------------------------------------------------------------------------------------- stadio sim
def sim_one(f, out_dir, rois):
    o = os.path.join(out_dir, os.path.basename(f))
    if os.path.exists(o):
        return "skip"
    A, roi = os.path.basename(f).split("_")[:2]
    rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]
    upp = float(rw.um_per_px); side = int(rw.side_px)
    t = pd.read_parquet(f)
    geoms = shapely.from_wkb(t["wkb"].to_numpy())
    lab, ncf = rasterize(geoms, t["idx"].to_numpy(), (side, side), upp)
    m = measure(lab, upp)
    m["n_conflict"] = ncf
    m.to_parquet(o + ".tmp", index=False); os.replace(o + ".tmp", o)
    return "ok"


# ---------------------------------------------------------------------------------------------------- stadio real
def real_one(A, roi, out_dir):
    seg = SHAPE_SEG[A]
    o = os.path.join(out_dir, f"{A}_{roi}_{seg}_shape.parquet")
    if os.path.exists(o):
        return
    gen = pd.read_parquet(f"{R3}/real/{A}_{roi}_{seg}_gen.parquet")
    lab = np.load(f"{MASKS}/{A}_{roi}_{seg}_native.npz")["labels"].astype(np.int64)
    upp = float(gen["um_per_px"].iloc[0])
    keep = np.zeros(int(lab.max()) + 1, bool); keep[gen["label"].to_numpy().astype(np.int64)] = True
    lab[~keep[lab]] = 0
    m = measure(lab, upp)
    pr = contact_pairs(lab, 2)
    g = gen.set_index("label").loc[m["label"]]
    chk = dict(archetype=A, roi=roi, seg=seg, n_gen=len(gen), n_lab=len(m),
               dx=float(np.max(np.abs(m["cx"].to_numpy() - g["x"].to_numpy()))),
               dy=float(np.max(np.abs(m["cy"].to_numpy() - g["y"].to_numpy()))),
               darea=float(np.max(np.abs(m["area"].to_numpy() - g["area_nuc"].to_numpy()))))
    pr.to_parquet(os.path.join(out_dir, f"{A}_{roi}_{seg}_pairs.parquet"), index=False)
    pd.DataFrame([chk]).to_csv(os.path.join(out_dir, f"{A}_{roi}_{seg}_cop.csv"), index=False)
    m.to_parquet(o + ".tmp", index=False); os.replace(o + ".tmp", o)


# ---------------------------------------------------------------------------------------------------- stadio cp2
def read_valid(A, roi):
    from PIL import Image
    v = np.asarray(Image.open(f"results/R2/tissue_masks/{A}_{roi}_valid_ds4.png"))
    if v.ndim == 3:
        v = v[..., 0]
    return v > (127 if v.dtype == np.uint8 else 0.5)


def cp2_one(A, roi, seg, out_dir, rois):
    o = os.path.join(out_dir, f"{A}_{roi}_{seg}_cp2.parquet")
    if os.path.exists(o):
        return "skip"
    rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]
    upp = float(rw.um_per_px); side_um = float(rw.side_px) * upp
    lab = np.load(f"{MASKS}/{A}_{roi}_{seg}_native.npz")["labels"].astype(np.int64)
    v = read_valid(A, roi); pxv = side_um / v.shape[1]
    dist = ndi.distance_transform_edt(v) * pxv               # distanza (µm) dal tessuto non valido, griglia ds4
    r, c = np.nonzero(lab); l = lab[r, c]; M = int(lab.max()) + 1
    n = np.bincount(l, minlength=M); ids = np.nonzero(n)[0]
    cx = np.bincount(l, c.astype(float), M)[ids] / n[ids] * upp; cy = np.bincount(l, r.astype(float), M)[ids] / n[ids] * upp
    ci = np.clip((cx / pxv).astype(int), 0, v.shape[1] - 1); ri = np.clip((cy / pxv).astype(int), 0, v.shape[0] - 1)
    inside = v[ri, ci]; ev = dist[ri, ci] >= CP2_MARGIN_UM
    keep = np.zeros(M, bool); keep[ids[inside]] = True
    lab[~keep[lab]] = 0
    m = measure(lab, upp, frac_out=False)
    evl = set(ids[inside & ev].tolist())
    m["ev"] = m["label"].isin(evl)
    m["archetype"] = A; m["roi"] = roi; m["seg"] = seg
    m.to_parquet(o + ".tmp", index=False); os.replace(o + ".tmp", o)
    return "ok"


# ---------------------------------------------------------------------------------------------------- test C-R
def brute(lab, d):
    """Arbitro a cicli espliciti (indipendente dai filtri): n_bnd e n_c<d> per etichetta."""
    H, W = lab.shape; nb = {}; nc = {}
    for i in range(H):
        for j in range(W):
            v = lab[i, j]
            if v == 0:
                continue
            isb = False
            for (a, b) in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
                w = lab[a, b] if 0 <= a < H and 0 <= b < W else 0
                if w != v:
                    isb = True
            if not isb:
                continue
            nb[v] = nb.get(v, 0) + 1
            hit = False
            for a in range(max(0, i - d), min(H, i + d + 1)):
                for b in range(max(0, j - d), min(W, j + d + 1)):
                    w = lab[a, b]
                    if w > 0 and w != v:
                        hit = True
            if hit:
                nc[v] = nc.get(v, 0) + 1
    return nb, nc


def run_tests():
    upp = 0.2737612671531527; res = []
    def ck(name, ok, val):
        res.append((name, bool(ok), val)); print(f"{'PASS' if ok else 'FAIL'} {name} {val}")
    # T1 disco r = 5 µm, centro non intero: area, centroide, solidita'
    x0, y0, R = 30.123, 27.871, 5.0
    disc = shapely.Point(x0, y0).buffer(R, quad_segs=256)
    lab, _ = rasterize([disc], [1], (220, 220), upp); m = measure(lab, upp)
    ck("T1 area disco |dA| <= P*upp", abs(m.area[0] - np.pi * R ** 2) <= 2 * np.pi * R * upp, round(float(m.area[0] - np.pi * R ** 2), 4))
    ck("T1 centroide disco <= 0.25 px", max(abs(m.cx[0] - x0), abs(m.cy[0] - y0)) <= 0.25 * upp, round(float(max(abs(m.cx[0] - x0), abs(m.cy[0] - y0)) / upp), 4))
    ck("T1 solidita' disco >= 0.96", m.sol[0] >= 0.96, round(float(m.sol[0]), 4))
    ck("T1 nessun contatto (un solo nucleo)", m.n_c2[0] == 0, int(m.n_c2[0]))
    # T2 due semidischi con corda verticale comune
    h1 = disc.intersection(shapely.box(0, 0, x0, 100)); h2 = disc.intersection(shapely.box(x0, 0, 100, 100))
    lab, ncf = rasterize([h1, h2], [1, 2], (220, 220), upp); m = measure(lab, upp)
    f = (m.n_c2 / m.n_bnd).to_numpy()
    ck("T2 frazione di perimetro a contatto in [0.38, 0.50]", np.all((f >= 0.38) & (f <= 0.50)), np.round(f, 3).tolist())
    ck("T2 conflitti di pixel <= numero di righe della corda", ncf <= int(2 * R / upp) + 2, ncf)
    # T3 quadrati separati da g colonne vuote: contatto con delta = 2 sse g <= 1; con delta = 1 sse g = 0
    for g in (0, 1, 2, 3):
        L = np.zeros((40, 60), np.int64); L[10:30, 5:25] = 1; L[10:30, 25 + g:45 + g] = 2
        m = measure(L, upp)
        ck(f"T3 g={g} delta=2", (m.n_c2[0] > 0) == (g <= 1), int(m.n_c2[0]))
        ck(f"T3 g={g} delta=1", (m.n_c1[0] > 0) == (g == 0), int(m.n_c1[0]))
    # T4 forma a L (3 quadrati di lato a): solidita' continua 6/7
    a = 5.0
    Lp = shapely.Polygon([(10, 10), (10 + 2 * a, 10), (10 + 2 * a, 10 + a), (10 + a, 10 + a), (10 + a, 10 + 2 * a), (10, 10 + 2 * a)])
    lab, _ = rasterize([Lp], [1], (120, 120), upp); m = measure(lab, upp)
    ck("T4 solidita' L = 6/7 +- 0.03", abs(m.sol[0] - 6 / 7) <= 0.03, round(float(m.sol[0]), 4))
    # T5 arbitro a cicli su 3 immagini casuali (dischi dipinti): n_bnd, n_c1, n_c2 identici
    rng = np.random.default_rng(20261015)
    for t in range(3):
        L = np.zeros((90, 90), np.int64); yy, xx = np.mgrid[0:90, 0:90]
        for k in range(1, 26):
            cy_, cx_, rr_ = rng.uniform(0, 90), rng.uniform(0, 90), rng.uniform(3, 9)
            L[(yy - cy_) ** 2 + (xx - cx_) ** 2 <= rr_ ** 2] = k
        m = measure(L, upp).set_index("label")
        for d in (1, 2):
            nb, nc = brute(L, d)
            okb = all(m.loc[k, "n_bnd"] == nb.get(k, 0) for k in m.index)
            okc = all(m.loc[k, f"n_c{d}"] == nc.get(k, 0) for k in m.index)
            ck(f"T5 img{t} delta={d} arbitro a cicli (n_bnd, n_c)", okb and okc, int(sum(nc.values())))
    # T6 frac_out: due rettangoli adiacenti di larghezze 10 e 30 px: Voronoi dei centroidi a meta' fra i centroidi
    L = np.zeros((20, 50), np.int64); L[:, 5:15] = 1; L[:, 15:45] = 2
    m = measure(L, upp).set_index("label")
    # centroidi a x = 9.5 e 29.5 px: bisettrice a 19.5 -> pixel 15..19 del rettangolo 2 (5 colonne su 30) piu' vicini a 1
    ck("T6 frac_out analitico (0, 5/30)", abs(m.loc[1, "frac_out"]) < 1e-12 and abs(m.loc[2, "frac_out"] - 5 / 30) < 1e-12,
       [round(float(m.loc[1, "frac_out"]), 4), round(float(m.loc[2, "frac_out"]), 4)])
    # T7 raster di un quadrato: numero di punti interi dentro
    sq = shapely.box(3.3 * upp, 4.6 * upp, 13.3 * upp, 20.6 * upp)
    lab, _ = rasterize([sq], [1], (40, 40), upp)
    ck("T7 quadrato: 10 x 16 pixel", int((lab == 1).sum()) == 160, int((lab == 1).sum()))
    # T8 coppie di contatto: due quadrati a contatto + uno lontano
    L = np.zeros((40, 80), np.int64); L[5:15, 5:15] = 1; L[5:15, 16:26] = 2; L[5:15, 40:50] = 3
    pr = contact_pairs(L, 2)
    ck("T8 coppie di contatto = {(1, 2)}", len(pr) == 1 and pr.label_i[0] == 1 and pr.label_j[0] == 2, pr.values.tolist())
    nf = sum(1 for r in res if not r[1])
    print(f"RESULT {len(res) - nf}/{len(res)} PASS (mutante {MUT})")
    return nf


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "test":
        sys.exit(1 if run_tests() else 0)
    elif a[0] == "sim":
        from multiprocessing import Pool
        os.makedirs(a[2], exist_ok=True); rois = rois_table()
        fs = sorted(glob.glob(os.path.join(a[1], "*.parquet"))); nc = int(a[3]) if len(a) > 3 else 30
        with Pool(nc) as p:
            r = p.starmap(sim_one, [(f, a[2], rois) for f in fs])
        print(f"[sim] {len(fs)} file, {sum(x == 'ok' for x in r)} nuovi")
    elif a[0] == "real":
        os.makedirs(a[1], exist_ok=True); rois = rois_table()
        if len(a) > 3:                                       # un solo ROI (sviluppo): real <out_dir> <A> <roi>
            real_one(a[2], a[3], a[1]); sys.exit(0)
        for _, rw in rois.iterrows():
            real_one(rw.archetype, rw.roi_id, a[1]); print(rw.archetype, rw.roi_id, "done", flush=True)
    elif a[0] == "cp2":
        from multiprocessing import Pool
        os.makedirs(a[1], exist_ok=True); rois = rois_table(); nc = int(a[2]) if len(a) > 2 else 20
        jobs = [(rw.archetype, rw.roi_id, s, a[1], rois) for _, rw in rois.iterrows() for s in ("cellpose_rgb", "stardist_he")]
        with Pool(nc) as p:
            r = p.starmap(cp2_one, jobs)
        print(f"[cp2] {len(jobs)} unita', {sum(x == 'ok' for x in r)} nuove")
    else:
        raise SystemExit("modo sconosciuto")
