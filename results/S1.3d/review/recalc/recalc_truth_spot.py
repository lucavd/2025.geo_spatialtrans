# results/S1.3d/review/recalc/recalc_truth_spot.py — revisione «ricalcolo» S1.3d: controllo a campione della verità dell'arbitro con codice proprio.
# Per i 24 nulli di K: (a) area esatta (Fraction, shoelace) di ogni cella contro truth/<key>_o1_cells.parquet (a);
# (b) identità del json: union = sum_a - sum_ov; symd = (union - inside) + (A_reg - inside); A_reg = area esatta della regione;
# (c) coppie violate da V0 (res viol o1_V0) con entrambe le celle convesse (1 parte, 1 anello): intersezione esatta Sutherland-Hodgman in Fraction
#     contro truth/<key>_o1_pairs.parquet (a). Uso: .venv/bin/python recalc_truth_spot.py [n_proc]
import sys, json, os
from fractions import Fraction as Fr
from multiprocessing import Pool
import pyarrow.parquet as pq, numpy as np
D = "/mnt/micron/geo_spatialtrans/S1.3d"
def rings(path):
    t = pq.read_table(path).to_pandas().sort_values(["cell", "part", "ring", "k"])
    out = {}
    for (c, p, r), d in t.groupby(["cell", "part", "ring"], sort=True):
        out.setdefault(int(c), {}).setdefault(int(p), []).append([(Fr(float(x)), Fr(float(y))) for x, y in zip(d.x.to_numpy(), d.y.to_numpy())])
    return out
def sh2(R):  # 2 * area con segno
    n = len(R); return sum(R[q][0] * R[(q + 1) % n][1] - R[(q + 1) % n][0] * R[q][1] for q in range(n))
def area(cell):
    return sum(abs(sh2(rr[0])) - sum(abs(sh2(h)) for h in rr[1:]) for rr in cell.values()) / 2
def convex_ccw(cell):
    if len(cell) != 1 or len(list(cell.values())[0]) != 1: return None
    R = list(cell.values())[0][0]
    if sh2(R) < 0: R = R[::-1]
    n = len(R)
    for q in range(n):
        a, b, c = R[q], R[(q + 1) % n], R[(q + 2) % n]
        if (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]) < 0: return None
    return R
def clip(P, C):  # Sutherland-Hodgman: P qualsiasi (qui convesso), C convesso CCW, aritmetica esatta
    out = P
    n = len(C)
    for q in range(n):
        a, b = C[q], C[(q + 1) % n]
        side = lambda p: (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
        inp = out; out = []
        if not inp: break
        for r in range(len(inp)):
            cur, prv = inp[r], inp[r - 1]; sc, sp = side(cur), side(prv)
            if sc >= 0:
                if sp < 0: t = sp / (sp - sc); out.append((prv[0] + t * (cur[0] - prv[0]), prv[1] + t * (cur[1] - prv[1])))
                out.append(cur)
            elif sp >= 0:
                t = sp / (sp - sc); out.append((prv[0] + t * (cur[0] - prv[0]), prv[1] + t * (cur[1] - prv[1])))
    return out
def one(key):
    import subprocess
    cells = rings(f"{D}/geom/{key}_cells.parquet"); reg = rings(f"{D}/geom/{key}_region.parquet")
    tc = pq.read_table(f"{D}/truth/{key}_o1_cells.parquet").to_pandas(); tp = pq.read_table(f"{D}/truth/{key}_o1_pairs.parquet").to_pandas()
    js = json.load(open(f"{D}/truth/{key}_o1.json"))["o1"]
    ex = {c: area(g) for c, g in cells.items()}
    a_t = dict(zip(tc.cell.astype(int), tc.a))
    dc = [abs(float(ex[c]) - a_t[c]) for c in ex]; rel = [abs(float(ex[c]) - a_t[c]) / max(abs(float(ex[c])), 1e-300) for c in ex]
    A_reg = float(area(list(reg.values())[0]))
    idU = abs(js["union"] - (js["sum_a"] - js["sum_ov"])); idS = abs(js["symd"] - ((js["union"] - js["inside"]) + (js["A_reg"] - js["inside"])))
    res = dict(key=key, n_cells=len(ex), n_tc=len(tc), cell_max_abs=max(dc), cell_max_rel=max(rel), cells_float_equal=sum(d == 0 for d in dc),
               A_reg_rel=abs(A_reg - js["A_reg"]) / A_reg, sum_a_rel=abs(float(sum(ex.values())) - js["sum_a"]) / js["sum_a"], id_union=idU, id_symd=idS)
    vp = os.path.join(D, "review_recalc", "v0viol", key + ".csv")
    pairs = []
    if os.path.exists(vp):
        for line in open(vp).read().split("\n")[1:]:
            if line.strip(): i, j = map(int, line.replace('"', '').split(",")); pairs.append((i, j))
    tpd = {(int(i), int(j)): a for i, j, a in zip(tp.i, tp.j, tp.a)}
    nconv = 0; worst = 0.0; nchk = 0
    for i, j in pairs:
        A = convex_ccw(cells[i]); B = convex_ccw(cells[j])
        if A is None or B is None: continue
        nconv += 1; I = clip(A, B); a = float(abs(sh2(I)) / 2) if len(I) >= 3 else 0.0
        t = tpd[(min(i, j), max(i, j))]; worst = max(worst, abs(a - t)); nchk += 1
    res.update(v0_viol_pairs=len(pairs), convex_checked=nchk, convex_max_abs_diff=worst)
    return res
if __name__ == "__main__":
    keys = [k.strip() for k in open(f"{D}/keys_K.txt") if k.strip()]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 2) as p: R = p.map(one, keys)
    import pandas as pd
    df = pd.DataFrame(R); df.to_csv(f"{D}/review_recalc/truth_spot.csv", index=False)
    pd.set_option("display.width", 250); print(df.to_string())
    print("TOT celle", df.n_cells.sum(), "max |exact - truth| celle", df.cell_max_abs.max(), "max rel", df.cell_max_rel.max(), "| identità union max", df.id_union.max(), "symd max", df.id_symd.max(),
          "| coppie V0 violate", df.v0_viol_pairs.sum(), "convesse controllate", df.convex_checked.sum(), "max |S-H esatto - truth|", df.convex_max_abs_diff.max())
