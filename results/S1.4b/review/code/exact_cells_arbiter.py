#!/usr/bin/env python
"""review/code/exact_cells_arbiter.py — revisore "code": arbitro ESATTO (non campionario) di C-L2 sui 30 insiemi reali.
Cella_i = box ritagliato dai semipiani dei K vicini euclidei; certificato: per ogni vertice v, min_j pow_j(v) >= pow_i(v) - eps
su TUTTI i j (query a raggio sqrt(pow_i(v) + W) nel KD-tree dei generatori); se fallisce, si aggiungono i violatori e si
ripete. Una cella convessa che soddisfa tutti i semipiani ai vertici li soddisfa ovunque. Confronto con le celle di R
(griglia 1e-6): |dA|, differenza simmetrica, celle vuote."""
import glob, os, numpy as np, pandas as pd, shapely
from scipy.spatial import cKDTree
from multiprocessing import Pool
def clip(P, a, b):
    if len(P) == 0: return P
    s = P @ a - b; ins = s <= 0
    if ins.all(): return P
    if not ins.any(): return P[:0]
    out = []; m = len(P)
    for i in range(m):
        j = (i + 1) % m
        if ins[i]: out.append(P[i])
        if ins[i] != ins[j]: t = s[i] / (s[i] - s[j]); out.append(P[i] + t * (P[j] - P[i]))
    return np.array(out)
def run(f):
    d = pd.read_parquet(f); reg = shapely.from_wkb(open(f.replace(".parquet", "_region.wkb"), "rb").read()); cells = shapely.from_wkb(d.wkb.to_numpy())
    x, y, w = d.x.to_numpy(), d.y.to_numpy(), d.w.to_numpy(); n = len(x); mx, my = x.mean(), y.mean(); xc, yc = x - mx, y - my; q = xc**2 + yc**2 - w; W = w.max()
    bb = reg.bounds; B = np.array([[bb[0] - 5, bb[1] - 5], [bb[2] + 5, bb[1] - 5], [bb[2] + 5, bb[3] + 5], [bb[0] - 5, bb[3] + 5]]) - [mx, my]
    tree = cKDTree(np.column_stack([xc, yc])); _, nn = tree.query(np.column_stack([xc, yc]), k=min(41, n))
    polys = []; n_extra = 0
    for i in range(n):
        js = [j for j in nn[i] if j != i]; done = set(js)
        while True:
            P = B.copy()
            for j in js:
                P = clip(P, 2 * np.array([xc[j] - xc[i], yc[j] - yc[i]]), q[j] - q[i])
                if len(P) == 0: break
            if len(P) == 0: break
            pi = (P[:, 0] - xc[i])**2 + (P[:, 1] - yc[i])**2 - w[i]
            viol = set()
            for v, pv in zip(P, pi):
                for j in tree.query_ball_point(v, np.sqrt(max(pv + W, 0)) + 1e-9):
                    if j != i and (v[0] - xc[j])**2 + (v[1] - yc[j])**2 - w[j] < pv - 1e-9: viol.add(j)
            viol -= done
            if not viol: break
            n_extra += len(viol); js += sorted(viol); done |= viol
        polys.append(shapely.Polygon(P + [mx, my]) if len(P) >= 3 else shapely.Polygon())
    ex = shapely.intersection(np.array(polys, dtype=object), reg)
    aR, aE = shapely.area(cells), shapely.area(ex); sym = shapely.area(shapely.symmetric_difference(cells, ex)); per = shapely.length(ex)
    return dict(set=os.path.basename(f)[:-8], n=n, n_extra_halfplanes=n_extra, empty_R=int((aR <= 0).sum()), empty_exact=int((aE <= 0).sum()),
                max_abs_dA=float(np.max(np.abs(aR - aE))), max_dA_over_perim=float(np.max(np.abs(aR - aE) / np.maximum(per, 1e-12))),
                max_symdiff=float(sym.max()), n_symdiff_gt_perim_1e6=int((sym > per * 1e-6 + 1e-12).sum()), sum_rel=float(abs(aR.sum() - aE.sum()) / aE.sum()))
if __name__ == "__main__":
    fs = sorted(glob.glob("/mnt/micron/geo_spatialtrans/S1.4b/arb_in/*.parquet"))
    with Pool(30) as p: R = p.map(run, fs)
    r = pd.DataFrame(R); r.to_csv("results/S1.4b/review/code/exact_cells_arbiter.csv", index=False); pd.set_option("display.width", 220); print(r.to_string())
