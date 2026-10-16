#!/usr/bin/env python
"""review/code/02_obs_recompute.py — ricalcolo indipendente dello stadio obs (B-1, B-2, B-7, B-8, S-1, S-2, S-3 per file).

Dai soli prodotti del raster (sim_out, real_shape) e dalle celle reali esportate (01_export.R): territori osservati
ricostruiti con shapely (Voronoi dei centroidi raster ritagliato alla regione del ROI), cerchi di pari area con cKDTree,
KS con codice proprio. Confronta con S1.4c_obs_raw.csv. Uso (root del repo):
  .venv/bin/python results/S1.4c/review/code/02_obs_recompute.py [ncores]
"""
import os, sys, glob
import numpy as np, pandas as pd, shapely
from shapely import wkt
from scipy.spatial import cKDTree
from multiprocessing import Pool

S = "/mnt/micron/geo_spatialtrans/S1.4c"; SC = "/mnt/micron/geo_spatialtrans/S1.4c_review/code"; EXP = SC + "/exp"
SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}

def ks(a, b):
    a = np.sort(np.asarray(a, float)); b = np.sort(np.asarray(b, float)); v = np.unique(np.concatenate([a, b]))
    return float(np.max(np.abs(np.searchsorted(a, v, "right") / len(a) - np.searchsorted(b, v, "right") / len(b))))

def circle_hits(x, y, r):
    """hit[i] = il cerchio i interseca (R_i + R_j > d) il cerchio di un'altra cella."""
    t = cKDTree(np.column_stack([x, y])); hit = np.zeros(len(x), bool); rmax = r.max()
    pairs = t.query_pairs(2 * rmax, output_type="ndarray")
    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        d = np.hypot(x[i] - x[j], y[i] - y[j]); bad = r[i] + r[j] > d
        hit[i[bad]] = True; hit[j[bad]] = True
    return hit

def voronoi_obs(cx, cy, region):
    """a_t, d (centroide -> bordo del proprio territorio) per ogni punto; NaN se fuori dalla regione."""
    n = len(cx); a_t = np.full(n, np.nan); d = np.full(n, np.nan)
    ok = np.isfinite(cx) & np.isfinite(cy)
    pts = shapely.points(cx[ok], cy[ok]); inside = shapely.intersects(pts, region)
    id_ = np.nonzero(ok)[0][inside]; P = shapely.points(cx[id_], cy[id_])
    vor = shapely.voronoi_polygons(shapely.multipoints(P), extend_to=shapely.buffer(shapely.envelope(region), 1000.0))
    cells = np.asarray(shapely.get_parts(vor))
    # assegna ogni cella di Voronoi al suo generatore (il punto contenuto)
    tree = shapely.STRtree(cells); pi, ci = tree.query(P, predicate="intersects")
    own = np.full(len(P), -1)
    for p_, c_ in zip(pi, ci):
        if own[p_] < 0 or shapely.contains(cells[c_], P[p_]): own[p_] = c_
    T = shapely.intersection(cells[own], region)
    a_t[id_] = shapely.area(T); d[id_] = shapely.distance(P, shapely.boundary(T))
    return a_t, d, int(np.sum(~inside))

_REG = None
def regions():
    global _REG
    if _REG is None:
        r = pd.read_csv(EXP + "/regions.csv"); _REG = {(a, b): wkt.loads(w) for a, b, w in zip(r.archetype, r.roi, r.wkt)}
    return _REG

def real_ref(A, roi):
    cc = pd.read_parquet(f"{EXP}/cells_{A}_{roi}.parquet"); ev = cc.ev_prep.to_numpy()
    rs = pd.read_parquet(f"{S}/real_shape/{A}_{roi}_{SEG[A]}_shape.parquet").set_index("label")
    lab = cc.shape_label.to_numpy()[ev]; z = rs.loc[lab]
    pr = pd.read_csv(EXP + "/prep.csv"); pr = pr[(pr.archetype == A) & (pr.roi == roi)].iloc[0]
    return dict(cc=cc, ev=ev, nc=(cc.a_n / cc.a_t).to_numpy()[ev], an=cc.a_n.to_numpy()[ev], s1=float(np.mean(z.n_c2 > 0)),
                s2=(z.n_c2 / z.n_bnd).to_numpy(), s3=z.sol.to_numpy(), b7=float(pr.b7_real), b8=float(pr.b8_real))

def one(f):
    b = os.path.basename(f)[:-8]; A, roi, fold, s, V = b.split("_")
    R = real_ref(A, roi); n = len(R["cc"]); ev = R["ev"]
    m = pd.read_parquet(f)
    def v(col):
        z = np.full(n, np.nan); z[m.label.to_numpy() - 1] = m[col].to_numpy(); return z
    cx, cy, area, sol, nb, c2 = v("cx"), v("cy"), v("area"), v("sol"), v("n_bnd"), v("n_c2")
    a_t, d, nexc = voronoi_obs(cx, cy, regions()[(A, roi)])
    e = ev & np.isfinite(cx) & np.isfinite(a_t)
    ok = np.isfinite(cx) & np.isfinite(area) & (area > 0); idx = np.nonzero(ok)[0]
    hit = np.zeros(n, bool); hit[idx] = circle_hits(cx[idx], cy[idx], np.sqrt(area[idx] / np.pi))
    b7 = float(np.mean(np.sqrt(area[e] / np.pi) > d[e])); b8 = float(np.mean(hit[ev & ok]))
    return dict(archetype=A, roi=roi, fold=int(fold[4:]), seed_idx=int(s[1:]), variant={"Kc": "K-c", "K1": "K-1", "K2": "K-2"}[V],
                n_eval=int(e.sum()), n_excluded=nexc, b7=b7, b8=b8,
                B1=ks(area[e] / a_t[e], R["nc"]), B2=ks(area[e], R["an"]), B7=abs(b7 - R["b7"]), B8=abs(b8 - R["b8"]),
                S1=abs(float(np.mean(c2[e] > 0)) - R["s1"]), S2=ks(c2[e] / nb[e], R["s2"]), S3=ks(sol[e], R["s3"]))

if __name__ == "__main__":
    nc = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    fs = sorted(glob.glob(f"{S}/sim_out/*_s[01][0-9]_K[c12].parquet"))
    fs = [f for f in fs if "_s00_" not in f]
    with Pool(nc) as p:
        rows = p.map(one, fs, chunksize=1)
    out = pd.DataFrame(rows); out.to_csv(SC + "/02_obs_recomputed.csv", index=False)
    ref = pd.read_csv("results/S1.4c/S1.4c_obs_raw.csv")
    k = ["archetype", "roi", "fold", "seed_idx", "variant"]
    j = out.merge(ref, on=k, suffixes=("_rv", "_rep"))
    lines = [f"file ricalcolati {len(out)}, appaiati con obs_raw {len(j)}"]
    for col in ["n_eval", "n_excluded", "b7", "b8", "B1", "B2", "B7", "B8", "S1", "S2", "S3"]:
        dd = np.abs(j[col + "_rv"] - j[col + "_rep"])
        lines.append(f"{col}: max |diff| = {dd.max():.3e}; n(diff > 1e-9) = {(dd > 1e-9).sum()}")
    open(SC + "/02_obs_recompute_summary.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
