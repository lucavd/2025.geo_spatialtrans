"""tools/S1.5d_unit_metrics.py - S1.5d: metriche per unita' simulata. Stesse funzioni dei riferimenti veri: tools/S1.5c_metrics.py (nn_rel, markconn, cv_local4,
P-5...P-7) e tools/S1.5d_metrics.py (H_rel, xcv2, nn_cv, frac_gt3: P-11...P-13, D-desc-2), con le celle interne della regola di R3 scritte dal banco; P-14:
90deg percentile dell'area del nucleo delle celle interne per tipo principale; D-desc-3 (distanza dal vicino / somma dei raggi). Ingresso:
/mnt/micron/geo_spatialtrans/S1.5d/cells/<braccio>/*.parquet (tools/S1.5d_bench.R gen). Uscita: un csv per unita' in .../S1.5d/umet/<braccio>/
(rilanciabile) e results/S1.5d/S1.5d_unit_metrics.csv. Uso: .venv/bin/python tools/S1.5d_unit_metrics.py [worker] [bracci separati da virgole]"""
import sys, os, glob, importlib.util, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial import cKDTree
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
m5c = mod("m5c", "tools/S1.5c_metrics.py"); m5d = mod("m5d", "tools/S1.5d_metrics.py")
OUT = "/mnt/micron/geo_spatialtrans/S1.5d"; RES = "results/S1.5d"
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
GEOM = {}
def one(f):
    arm = f.split("/")[-2]; A, roi, sd = os.path.basename(f)[:-8].split("_"); o = f"{OUT}/umet/{arm}/{A}_{roi}_{sd}.csv"
    if os.path.exists(o): return o
    e = pd.read_parquet(f); I = e.interior.values.astype(bool)
    rows = []
    for k, v in m5c.roi_metrics(e.x.values, e.y.values, e.cell_type.values, I, MAIN[A]).items():
        rows.append(dict(metric=k[0], cell_type=k[1], bin=k[2] if len(k) > 2 else "", value=v) if isinstance(k, tuple) else dict(metric=k, cell_type="", bin="", value=v))
    g = GEOM.get((A, roi)) or GEOM.setdefault((A, roi), m5d.RoiGeom(A, roi))
    for k, v in m5d.position_metrics(e.x.values, e.y.values, I, e.a_t.values, g).items(): rows.append(dict(metric=k, cell_type="", bin="", value=v))
    for t in MAIN[A]:
        k = I & (e.cell_type.values == t)
        if k.sum(): rows.append(dict(metric="q90_an", cell_type=t, bin="", value=float(np.quantile(e.a_n.values[k], 0.9))))
    P = e[["x", "y"]].values; dd, ii = cKDTree(P).query(P, k=2); r = np.sqrt(e.a_n.values / np.pi); gr = dd[:, 1] / (r + r[ii[:, 1]])
    rows += [dict(metric="gapratio_q05", cell_type="", bin="", value=np.quantile(gr[I], 0.05)), dict(metric="gapratio_lt1", cell_type="", bin="", value=np.mean(gr[I] < 1))]
    d = pd.DataFrame(rows).assign(arm=arm, archetype=A, roi_id=roi, seed=int(sd[1:])); os.makedirs(os.path.dirname(o), exist_ok=True)
    d.to_csv(o + ".tmp", index=False); os.replace(o + ".tmp", o); return o
if __name__ == "__main__":
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    arms = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    fs = sorted(glob.glob(f"{OUT}/cells/*/*.parquet")); fs = [f for f in fs if arms is None or f.split("/")[-2] in arms]
    with Pool(nw) as p: p.map(one, fs, chunksize=4)
    allf = sorted(glob.glob(f"{OUT}/umet/*/*.csv"))
    pd.concat([pd.read_csv(o, keep_default_na=False, na_values=[""]) for o in allf]).to_csv(f"{RES}/S1.5d_unit_metrics.csv", index=False); print(len(allf), "unita'")
