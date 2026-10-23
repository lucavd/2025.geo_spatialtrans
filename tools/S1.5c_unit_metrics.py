"""tools/S1.5c_unit_metrics.py - S1.5c: metriche spaziali per unita' simulata (nn_rel, markconn, cv_local4 con tools/S1.5c_metrics.py, le stesse
funzioni dei riferimenti reali) e descrittiva D-desc-3 (distanza dal vicino / somma dei raggi). Ingresso: /mnt/micron/geo_spatialtrans/S1.5c/cells/<braccio>/*.parquet
(scritti da tools/S1.5c_bench.R gen). Uscita: un csv per unita' in .../S1.5c/umet/<braccio>/ (rilanciabile) e results/S1.5c/S1.5c_unit_metrics.csv.
Uso: .venv/bin/python tools/S1.5c_unit_metrics.py [worker]"""
import sys, os, glob, importlib.util, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial import cKDTree
s = importlib.util.spec_from_file_location("met", "tools/S1.5c_metrics.py"); met = importlib.util.module_from_spec(s); s.loader.exec_module(met)
OUT = "/mnt/micron/geo_spatialtrans/S1.5c"; RES = "results/S1.5c"
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
def one(f):
    arm = f.split("/")[-2]; A, roi, sd = os.path.basename(f)[:-8].split("_")
    o = f"{OUT}/umet/{arm}/{A}_{roi}_{sd}.csv"
    if os.path.exists(o): return o
    e = pd.read_parquet(f); m = met.roi_metrics(e.x.values, e.y.values, e.cell_type.values, e.interior.values, MAIN[A])
    rows = []
    for k, v in m.items():
        if isinstance(k, tuple): rows.append(dict(metric=k[0], cell_type=k[1], bin=k[2] if len(k) > 2 else "", value=v))
        else: rows.append(dict(metric=k, cell_type="", bin="", value=v))
    P = e[["x", "y"]].values; dd, ii = cKDTree(P).query(P, k=2); r = np.sqrt(e.a_n.values / np.pi); g = dd[:, 1] / (r + r[ii[:, 1]]); I = e.interior.values
    rows += [dict(metric="gapratio_q05", cell_type="", bin="", value=np.quantile(g[I], 0.05)), dict(metric="gapratio_lt1", cell_type="", bin="", value=np.mean(g[I] < 1))]
    d = pd.DataFrame(rows).assign(arm=arm, archetype=A, roi_id=roi, seed=int(sd[1:])); os.makedirs(os.path.dirname(o), exist_ok=True)
    d.to_csv(o + ".tmp", index=False); os.replace(o + ".tmp", o); return o
if __name__ == "__main__":
    fs = sorted(glob.glob(f"{OUT}/cells/*/*.parquet")); nw = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    with Pool(nw) as p: outs = p.map(one, fs, chunksize=4)
    pd.concat([pd.read_csv(o) for o in outs]).to_csv(f"{RES}/S1.5c_unit_metrics.csv", index=False); print(len(outs), "unita'")
