"""tools/S1.5d_real_refs.py - S1.5d: valori reali per ROI delle metriche della disposizione (tools/S1.5d_metrics.py), calcolati PRIMA di ogni
simulazione; con gli stessi pixel e quadrati, i valori dei bracci U e N-a di S1.5c (10 semi), per gli errori standard tipici.
Nuclei veri con territorio (etichette R6 PH-1, insieme di S1.5b/S1.5c). Uso: .venv/bin/python tools/S1.5d_real_refs.py [n_proc]
Uscite: results/S1.5d/S1.5d_real_refs.csv (ROI x metrica), S1.5d_real_intervals.csv ([min, max] dei 5 ROI), S1.5d_ref_arms.csv (U, N-a)."""
import sys, importlib.util, numpy as np, pandas as pd
from multiprocessing import Pool
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py"); met = mod("met", "tools/S1.5d_metrics.py")
OUT5C = "/mnt/micron/geo_spatialtrans/S1.5c/cells"
def one(args):
    A, r = args; g = met.RoiGeom(A, r); _, d = cal.load(A); e = d[d.roi_id == r]
    rows = [dict(archetype=A, roi_id=r, source="real", seed=0, **met.position_metrics(e.x.values, e.y.values, e.interior.values, e.a_t.values, g))]
    for arm in ["U", "N-a"]:
        for s in range(1, 11):
            c = pd.read_parquet(f"{OUT5C}/{arm}/{A}_{r}_s{s:02d}.parquet")
            rows.append(dict(archetype=A, roi_id=r, source=arm, seed=s, **met.position_metrics(c.x.values, c.y.values, c.interior.values, c.a_t.values, g)))
    return rows
if __name__ == "__main__":
    R = met.rois_table(); jobs = list(zip(R.archetype, R.roi_id))
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as p: out = [x for o in p.map(one, jobs) for x in o]
    D = pd.DataFrame(out); real = D[D.source == "real"].drop(columns=["source", "seed"])
    real.to_csv("results/S1.5d/S1.5d_real_refs.csv", index=False)
    L = real.melt(id_vars=["archetype", "roi_id"], var_name="metric")
    L.groupby(["archetype", "metric"]).value.agg(["min", "max", "mean", "count"]).reset_index().to_csv("results/S1.5d/S1.5d_real_intervals.csv", index=False)
    D[D.source != "real"].to_csv("results/S1.5d/S1.5d_ref_arms.csv", index=False); print("ok", len(D))
