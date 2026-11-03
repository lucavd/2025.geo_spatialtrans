"""tools/R5_ph3_reach.py — R5, post hoc PH-3 (dopo la revisione, rilievi F07, C12): intervallo raggiungibile dalla statistica di R4 PH-1 con il nucleo
dei bordi (K2 mediano della sezione) al variare di r in [0.05, 50], per ogni ROI; dice se l'inversione fallisce perche' l'osservato e' fuori portata.
Uso: .venv/bin/python tools/R5_ph3_reach.py <outdir> <workers>. Scrive <outdir>/tables/R5_PH3_reach.csv."""
import sys, os, json
import numpy as np, pandas as pd
from multiprocessing import Pool
OUT, NW = sys.argv[1], int(sys.argv[2])
sys.argv = ["R5_nc.py", OUT, f"{OUT}/kernels.json"]
sys.path.insert(0, "tools")
import importlib.util
spec = importlib.util.spec_from_file_location("nc", "tools/R5_nc.py"); NC = importlib.util.module_from_spec(spec); spec.loader.exec_module(NC)
import R5_common as C
def job(a):
    tag, Rw, kern = a
    Dd, IX, IXf, Dobs, am, acell, S, U, base = NC.prep(tag, Rw)
    kt = (kern["sigma_r"], kern["sigma_c"], kern["w"], kern["ell"]); di = NC.d_image(Dd, kern)
    st = {}
    for r in (0.05, 1.0, 50.0):
        Y = C.nc_simulate(S, r, Dd["pos"], kt, di, Dd["c"], U); st[r] = float(np.mean([C.star_fast(y, IX) for y in Y]))
    return {"tag": tag, "archetype": base["archetype"], "obs_star": base["obs_star"], "stat_r0.05": st[0.05], "stat_r1": st[1.0], "stat_r50": st[50.0],
            "reachable": bool(min(st[0.05], st[50.0]) <= base["obs_star"] <= max(st[0.05], st[50.0]))}
if __name__ == "__main__":
    K = json.load(open(f"{OUT}/kernels.json")); rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    jobs = [(f"{Rw.archetype}_{Rw.roi_id}", Rw, K[NC.KEY[Rw.archetype]]["K2"]) for _, Rw in rois.iterrows()]
    with Pool(NW) as p: D = pd.DataFrame(p.map(job, jobs, chunksize=1))
    D.to_csv(f"{OUT}/tables/R5_PH3_reach.csv", index=False); print(D.round(3).to_string())
