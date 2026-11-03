"""tools/R5_ph2_synth.py — R5, controprova di PH-2 (post hoc): su conteggi sintetici generati sulla geometria vera di una ROI
(densita' per cellula osservate come verita', rapporto r_vero, sfocatura gaussiana sigma_vero, scarto della sezione, parita', Poisson; 1 replica),
la stessa procedura di PH-2 (r invertito sulla statistica di R4 PH-1 per ogni sigma della griglia, statistica lontana dal bordo prevista)
ritrova sigma_vero?  Uso: .venv/bin/python tools/R5_ph2_synth.py <outdir> <workers> -> <outdir>/ph2_synth/*.json"""
import sys, os, json, traceback
import numpy as np, pandas as pd
from multiprocessing import Pool
OUT, NW = sys.argv[1], int(sys.argv[2])
sys.argv = ["R5_nc.py", OUT, f"{OUT}/kernels.json"]
sys.path.insert(0, "tools")
import importlib.util
spec = importlib.util.spec_from_file_location("nc", "tools/R5_nc.py"); NC = importlib.util.module_from_spec(spec); spec.loader.exec_module(NC)
import R5_common as C
GRID = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
CASES = [("A2_r1", 0.75), ("A5_r1", 1.5)]
S_TRUE = [0.7, 2.0]

def job(a):
    tag, Rw, r_true, s_true, s_fit, kern = a
    f = f"{OUT}/ph2_synth/{tag}__true{s_true:.1f}__fit{s_fit:.1f}.json"
    if os.path.exists(f): return f
    try:
        Dd, IX, IXf, Dobs, am, acell, S, U, base = NC.prep(tag, Rw)
        di = NC.d_image(Dd, kern)
        kt_true = (s_true, s_true, 0.0, 30.0)
        Ug = np.random.default_rng([C.SEED, 999, int(Rw.name), int(s_true * 10)]).uniform(size=(1, len(Dd["u"])))
        y = C.nc_simulate(S, r_true, Dd["pos"], kt_true, di, Dd["c"], Ug)[0]
        nuc_g = C.labels_to_grid(Dd["mask"], Dd["upp"]); exp_g = C.expand_labels_grid(nuc_g, NC.D_EXP)
        Dy = C.observed_cell_density(Dd["L"], Dd["keep"], y); ay = float(y[am].mean()) if am.any() else 0.0
        Sy = C.nc_setup(nuc_g, exp_g, Dy, float(np.median(list(Dy.values()))), ay)
        o_star, o_far = C.star_fast(y, IX), C.star_fast(y, IXf)
        kt = (s_fit, s_fit, 0.0, 30.0)
        r, Sr = C.invert_r_refined(o_star, Sy, Dd["pos"], kt, di, Dd["c"], U, IX, Dd["L"], Dd["keep"], Dy, ay, am)
        pf = float("nan")
        if np.isfinite(r):
            Y = C.nc_simulate(Sr, r, Dd["pos"], kt, di, Dd["c"], U); pf = float(np.mean([C.star_fast(v, IXf) for v in Y]))
        res = {"tag": tag, "r_true": r_true, "sigma_true": s_true, "sigma_in": s_fit, "r": r, "pred_far": pf, "obs_star": o_star, "obs_far": o_far}
    except Exception:
        res = {"tag": tag, "sigma_true": s_true, "sigma_in": s_fit, "error": traceback.format_exc()}
    json.dump(res, open(f + ".tmp", "w")); os.replace(f + ".tmp", f); return f

if __name__ == "__main__":
    os.makedirs(f"{OUT}/ph2_synth", exist_ok=True)
    K = json.load(open(f"{OUT}/kernels.json")); rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    jobs = []
    for tag, r_true in CASES:
        Rw = rois[(rois.archetype + "_" + rois.roi_id) == tag].iloc[0]
        for st in S_TRUE:
            for sf in GRID: jobs.append((tag, Rw, r_true, st, sf, K[NC.KEY[Rw.archetype]]["K2"]))
    with Pool(NW) as p:
        for i, f in enumerate(p.imap(job, jobs, chunksize=6), 1): print(i, len(jobs), f, flush=True)
    print("done", flush=True)
