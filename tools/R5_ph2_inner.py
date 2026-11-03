"""tools/R5_ph2_inner.py — R5, post hoc PH-2 (dichiarato dopo CP-1): larghezza della sfocatura DENTRO il tessuto stimata dal contrasto nucleo/corona.
Per ogni ROI e per sigma_in in una griglia (gaussiana sola, w = 0, scarto della sezione), si inverte r sul rapporto di P-5 e si prevede il rapporto
senza i bin entro 2 µm dal bordo del nucleo (CP-1). sigma_in* = valore per cui previsto = osservato (interpolazione lineare sulla griglia).
Uso: .venv/bin/python tools/R5_ph2_inner.py <outdir> <workers>. Scrive <outdir>/ph2/<tag>__s<sigma>.json."""
import sys, os, json, time, traceback
import numpy as np, pandas as pd
from multiprocessing import Pool
OUT, NW = sys.argv[1], int(sys.argv[2])
sys.argv = ["R5_nc.py", OUT, f"{OUT}/kernels.json"]
sys.path.insert(0, "tools")
import importlib.util
spec = importlib.util.spec_from_file_location("nc", "tools/R5_nc.py"); NC = importlib.util.module_from_spec(spec); spec.loader.exec_module(NC)
import R5_common as C
GRID = [float(x) for x in os.environ.get("R5_PH2_GRID", "0.5,1,1.5,2,3,4").split(",")]
D0 = os.environ.get("R5_PH2_DELTA0") == "1"            # sensibilita': scarto della sezione posto a 0
TAIL = os.environ.get("R5_PH2_TAIL") == "1"            # variante: coda della sezione (w, ell di K2) mantenuta, varia solo sigma
SUBDIR = ("ph2_d0" if D0 else "ph2") + ("_tail" if TAIL else "")
TAGS = [t for t in os.environ.get("R5_PH2_TAGS", "").split(",") if t]

def job(a):
    tag, Rw, s, kern = a
    f = f"{OUT}/{SUBDIR}/{tag}__s{s:.2f}.json" if s < 0.5 else f"{OUT}/{SUBDIR}/{tag}__s{s:.1f}.json"
    if os.path.exists(f): return f
    try:
        Dd, IX, IXf, Dobs, am, acell, S, U, base = NC.prep(tag, Rw)
        k = dict(kern); k.update({"sigma_r": s, "sigma_c": s, "w": kern["w"] if TAIL else 0.0})
        if D0: k.update({"dr": 0.0, "dc": 0.0})
        kt = (s, s, k["w"], kern["ell"]) if TAIL else (s, s, 0.0, 30.0); di = NC.d_image(Dd, k)
        r, Sr = C.invert_r_refined(base["obs_star"], S, Dd["pos"], kt, di, Dd["c"], U, IX, Dd["L"], Dd["keep"], Dobs, acell, am)
        pf = float("nan")
        if np.isfinite(r):
            Y = C.nc_simulate(Sr, r, Dd["pos"], kt, di, Dd["c"], U); pf = float(np.mean([C.star_fast(y, IXf) for y in Y]))
        res = {"tag": tag, "archetype": base["archetype"], "sigma_in": s, "r": r, "pred_far": pf, "obs_star": base["obs_star"], "obs_far": base["obs_star_far"], "delta0": D0, "tail_w": k["w"], "tail_ell": kern["ell"] if TAIL else None}
    except Exception:
        res = {"tag": tag, "sigma_in": s, "error": traceback.format_exc()}
    json.dump(res, open(f + ".tmp", "w")); os.replace(f + ".tmp", f); return f

if __name__ == "__main__":
    os.makedirs(f"{OUT}/{SUBDIR}", exist_ok=True)
    K = json.load(open(f"{OUT}/kernels.json")); rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    jobs = [(f"{Rw.archetype}_{Rw.roi_id}", Rw, s, K[NC.KEY[Rw.archetype]]["K2"]) for _, Rw in rois.iterrows() for s in GRID if not TAGS or f"{Rw.archetype}_{Rw.roi_id}" in TAGS]
    with Pool(NW) as p:
        for i, f in enumerate(p.imap(job, jobs, chunksize=len(GRID)), 1):
            if i % 6 == 0: print(i, len(jobs), flush=True)
    print("done", flush=True)
