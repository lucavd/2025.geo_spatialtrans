"""tools/S2.0_inner.py — S2.0, M2: controllo interno (P-3) di ogni regola nelle 30 ROI. Come la CP-1 di R5 (tools/R5_nc.py: r invertito sulla
statistica di R4 PH-1, statistica senza i bin entro 2 µm dal bordo del nucleo prevista), con la sfocatura della regola (parametri mediani della
sezione, <params.json>) al posto del nucleo di R5; per C2 e C3 la maschera del tessuto della ROI. 5 repliche, uniformi con seme 20261103.
In piu' 'C0_r5u' = C0 con le uniformi di R5 (deve riprodurre la CP-1 di R5). Uso: .venv/bin/python tools/S2.0_inner.py <outdir> <params.json>
-> <outdir>/inner/<tag>__<rule>.json (processi dal file <outdir>/WORKERS)."""
import sys, os, json, time, traceback
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
sys.argv_saved = list(sys.argv)
sys.argv = ["R5_nc.py", R5RUN, f"{R5RUN}/kernels.json"]
import importlib.util
_spec = importlib.util.spec_from_file_location("r5nc", "tools/R5_nc.py"); NC = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(NC)
sys.argv = sys.argv_saved
import R5_common as C
import S2_0_common as S
import R5_tiles as T5
from R2_common import bubble_mask, downsample
from R1_verify import read_region

def roi_tissue(Rw, Dd):
    """Tessuto H&E della ROI (maschera di R5: R2 unita a densita' ottica > 0.15, bolle escluse) sulla griglia di 0.5 µm del campo F."""
    p, H, W, upp_img = T5.image_info(NC.DS[Rw.archetype])
    rgb = read_region(p, int(Rw.r0), int(Rw.r1), int(Rw.c0), int(Rw.c1))
    assert rgb.shape[:2] == Dd["mask"].shape, (rgb.shape, Dd["mask"].shape)
    rgb2 = downsample(rgb, 2); u2 = Dd["upp"] * 2
    tm = T5.tissue_mask_r5(rgb2, u2) & ~bubble_mask(rgb2, u2)
    H0, W0 = Dd["mask"].shape; step = S.SUB_UM
    n0, n1 = int(np.floor(H0 * Dd["upp"] / step)), int(np.floor(W0 * Dd["upp"] / step))
    ir = np.minimum(((np.arange(n0) + 0.5) * step / Dd["upp"] / 2).astype(int), tm.shape[0] - 1)
    ic = np.minimum(((np.arange(n1) + 0.5) * step / Dd["upp"] / 2).astype(int), tm.shape[1] - 1)
    return tm[np.ix_(ir, ic)].astype(float)

def rule_prm(PJ, A, rule):
    return dict(PJ[NC.KEY[A]][rule])

OUT, PJF = sys.argv[1], sys.argv[2]
RULES = ("C0", "C0_r5u", "C1", "C2", "C3")

def job(args):
    tag, Rw, rule = args
    fo = f"{OUT}/inner/{tag}__{rule}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        PJ = json.load(open(PJF)); rl = "C0" if rule == "C0_r5u" else rule
        prm = rule_prm(PJ, Rw.archetype, rl)
        Dd, IX, IXf, Dobs, am, acell, S_, U5, base = NC.prep(tag, Rw)
        U = U5 if rule == "C0_r5u" else np.random.default_rng([S.SEED, int(Rw.name)]).uniform(size=U5.shape)
        di = NC.d_image(Dd, prm)
        fm = roi_tissue(Rw, Dd) if rl in ("C2", "C3") else None
        blur = S.make_blur(rl, prm, Dd["pos"], di, fm)
        r, Sr = S.invert_r_refined_b(base["obs_star"], S_, blur, Dd["c"], U, IX, Dd["L"], Dd["keep"], Dobs, acell, am)
        pf = float("nan")
        if np.isfinite(r):
            Y = S.nc_simulate_b(Sr, r, blur, Dd["c"], U); pf = float(np.mean([C.star_fast(y, IXf) for y in Y]))
        res = {"tag": tag, "archetype": Rw.archetype, "rule": rule, "r": r, "pred_far": pf, "obs_star": base["obs_star"],
               "obs_far": base["obs_star_far"], "ratio_far": pf / base["obs_star_far"] if np.isfinite(pf) else float("nan"),
               "prm": prm, "d_img": list(map(float, di)), "tissue_frac": float(fm.mean()) if fm is not None else None, "sec": round(time.time() - t0, 1)}
    except Exception:
        res = {"tag": tag, "rule": rule, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

if __name__ == "__main__":
    os.makedirs(f"{OUT}/inner", exist_ok=True)
    rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    only = [t for t in os.environ.get("S20_TAGS", "").split(",") if t]
    jobs = [(f"{Rw.archetype}_{Rw.roi_id}", Rw, rl) for _, Rw in rois.iterrows() for rl in RULES
            if (not only or f"{Rw.archetype}_{Rw.roi_id}" in only)]
    jobs = [j for j in jobs if not os.path.exists(f"{OUT}/inner/{j[0]}__{j[2]}.json")]
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 8
    print(f"inner jobs {len(jobs)} workers {nw}", flush=True); t0 = time.time()
    with Pool(nw) as p:
        for i, f in enumerate(p.imap(job, jobs, chunksize=len(RULES)), 1):
            print(f"{i}/{len(jobs)} {time.time() - t0:.0f}s {os.path.basename(f)}", flush=True)
    print("done", flush=True)
