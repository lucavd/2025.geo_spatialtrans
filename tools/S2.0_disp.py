"""tools/S2.0_disp.py — S2.0, M3 (P-4): prototipo del modello del documento sulle 30 ROI e indice di dispersione dentro il nucleo.
Nuclei e impronte di R4 (tools/R5_nc.py prep), totali osservati per cellula, r invertito con la regola scelta (se l'inversione non riesce:
statistica osservata pesata per cellula, DEC-034, con un passo di correzione delle densita'), sfocatura della regola scelta, scarto della sezione,
parita' della ROI, Poisson; 5 repliche con seme 20261103. Bin profondi: centro a >= 1 µm dentro la maschera del nucleo.
Uso: .venv/bin/python tools/S2.0_disp.py <outdir> <params.json> <rule> -> <outdir>/disp/<tag>.json"""
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

import pyarrow.parquet as pq
OUT, PJF, RULE = sys.argv[1], sys.argv[2], sys.argv[3]
NREP = 5
STRIPES = os.environ.get("S20_STRIPES") == "1"          # post hoc PH-2 (dichiarata dopo P-4)
NBR = os.environ.get("S20_NBR") == "1"                  # post hoc PH-3 (dichiarata dopo PH-2)
DDIR = ("disp_nbr" if NBR else "disp_stripes") if STRIPES else "disp"

def roi_rowscols(Rw):
    sq = f"{NC.D}/{NC.DS[Rw.archetype]}/binned_outputs_x/binned_outputs/square_002um"
    tp = pq.read_table(f"{sq}/spatial/tissue_positions.parquet").to_pandas()
    box = tp[(tp.pxl_row_in_fullres >= Rw.r0) & (tp.pxl_row_in_fullres < Rw.r1) & (tp.pxl_col_in_fullres >= Rw.c0) & (tp.pxl_col_in_fullres < Rw.c1) & (tp.in_tissue == 1)]
    return box.array_row.to_numpy(), box.array_col.to_numpy()

def cv8(u, rows, cols):
    """CV dei bin di 8 µm interamente presenti (16 bin di 2 µm nella ROI e nel tessuto)."""
    k = pd.DataFrame({"g": (rows // 4) * 100000 + (cols // 4), "u": u})
    g = k.groupby("g").u.agg(["sum", "size"]); g = g[g["size"] == 16]["sum"]
    return float(g.std() / g.mean()) if len(g) > 1 and g.mean() > 0 else float("nan"), int(len(g))

def nbr_corr(u, lab, c, keep, rows, cols, min_bins=4):
    """Post hoc PH-3: correlazione dei residui standardizzati (y - m)/sqrt(m) fra bin dello stesso nucleo profondo a distanza 1 e 2
    lungo le righe (h) e le colonne (v); m come in psi_index."""
    keep = np.asarray(keep); y = np.asarray(u, float); cc = np.asarray(c, float); lab = np.asarray(lab)
    sel = np.isin(lab, keep) & (lab > 0)
    z = np.full(len(y), np.nan)
    l = lab[sel]; uniq, inv, n = np.unique(l, return_inverse=True, return_counts=True)
    sy = np.bincount(inv, weights=y[sel]); sc = np.bincount(inv, weights=cc[sel])
    ok = ((n >= min_bins) & (sy > 0))[inv]
    m = (sy / n)[inv] * cc[sel] / (sc / n)[inv]
    zz = np.where(ok, (y[sel] - m) / np.sqrt(np.where(ok, m, 1.0)), np.nan); z[sel] = zz
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(rows, cols))}
    out = {}
    for nm, dr, dc in (("h1", 0, 1), ("h2", 0, 2), ("v1", 1, 0), ("v2", 2, 0)):
        ia, ib = [], []
        for i, (a, b) in enumerate(zip(rows, cols)):
            j = key.get((int(a) + dr, int(b) + dc))
            if j is not None and lab[i] > 0 and lab[i] == lab[j] and np.isfinite(z[i]) and np.isfinite(z[j]): ia.append(i); ib.append(j)
        out[nm] = float(np.corrcoef(z[ia], z[ib])[0, 1]) if len(ia) > 10 else float("nan"); out[nm + "_n"] = len(ia)
    return out

def job(args):
    tag, Rw = args
    fo = f"{OUT}/{DDIR}/{tag}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        PJ = json.load(open(PJF)); prm = rule_prm(PJ, Rw.archetype, RULE)
        Dd, IX, IXf, Dobs, am, acell, S_, U5, base = NC.prep(tag, Rw)
        if STRIPES:      # post hoc PH-2: fattori di cattura per riga e colonna stimati sulla sezione (M1, adattamento completo con la parita')
            rows0, cols0 = roi_rowscols(Rw)
            Fz = np.load(f"{OUT}/stripes/{NC.KEY[Rw.archetype]}_factors.npz")
            cf = np.exp(np.nan_to_num(Fz["alpha"])[rows0] + np.nan_to_num(Fz["beta"])[cols0]); Dd = dict(Dd); Dd["c"] = cf / cf.mean()
        Ui = np.random.default_rng([S.SEED, int(Rw.name)]).uniform(size=U5.shape)
        di = NC.d_image(Dd, prm)
        fm = roi_tissue(Rw, Dd) if RULE in ("C2", "C3") else None
        blur = S.make_blur(RULE, prm, Dd["pos"], di, fm)
        r, Sr = S.invert_r_refined_b(base["obs_star"], S_, blur, Dd["c"], Ui, IX, Dd["L"], Dd["keep"], Dobs, acell, am)
        r_src = "inverted"
        if not np.isfinite(r):
            r = base["obs_star"]; r_src = "observed_DEC034"
            Sr = S.refine_b(S_, r, blur, Dd["c"], Dd["L"], Dd["keep"], Dobs, acell, am)
        Ud = np.random.default_rng([S.SEED, 1000 + int(Rw.name)]).uniform(size=(NREP, len(Dd["u"])))
        Y = S.nc_simulate_b(Sr, r, blur, Dd["c"], Ud)
        lam = Dd["c"] * np.maximum(blur(C.nc_field(Sr, r)), 1e-12)
        lab = S.deep_labels(Dd["mask"], Dd["rr"], Dd["cc"], Dd["upp"], 1.0)
        po, mo, no = S.psi_index(Dd["u"], lab, Dd["c"], Dd["keep"])
        sims = [S.psi_index(y, lab, Dd["c"], Dd["keep"]) for y in Y]
        ps_ = float(np.mean([s[0] for s in sims]))
        rows, cols = roi_rowscols(Rw); assert len(rows) == len(Dd["u"])
        cvo, n8 = cv8(Dd["u"], rows, cols); cvs = float(np.mean([cv8(y, rows, cols)[0] for y in Y]))
        res = {"tag": tag, "archetype": Rw.archetype, "rule": RULE, "r": float(r), "r_source": r_src, "psi_obs": po, "psi_sim": ps_,
               "psi_sim_reps": [s[0] for s in sims], "ratio": po / ps_, "mbar_obs": mo, "v": (po - ps_) / mo, "n_nuclei": no,
               "zero_frac_obs": float((Dd["u"] == 0).mean()), "zero_frac_sim": float(np.mean([(y == 0).mean() for y in Y])),
               "cv8_obs": cvo, "cv8_sim": cvs, "n_bins8": n8, "mean_lam": float(lam.mean()), "mean_obs": float(Dd["u"].mean()),
               "sec": round(time.time() - t0, 1)}
        if NBR:
            no_ = nbr_corr(Dd["u"], lab, Dd["c"], Dd["keep"], rows, cols)
            ns_ = [nbr_corr(y, lab, Dd["c"], Dd["keep"], rows, cols) for y in Y]
            for k_ in ("h1", "h2", "v1", "v2"):
                res[f"rho_{k_}_obs"] = no_[k_]; res[f"rho_{k_}_sim"] = float(np.mean([s[k_] for s in ns_])); res[f"n_{k_}"] = no_[k_ + "_n"]
    except Exception:
        res = {"tag": tag, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

if __name__ == "__main__":
    os.makedirs(f"{OUT}/{DDIR}", exist_ok=True)
    rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    jobs = [(f"{Rw.archetype}_{Rw.roi_id}", Rw) for _, Rw in rois.iterrows()]
    jobs = [j for j in jobs if not os.path.exists(f"{OUT}/{DDIR}/{j[0]}.json")]
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 8
    print(f"disp jobs {len(jobs)} workers {nw} rule {RULE}", flush=True); t0 = time.time()
    with Pool(nw) as p:
        for i, f in enumerate(p.imap_unordered(job, jobs, chunksize=1), 1):
            print(f"{i}/{len(jobs)} {time.time() - t0:.0f}s {os.path.basename(f)}", flush=True)
    print("done", flush=True)
