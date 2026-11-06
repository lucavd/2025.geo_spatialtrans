"""tools/S2.4_explore.py — S2.4, esplorazioni dichiarate prima della pre-registrazione (docs/steps/S2.4.md, «Esplorazioni»).
E-1 real  : 30 aree vere: strisce con lo stimatore di S2.0 + PH-5 (fattori per riga e colonna, parte fine, parita'), fasce profonde dei
            lumi (8-16, 16-32, 32-64 um), frazione di quadratini vuoti e UMI medio in tessuto. Solo dati veri.
E-2 sim   : le stesse statistiche sul pannello di S2.3 (variante main: niente strisce, niente seconda coda, fondo di Poisson).
E-3 glass : 5 sezioni: vetro osservato / previsto con la seconda coda di DEC-051 (valori senza correzioni), per fasce fino a 512 um e
            nella zona del fondo lontano (>= 300 um); imbottitura della FFT 300 (come S2.3) e 1500 quadratini; fondo b rifatto a coda fissa.
E-4 mfac  : fattore di sovra-generazione max(c)/media(c) del generatore delle strisce (design sezione 7) su area di 1 mm2 con margine,
            tessera di 500 um con bordo di 280 um, area di cattura intera (3350 x 3350 quadratini).
Uso: .venv/bin/python tools/S2.4_explore.py real|sim|glass|mfac <outdir> [args]"""
import sys, os, json, time, glob, traceback
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
import importlib.util
def _load(name, path):
    sp = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m

DEEP_BANDS = [(2.0, 4.0), (4.0, 8.0), (8.0, 16.0), (16.0, 32.0), (32.0, 64.0)]
MIN_HALF_AREA = 100          # bin di tessuto per meta' riga (sezioni: 200; un'area di 1 mm2 ha al piu' ~500 bin per riga)

def stripe_summary(PH, U, T):
    """Stimatore di S2.0 + PH-5 (maschera senza righe e colonne a zero, 10 iterazioni, parte fine W = 21) su un'area."""
    PH.MIN_HALF = MIN_HALF_AREA
    T2, zr, zc = PH.fix_mask(U, T)
    o = PH.summarize(U, T2, 10)
    out = {"zero_rows_removed": zr, "zero_cols_removed": zc}
    for ax, v in o.items():
        for k, x in v.items(): out[f"stripes_{ax}_{k}"] = x
    return out

def area_matrix(i, j, u, tis):
    # la matrice parte da una riga e da una colonna pari della griglia di Space Ranger: la parita' dello stimatore (indice locale)
    # coincide cosi' con quella di array_row / array_col (prima versione: indice relativo all'area, parita' invertita se i0 dispari)
    i0, j0 = i.min() - i.min() % 2, j.min() - j.min() % 2; nr, nc = i.max() - i0 + 1, j.max() - j0 + 1
    U = np.zeros((nr, nc)); T = np.zeros((nr, nc), bool)
    U[i - i0, j - j0] = u; T[i - i0, j - j0] = tis
    return U, T, int(i0), int(j0)

def job_real(Rw, OUT):
    tag = f"{Rw.archetype}_{Rw.roi_id}"; fo = f"{OUT}/real/{tag}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        import shapely
        RL = _load("s23real", "tools/S2.3_real.py"); RL.BANDS = DEEP_BANDS
        PH = _load("ph5", "tools/S2.0_ph5_stripes_hp.py")
        A, upp = Rw.archetype, Rw.um_per_px
        PJ = json.load(open(f"{RL.S20}/params_flat.json")); prm = PJ[RL.NC.KEY[A]]["C2"]
        Dd = RL.NC.roi_data(Rw); di = RL.NC.d_image(Dd, prm)
        rowsA, colsA, posA = RL.box_bins(Rw, False)
        z = np.load(f"{RL.R5RUN}/section/{RL.NC.KEY[A]}.npz"); uA = z["U"][rowsA, colsA].astype(float)
        posA_al = posA - di[None, :]; side = float(Rw.side_px * upp)
        R = RL.region_geom(tag); inreg = shapely.contains_xy(R, posA_al[:, 1], posA_al[:, 0])
        fr = (posA_al[:, 0] >= RL.FRAME_UM) & (posA_al[:, 0] <= side - RL.FRAME_UM) & (posA_al[:, 1] >= RL.FRAME_UM) & (posA_al[:, 1] <= side - RL.FRAME_UM)
        P = inreg & fr
        res = {"tag": tag, "archetype": A, "roi_id": Rw.roi_id, "side_um": side}
        U, T, i0, j0 = area_matrix(rowsA, colsA, uA, inreg)
        res["array_row0"] = i0; res["array_col0"] = j0
        res.update(stripe_summary(PH, U, T))
        res["n_tissue_bins"] = int(P.sum()); res["mean_umi_tissue"] = float(uA[P].mean()); res["frac_zero_tissue"] = float((uA[P] == 0).mean())
        res.update({f"lumen_{k}": v for k, v in RL.lumen_stats(R, posA_al, uA, inreg, side).items()})
        res["sec"] = round(time.time() - t0, 1)
    except Exception:
        res = {"tag": tag, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

def job_sim(f, Rw, OUT):
    tag_k = os.path.basename(f)[:-5]; fo = f"{OUT}/sim/{tag_k}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        import shapely
        RL = _load("s23real", "tools/S2.3_real.py"); RL.BANDS = DEEP_BANDS
        PH = _load("ph5", "tools/S2.0_ph5_stripes_hp.py")
        meta = json.load(open(f)); tag = meta["tag"]; g = meta["grid"]
        side = float(Rw.side_px * Rw.um_per_px)
        tot = np.fromfile(f[:-5] + ".bin", dtype=np.int32).reshape(g["n_row"], g["n_col"]).astype(float)
        i, j = np.divmod(np.arange(g["n_row"] * g["n_col"]), g["n_col"])
        x = g["x0"] + (j + 0.5) * g["bin_um"]; y = g["y0"] + (i + 0.5) * g["bin_um"]
        R = RL.region_geom(tag); inreg_all = shapely.contains_xy(R, x, y)
        res = {"tag": tag, "archetype": meta["archetype"], "k": meta["k"], "variant": meta["variant"], "side_um": side}
        res.update(stripe_summary(PH, tot, inreg_all.reshape(tot.shape)))
        box = (x >= 0) & (x < side) & (y >= 0) & (y < side)
        x, y, u, inreg = x[box], y[box], tot.ravel()[box], inreg_all[box]
        fr = (y >= RL.FRAME_UM) & (y <= side - RL.FRAME_UM) & (x >= RL.FRAME_UM) & (x <= side - RL.FRAME_UM)
        P = inreg & fr
        res["n_tissue_bins"] = int(P.sum()); res["mean_umi_tissue"] = float(u[P].mean()); res["frac_zero_tissue"] = float((u[P] == 0).mean())
        res.update({f"lumen_{k}": v for k, v in RL.lumen_stats(R, np.c_[y, x], u, inreg, side).items()})
        res["sec"] = round(time.time() - t0, 1)
    except Exception:
        res = {"file": f, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

# ---------------- E-3: vetro delle sezioni ----------------
DEC051 = "results/S2.3/S2.3_R2_fits_PH1_nocorr.csv"
GL_BANDS = [(16, 32), (32, 64), (64, 128), (128, 256), (256, 512), (300, 10000), (512, 10000)]

def job_glass(a):
    KEY, PADB, OUT = a
    fo = f"{OUT}/glass/{KEY}_pad{PADB}.json"
    if os.path.exists(fo): return fo
    from scipy import ndimage, fft as sfft, optimize
    GL = _load("s23glass", "tools/S2.3_glass.py"); C = GL.C
    t0 = time.time()
    z = np.load(f"{GL.R5RUN}/section/{KEY}.npz"); U = z["U"].astype(float); IT = z["in_tissue"].astype(bool)
    prm = json.load(open(f"{GL.S20}/params_flat.json"))[KEY]["C2"]
    b = float(pd.read_csv("results/R5/R5_kernel_descriptive.csv").set_index("key").loc[KEY, "b_far"])
    fits = pd.read_csv(DEC051); fa = fits[(fits.key == KEY) & (fits.candidate == "A_second_tail")].iloc[0]
    w2, l2 = float(fa.w2), float(fa.ell2)
    Up = np.pad(U, PADB); ITp = np.pad(IT, PADB); shape = Up.shape
    inner = np.zeros(shape, bool); inner[PADB + GL.BORDER:PADB + U.shape[0] - GL.BORDER, PADB + GL.BORDER:PADB + U.shape[1] - GL.BORDER] = True
    din = ndimage.distance_transform_edt(ITp) * GL.BIN; pi = np.where(ITp, np.exp(-din / prm["D"]), 0.0)
    dg = ndimage.distance_transform_edt(~ITp) * GL.BIN
    K = GL.kft(shape); kr, kc, box = K
    T = np.where(ITp, np.maximum(Up - b, 0.0), 0.0)
    Pm, Tf = GL.model(T, pi, prm, shape, K)
    E2 = sfft.irfft2(Tf * C.exp_ft(kr, kc, l2) * box, s=shape, workers=1)
    glass = (~ITp) & inner
    lam_core = (1 - w2) * Pm; lam_tail = w2 * E2
    out = {"key": KEY, "pad_bins": PADB, "b_far": b, "w2": w2, "ell2": l2, "max_dist_um": float(dg[glass].max())}
    # massa della seconda coda: frazione che esce dall'area di cattura (fuori dal riquadro dei dati, prima dell'imbottitura)
    capt = np.zeros(shape, bool); capt[PADB:PADB + U.shape[0], PADB:PADB + U.shape[1]] = True
    out["tail_mass_inside_capture_frac"] = float(E2[capt].sum() / max(T.sum(), 1e-9))
    rows = []
    for lo, hi in GL_BANDS:
        m = glass & (dg >= lo) & (dg < hi)
        if m.sum() == 0: continue
        rows.append({"band": f"{lo}-{hi}", "n_bins": int(m.sum()), "obs_mean": float(Up[m].mean()), "b": b,
                     "core_mean": float(lam_core[m].mean()), "tail_mean": float(lam_tail[m].mean()),
                     "ratio_with_tail": float(Up[m].sum() / (b * m.sum() + lam_core[m].sum() + lam_tail[m].sum())),
                     "ratio_no_tail": float(Up[m].sum() / (b * m.sum() + Pm[m].sum()))})
    out["bands"] = rows
    # fondo rifatto a coda fissa (Poisson sul vetro a >= 16 um, tutta la distanza)
    fm = glass & (dg >= 16); y = Up[fm]; kk = (lam_core + lam_tail)[fm]
    g = lambda bb: float(np.sum(y / (bb + kk)) - fm.sum())
    try: out["b_refit_tail_fixed"] = float(optimize.brentq(g, 1e-6, 5.0))
    except Exception as e: out["b_refit_tail_fixed"] = None; out["b_refit_err"] = str(e)
    # coda rifatta con questa imbottitura, b fisso (come S2.3 PH-1) e con b libero; vetro a 16-500 um come S2.3
    fit = glass & (dg >= 16) & (dg < 500); y = Up[fit]; core_g = Pm[fit]
    def nll(lam): return float(np.sum(lam - y * np.log(np.maximum(lam, 1e-12))))
    def E2f(l): return sfft.irfft2(Tf * C.exp_ft(kr, kc, l) * box, s=shape, workers=1)[fit]
    def fA(th): ww, ll = 1 / (1 + np.exp(-th[0])), np.exp(th[1]); return nll(b + (1 - ww) * core_g + ww * E2f(ll))
    def fA3(th): ww, ll, bb = 1 / (1 + np.exp(-th[0])), np.exp(th[1]), np.exp(th[2]); return nll(bb + (1 - ww) * core_g + ww * E2f(ll))
    best = None
    for st in ((np.log(w2 / (1 - w2)), np.log(l2)), (-3, np.log(80)), (-4, np.log(300))):
        rr = optimize.minimize(fA, st, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 0.01, "maxiter": 400})
        if best is None or rr.fun < best.fun: best = rr
    out["refit_bfixed"] = {"w2": float(1 / (1 + np.exp(-best.x[0]))), "ell2": float(np.exp(best.x[1])), "nll": float(best.fun)}
    best3 = None
    for st in ((np.log(w2 / (1 - w2)), np.log(l2), np.log(b)), (-3, np.log(150), np.log(b / 2))):
        rr = optimize.minimize(fA3, st, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 0.01, "maxiter": 600})
        if best3 is None or rr.fun < best3.fun: best3 = rr
    out["refit_bfree"] = {"w2": float(1 / (1 + np.exp(-best3.x[0]))), "ell2": float(np.exp(best3.x[1])), "b": float(np.exp(best3.x[2])), "nll": float(best3.fun)}
    out["sec"] = round(time.time() - t0, 1)
    json.dump(out, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

# ---------------- E-4: fattore di sovra-generazione ----------------
# parametri delle strisce: design dello Step 2 sezione 9 (B-397...B-401): parita' righe, colonne; SD alpha, beta; passo 1 alpha, beta
STRIPES = {"A1": (1.417, 1.750, 0.193, 0.185, -0.30, -0.27), "A2": (1.578, 1.724, 0.221, 0.191, -0.31, -0.21),
           "A4": (0.996, 1.105, 0.211, 0.251, -0.38, -0.46), "A5": (1.201, 1.231, 0.242, 0.180, -0.31, -0.29),
           "A6": (0.947, 0.967, 0.183, 0.188, -0.17, -0.21)}
def ma1(rng, n, sd, rho):
    th = (1 - np.sqrt(1 - 4 * rho ** 2)) / (2 * rho) if rho != 0 else 0.0
    e = rng.normal(0, 1, n + 1); return sd * (e[1:] + th * e[:-1]) / np.sqrt(1 + th ** 2)
def job_mfac(a):
    A, rep, OUT = a
    pr, pc, sa, sb, ra, rb = STRIPES[A]; rng = np.random.default_rng([20261114, 40, ["A1", "A2", "A4", "A5", "A6"].index(A), rep])
    N = 3350; al = ma1(rng, N, sa, ra); be = ma1(rng, N, sb, rb)
    ea, eb = np.exp(al), np.exp(be); par_r = np.where(np.arange(N) % 2 == 1, pr, 1.0); par_c = np.where(np.arange(N) % 2 == 1, pc, 1.0)
    fr, fc = ea * par_r, eb * par_c               # c = fr_i * fc_j (separabile)
    def mx(i0, i1, j0, j1, norm): return float(fr[i0:i1].max() * fc[j0:j1].max() / norm)
    def mean_(i0, i1, j0, j1): return float(fr[i0:i1].mean() * fc[j0:j1].mean())
    mg = mean_(0, N, 0, N); out = {"A": A, "rep": rep, "M_global": mx(0, N, 0, N, mg)}
    # area di 1 mm2 + margine 280 um = 780 quadratini; normalizzazione sull'area (tessuto ~ area intera)
    areas = [mx(i, i + 780, j, j + 780, mean_(i + 140, i + 640, j + 140, j + 640)) for i in range(0, N - 780, 390) for j in range(0, N - 780, 390)]
    tiles = [mx(max(0, i - 140), min(N, i + 250 + 140), max(0, j - 140), min(N, j + 250 + 140), mg) for i in range(0, N, 250) for j in range(0, N, 250)]
    out.update({"M_area_median": float(np.median(areas)), "M_area_max": float(np.max(areas)), "M_tile_mean": float(np.mean(tiles)),
                "M_tile_max": float(np.max(tiles)), "M_tile_min": float(np.min(tiles))})
    return out

if __name__ == "__main__":
    mode, OUT = sys.argv[1], sys.argv[2]
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 4
    rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois["tag"] = rois.archetype + "_" + rois.roi_id
    t0 = time.time()
    if mode == "real":
        os.makedirs(f"{OUT}/real", exist_ok=True)
        jobs = [(Rw, OUT) for _, Rw in rois.iterrows()]
        with Pool(min(nw, len(jobs))) as p: p.starmap(job_real, jobs, chunksize=1)
    elif mode == "sim":
        PD, VAR = sys.argv[3], sys.argv[4]; os.makedirs(f"{OUT}/sim", exist_ok=True)
        files = sorted(glob.glob(f"{PD}/{VAR}/*_k*.json"))
        jobs = [(f, rois.set_index("tag").loc[json.load(open(f))["tag"]], OUT) for f in files]
        with Pool(min(nw, len(jobs))) as p: p.starmap(job_sim, jobs, chunksize=1)
    elif mode == "glass":
        os.makedirs(f"{OUT}/glass", exist_ok=True)
        jobs = [(k, pb, OUT) for pb in (300, 1500) for k in ("A1", "A23", "A4", "A5", "A6")]
        with Pool(min(nw, len(jobs))) as p: p.map(job_glass, jobs, chunksize=1)
    elif mode == "mfac":
        jobs = [(A, r, OUT) for A in STRIPES for r in range(20)]
        with Pool(min(nw, len(jobs))) as p: res = p.map(job_mfac, jobs, chunksize=1)
        pd.DataFrame(res).to_csv(f"{OUT}/E4_mfac.csv", index=False)
    print(f"{mode} done {time.time() - t0:.0f}s", flush=True)
