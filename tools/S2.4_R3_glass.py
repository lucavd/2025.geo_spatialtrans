"""tools/S2.4_R3_glass.py — S2.4, regola R-3 e previsione P-5 (docs/steps/S2.4.md).
Per sezione: sorgente T = max(u - b_far, 0) sul tessuto (in_tissue), come S2.3 senza correzioni; modello del vetro
lambda = b + (1 - w2) * core + w2 * E(ell2) * T, con core = regola di S2.0 (C2, parametri della sezione) integrata sul quadratino;
imbottitura della trasformata PAD = 1500 quadratini (3 mm) per lato (E-3: con 300 la coda rientra dal lato opposto).
Adattamento di Poisson a tre parametri sul vetro a >= 16 um dal tessuto, a qualunque distanza (esclusi i 50 quadratini al bordo
dell'area di cattura): verosimiglianza profilo su ell2 (griglia log da 10 a 5000 um, 41 punti, poi Brent nell'intervallo del minimo),
con (b, w2) ottimizzati per ogni ell2 (L-BFGS-B con gradiente esatto, b >= 0, 0 <= w2 <= 0.5; il problema in (b, w2) e' convesso).
P-5: osservato / previsto per fasce 16-32, 32-64, 64-128, 128-256, 256-512, >= 512 um; bootstrap delle tessere di 300 um (B = 2000,
seme 20261115). Uscite: <outdir>/R3_<KEY>.json; results/S2.4/S2.4_R3_fits.csv, S2.4_P5_glass.csv (con `summary`).
Uso: .venv/bin/python tools/S2.4_R3_glass.py fit <outdir> [KEY ...] | summary <outdir>"""
import sys, os, json, time
import numpy as np, pandas as pd
from scipy import ndimage, fft as sfft, optimize
sys.path.insert(0, "tools")
import importlib.util
_sp = importlib.util.spec_from_file_location("s23glass", "tools/S2.3_glass.py"); GL = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(GL)
C = GL.C
KEYS = ["A1", "A23", "A4", "A5", "A6"]
PADB = 1500; DMIN = 16.0; SEED_BOOT = 20261115; B = 2000; TILE = 150
BANDS = [(16, 32), (32, 64), (64, 128), (128, 256), (256, 512), (512, 100000)]
ELL_GRID = np.exp(np.linspace(np.log(10.0), np.log(5000.0), 41))

def fit_bw(y, core, e2, x0):
    """min sum(lam - y log lam), lam = b + core + w (e2 - core); vincoli b >= 0, 0 <= w <= 0.5."""
    dlt = e2 - core
    def f(th):
        b, w = th; lam = np.maximum(b + core + w * dlt, 1e-12); r = 1.0 - y / lam
        return float(np.sum(lam - y * np.log(lam))), np.array([np.sum(r), np.sum(r * dlt)])
    rr = optimize.minimize(f, x0, jac=True, method="L-BFGS-B", bounds=[(0.0, 50.0), (0.0, 0.5)], options={"ftol": 1e-14, "gtol": 1e-8, "maxiter": 500})
    return rr.fun, rr.x

def run(KEY, OUT):
    fo = f"{OUT}/R3_{KEY}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    z = np.load(f"{GL.R5RUN}/section/{KEY}.npz"); U = z["U"].astype(float); IT = z["in_tissue"].astype(bool)
    prm = json.load(open(f"{GL.S20}/params_flat.json"))[KEY]["C2"]
    bfar = float(pd.read_csv("results/R5/R5_kernel_descriptive.csv").set_index("key").loc[KEY, "b_far"])
    f51 = pd.read_csv("results/S2.3/S2.3_R2_fits_PH1_nocorr.csv"); f51 = f51[(f51.key == KEY) & (f51.candidate == "A_second_tail")].iloc[0]
    Up = np.pad(U, PADB); ITp = np.pad(IT, PADB); shape = Up.shape
    inner = np.zeros(shape, bool); inner[PADB + GL.BORDER:PADB + U.shape[0] - GL.BORDER, PADB + GL.BORDER:PADB + U.shape[1] - GL.BORDER] = True
    din = ndimage.distance_transform_edt(ITp) * GL.BIN; pi = np.where(ITp, np.exp(-din / prm["D"]), 0.0); del din
    dg = ndimage.distance_transform_edt(~ITp) * GL.BIN
    K = GL.kft(shape); kr, kc, box = K
    T = np.where(ITp, np.maximum(Up - bfar, 0.0), 0.0)
    Pm, Tf = GL.model(T, pi, prm, shape, K); del pi
    glass = (~ITp) & inner
    fm = glass & (dg >= DMIN)
    y = Up[fm]; core = Pm[fm]; dist = dg[fm]
    def E2(l): return sfft.irfft2(Tf * C.exp_ft(kr, kc, l) * box, s=shape, workers=1)
    prof = []; x0 = np.array([bfar, float(f51.w2)])
    for l in ELL_GRID:
        e2 = E2(l)[fm]; nll, x = fit_bw(y, core, e2, x0); prof.append((float(l), nll, float(x[0]), float(x[1])))
    prof = sorted(prof); k = int(np.argmin([p[1] for p in prof]))
    lo = prof[max(k - 1, 0)][0]; hi = prof[min(k + 1, len(prof) - 1)][0]
    cache = {}
    def g(logl):
        l = float(np.exp(logl)); e2 = E2(l)[fm]; nll, x = fit_bw(y, core, e2, np.array([prof[k][2], max(prof[k][3], 1e-4)])); cache[logl] = (nll, x); return nll
    if lo < hi:
        rr = optimize.minimize_scalar(g, bounds=(np.log(lo), np.log(hi)), method="bounded", options={"xatol": 1e-4})
        lbest = float(np.exp(rr.x)); nbest, xb = cache[rr.x]
    else:
        lbest, nbest, xb = prof[k][0], prof[k][1], np.array([prof[k][2], prof[k][3]])
    bh, wh = float(xb[0]), float(xb[1])
    # stesso adattamento con i valori di DEC-051 e b_far (riferimento) e nll del solo modello di S2.0
    e51 = E2(float(f51.ell2))
    nll51 = float(np.sum((bfar + (1 - f51.w2) * core + f51.w2 * e51[fm]) - y * np.log(bfar + (1 - f51.w2) * core + f51.w2 * e51[fm])))
    lam = bh + (1 - wh) * Pm + wh * E2(lbest)
    tiles = ((np.arange(shape[0])[:, None] - PADB) // TILE) * 100000 + ((np.arange(shape[1])[None, :] - PADB) // TILE)
    rng = np.random.default_rng([SEED_BOOT, KEYS.index(KEY)]); rows = []
    capt = np.zeros(shape, bool); capt[PADB:PADB + U.shape[0], PADB:PADB + U.shape[1]] = True
    for a, b_ in BANDS:
        m = glass & (dg >= a) & (dg < b_)
        if m.sum() == 0: continue
        tid = tiles[m]; uo = Up[m]; up = lam[m]; ut, inv = np.unique(tid, return_inverse=True)
        so = np.bincount(inv, weights=uo); sp = np.bincount(inv, weights=up)
        bs = np.array([so[kk].sum() / sp[kk].sum() for kk in (rng.integers(0, len(ut), len(ut)) for _ in range(B))])
        rows.append({"key": KEY, "band": f"{a}-{b_}", "n_bins": int(m.sum()), "n_tiles": int(len(ut)), "obs_mean": float(uo.mean()),
                     "pred_mean": float(up.mean()), "ratio": float(uo.sum() / up.sum()), "lo": float(np.percentile(bs, 2.5)), "hi": float(np.percentile(bs, 97.5)),
                     "se_boot": float(bs.std(ddof=1))})
    tis = ITp & inner
    out = {"key": KEY, "pad_bins": PADB, "dmin_um": DMIN, "n_glass_fit": int(fm.sum()), "max_dist_um": float(dist.max()), "b_far": bfar,
           "w2_dec051": float(f51.w2), "ell2_dec051": float(f51.ell2), "nll_dec051_bfar": nll51,
           "b": bh, "w2": wh, "ell2": lbest, "nll": float(nbest), "profile": prof,
           "at_bound_b0": bool(bh <= 1e-8), "at_bound_w2": bool(wh >= 0.5 - 1e-8 or wh <= 1e-8),
           "ell_at_grid_edge": bool(k in (0, len(prof) - 1)),
           "tissue_obs_over_pred": float(Up[tis].sum() / lam[tis].sum()),
           "tail_mass_in_capture": float((wh * E2(lbest))[capt].sum() / max(wh * T.sum(), 1e-12)) if wh > 0 else None,
           "bands": rows, "sec": round(time.time() - t0, 1)}
    json.dump(out, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

def summary(OUT):
    R = [json.load(open(f"{OUT}/R3_{k}.json")) for k in KEYS if os.path.exists(f"{OUT}/R3_{k}.json")]
    F = pd.DataFrame([{k: v for k, v in r.items() if k not in ("profile", "bands")} for r in R])
    Bd = pd.DataFrame([b for r in R for b in r["bands"]])
    Bd["verdict"] = np.where((Bd.lo >= 0.8) & (Bd.hi <= 1.25), "confermata", np.where((Bd.hi < 0.8) | (Bd.lo > 1.25), "smentita", "indeterminata"))
    os.makedirs("results/S2.4", exist_ok=True)
    F.to_csv("results/S2.4/S2.4_R3_fits.csv", index=False); Bd.to_csv("results/S2.4/S2.4_P5_glass.csv", index=False)
    pd.set_option("display.width", 250)
    print(F[["key", "b_far", "b", "w2_dec051", "w2", "ell2_dec051", "ell2", "nll_dec051_bfar", "nll", "at_bound_b0", "ell_at_grid_edge", "tissue_obs_over_pred", "max_dist_um", "sec"]].round(4).to_string(index=False))
    print(Bd[["key", "band", "n_tiles", "ratio", "lo", "hi", "verdict"]].round(3).to_string(index=False))

if __name__ == "__main__":
    mode, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
    if mode == "fit":
        from multiprocessing import Pool
        keys = sys.argv[3:] or KEYS
        nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 3
        with Pool(min(nw, len(keys))) as p: p.starmap(run, [(k, OUT) for k in keys], chunksize=1)
        summary(OUT)
    else:
        summary(OUT)
