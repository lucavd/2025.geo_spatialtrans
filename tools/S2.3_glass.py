"""tools/S2.3_glass.py — S2.3, previsione P-4 e regola R-2 (BL-206): vetro osservato / previsto per fasce di distanza dal tessuto sulle
cinque sezioni intere. Definizioni: docs/steps/S2.3.md («Fondo sul vetro vicino al tessuto», R-2).
Sorgente T = max(u - b, 0) sui quadratini di tessuto (in_tissue di Space Ranger), b = fondo del vetro lontano (R5, B-342...B-366);
previsione P = b + [(1 - w)(G(sigma_e) * (T pi) + G(sigma_in) * (T (1 - pi))) + w E(ell) * T] integrata sul quadratino, griglia di 2 um,
pi = exp(-d / D) con d dalla trasformata di distanza di in_tissue; due correzioni moltiplicative di T sul tessuto (limitate a [0.2, 5]).
Fasce 16-32, 32-64, 64-128, 128-256 um dal tessuto; esclusi i 50 quadratini al bordo dell'area di cattura; bootstrap delle tessere di 300 um.
R-2: due candidate a due parametri adattate per sezione (Poisson sul vetro a 16-500 um): A seconda coda (w2, ell2), B b1 exp(-d / Lb).
Uso: .venv/bin/python tools/S2.3_glass.py <outdir>  -> <outdir>/glass_<KEY>.json, results/S2.3/S2.3_P4_glass.csv, S2.3_R2_fits.csv"""
import sys, os, json, time
import numpy as np, pandas as pd
from scipy import ndimage, fft as sfft, optimize
sys.path.insert(0, "tools")
import R5_common as C
import S2_3_common as S
R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"; S20 = "/mnt/micron/geo_spatialtrans/S2.0/run1_wfix"
KEYS = ["A1", "A23", "A4", "A5", "A6"]
BANDS = [(16, 32), (32, 64), (64, 128), (128, 256)]
BORDER = 50; TILE = 150; PAD = 300; BIN = 2.0; SIG_IN = 0.7
NOCORR = os.environ.get("S23_GLASS_NOCORR") == "1"      # post hoc PH-1 (dichiarata dopo P-4): le correzioni gonfiano T (vedi docs/steps/S2.3.md)
SUFFIX = "_PH1_nocorr" if NOCORR else ""

def kft(shape):
    kr, kc = C.freq_grid(shape, BIN)
    return kr, kc, C.box_ft(kr, kc, BIN)

def conv(X, Kf, shape):
    return sfft.irfft2(sfft.rfft2(X, workers=1) * Kf, s=shape, workers=1)

def model(T, pi, prm, shape, K):
    kr, kc, box = K
    Ge = C.gauss_ft(kr, kc, prm["sigma_e"], prm["sigma_e"]); Gi = C.gauss_ft(kr, kc, SIG_IN, SIG_IN); E = C.exp_ft(kr, kc, prm["ell"])
    Tf = sfft.rfft2(T, workers=1); Tpf = sfft.rfft2(T * pi, workers=1)
    F = (1 - prm["w"]) * (Ge * Tpf + Gi * (Tf - Tpf)) + prm["w"] * E * Tf
    return sfft.irfft2(F * box, s=shape, workers=1), Tf

def run(KEY, OUT):
    t0 = time.time()
    z = np.load(f"{R5RUN}/section/{KEY}.npz"); U = z["U"].astype(float); IT = z["in_tissue"]
    prm = json.load(open(f"{S20}/params_flat.json"))[KEY]["C2"]
    b = float(pd.read_csv("results/R5/R5_kernel_descriptive.csv").set_index("key").loc[KEY, "b_far"])
    Up = np.pad(U, PAD); ITp = np.pad(IT, PAD); shape = Up.shape
    inner = np.zeros(shape, bool); inner[PAD + BORDER:PAD + U.shape[0] - BORDER, PAD + BORDER:PAD + U.shape[1] - BORDER] = True
    din = ndimage.distance_transform_edt(ITp) * BIN; pi = np.where(ITp, np.exp(-din / prm["D"]), 0.0)
    dg = ndimage.distance_transform_edt(~ITp) * BIN
    K = kft(shape)
    T = np.where(ITp, np.maximum(Up - b, 0.0), 0.0)
    for _ in range(0 if NOCORR else 2):          # post hoc PH-1: senza le correzioni (S23_GLASS_NOCORR=1)
        Pm, _ = model(T, pi, prm, shape, K)
        ratio = np.where(ITp & (Pm > 1e-9), (Up - b) / np.maximum(Pm, 1e-9), 1.0)
        T = np.where(ITp, T * np.clip(ratio, 0.2, 5.0), 0.0)
    Pm, Tf = model(T, pi, prm, shape, K)
    P = b + Pm
    glass = (~ITp) & inner
    tiles = ((np.arange(shape[0])[:, None] - PAD) // TILE) * 1000 + ((np.arange(shape[1])[None, :] - PAD) // TILE)
    out = {"key": KEY, "b_far": b, "prm": prm, "tissue_obs_over_pred": float(Up[ITp & inner].sum() / P[ITp & inner].sum())}
    rng = np.random.default_rng([S.SEED_BOOT, KEYS.index(KEY)])
    rows = []
    for lo, hi in BANDS:
        m = glass & (dg >= lo) & (dg < hi)
        tid = tiles[m]; uo = Up[m]; up = P[m]
        ut, inv = np.unique(tid, return_inverse=True)
        so = np.bincount(inv, weights=uo); sp = np.bincount(inv, weights=up)
        bs = []
        for _ in range(S.B):
            k = rng.integers(0, len(ut), len(ut)); bs.append(so[k].sum() / sp[k].sum())
        rows.append({"key": KEY, "band": f"{lo}-{hi}", "n_bins": int(m.sum()), "n_tiles": int(len(ut)), "obs_mean": float(uo.mean()),
                     "pred_mean": float(up.mean()), "ratio": float(uo.sum() / up.sum()), "lo": float(np.percentile(bs, 2.5)), "hi": float(np.percentile(bs, 97.5))})
    # R-2: candidate adattate sul vetro a 16-500 um (verosimiglianza di Poisson)
    fitm = glass & (dg >= 16) & (dg < 500)
    y = Up[fitm]
    kr, kc, box = K
    core_glass = Pm[fitm]
    def nll(lam): return float(np.sum(lam - y * np.log(np.maximum(lam, 1e-12))))
    nll0 = nll(P[fitm])
    def fA(th):
        w2, l2 = 1 / (1 + np.exp(-th[0])), np.exp(th[1])
        E2 = sfft.irfft2(Tf * C.exp_ft(kr, kc, l2) * box, s=shape, workers=1)[fitm]
        return nll(b + (1 - w2) * core_glass + w2 * E2)
    def fB(th):
        b1, lb = np.exp(th[0]), np.exp(th[1])
        return nll(P[fitm] + b1 * np.exp(-dg[fitm] / lb))
    fits = {}
    for nm, f, starts in (("A_second_tail", fA, [(-3, np.log(80)), (-2, np.log(150)), (-4, np.log(300))]),
                          ("B_glass_term", fB, [(np.log(0.05), np.log(50)), (np.log(0.2), np.log(100)), (np.log(0.02), np.log(200))])):
        best = None
        for s0 in starts:
            r = optimize.minimize(f, np.array(s0, float), method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 400})
            if best is None or r.fun < best.fun: best = r
        th = best.x
        par = {"w2": float(1 / (1 + np.exp(-th[0]))), "ell2": float(np.exp(th[1]))} if nm.startswith("A") else {"b1": float(np.exp(th[0])), "Lb": float(np.exp(th[1]))}
        fits[nm] = {"nll": float(best.fun), "aic": 2 * 2 + 2 * float(best.fun), **par}
    out.update({"nll_design": nll0, "fits": fits, "sec": round(time.time() - t0, 1)})
    json.dump({"summary": out, "bands": rows}, open(f"{OUT}/glass_{KEY}{SUFFIX}.json", "w"))
    return rows, out

if __name__ == "__main__":
    OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
    allr, fr = [], []
    for KEY in KEYS:
        rows, out = run(KEY, OUT); allr += rows
        for nm, f in out["fits"].items(): fr.append({"key": KEY, "candidate": nm, "nll_design": out["nll_design"], "tissue_obs_over_pred": out["tissue_obs_over_pred"], **f})
        print(KEY, f"{out['sec']} s", [f"{r['band']}: {r['ratio']:.3f} [{r['lo']:.3f}, {r['hi']:.3f}]" for r in rows], flush=True)
    os.makedirs("results/S2.3", exist_ok=True)
    pd.DataFrame(allr).to_csv(f"results/S2.3/S2.3_P4_glass{SUFFIX}.csv", index=False); pd.DataFrame(fr).to_csv(f"results/S2.3/S2.3_R2_fits{SUFFIX}.csv", index=False)
