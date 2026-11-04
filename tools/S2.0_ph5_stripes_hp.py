"""tools/S2.0_ph5_stripes_hp.py — S2.0, post hoc PH-5 (dalla revisione G-5 a, dubbi D-1 e D-2 del revisore delle strisce).
(1) Righe e colonne con >= 20 bin di tessuto e 0 UMI escono dalla maschera del tessuto prima della stima (il fattore sentinella log(1e-12)
abbassava il campo dei vicini e la media della parita'). (2) La parte lenta dei fattori non e' identificata (deriva con le iterazioni):
per il simulatore si usa la parte fine, scarto dalla media mobile centrata su W = 21 righe o colonne (media dei soli valori finiti).
Per ogni sezione, a 10 e 20 iterazioni: rapporti dispari/pari, SD riproducibile, correlazione fra meta' e passo 1 della parte fine.
(3) Calibrazione dello stesso procedimento su sezioni sintetiche (generatore di G-1, 800 x 800, struttura del tessuto 0.2 e 0.4) con
fattori noti di SD 0.2: indipendenti o media mobile di ordine 1 con passo 1 -0.2 e -0.3, 4 semi.
Uscite: results/S2.0/S2.0_PH5_sections.csv, results/S2.0/S2.0_PH5_calib.csv. Uso: .venv/bin/python tools/S2.0_ph5_stripes_hp.py <nproc>"""
import sys, os, importlib.util
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
import S2_0_common as S, R5_common as C
R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"; W = 21; MIN_HALF = 200
def hp(f):
    f = np.asarray(f, float); ok = np.isfinite(f); k = np.ones(W)
    num = np.convolve(np.where(ok, f, 0.0), k, "same"); den = np.convolve(ok.astype(float), k, "same")
    with np.errstate(invalid="ignore", divide="ignore"): m = num / den
    return np.where(ok & (den >= W // 2), f - m, np.nan)
def lag1(r1, r2):
    a = np.isfinite(r1) & np.isfinite(r2); b = np.isfinite(r1[:-1]) & np.isfinite(r2[1:])
    c0 = np.cov(r1[a], r2[a])[0, 1]; c1 = np.cov(r1[:-1][b], r2[1:][b])[0, 1]
    c1b = np.cov(r1[1:][np.isfinite(r1[1:]) & np.isfinite(r2[:-1])], r2[:-1][np.isfinite(r1[1:]) & np.isfinite(r2[:-1])])[0, 1]
    return float(0.5 * (c1 + c1b) / c0) if c0 > 0 else np.nan
def fix_mask(U, T, min_bins=20):
    T = T.astype(bool).copy(); Y = U * T
    zr = (T.sum(1) >= min_bins) & (Y.sum(1) == 0); zc = (T.sum(0) >= min_bins) & (Y.sum(0) == 0)
    T[zr, :] = False; T[:, zc] = False; return T, int(zr.sum()), int(zc.sum())
def summarize(U, T, n_iter):
    al, be, _, _ = S.stripe_factors(U, T, n_iter=n_iter); out = {}
    for ax, f, nm in ((0, al, "rows"), (1, be, "cols")):
        idx = np.arange(len(f)); mu, _ = S.periodic_split(f, idx)
        h1, h2, n1, n2 = S.half_split_factors(U, T, ax, n_iter=n_iter); ok = (n1 >= MIN_HALF) & (n2 >= MIN_HALF)
        _, r1 = S.periodic_split(np.where(ok, h1, np.nan), idx); _, r2 = S.periodic_split(np.where(ok, h2, np.nan), idx)
        q1, q2 = hp(r1), hp(r2); g = np.isfinite(q1) & np.isfinite(q2)
        out[nm] = {"odd_over_even": float(np.exp(mu[1] - mu[0])), "sd_repro_all": S.repro_sd(r1, r2), "sd_repro_hp": S.repro_sd(q1, q2),
                   "corr_halves_hp": float(np.corrcoef(q1[g], q2[g])[0, 1]), "lag1_repro_hp": lag1(q1, q2), "n_hp": int(g.sum())}
    return out
def job_sec(a):
    k, n_iter = a
    z = np.load(f"{R5RUN}/section/{k}.npz"); U = z["U"].astype(float); T0 = z["in_tissue"].astype(bool)
    T, zr, zc = fix_mask(U, T0); o = summarize(U, T, n_iter)
    return [{"key": k, "n_iter": n_iter, "axis": ax, "zero_rows_removed": zr, "zero_cols_removed": zc, **v} for ax, v in o.items()]
def job_syn(a):
    amp, rho, seed = a
    sp = importlib.util.spec_from_file_location("t", "tools/S2.0_test.py"); T_ = importlib.util.module_from_spec(sp); sp.loader.exec_module(T_)
    from scipy import ndimage
    n = 800; rng = np.random.default_rng([S.SEED, 950 + seed, int(100 * amp), int(-100 * rho)])
    yy, xx = np.mgrid[0:n, 0:n]; Tm = np.zeros((n, n), bool)
    for cy, cx, r in ((250, 250, 220), (560, 520, 260), (200, 620, 150), (640, 180, 140)): Tm |= (yy - cy) ** 2 + (xx - cx) ** 2 < r ** 2
    z = ndimage.gaussian_filter(rng.normal(size=(n, n)), 10); z /= z.std(); m = 5.0 * np.exp(amp * z)
    def fac():
        if rho == 0: return rng.normal(0, 0.2, n)
        th = (1 - np.sqrt(1 - 4 * rho ** 2)) / (2 * rho); e = rng.normal(0, 1, n + 1); return 0.2 * (e[1:] + th * e[:-1]) / np.sqrt(1 + th ** 2)
    al, be = fac(), fac(); c = (T_.PAR / T_.PAR.mean())[C.parity_class(yy, xx)]
    U = rng.poisson(m * c * np.exp(al[:, None] + be[None, :]) * Tm).astype(float)
    o = summarize(U, Tm, 10); rows = []
    for ax, f in (("rows", al), ("cols", be)):
        q = hp(f); g = np.isfinite(q[:-1]) & np.isfinite(q[1:])
        rows.append({"amp": amp, "rho_true": rho, "seed": seed, "axis": ax, "sd_true": 0.2, "sd_true_hp_sample": float(np.nanstd(q)),
                     "lag1_true_hp_sample": float(np.corrcoef(q[:-1][g], q[1:][g])[0, 1]), **o[ax]})
    return rows
if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    secs = [(k, n) for k in ("A1", "A23", "A4", "A5", "A6") for n in (10, 20)]
    syn = [(amp, rho, s) for amp in (0.2, 0.4) for rho in (0.0, -0.2, -0.3) for s in range(4)]
    with Pool(nproc) as p:
        A = sum(p.map(job_sec, secs, chunksize=1), []); Bc = sum(p.map(job_syn, syn, chunksize=1), [])
    A = pd.DataFrame(A); A.to_csv("results/S2.0/S2.0_PH5_sections.csv", index=False)
    Bc = pd.DataFrame(Bc); Bc.to_csv("results/S2.0/S2.0_PH5_calib.csv", index=False)
    print(A.round(4).to_string(index=False))
    print(Bc.groupby(["amp", "rho_true"])[["sd_true_hp_sample", "sd_repro_hp", "lag1_true_hp_sample", "lag1_repro_hp"]].mean().round(3).to_string())
