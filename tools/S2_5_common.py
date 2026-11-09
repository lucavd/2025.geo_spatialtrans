"""tools/S2_5_common.py — S2.5: funzioni comuni al lato vero (tools/S2.5_real.py) e al lato simulato (tools/S2.5_stats.py).
Pure, senza stato; definizioni in docs/steps/S2.5.md."""
import numpy as np, pandas as pd
from numpy.fft import rfft2, irfft2
B = 2000
KAP = "results/R7/coverage/R7_kappa_final_v4.csv"; CVF = "results/R7/coverage/R7_cv_factor.csv"

def stats(u, g):
    """UMI medio, frazione di vuoti, geni rilevati medi su un insieme di quadratini o blocchi."""
    u = np.asarray(u, float); g = np.asarray(g, float)
    return {"n": int(len(u)), "mean_umi": float(u.mean()) if len(u) else np.nan, "frac_empty": float((u == 0).mean()) if len(u) else np.nan,
            "mean_genes": float(g.mean()) if len(g) else np.nan}

def aniso(U, T, rmax=4, min_pairs=100):
    """Anisotropia del campo log(1 + UMI) dei blocchi di 8 um in tessuto (T): autocorrelazione con maschera (FFT con imbottitura), tensore
    dei momenti secondi sum rho(d) d d^T / sum rho(d) sulle distanze 0 < |d| <= rmax blocchi con rho > 0 e almeno min_pairs coppie;
    indice = 1 - lambda_min / lambda_max (0 = isotropo), angolo dell'asse principale in gradi [0, 180) (0 = lungo le colonne, asse x)."""
    V = np.log1p(np.asarray(U, float)); T = np.asarray(T, bool)
    if T.sum() < 2 * min_pairs: return np.nan, np.nan
    x = np.where(T, V - V[T].mean(), 0.0); m = T.astype(float)
    sh = (2 * V.shape[0], 2 * V.shape[1])
    F = rfft2(x, sh); Fm = rfft2(m, sh)
    num = irfft2(F * np.conj(F), sh); den = irfft2(Fm * np.conj(Fm), sh)
    var = (x ** 2).sum() / m.sum()
    S = np.zeros((2, 2)); w = 0.0
    for dy in range(-rmax, rmax + 1):
        for dx in range(-rmax, rmax + 1):
            if (dy, dx) == (0, 0) or dy * dy + dx * dx > rmax * rmax: continue
            n = den[dy % sh[0], dx % sh[1]]
            if n < min_pairs: continue
            rho = num[dy % sh[0], dx % sh[1]] / n / var
            if rho <= 0: continue
            d = np.array([dy, dx], float); S += rho * np.outer(d, d); w += rho
    if w <= 0: return np.nan, np.nan
    ev, vec = np.linalg.eigh(S / w)
    return float(1 - ev[0] / ev[1]), float(np.degrees(np.arctan2(vec[0, 1], vec[1, 1])) % 180)

def corrected_interval(v, A, rng, side="real", q="bin"):
    """Intervallo di R7 (DEC-060, DEC-061): bootstrap percentile della media sulle aree (B = 2000), semi-ampiezze per kappa_boot del tessuto
    (lato vero) o della riga sim_5 (lato simulato), per il fattore dell'asimmetria se CV > 0.45 (solo lato vero); limite inferiore >= 0."""
    v = np.asarray(v, float); est = float(v.mean())
    bs = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(B)])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    K = pd.read_csv(KAP).set_index(["archetype", "unit"])
    k = float(K.loc[("all", "sim_5") if side == "simulated" else (A, "area"), f"kappa_boot_{q}"])
    cv = float(v.std(ddof=1) / abs(est)) if est != 0 else np.nan; f = 1.0
    if side == "real" and np.isfinite(cv) and cv > 0.45:
        C = pd.read_csv(CVF); f = float(C.cv_factor.iloc[-1]) if cv > 1.0 else float(np.interp(cv, C.cv.values, C.cv_factor.values))
    lo_c, hi_c = est - k * f * (est - lo), est + k * f * (hi - est)
    return {"est": est, "lo_boot": float(lo), "hi_boot": float(hi), "kappa": k, "cv": cv, "cv_factor": f, "lo": max(0.0, lo_c), "hi": hi_c, "n": int(len(v))}

def verdict(s_lo, s_hi, r_lo, r_hi):
    """Tre esiti del protocollo: intervallo del simulato tutto dentro l'intervallo di rilevanza, tutto fuori, a cavallo."""
    if not all(np.isfinite([s_lo, s_hi, r_lo, r_hi])): return "non valutabile"
    if s_lo >= r_lo and s_hi <= r_hi: return "confermata"
    if s_hi < r_lo or s_lo > r_hi: return "smentita"
    return "indeterminata"
