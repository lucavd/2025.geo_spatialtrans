"""tools/S1.5d_explore_sizeform.py - S1.5d, esplorazione dichiarata: forme del modello di taglia valutate sui territori VERI (nessuna simulazione).
Per tessuto e tipo con almeno 200 celle interne; pieghe lascia-fuori-un-ROI (parametri stimati sulle altre 4 ROI, valutati sulla ROI esclusa).
M0 = attuale (log A_N = a + b log A_T + sigma e, sigma costante); M1 = sigma che dipende da log A_T (sigma = exp(c0 + c1 (log A_T - media)),
stima per massima verosimiglianza); M2 = logit(A_N / A_T) = a + b log A_T + sigma e (rapporto limitato a 1). Per ogni forma: frazione attesa
di N/C >= 0.95 (con A_N > A_T contato come >= 0.95), area mediana del nucleo attesa (mediana della miscela, per simulazione con 20 estrazioni),
contro il vero della ROI esclusa. Uscita: results/S1.5d/S1.5d_explore_sizeform.csv"""
import importlib.util, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py"); rng = np.random.default_rng(20261025)
def fit_m1(x, y):
    xm = x.mean(); X = np.c_[np.ones_like(x), x]
    def nll(p):
        mu = X @ p[:2]; s = np.exp(p[2] + p[3] * (x - xm)); return np.sum(np.log(s) + 0.5 * ((y - mu) / s) ** 2)
    b = np.polyfit(x, y, 1); p0 = [b[1], b[0], np.log((y - np.polyval(b, x)).std()), 0.0]
    p = minimize(nll, p0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-6, fatol=1e-8)).x
    return lambda xx: (p[0] + p[1] * xx, np.exp(p[2] + p[3] * (xx - xm))), p[3]
rows = []
for A in cal.ARCH:
    _, d = cal.load(A); e = d[d.interior.astype(bool) & ~d.type.isin(list(cal.UNTYPED))]
    for t, g in e.groupby("type"):
        if len(g) < 200: continue
        for r in sorted(g.roi_id.unique()):
            tr, te = g[g.roi_id != r], g[g.roi_id == r]
            if len(te) < 30 or len(tr) < 100: continue
            x, y = np.log(tr.a_t.values), np.log(tr.a_n.values); xt = np.log(te.a_t.values)
            b = np.polyfit(x, y, 1); s0 = (y - np.polyval(b, x)).std(ddof=2)
            f1, c1 = fit_m1(x, y)
            z = np.log(np.clip(tr.a_n.values / tr.a_t.values, 1e-4, 0.999)); lz = np.log(np.clip(tr.a_n / tr.a_t, 1e-4, 0.999))
            lg = np.log(np.clip(tr.a_n.values / tr.a_t.values, 1e-4, 0.999) / (1 - np.clip(tr.a_n.values / tr.a_t.values, 1e-4, 0.999)))
            b2 = np.polyfit(x, lg, 1); s2 = (lg - np.polyval(b2, x)).std(ddof=2)
            E = rng.standard_normal((20, len(xt)))
            an0 = np.exp(np.polyval(b, xt) + s0 * E); mu1, sg1 = f1(xt); an1 = np.exp(mu1 + sg1 * E)
            q2 = 1 / (1 + np.exp(-(np.polyval(b2, xt) + s2 * E))); an2 = q2 * te.a_t.values
            at = te.a_t.values
            o = dict(archetype=A, cell_type=t, roi_id=r, n=len(te), slope_logsigma=c1,
                     real_nc095=np.mean(te.a_n / te.a_t >= 0.95), real_med_an=np.median(te.a_n),
                     M0_nc095=np.mean(an0 / at >= 0.95), M1_nc095=np.mean(an1 / at >= 0.95), M2_nc095=np.mean(an2 / at >= 0.95),
                     M0_med_an=np.median(an0), M1_med_an=np.median(an1), M2_med_an=np.median(an2),
                     real_q90_an=np.quantile(te.a_n, 0.9), M0_q90_an=np.quantile(an0, 0.9), M1_q90_an=np.quantile(an1, 0.9), M2_q90_an=np.quantile(an2, 0.9),
                     real_q10_an=np.quantile(te.a_n, 0.1), M0_q10_an=np.quantile(an0, 0.1), M1_q10_an=np.quantile(an1, 0.1), M2_q10_an=np.quantile(an2, 0.1))
            rows.append(o)
pd.DataFrame(rows).to_csv("results/S1.5d/S1.5d_explore_sizeform.csv", index=False); print("ok", len(rows))
