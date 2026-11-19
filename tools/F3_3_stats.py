"""tools/F3_3_stats.py - F3.3: statistiche pre-registrate (file dello step, sezione «Attese pre-registrate»).
Unita' indipendente: l'area (nel linfonodo il quarto di 500 um dentro campione e l'area di 500 um fuori, DEC-093 R-9).
- regola V (giudizio di uno scarto): stima +- w, w = t(0.975, df) * errore standard, contro il margine +-h; tre esiti.
- spostamento S: media su celle (tessuto x misura) e insiemi del livello di |media vera dell'insieme - media vera dentro campione| / h.
- degrado G del simulato: media sulle celle non smentite dentro campione di (|S - vero fuori| - |S - vero dentro|) / h.
- correzione per il caso: S* = S - mediana di S sotto nessuno spostamento (L1: permutazione delle aree fra dentro e fuori nella stessa
  sezione; L2: modello parametrico t con la deviazione standard fra aree combinata); intervallo [S - q97.5, S - q2.5] della nulla.
- componenti di varianza fra pazienti e fra aree (stima dei momenti, ANOVA a una via), intervallo da bootstrap parametrico."""
import numpy as np, pandas as pd
from scipy import stats

def vrule(g, w, h):
    if not np.isfinite(g) or not np.isfinite(w): return "non valutabile"
    if abs(g) + w <= h: return "confermata"
    if abs(g) - w > h: return "smentita"
    return "indeterminata"

def pooled_sd(groups):
    """deviazione standard combinata fra aree dentro gli insiemi con almeno 2 aree; df = somma (n - 1)."""
    ss, df = 0.0, 0
    for g in groups:
        g = np.asarray(g, float)
        if len(g) >= 2: ss += ((g - g.mean()) ** 2).sum(); df += len(g) - 1
    return (np.sqrt(ss / df) if df > 0 else np.nan), df

def stat_S(cells):
    v = [abs(np.mean(o) - np.mean(c["in"])) / c["h"] for c in cells for o in c["out"].values()]
    return float(np.mean(v)) if v else np.nan

def stat_G(cells):
    v = [(abs(c["sim"] - np.mean(o)) - abs(c["sim"] - np.mean(c["in"]))) / c["h"] for c in cells if c["good"] for o in c["out"].values()]
    return float(np.mean(v)) if v else np.nan

def null_draw(cells, rng, mode):
    """una realizzazione sotto nessuno spostamento. mode 'perm': aree rimescolate fra dentro e fuori (stessa sezione);
    'param': medie fuori e dentro da t con la deviazione standard combinata della cella."""
    out = []
    for c in cells:
        if mode == "perm":
            pool = np.concatenate([np.asarray(c["in"], float)] + [np.asarray(o, float) for o in c["out"].values()]); p = rng.permutation(pool)
            k = len(c["in"]); new = {"in": p[:k]}; o2 = {}; i = k
            for nm, o in c["out"].items(): o2[nm] = p[i:i + len(o)]; i += len(o)
        else:
            s, df = c["sd"], c["df"]; m0 = float(np.mean(c["in"]))
            new = {"in": np.array([m0 + s / np.sqrt(len(c["in"])) * rng.standard_t(df)])}
            o2 = {nm: np.array([m0 + s / np.sqrt(len(o)) * rng.standard_t(df)]) for nm, o in c["out"].items()}
        out.append({**c, "in": new["in"], "out": o2})
    return out

def corrected(cells, fun, mode, B, seed):
    obs = fun(cells); rng = np.random.default_rng(seed)
    nul = np.array([fun(null_draw(cells, rng, mode)) for _ in range(B)])
    q = np.nanpercentile(nul, [2.5, 50, 97.5])
    return {"obs": obs, "null_median": q[1], "est": obs - q[1], "lo": obs - q[2], "hi": obs - q[0]}

def verdict_band(lo, hi, band):
    if lo >= band[0] and hi <= band[1]: return "confermata"
    if hi < band[0] or lo > band[1]: return "smentita"
    return "indeterminata"

def var_components(groups, B, seed):
    """ANOVA a una via (stima dei momenti); rho = sigma_pazienti / sigma_aree; intervallo percentile da bootstrap parametrico."""
    def est(gs):
        n = np.array([len(g) for g in gs], float); k = len(gs); N = n.sum(); mu = [np.mean(g) for g in gs]; gm = np.concatenate(gs).mean()
        msw = sum(((np.asarray(g) - np.mean(g)) ** 2).sum() for g in gs) / (N - k)
        msb = sum(ni * (m - gm) ** 2 for ni, m in zip(n, mu)) / (k - 1); n0 = (N - (n ** 2).sum() / N) / (k - 1)
        sp2 = max((msb - msw) / n0, 0.0); return np.sqrt(sp2), np.sqrt(msw)
    sp, sa = est(groups); rho = sp / sa if sa > 0 else np.nan; rng = np.random.default_rng(seed); bs = []
    for _ in range(B):
        g2 = [rng.normal(0, sp) + rng.normal(0, sa, len(g)) for g in groups]; a, b = est(g2); bs.append(a / b if b > 0 else np.nan)
    lo, hi = np.nanpercentile(bs, [2.5, 97.5]); return {"sigma_patient": sp, "sigma_area": sa, "rho": rho, "lo": lo, "hi": hi}
