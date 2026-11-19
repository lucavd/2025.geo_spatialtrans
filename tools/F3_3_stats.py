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
    out = []; perms = {}
    for c in cells:
        if mode == "perm":
            pool = np.concatenate([np.asarray(c["in"], float)] + [np.asarray(o, float) for o in c["out"].values()])
            g = c.get("group")
            if g is not None:      # stessa permutazione delle aree per tutte le misure dello stesso tessuto (le aree sono le unita')
                if g not in perms: perms[g] = rng.permutation(len(pool))
                assert len(perms[g]) == len(pool); p = pool[perms[g]]
            else: p = rng.permutation(pool)
            k = len(c["in"]); new = {"in": p[:k]}; o2 = {}; i = k
            for nm, o in c["out"].items(): o2[nm] = p[i:i + len(o)]; i += len(o)
        elif mode == "resid":      # sensibilita' (non pre-registrata): residui per area ricampionati insieme per tutte le misure del tessuto
            g = c["group"]; Rm = c["resid"]; m0 = float(np.mean(c["in"]))
            if g not in perms: perms[g] = {nm: rng.integers(0, Rm.shape[0], n) for nm, n in [("in", len(c["in"]))] + [(k, len(o)) for k, o in c["out"].items()]}
            ix = perms[g]; new = {"in": np.full(len(c["in"]), m0 + Rm[ix["in"], c["mi"]].mean())}
            o2 = {nm: np.full(len(o), m0 + Rm[ix[nm], c["mi"]].mean()) for nm, o in c["out"].items()}
        else:
            s, df = c["sd"], c["df"]; m0 = float(np.mean(c["in"]))
            new = {"in": np.full(len(c["in"]), m0 + s / np.sqrt(len(c["in"])) * rng.standard_t(df))}       # numerosita' conservate (servono a Q_S, Q_G)
            o2 = {nm: np.full(len(o), m0 + s / np.sqrt(len(o)) * rng.standard_t(df)) for nm, o in c["out"].items()}
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


# ---- correzione del 2026-11-19 (cancello G-1 f fallito: S* e G* sottostimano gli spostamenti grandi) ----
def _var_mean(c, n): return c["sd"] ** 2 / n
def stat_QS(cells):
    """spostamento quadratico vero medio, (D^2 - Var D) / h^2, stimatore corretto della media di delta^2 / h^2."""
    v = [((np.mean(o) - np.mean(c["in"])) ** 2 - _var_mean(c, len(o)) - _var_mean(c, len(c["in"]))) / c["h"] ** 2
         for c in cells for o in c["out"].values()]
    return float(np.mean(v)) if v else np.nan
def stat_QG(cells):
    """aumento medio del quadrato dello scarto vero del simulato, fuori meno dentro, sulle celle non smentite dentro campione."""
    v = []
    for c in cells:
        if not c["good"]: continue
        gi = (c["sim"] - np.mean(c["in"])) ** 2 - _var_mean(c, len(c["in"])) - c.get("se_sim", 0.0) ** 2
        for o in c["out"].values():
            go = (c["sim"] - np.mean(o)) ** 2 - _var_mean(c, len(o)) - c.get("se_sim", 0.0) ** 2
            v.append((go - gi) / c["h"] ** 2)
    return float(np.mean(v)) if v else np.nan
def ssqrt(x): return float(np.sign(x) * np.sqrt(abs(x)))
def boot_draw(cells, rng):
    """medie ridisegnate intorno a quelle osservate (modello t con la deviazione standard combinata), per l'intervallo."""
    out = []
    for c in cells:
        s, df = c["sd"], c["df"]
        new_in = np.full(len(c["in"]), np.mean(c["in"]) + s / np.sqrt(len(c["in"])) * rng.standard_t(df))
        o2 = {nm: np.full(len(o), np.mean(o) + s / np.sqrt(len(o)) * rng.standard_t(df)) for nm, o in c["out"].items()}
        out.append({**c, "in": new_in, "out": o2})
    return out
def q_estimate(cells, fun, B, seed, null_mode=None):
    """stima in unita' di h (radice con segno del quadrato medio), intervallo percentile dal bootstrap parametrico; facoltativa la
    distribuzione sotto nessuno spostamento (permutazione o modello t) come controllo."""
    obs = fun(cells); rng = np.random.default_rng(seed)
    bs = np.array([fun(boot_draw(cells, rng)) for _ in range(B)]); lo, hi, med = np.nanpercentile(bs, [2.5, 97.5, 50])
    b = med - obs; lo, hi = lo - b, hi - b      # il ridisegno aggiunge varianza: si toglie lo scarto della mediana (percentile corretto)
    r = {"obs_q": obs, "est": ssqrt(obs), "lo": ssqrt(lo), "hi": ssqrt(hi), "boot_bias_q": b}
    if null_mode:
        rng2 = np.random.default_rng(seed + 1); nul = np.array([fun(null_draw(cells, rng2, null_mode)) for _ in range(B)])
        a, b = np.nanpercentile(nul, [2.5, 97.5]); r.update(null_lo=ssqrt(a), null_hi=ssqrt(b))
    return r
