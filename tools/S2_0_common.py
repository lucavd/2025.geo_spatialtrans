"""tools/S2_0_common.py — S2.0: funzioni pure delle tre misure (testate in tools/S2.0_test.py).
Pre-registrazione: docs/steps/S2.0.md (cancelli afef06b, ok di Luca a8690b3).

M2, regole d'uscita nel modello in avanti dei bordi di R5 (sorgente T = f * rho sulla sottogriglia di 0.5 µm, scarto delta, a e b profilati,
fattori di parita', Poisson). Sfocatura interna SIG_IN = 0.7 µm (PH-5 di R5). Regole:
  C0  (1 - w) G(sigma) + w E(ell) ovunque (nucleo di R5)
  C1  due popolazioni: (1 - w)[(1 - q) G(SIG_IN) + q G(sigma_m)] + w E(ell) ovunque  (p = q (1 - w) = massa della parte larga)
  C2  zona di bordo: sorgenti pesate pi(d) = exp(-d / D) (d = distanza dal bordo, dentro il tessuto) sul nucleo largo G(sigma_e),
      1 - pi sul nucleo stretto G(SIG_IN); coda w E(ell) per tutte
  C3  trasporto sul vetro: F = (1 - w) G(SIG_IN) * T + w E(ell) * T; la parte di F fuori dal tessuto si sposta ancora con G(sigma_g)
M1, strisce: y_ij ~ m_ij exp(alpha_i + beta_j), m campo del tessuto lisciato senza la riga (colonna) stessa.
M3, indice di dispersione dei bin profondi di ogni nucleo attorno alla media del nucleo corretta per la parita'."""
import os
import numpy as np
from scipy import fft as sfft, optimize, ndimage
from scipy import stats as _stats
import R5_common as C

DT = np.float64 if os.environ.get("S20_F64") == "1" else np.float32     # doppia precisione solo per i test di identita' (G-1)
CT = np.complex128 if DT is np.float64 else np.complex64
SEED = 20261103
B = 2000
SIG_IN = 0.7
SUB, SUB_UM, BIN_UM = C.SUB, C.SUB_UM, C.BIN_UM
RULES = {"C0": ("sigma", "w", "ell", "dr", "dc"), "C1": ("q", "sigma_m", "w", "ell", "dr", "dc"),
         "C2": ("sigma_e", "D", "w", "ell", "dr", "dc"), "C3": ("sigma_g", "w", "ell", "dr", "dc")}
K_AIC = {"C0": 7, "C1": 8, "C2": 8, "C3": 7}     # parametri liberi piu' a e b
LOGB = {"sigma": (np.log(0.2), np.log(10.0)), "D": (np.log(0.5), np.log(1e4)), "sigma_g": (np.log(0.05), np.log(30.0)),
        "ell": (np.log(5.0), np.log(200.0))}
DLT = (-6.0, 6.0)

# ------------------------------------------------------------------ M2: regole d'uscita (finestre di bordo)
def region_extra(P, f_sub):
    """Sorgente in spazio reale, maschera del tessuto e distanza dal bordo (dentro il tessuto), sulla griglia allargata P['pshape']."""
    ps = P["pshape"]; n0, n1 = f_sub.shape
    T = np.zeros(ps, DT); T[:n0, :n1] = C.source_sub(f_sub, P["rho"])
    f = np.zeros(ps, DT); f[:n0, :n1] = f_sub
    din = np.zeros(ps, DT); din[:n0, :n1] = ndimage.distance_transform_edt(f_sub > 0.5) * SUB_UM
    return {"T": T, "f": f, "din": din}

def _G(kr, kc, s):
    return C.gauss_ft(kr, kc, s, s)

def _irfft(F, ps):
    return sfft.irfft2(F, s=ps, workers=1)

def _rfft(x, ps):
    return sfft.rfft2(np.asarray(x, DT), s=ps, workers=1)

def kernel_lin(rule, prm, kr, kc):
    """Nucleo delle regole lineari (C0, C1), nello stesso ordine di operazioni di R5 kernel_ft (identita' esatta C1(q=1) = C0 = R5)."""
    w, ell = prm["w"], prm["ell"]
    if rule == "C0":
        return C.kernel_ft(kr, kc, prm["sigma"], prm["sigma"], w, ell)
    q = prm["q"]
    g = (1.0 - q) * _G(kr, kc, SIG_IN) + q * _G(kr, kc, prm["sigma_m"])
    return g if w <= 0 else (1.0 - w) * g + w * C.exp_ft(kr, kc, ell)

def field_ft(rule, prm, Tft, kr, kc, ps, X=None):
    """Trasformata del campo dei conteggi (prima di bin, scarto e parita') secondo la regola."""
    kr, kc = kr.astype(DT, copy=False), kc.astype(DT, copy=False)
    if rule in ("C0", "C1"):
        return Tft * kernel_lin(rule, prm, kr, kc)
    w, ell = prm["w"], prm["ell"]
    E = C.exp_ft(kr, kc, ell) if w > 0 else 0.0
    if rule == "C2":
        pi = np.where(X["f"] > 0.5, np.exp(-X["din"] / prm["D"]), 0.0).astype(DT)
        Te = _rfft(X["T"] * pi, ps); Ti = Tft - Te
        return (1.0 - w) * (Te * _G(kr, kc, prm["sigma_e"]) + Ti * _G(kr, kc, SIG_IN)) + w * E * Tft
    if rule == "C3":
        Fin = _irfft(Tft * ((1.0 - w) * _G(kr, kc, SIG_IN) + w * E), ps)
        out = Fin * (1.0 - X["f"])
        Lf = Fin * X["f"] + _irfft(_rfft(out, ps) * _G(kr, kc, prm["sigma_g"]), ps)
        return _rfft(Lf, ps)
    raise ValueError(rule)

def _sample(F, ps, region_bins):
    h = _irfft(F, ps)
    j0 = SUB // 2 - 1
    return h[j0: region_bins[0] * SUB: SUB, j0: region_bins[1] * SUB: SUB].astype(np.float64)

def sample_bins(Fft, kr, kc, ps, region_bins, dr, dc):
    """Come R5 expected_bins: integrazione sul bin, scarto delta, campionamento al centro del bin."""
    kr, kc = kr.astype(DT, copy=False), kc.astype(DT, copy=False)
    return _sample(Fft * (C.box_ft(kr, kc) * C.shift_ft(kr, kc, dr - SUB_UM / 2, dc - SUB_UM / 2)).astype(CT), ps, region_bins)

def g_rule(rule, prm, P, X=None):
    kr, kc, ps = P["kr"], P["kc"], P["pshape"]
    if rule in ("C0", "C1"):        # stesso ordine di R5 expected_bins: Tft * (nucleo * bin * scarto)
        k_r, k_c = kr.astype(DT, copy=False), kc.astype(DT, copy=False)
        K = (kernel_lin(rule, prm, k_r, k_c) * C.box_ft(k_r, k_c) * C.shift_ft(k_r, k_c, prm["dr"] - SUB_UM / 2, prm["dc"] - SUB_UM / 2)).astype(CT)
        return _sample(K * P["Tft"], ps, P["region_bins"])     # K * Tft: e' l'ordine che NumPy usa in R5 (elisione del temporaneo); conta a 32 bit
    return sample_bins(field_ft(rule, prm, P["Tft"], kr, kc, ps, X), kr, kc, ps, P["region_bins"], prm["dr"], prm["dc"])

def to_prm(rule, th):
    th = np.asarray(th, float)
    if rule == "C0":
        return {"sigma": float(np.exp(th[0])), "w": float(np.clip(th[1], 0, 1)), "ell": float(np.exp(th[2])), "dr": float(th[3]), "dc": float(th[4])}
    if rule == "C1":
        return {"q": float(np.clip(th[0], 0, 1)), "sigma_m": float(np.exp(th[1])), "w": float(np.clip(th[2], 0, 1)), "ell": float(np.exp(th[3])),
                "dr": float(th[4]), "dc": float(th[5])}
    if rule == "C2":
        return {"sigma_e": float(np.exp(th[0])), "D": float(np.exp(th[1])), "w": float(np.clip(th[2], 0, 1)), "ell": float(np.exp(th[3])),
                "dr": float(th[4]), "dc": float(th[5])}
    if rule == "C3":
        return {"sigma_g": float(np.exp(th[0])), "w": float(np.clip(th[1], 0, 1)), "ell": float(np.exp(th[2])), "dr": float(th[3]), "dc": float(th[4])}
    raise ValueError(rule)

def to_theta(rule, prm):
    lg = lambda k, b: float(np.clip(np.log(prm[k]), *b))
    w = float(np.clip(prm["w"], 0, 1)); ell = lg("ell", LOGB["ell"]); d = [float(np.clip(prm["dr"], *DLT)), float(np.clip(prm["dc"], *DLT))]
    if rule == "C0":
        return [lg("sigma", LOGB["sigma"]), w, ell] + d
    if rule == "C1":
        return [float(np.clip(prm["q"], 0, 1)), lg("sigma_m", LOGB["sigma"]), w, ell] + d
    if rule == "C2":
        return [lg("sigma_e", LOGB["sigma"]), lg("D", LOGB["D"]), w, ell] + d
    return [lg("sigma_g", LOGB["sigma_g"]), w, ell] + d

def bounds(rule):
    w, ell = (0.0, 1.0), LOGB["ell"]
    return {"C0": [LOGB["sigma"], w, ell, DLT, DLT], "C1": [(0.0, 1.0), LOGB["sigma"], w, ell, DLT, DLT],
            "C2": [LOGB["sigma"], LOGB["D"], w, ell, DLT, DLT], "C3": [LOGB["sigma_g"], w, ell, DLT, DLT]}[rule]

def fit_rule(rule, P, X, y, sel, starts, n_opt=3, polish=True):
    """Massima verosimiglianza di Poisson della regola sui bin sel (a, b profilati come R5). starts: lista di dict di parametri.
    Le partenze si valutano tutte; si ottimizza dalle n_opt migliori e da quelle marcate 'keep' (equivalenti a C0 o a PH-5)."""
    yv, cv = y[sel].astype(float), P["c"][sel].astype(float)
    cache = {}
    def obj(th):
        key = tuple(np.round(th, 12))
        if key not in cache:
            g = g_rule(rule, to_prm(rule, th), P, X)[sel]
            a, b, f = C.fit_ab(yv, g, cv)
            cache[key] = (f, a, b)
        return cache[key][0]
    th0 = [to_theta(rule, s) for s in starts]
    f0 = [obj(t) for t in th0]
    order = list(np.argsort(f0))
    pick = order[:n_opt] + [i for i, s in enumerate(starts) if s.get("keep") and i not in order[:n_opt]]
    best = None
    for i in pick:
        res = optimize.minimize(obj, th0[i], method="L-BFGS-B", bounds=bounds(rule),
                                options={"maxiter": 500, "ftol": 1e-9, "gtol": 1e-5, "eps": 1e-3})
        if best is None or res.fun < best.fun:
            best = res
    if polish:
        for meth, opts in (("Nelder-Mead", {"xatol": 1e-4, "fatol": 1e-5, "maxfev": 400}),
                           ("L-BFGS-B", {"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7, "eps": 1e-4})):
            res = optimize.minimize(obj, best.x, method=meth, bounds=bounds(rule), options=opts)
            if res.fun < best.fun - 1e-9:
                best = res
    obj(best.x); f, a, b = cache[tuple(np.round(best.x, 12))]
    out = {"rule": rule, **to_prm(rule, best.x), "a": float(a), "b": float(b), "nll": float(f), "n_bins": int(sel.sum()),
           "aic": float(2 * K_AIC[rule] + 2 * f), "n_starts": len(starts), "n_opt": len(pick), "nfev": len(cache),
           "start_nll_min": float(min(f0))}
    if rule == "C1":
        out["p"] = out["q"] * (1.0 - out["w"])
    return out

def nll_at(rule, prm, P, X, y, sel):
    g = g_rule(rule, prm, P, X)[sel]
    return C.fit_ab(y[sel].astype(float), g, P["c"][sel].astype(float))

def predict_rule(rule, prm, P, X):
    return P["c"] * (prm["a"] * np.maximum(g_rule(rule, prm, P, X), 0.0) + prm["b"])

# ------------------------------------------------------------------ M2: controllo interno (ROI, riferimento dell'immagine)
def make_blur(rule, prm, pos_um, d_img, fmask=None, pad_um=300.0, step=SUB_UM):
    """Chiusura F -> conteggi attesi per bin alle posizioni pos_um (come R5 blur_at: riflessione ai bordi, integrazione sul bin,
    scarto d_img, interpolazione bilineare). fmask: tessuto sulla stessa griglia di F (serve a C2 e C3)."""
    p = int(round(pad_um / step))
    cache = {}
    def setup(shape):
        if "ps" in cache: return
        ps = (shape[0] + 2 * p, shape[1] + 2 * p); kr, kc = C.freq_grid(ps, step)
        w, ell = prm["w"], prm["ell"]
        E = C.exp_ft(kr, kc, ell) if w > 0 else 0.0
        cache.update(ps=ps, kr=kr, kc=kc, E=E, box=C.box_ft(kr, kc, BIN_UM))
        if rule in ("C2", "C3"):
            fm = np.pad(np.asarray(fmask, np.float64), p, mode="reflect")
            cache["f"] = fm
            if rule == "C2":
                din = ndimage.distance_transform_edt(fm > 0.5) * step
                cache["pi"] = np.where(fm > 0.5, np.exp(-din / prm["D"]), 0.0)
    def blur(F):
        setup(F.shape)
        G = np.pad(np.asarray(F, np.float64), p, mode="reflect"); ps, kr, kc, E = cache["ps"], cache["kr"], cache["kc"], cache["E"]
        w = prm["w"]
        if rule == "C0":
            Ft = sfft.rfft2(G, workers=1) * C.kernel_ft(kr, kc, prm["sigma"], prm["sigma"], w, prm["ell"])
        elif rule == "C1":
            q = prm["q"]
            Ft = sfft.rfft2(G, workers=1) * ((1.0 - w) * ((1.0 - q) * _G(kr, kc, SIG_IN) + q * _G(kr, kc, prm["sigma_m"])) + w * E)
        elif rule == "C2":
            T = sfft.rfft2(G, workers=1); Te = sfft.rfft2(G * cache["pi"], workers=1)
            Ft = (1.0 - w) * (Te * _G(kr, kc, prm["sigma_e"]) + (T - Te) * _G(kr, kc, SIG_IN)) + w * E * T
        elif rule == "C3":
            Fin = sfft.irfft2(sfft.rfft2(G, workers=1) * ((1.0 - w) * _G(kr, kc, SIG_IN) + w * E), s=ps, workers=1)
            Lf = Fin * cache["f"] + sfft.irfft2(sfft.rfft2(Fin * (1.0 - cache["f"]), workers=1) * _G(kr, kc, prm["sigma_g"]), s=ps, workers=1)
            Ft = sfft.rfft2(Lf, workers=1)
        h = sfft.irfft2(Ft * cache["box"], s=ps, workers=1)
        q_ = (np.asarray(pos_um) - np.asarray(d_img)[None, :]) / step - 0.5 + p
        return ndimage.map_coordinates(h, [q_[:, 0], q_[:, 1]], order=1, mode="nearest")
    return blur

def nc_simulate_b(S, r, blur, c, U):
    lam = c * np.maximum(blur(C.nc_field(S, r)), 1e-12)
    return _stats.poisson.ppf(U, lam[None, :])

def invert_r_b(target, S, blur, c, U, IX, lo=0.05, hi=50.0, iters=22):
    def stat(r):
        Y = nc_simulate_b(S, r, blur, c, U)
        return float(np.mean([C.star_fast(y, IX) for y in Y]))
    a, b = np.log(lo), np.log(hi)
    fa, fb = stat(lo) - target, stat(hi) - target
    if fa * fb > 0:
        return float("nan"), fa, fb
    for _ in range(iters):
        m = 0.5 * (a + b); fm = stat(np.exp(m)) - target
        if fa * fm <= 0: b, fb = m, fm
        else: a, fa = m, fm
    return float(np.exp(0.5 * (a + b))), fa, fb

def refine_b(S, r, blur, c, L, keep, Dobs, acell_obs, acell_mask):
    lam = c * np.maximum(blur(C.nc_field(S, r)), 1e-12)
    Dsim = C.observed_cell_density(L, keep, lam)
    Q = dict(S); D = S["D"].copy()
    for i, l in enumerate(S["labs"]):
        l = int(l)
        if l in Dobs and l in Dsim and Dsim[l] > 0:
            D[i] = D[i] * Dobs[l] / Dsim[l]
    Q["D"] = D
    a_sim = float(lam[acell_mask].mean()) if acell_mask.any() else 0.0
    Q["acell"] = S["acell"] * (acell_obs / a_sim) if a_sim > 0 else S["acell"]
    return Q

def invert_r_refined_b(target, S, blur, c, U, IX, L, keep, Dobs, acell_obs, acell_mask, n_refine=2):
    """Come R5 invert_r_refined, con la sfocatura della regola."""
    r, fa, fb = invert_r_b(target, S, blur, c, U, IX)
    for _ in range(n_refine):
        if not np.isfinite(r): break
        S = refine_b(S, r, blur, c, L, keep, Dobs, acell_obs, acell_mask)
        r, fa, fb = invert_r_b(target, S, blur, c, U, IX)
    return r, S

# ------------------------------------------------------------------ M1: strisce
def _kernel1d(sigma_bins, drop_center):
    R = int(np.ceil(4 * sigma_bins)); x = np.arange(-R, R + 1)
    k = np.exp(-0.5 * (x / sigma_bins) ** 2)
    if drop_center: k[R] = 0.0
    return k / k.sum()

def _smooth_excl(Z, M, sigma_bins, axis_excl):
    """Convoluzione normalizzata di Z sui bin M con gaussiana separabile; lungo axis_excl il peso del bin centrale e' 0
    (la riga, o la colonna, stessa non entra nel proprio campo del tessuto)."""
    kr = _kernel1d(sigma_bins, axis_excl == 0); kc = _kernel1d(sigma_bins, axis_excl == 1)
    num = ndimage.correlate1d(ndimage.correlate1d(Z * M, kr, axis=0, mode="constant"), kc, axis=1, mode="constant")
    den = ndimage.correlate1d(ndimage.correlate1d(M.astype(float), kr, axis=0, mode="constant"), kc, axis=1, mode="constant")
    return np.where(den > 1e-3, num / np.maximum(den, 1e-3), 0.0)

def stripe_factors(U, T, sigma_um=10.0, n_iter=10, min_bins=20):
    """Fattori moltiplicativi per riga (alpha) e colonna (beta) sui bin del tessuto T: y ~ m exp(alpha_i + beta_j), con m lisciato
    senza la riga (per alpha) o la colonna (per beta). Ritorna alpha, beta (NaN dove meno di min_bins bin), m_row, m_col."""
    s = sigma_um / BIN_UM; Y = U.astype(float); M = T.astype(bool)
    nr, nc = Y.shape
    al, be = np.zeros(nr), np.zeros(nc)
    nrow, ncol = M.sum(1), M.sum(0)
    okr, okc = nrow >= min_bins, ncol >= min_bins
    for _ in range(n_iter):
        Z = Y / np.exp(al[:, None] + be[None, :])
        mr = _smooth_excl(Z, M, s, 0)
        num = (Y * M).sum(1); den = (mr * np.exp(be)[None, :] * M).sum(1)
        al = np.where(okr & (den > 0), np.log(np.maximum(num, 1e-12) / np.maximum(den, 1e-12)), 0.0)
        Z = Y / np.exp(al[:, None] + be[None, :])
        mc = _smooth_excl(Z, M, s, 1)
        num = (Y * M).sum(0); den = (mc * np.exp(al)[:, None] * M).sum(0)
        be = np.where(okc & (den > 0), np.log(np.maximum(num, 1e-12) / np.maximum(den, 1e-12)), 0.0)
        wa, wb = nrow * okr, ncol * okc
        al = al - (al * wa).sum() / wa.sum(); be = be - (be * wb).sum() / wb.sum()
    Z = Y / np.exp(al[:, None] + be[None, :])
    mr, mc = _smooth_excl(Z, M, s, 0), _smooth_excl(Z, M, s, 1)
    return np.where(okr, al, np.nan), np.where(okc, be, np.nan), mr, mc

def half_masks(T, axis):
    """Due meta' disgiunte dei bin del tessuto di ogni riga (axis=0: sinistra/destra) o colonna (axis=1: sopra/sotto),
    divise alla mediana delle posizioni dei bin di tessuto della riga (colonna)."""
    M = T.astype(bool)
    if axis == 1:
        A, Bm = half_masks(M.T, 0); return A.T, Bm.T
    idx = np.arange(M.shape[1])[None, :].repeat(M.shape[0], 0)
    cnt = M.sum(1)
    rank = np.cumsum(M, axis=1) * M
    first = M & (rank >= 1) & (rank <= (cnt // 2)[:, None])
    second = M & (rank > (cnt // 2)[:, None])
    return first, second

def half_factors(U, M_half, m_excl, other, axis):
    """Fattore di ogni riga (axis=0) o colonna (axis=1) dai soli bin di M_half, con campo m e fattori dell'altro asse fissati."""
    Y = U.astype(float)
    if axis == 0:
        num = (Y * M_half).sum(1); den = (m_excl * np.exp(np.nan_to_num(other))[None, :] * M_half).sum(1); n = M_half.sum(1)
    else:
        num = (Y * M_half).sum(0); den = (m_excl * np.exp(np.nan_to_num(other))[:, None] * M_half).sum(0); n = M_half.sum(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.log(num / den)
    return np.where((den > 0) & (num > 0), f, np.nan), n

def periodic_split(f, idx):
    """Parte periodica (media per parita' dell'indice) e resto non periodico."""
    f = np.asarray(f, float); par = np.asarray(idx) % 2
    mu = np.array([np.nanmean(f[par == k]) for k in (0, 1)])
    return mu, f - mu[par]

def repro_sd(a, b):
    """Deviazione standard riproducibile: sqrt(max(cov(a, b), 0)) sulle coppie finite."""
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3: return float("nan")
    cv = np.cov(a[ok], b[ok])[0, 1]
    return float(np.sqrt(max(cv, 0.0)))

def block_boot_sd(a, b, block=25, B=B, seed=SEED):
    """Bootstrap per blocchi consecutivi di indici della deviazione standard riproducibile."""
    ok = np.isfinite(a) & np.isfinite(b); idx = np.flatnonzero(ok)
    blocks = [idx[(idx >= s) & (idx < s + block)] for s in range(0, len(a), block)]
    blocks = [x for x in blocks if len(x)]
    rng = np.random.default_rng(seed); reps = np.empty(B)
    for k in range(B):
        pick = rng.integers(0, len(blocks), len(blocks)); ii = np.concatenate([blocks[j] for j in pick])
        cv = np.cov(a[ii], b[ii])[0, 1]; reps[k] = np.sqrt(max(cv, 0.0))
    return float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5)), reps

def tile_parity(U, T, tb=150, min_frac=0.5, min_half=1000):
    """Rapporti colonne dispari/pari e righe dispari/pari per tessera, stimati due volte (meta' sinistra/destra per le colonne,
    sopra/sotto per le righe) sui bin del tessuto. Ritorna un array di record (tessera, lc1, lc2, lr1, lr2) in logaritmo."""
    Y = U.astype(float); M = T.astype(bool); nr, nc = Y.shape
    rr, cc = np.indices((tb, tb)); out = []
    for i in range(nr // tb):
        for j in range(nc // tb):
            y = Y[i*tb:(i+1)*tb, j*tb:(j+1)*tb]; m = M[i*tb:(i+1)*tb, j*tb:(j+1)*tb]
            if m.mean() < min_frac: continue
            rec = [i, j]
            ok = True
            for half_ax in (1, 0):            # colonne: meta' sinistra/destra; righe: meta' sopra/sotto
                for h in (0, 1):
                    hm = (cc < tb // 2) if half_ax == 1 else (rr < tb // 2)
                    hm = hm if h == 0 else ~hm
                    sel = m & hm
                    if sel.sum() < min_half: ok = False; break
                    par = (cc % 2) if half_ax == 1 else (rr % 2)
                    o, e = y[sel & (par == 1)].mean(), y[sel & (par == 0)].mean()
                    if not (o > 0 and e > 0): ok = False; break
                    rec.append(np.log(o / e))
                if not ok: break
            if ok: out.append(rec)
    return np.array(out, float).reshape(-1, 6)

def boot_rows(stat, n, B=B, seed=SEED):
    rng = np.random.default_rng(seed); reps = np.array([stat(rng.integers(0, n, n)) for _ in range(B)])
    return float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5)), reps

# ------------------------------------------------------------------ M3: dispersione dentro il nucleo
def deep_labels(mask, rr, cc, upp, min_um=1.0):
    """Etichetta del nucleo per i bin il cui centro sta a >= min_um dentro la maschera (0 altrimenti)."""
    d_in = ndimage.distance_transform_edt(mask > 0) * upp
    H, W = mask.shape
    ins = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
    lab = np.zeros(len(rr), np.int64)
    r_, c_ = rr[ins], cc[ins]
    lab[ins] = np.where(d_in[r_, c_] >= min_um, mask[r_, c_], 0)
    return lab

def psi_index(u, lab, c, keep, min_bins=4):
    """Indice di dispersione dei bin profondi: per nucleo X2 = sum (y - m)^2 / m con m = media del nucleo * c / media(c);
    psi = sum X2 / sum (n - 1). Ritorna psi, m medio (per bin), numero di nuclei."""
    keep = np.asarray(keep); sel = np.isin(lab, keep) & (lab > 0)
    l, y, cc = lab[sel], np.asarray(u, float)[sel], np.asarray(c, float)[sel]
    uniq, inv, n = np.unique(l, return_inverse=True, return_counts=True)
    sy = np.bincount(inv, weights=y); sc = np.bincount(inv, weights=cc)
    ok_c = (n >= min_bins) & (sy > 0)
    okb = ok_c[inv]
    m = (sy / n)[inv] * cc / (sc / n)[inv]
    x2 = np.where(okb, (y - m) ** 2 / np.where(okb, m, 1.0), 0.0)
    X2 = np.bincount(inv, weights=x2)
    df = (n - 1) * ok_c
    psi = float(X2[ok_c].sum() / df.sum()) if df.sum() > 0 else float("nan")
    return psi, float(m[okb].mean()) if okb.any() else float("nan"), int(ok_c.sum())

def bin8(U2):
    """Somma 2 -> 8 µm: ogni bin di 8 µm somma i 16 bin di 2 µm (righe e colonne a blocchi di 4)."""
    n0, n1 = U2.shape[0] // 4, U2.shape[1] // 4
    return U2[:n0 * 4, :n1 * 4].reshape(n0, 4, n1, 4).sum(axis=(1, 3))


# ------------------------------------------------------------------ M2: partenze (come pre-registrate)
def starts_for(rule, k2, ph5=None):
    """Partenze dell'adattamento: soluzione C0 (K2 di R5) della finestra, griglia grossolana sui parametri nuovi, partenza stretta
    (PH-5 di R5 se c'e', altrimenti w, ell, scarto di C0). 'keep' = sempre ottimizzata (equivalente a C0 o al nucleo stretto)."""
    base = {k: float(k2[k]) for k in ("w", "ell", "dr", "dc")}
    nar = {k: float(ph5[k]) for k in ("w", "ell", "dr", "dc")} if ph5 else dict(base)
    s0 = float(k2["sigma"])
    if rule == "C1":
        return ([dict(base, q=1.0, sigma_m=s0, keep=True)] + [dict(base, q=q, sigma_m=s) for q in (0.2, 0.5, 0.8) for s in (2.0, 3.5, 5.0)]
                + [dict(nar, q=0.0, sigma_m=3.5, keep=True)])
    if rule == "C2":
        return ([dict(base, sigma_e=s0, D=1e4, keep=True)] + [dict(base, sigma_e=s0, D=D) for D in (2.0, 5.0, 15.0, 50.0)]
                + [dict(nar, sigma_e=s0, D=0.5, keep=True)])
    if rule == "C3":
        return [dict(base, sigma_g=s) for s in (2.0, 5.0, 10.0)] + [dict(nar, sigma_g=0.05, keep=True)]
    raise ValueError(rule)


def half_split_factors(U, T, ax, sigma_um=10.0, n_iter=10):
    """P-1 come pre-registrato: il fattore di ogni riga (ax=0) o colonna (ax=1) si stima due volte, con due adattamenti separati sulle due
    meta' disgiunte dei bin di tessuto (sinistra/destra per le righe, sopra/sotto per le colonne): nessun campo, nessun fattore condiviso.
    Ritorna (f1, f2, n1, n2): fattori delle due meta' (NaN dove mancano) e numero di bin di tessuto per meta'."""
    A, Bm = half_masks(T, ax)
    out = []
    for Mh in (A, Bm):
        al, be, _, _ = stripe_factors(U, Mh, sigma_um, n_iter)
        out.append((al if ax == 0 else be, Mh.sum(1) if ax == 0 else Mh.sum(0)))
    return out[0][0], out[1][0], out[0][1], out[1][1]
