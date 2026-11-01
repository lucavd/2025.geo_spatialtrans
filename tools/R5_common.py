"""tools/R5_common.py — R5: funzioni pure del modello in avanti dei bordi (testate in tools/R5_test.py).
Pre-registrazione: docs/steps/R5.md (cancelli c5c9fce, ok di Luca ceca77b).

Convenzioni. Coordinate in µm nel sistema della griglia (righe, colonne di array); bin di 2 µm, sottogriglia di 0.5 µm (4 x 4 per bin).
La regione di lavoro parte dallo spigolo di un bin. Il campo dei conteggi è il campo H&E sfocato e spostato di delta:
lambda(x) = c_p(x) * [a * (K * T)(x - delta) + b]."""
import numpy as np
from scipy import fft as sfft, optimize, ndimage

SEED = 20261031
B = 2000
BIN_UM = 2.0
SUB = 4
SUB_UM = BIN_UM / SUB
RHO_SIGMA_UM = 10.0          # lisciatura della densita' del tessuto (convoluzione normalizzata)
EDGE_IN_UM = 10.0            # bin di tessuto entro 10 um dal bordo entrano nella verosimiglianza
INTERIOR_UM = 15.0           # variante provata e scartata in G-1 (densita' dai soli bin a >= 15 µm dal bordo; R5_RHO_INTERIOR=1)
FAR_GLASS_UM = 200.0         # variante provata e scartata in G-1 (fondo fissato dal vetro lontano); sui dati veri resta come sensibilita' S-0
BOUNDS = {"sigma": (0.2, 10.0), "w": (0.0, 0.5), "ell": (5.0, 200.0), "delta": (-6.0, 6.0)}

# ------------------------------------------------------------------ Fourier
def freq_grid(shape, d=SUB_UM):
    """Frequenze (cicli/µm) per rfft2 su una griglia `shape` di passo d (precisione singola)."""
    kr = sfft.fftfreq(shape[0], d).astype(np.float32)[:, None]
    kc = sfft.rfftfreq(shape[1], d).astype(np.float32)[None, :]
    return kr, kc

def gauss_ft(kr, kc, sigma_r, sigma_c):
    return np.exp(-2.0 * np.pi ** 2 * (sigma_r ** 2 * kr ** 2 + sigma_c ** 2 * kc ** 2))

def exp_ft(kr, kc, ell):
    """FT del nucleo esponenziale 2D normalizzato E(r) = exp(-r/ell) / (2 pi ell^2)."""
    return (1.0 + (2.0 * np.pi * ell) ** 2 * (kr ** 2 + kc ** 2)) ** -1.5

def kernel_ft(kr, kc, sigma_r, sigma_c, w=0.0, ell=30.0):
    g = gauss_ft(kr, kc, sigma_r, sigma_c)
    return g if w <= 0 else (1.0 - w) * g + w * exp_ft(kr, kc, ell)

def box_ft(kr, kc, width=BIN_UM):
    """Integrazione sul quadratino (media su un quadrato di lato width)."""
    return np.sinc(width * kr) * np.sinc(width * kc)

def shift_ft(kr, kc, dr, dc):
    """f(x - d) <-> F(k) exp(-2 pi i k.d)."""
    return np.exp(-2j * np.pi * (kr * dr + kc * dc))

def pad_shape(region_sub, pad_sub):
    return (region_sub[0] + pad_sub, region_sub[1] + pad_sub)

def source_ft(T_sub, pshape):
    """rfft2 della sorgente sulla sottogriglia (precisione singola), con zeri fino a pshape (contro l'avvolgimento)."""
    return sfft.rfft2(np.asarray(T_sub, np.float32), s=pshape, workers=1)

def expected_bins(Tft, kr, kc, pshape, region_bins, sigma_r, sigma_c, w, ell, dr, dc):
    """(K * T)(x - delta) mediato su ogni bin, per i bin della regione (righe, colonne di bin).
    Il centro del bin cade fra i sottopixel SUB/2-1 e SUB/2: si campiona il sottopixel SUB/2-1 dopo uno spostamento di -SUB_UM/2."""
    F = Tft * (kernel_ft(kr, kc, sigma_r, sigma_c, w, ell) * box_ft(kr, kc) * shift_ft(kr, kc, dr - SUB_UM / 2, dc - SUB_UM / 2)).astype(np.complex64)
    h = sfft.irfft2(F, s=pshape, workers=1)
    j0 = SUB // 2 - 1
    return h[j0: region_bins[0] * SUB: SUB, j0: region_bins[1] * SUB: SUB].astype(np.float64)

# ------------------------------------------------------------------ sorgente
def bin_occupancy(f_sub):
    """Frazione di tessuto per bin dalla sottogriglia."""
    n0, n1 = f_sub.shape[0] // SUB, f_sub.shape[1] // SUB
    return f_sub[: n0 * SUB, : n1 * SUB].reshape(n0, SUB, n1, SUB).mean(axis=(1, 3))

def parity_class(rows, cols):
    return (np.asarray(rows) % 2) * 2 + (np.asarray(cols) % 2)

def parity_factors(U, cls, sel):
    """Fattori di cattura per le 4 classi di parita' (riga, colonna), media delle classi = 1."""
    m = np.array([U[sel & (cls == k)].mean() if np.any(sel & (cls == k)) else np.nan for k in range(4)])
    if np.any(~np.isfinite(m)) or np.nanmean(m) <= 0:
        return np.ones(4)
    return m / m.mean()

def rho_field(U_over_c, full_tissue, sigma_um=RHO_SIGMA_UM):
    """Densita' del tessuto per bin (UMI per bin) per convoluzione normalizzata sui bin interamente di tessuto."""
    s = sigma_um / BIN_UM
    m = full_tissue.astype(float)
    num = ndimage.gaussian_filter(U_over_c * m, s, mode="constant")
    den = ndimage.gaussian_filter(m, s, mode="constant")
    fallback = U_over_c[full_tissue].mean() if full_tissue.any() else 0.0
    return np.where(den > 1e-3, num / np.maximum(den, 1e-3), fallback)

def source_sub(f_sub, rho_bins):
    return f_sub * np.kron(rho_bins, np.ones((SUB, SUB)))

def edge_distance_bins(f_sub):
    """Distanza con segno (µm) del centro di ogni bin dal bordo del tessuto: positiva sul vetro, negativa nel tessuto."""
    t = f_sub > 0.5
    d_out = ndimage.distance_transform_edt(~t) * SUB_UM
    d_in = ndimage.distance_transform_edt(t) * SUB_UM
    d = np.where(t, -d_in, d_out)
    n0, n1 = f_sub.shape[0] // SUB, f_sub.shape[1] // SUB
    j0 = SUB // 2   # sottopixel subito dopo il centro del bin
    return d[j0: n0 * SUB: SUB, j0: n1 * SUB: SUB]

# ------------------------------------------------------------------ verosimiglianza
def fit_a(y, g, c, b, iters=60):
    """MLE di Poisson di a > 0 con b fissato (D-2): lambda = c (a g + b), Newton 1D su un problema convesso."""
    y = np.asarray(y, float); g = np.maximum(np.asarray(g, float), 0.0); c = np.asarray(c, float)
    a = max((y.sum() - (c * b).sum()) / max((c * g).sum(), 1e-12), 1e-6)
    nll = lambda a: float(np.sum(c * (a * g + b) - y * np.log(np.maximum(c * (a * g + b), 1e-300))))
    f = nll(a)
    for _ in range(iters):
        lam = np.maximum(c * (a * g + b), 1e-300)
        gr = np.sum(c * g * (1 - y / lam)); h = np.sum(y * (c * g) ** 2 / lam ** 2) + 1e-12
        t = 1.0
        while t > 1e-8:
            an = max(a - t * gr / h, 1e-12); fn = nll(an)
            if fn <= f + 1e-12: break
            t *= 0.5
        step = abs(an - a); a, f = an, fn
        if step <= 1e-13 * max(a, 1e-12): break
    return a, b, f

def fit_ab(y, g, c, a0=None, b0=None, iters=60):
    """MLE di Poisson di (a > 0, b >= 0) per lambda = c (a g + b): problema convesso, Newton proiettato."""
    y = np.asarray(y, float); g = np.asarray(g, float); c = np.asarray(c, float)
    g = np.maximum(g, 0.0)
    a = a0 if a0 else max(y.sum() / max((c * g).sum(), 1e-12), 1e-9)
    b = b0 if b0 else max(1e-3, 0.1 * y.mean())
    def nll(a, b):
        lam = c * (a * g + b)
        return float(np.sum(lam - y * np.log(np.maximum(lam, 1e-300))))
    f = nll(a, b)
    for _ in range(iters):
        lam = np.maximum(c * (a * g + b), 1e-300)
        r = 1.0 - y / lam
        ga, gb = np.sum(c * g * r), np.sum(c * r)
        w2 = y * c * c / lam ** 2
        haa, hab, hbb = np.sum(w2 * g * g), np.sum(w2 * g), np.sum(w2)
        H = np.array([[haa, hab], [hab, hbb]]) + 1e-12 * np.eye(2)
        try:
            step = np.linalg.solve(H, [ga, gb])
        except np.linalg.LinAlgError:
            break
        t = 1.0
        while t > 1e-8:
            an, bn = a - t * step[0], b - t * step[1]
            an, bn = max(an, 1e-12), max(bn, 1e-9)
            fn = nll(an, bn)
            if fn <= f + 1e-12:
                break
            t *= 0.5
        small = abs(an - a) <= 1e-13 * max(a, 1e-12) and abs(bn - b) <= 1e-13 * max(b, 1e-9)
        a, b, f = an, bn, fn
        if small:
            break
    return a, b, f

MODELS = {"K1": ("sigma", "dr", "dc"), "K2": ("sigma", "w", "ell", "dr", "dc"), "K3": ("sigma_r", "sigma_c", "w", "ell", "dr", "dc")}

def _to_params(model, th):
    p = dict(zip(MODELS[model], th))
    if model == "K3":
        sr, sc = np.exp(p["sigma_r"]), np.exp(p["sigma_c"])
    else:
        sr = sc = np.exp(p["sigma"])
    w = 0.0 if model == "K1" else BOUNDS["w"][1] / (1.0 + np.exp(-p["w"]))
    ell = 30.0 if model == "K1" else np.exp(p["ell"])
    return sr, sc, w, ell, p["dr"], p["dc"]

def _theta0(model, start):
    s = start
    if model == "K1":
        return [np.log(s["sigma"]), s["dr"], s["dc"]]
    lw = np.log(max(s["w"], 1e-4) / (BOUNDS["w"][1] - max(s["w"], 1e-4)))
    if model == "K2":
        return [np.log(s["sigma"]), lw, np.log(s["ell"]), s["dr"], s["dc"]]
    return [np.log(s.get("sigma_r", s["sigma"])), np.log(s.get("sigma_c", s["sigma"])), lw, np.log(s["ell"]), s["dr"], s["dc"]]

def _bounds(model):
    lg = (np.log(BOUNDS["sigma"][0]), np.log(BOUNDS["sigma"][1])); d = BOUNDS["delta"]
    lw = (-12.0, 12.0); le = (np.log(BOUNDS["ell"][0]), np.log(BOUNDS["ell"][1]))
    return {"K1": [lg, d, d], "K2": [lg, lw, le, d, d], "K3": [lg, lg, lw, le, d, d]}[model]

def fit_model(model, Tft, kr, kc, pshape, region_bins, y, c, sel, starts, b_fixed=None):
    """Massima verosimiglianza di Poisson sui bin `sel` (maschera della regione); a, b profilati. Ritorna il migliore fra gli avvii."""
    yv, cv = y[sel].astype(float), c[sel].astype(float)
    cache = {}
    def obj(th):
        key = tuple(np.round(th, 12))
        if key in cache:
            return cache[key][0]
        sr, sc, w, ell, dr, dc = _to_params(model, th)
        g = expected_bins(Tft, kr, kc, pshape, region_bins, sr, sc, w, ell, dr, dc)[sel]
        a, b, f = fit_ab(yv, g, cv) if b_fixed is None else fit_a(yv, g, cv, b_fixed)
        cache[key] = (f, a, b)
        return f
    best = None
    for st in starts:
        res = optimize.minimize(obj, _theta0(model, st), method="L-BFGS-B", bounds=_bounds(model),
                                options={"maxiter": 500, "ftol": 1e-9, "gtol": 1e-5, "eps": 1e-3})
        if best is None or res.fun < best.fun:
            best = res
    sr, sc, w, ell, dr, dc = _to_params(model, best.x)
    obj(best.x)
    f, a, b = cache[tuple(np.round(best.x, 12))]
    return {"model": model, "sigma_r": float(sr), "sigma_c": float(sc), "sigma": float(np.sqrt(sr * sc)), "w": float(w), "ell": float(ell),
            "dr": float(dr), "dc": float(dc), "a": float(a), "b": float(b), "nll": float(f), "converged": bool(best.success),
            "nit": int(best.nit), "nfev": int(best.nfev), "n_bins": int(sel.sum()), "b_fixed": b_fixed is not None}

def predict(fit, Tft, kr, kc, pshape, region_bins, c):
    g = expected_bins(Tft, kr, kc, pshape, region_bins, fit["sigma_r"], fit["sigma_c"], fit["w"], fit["ell"], fit["dr"], fit["dc"])
    return c * (fit["a"] * np.maximum(g, 0.0) + fit["b"])

# ------------------------------------------------------------------ incertezza e verdetti
def boot_median(values, B=B, rng=None):
    """Mediana delle unita' e intervallo percentile 2.5-97.5 dal bootstrap delle unita'."""
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    rng = rng or np.random.default_rng(SEED)
    idx = rng.integers(0, len(v), size=(B, len(v)))
    reps = np.median(v[idx], axis=1)
    return float(np.median(v)), float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5)), reps

def verdict(lo, hi, a, b):
    """Tre esiti del protocollo per un intervallo di rilevanza [a, b] (confermata se tutto dentro)."""
    if lo >= a and hi <= b:
        return "confermata"
    if hi < a or lo > b:
        return "smentita"
    return "indeterminata"

# ------------------------------------------------------------------ sezioni sintetiche (G-1)
def synth_tile(rng, sigma_r, sigma_c, w, ell, dr, dc, region_bins=(300, 300), a=1.0, b=0.3, r_nc=2.0, nb_size=5.0,
               parity=(0.6, 1.0, 0.9, 1.5), cell_density=1 / 120.0, mean_bin=20.0, geometry="blob"):
    """Regione sintetica con bordo curvo, un lumo e un frammento sul vetro; cellule (nucleo r_nc volte piu' denso del citoplasma,
    taglia lognormale); conteggi binomiali negativi del modello in avanti. Ritorna dict con f_sub (tessuto), T_true, y (bin), c (bin), rows, cols."""
    n0, n1 = region_bins[0] * SUB, region_bins[1] * SUB
    yy, xx = (np.mgrid[0:n0, 0:n1] + 0.5) * SUB_UM
    H0, H1 = n0 * SUB_UM, n1 * SUB_UM
    if geometry == "blob":       # tessuto a isola: bordi in tutte le direzioni dentro la tessera centrale
        tissue = (yy - H0 / 2) ** 2 + (xx - H1 / 2) ** 2 < 130.0 ** 2
        tissue &= (yy - H0 / 2 - 40) ** 2 + (xx - H1 / 2 + 30) ** 2 > 25.0 ** 2      # lumo
        tissue |= (yy - H0 / 2 + 60) ** 2 + (xx - H1 / 2 - 210) ** 2 < 15.0 ** 2     # frammento sul vetro
    else:                        # bordo curvo quasi orizzontale: solo la componente normale dello scarto e' determinata
        tissue = yy > H0 / 2 + 50 * np.sin(2 * np.pi * xx / 400.0)
        tissue &= (yy - (H0 / 2 + 120)) ** 2 + (xx - H1 / 2) ** 2 > 25.0 ** 2
        tissue |= (yy - (H0 / 2 - 110)) ** 2 + (xx - H1 / 3) ** 2 < 15.0 ** 2
    # centri uniformi sul tessuto dilatato del raggio della cellula, poi ritaglio: densita' attesa uniforme fino al bordo
    # (il modello assume T = f * rho con rho liscia; una rampa di densita' al bordo e' una violazione, studiata a parte)
    dil = ndimage.binary_dilation(tissue, iterations=int(np.ceil(5.5 / SUB_UM)))
    n_cells = rng.poisson(cell_density * dil.mean() * H0 * H1)
    T = np.zeros((n0, n1))
    cand = np.argwhere(dil)
    pick = cand[rng.integers(0, len(cand), size=n_cells)]
    L = rng.lognormal(0.0, 0.5, size=n_cells)
    rc, rn = 5.5 / SUB_UM, 3.0 / SUB_UM
    R = int(np.ceil(rc)) + 1
    oy, ox = np.mgrid[-R:R + 1, -R:R + 1]
    dist = np.hypot(oy, ox)
    cyto = (dist <= rc) & (dist > rn); nucl = dist <= rn
    stamp_base = cyto * 1.0 + nucl * r_nc
    for (py, px), l in zip(pick, L):
        y0, y1, x0, x1 = py - R, py + R + 1, px - R, px + R + 1
        sy0, sx0 = max(0, -y0), max(0, -x0)
        y0c, x0c, y1c, x1c = max(y0, 0), max(x0, 0), min(y1, n0), min(x1, n1)
        st = stamp_base[sy0: sy0 + (y1c - y0c), sx0: sx0 + (x1c - x0c)]
        T[y0c:y1c, x0c:x1c] += l * st
    T *= tissue
    T += 0.05 * tissue * T[tissue].mean()                                            # matrice
    T *= mean_bin / max(bin_occupancy(T).sum() / max(bin_occupancy(tissue.astype(float)).sum(), 1), 1e-9)
    pshape = pad_shape((n0, n1), 900)
    kr, kc = freq_grid(pshape)
    g = expected_bins(source_ft(T, pshape), kr, kc, pshape, region_bins, sigma_r, sigma_c, w, ell, dr, dc)
    rows, cols = np.mgrid[0:region_bins[0], 0:region_bins[1]]
    c = np.asarray(parity)[parity_class(rows, cols)]
    c = c / np.asarray(parity).mean()
    lam = c * (a * np.maximum(g, 0) + b)
    y = rng.negative_binomial(nb_size, nb_size / (nb_size + lam)) if nb_size else rng.poisson(lam)
    return {"f_sub": tissue.astype(float), "T_true": T, "y": y, "c_true": c, "lam": lam, "rows": rows, "cols": cols}

def prepare_region(y, f_sub, rows, cols, pad_sub=900):
    """Dati di una regione per il modello: parita' stimata, densita' lisciata, FFT della sorgente, distanze, occupazione."""
    occ = bin_occupancy(f_sub)
    dist = edge_distance_bins(f_sub)
    full = occ >= 1.0 - 1e-9
    cls = parity_class(rows, cols)
    pf = parity_factors(y.astype(float), cls, full & (dist <= -EDGE_IN_UM))
    c = pf[cls]
    import os as _os
    rho = rho_field(y / c, full & (dist <= -INTERIOR_UM) if _os.environ.get("R5_RHO_INTERIOR") == "1" else full)   # come pre-registrato
    T = source_sub(f_sub, rho)
    pshape = pad_shape(f_sub.shape, pad_sub)
    kr, kc = freq_grid(pshape)
    return {"occ": occ, "dist": dist, "c": c, "pf": pf, "rho": rho, "Tft": source_ft(T, pshape), "kr": kr, "kc": kc, "pshape": pshape,
            "region_bins": occ.shape}

def likelihood_mask(occ, dist, tile, glass_only=False):
    """Bin della tessera (slice) usati: occupazione < 1, piu' il tessuto entro EDGE_IN_UM dal bordo (salvo glass_only)."""
    m = np.zeros(occ.shape, bool)
    t = np.zeros(occ.shape, bool); t[tile] = True
    m = occ < 1.0 - 1e-9
    if not glass_only:
        m |= (occ >= 1.0 - 1e-9) & (dist >= -EDGE_IN_UM)
    return m & t

STARTS = {"K1": [{"sigma": 1.5, "dr": 0.0, "dc": 0.0}, {"sigma": 4.0, "dr": 0.0, "dc": 0.0}],
          "K2": [{"sigma": 1.5, "w": 0.05, "ell": 20.0, "dr": 0.0, "dc": 0.0}, {"sigma": 1.5, "w": 0.05, "ell": 80.0, "dr": 0.0, "dc": 0.0}]}

def rl_update(P, fit, y, f_sub):
    """Correzione della densita' del tessuto (Richardson-Lucy alla scala di RHO_SIGMA_UM): rho <- rho * G(y/c) / G(previsto/c) sui bin
    interamente di tessuto, cosi' la sorgente non eredita l'impoverimento del bordo dovuto alla diffusione stessa."""
    g = expected_bins(P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], fit["sigma_r"], fit["sigma_c"], fit["w"], fit["ell"], fit["dr"], fit["dc"])
    pred = fit["a"] * np.maximum(g, 0.0) + fit["b"]
    full = P["occ"] >= 1.0 - 1e-9
    s = RHO_SIGMA_UM / BIN_UM
    num = ndimage.gaussian_filter((y / P["c"]) * full, s, mode="constant")
    den = ndimage.gaussian_filter(pred * full, s, mode="constant")
    ratio = np.where(den > 1e-6, num / np.maximum(den, 1e-6), 1.0)
    rho = P["rho"] * np.clip(ratio, 0.2, 5.0)
    T = source_sub(f_sub, rho)
    Q = dict(P); Q["rho"] = rho; Q["Tft"] = source_ft(T, P["pshape"])
    return Q

def far_glass_background(y, f_sub, c, min_um=FAR_GLASS_UM):
    """Fondo per bin (prima della parita') dai bin di vetro a >= min_um dal tessuto della regione (solo per i test sintetici;
    sui dati veri il fondo si stima sulla sezione intera)."""
    d = edge_distance_bins(f_sub)
    m = d >= min_um
    return float((y[m] / c[m]).mean()) if m.any() else float("nan")

def fit_region(y, f_sub, rows, cols, tile, model="K2", n_rl=0, glass_only=False, starts=None, pad_sub=900, P=None, sel=None, b_fixed=None):
    """Adattamento completo di una regione: sorgente, verosimiglianza sulla tessera, n_rl passi di correzione della densita'."""
    P = P if P is not None else prepare_region(y, f_sub, rows, cols, pad_sub)
    sel = sel if sel is not None else likelihood_mask(P["occ"], P["dist"], tile, glass_only)
    fit = fit_model(model, P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, starts or STARTS[model], b_fixed)
    for _ in range(n_rl):
        P = rl_update(P, fit, y, f_sub)
        fit = fit_model(model, P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, [fit], b_fixed)
    fit["n_rl"] = n_rl
    return fit, P, sel

# ------------------------------------------------------------------ rapporto nucleo/corona prima della diffusione (BL-196)
from scipy import stats as _stats

def labels_to_grid(mask_native, upp, step=SUB_UM):
    """Etichette dei nuclei sulla griglia di step µm (vicino piu' prossimo) nel riferimento dell'immagine della ROI."""
    H, W = mask_native.shape
    n0, n1 = int(np.floor(H * upp / step)), int(np.floor(W * upp / step))
    ir = np.minimum(((np.arange(n0) + 0.5) * step / upp).astype(int), H - 1)
    ic = np.minimum(((np.arange(n1) + 0.5) * step / upp).astype(int), W - 1)
    return mask_native[np.ix_(ir, ic)]

def expand_labels_grid(nuc_grid, d_exp_um, step=SUB_UM):
    """Corona: ogni pixel fuori dai nuclei prende l'etichetta del nucleo piu' vicino se il bordo e' entro d_exp_um (come R6 'exp')."""
    dist, (ir, ic) = ndimage.distance_transform_edt(nuc_grid == 0, return_indices=True)
    lab = nuc_grid[ir, ic]
    lab = np.where(dist * step <= d_exp_um, lab, 0)
    return np.where(nuc_grid > 0, nuc_grid, lab)

def nc_setup(nuc_grid, exp_grid, cell_density, ref_density, acell_density):
    """Prepara le grandezze per cellula. cell_density: dict etichetta -> UMI per bin osservati (nucleo + corona);
    le etichette assenti prendono ref_density. Ritorna array per pixel usati da nc_field."""
    labs = np.unique(exp_grid[exp_grid > 0])
    idx = np.zeros(int(exp_grid.max()) + 1, np.int64); idx[labs] = np.arange(len(labs))
    in_cell = exp_grid > 0
    is_nuc = nuc_grid > 0
    pix = idx[exp_grid[in_cell]]
    A_n = np.bincount(pix[is_nuc[in_cell]], minlength=len(labs)).astype(float)
    A_r = np.bincount(pix[~is_nuc[in_cell]], minlength=len(labs)).astype(float)
    D = np.array([cell_density.get(int(l), ref_density) for l in labs], float)
    return {"in_cell": in_cell, "is_nuc_in": is_nuc[in_cell], "pix": pix, "A_n": A_n, "A_r": A_r, "D": D, "acell": float(acell_density),
            "shape": exp_grid.shape, "labs": labs}

def nc_field(S, r):
    """Campo vero (UMI per bin) con densita' r volte piu' alta nel nucleo che nella corona, stesso totale per cellula."""
    rho_c = S["D"] * (S["A_n"] + S["A_r"]) / np.maximum(S["A_r"] + r * S["A_n"], 1e-12)
    F = np.full(S["shape"], S["acell"], float)
    v = rho_c[S["pix"]] * np.where(S["is_nuc_in"], r, 1.0)
    F[S["in_cell"]] = v
    return F

def blur_at(F, pos_um, sigma_r, sigma_c, w, ell, d_img, pad_um=300.0, step=SUB_UM):
    """Campo F (griglia di step µm, riflesso ai bordi per pad_um) sfocato dal nucleo, mediato sul bin di 2 µm, spostato di d_img (µm,
    righe e colonne dell'immagine) e letto alle posizioni pos_um (n x 2, µm dall'origine della ROI) per interpolazione bilineare."""
    p = int(round(pad_um / step))
    G = np.pad(F, p, mode="reflect")
    pshape = pad_shape(G.shape, 0)
    kr, kc = freq_grid(pshape, step)
    Ft = sfft.rfft2(G, workers=1) * kernel_ft(kr, kc, sigma_r, sigma_c, w, ell) * box_ft(kr, kc)
    h = sfft.irfft2(Ft, s=pshape, workers=1)
    q = (np.asarray(pos_um) - np.asarray(d_img)[None, :]) / step - 0.5 + p
    return ndimage.map_coordinates(h, [q[:, 0], q[:, 1]], order=1, mode="nearest")

def per_cell_density(lab_nuc, lab_exp, keep, u, sel=None):
    """Densita' per cellula (UMI per bin) nel nucleo e nella corona dei nuclei tenuti (come R4 PH-1), sui bin sel."""
    sel = np.ones(len(u), bool) if sel is None else sel
    pos = {int(k): i for i, k in enumerate(keep)}
    keep_set = np.asarray(keep)
    out = {}
    for nm, s_, lab in (("nuc", np.isin(lab_nuc, keep_set) & sel, lab_nuc), ("ring", (lab_nuc == 0) & np.isin(lab_exp, keep_set) & sel, lab_exp)):
        ix = np.fromiter((pos[int(v)] for v in lab[s_]), np.int64, count=int(s_.sum()))
        out["u_" + nm] = np.bincount(ix, weights=u[s_], minlength=len(keep)); out["n_" + nm] = np.bincount(ix, minlength=len(keep))
    return out

def cellweighted(pc):
    ok = (pc["n_nuc"] > 0) & (pc["n_ring"] > 0)
    dn = pc["u_nuc"][ok] / pc["n_nuc"][ok]; dr = pc["u_ring"][ok] / pc["n_ring"][ok]
    return float(dn.mean() / dr.mean()) if ok.any() and dr.mean() > 0 else float("nan")

def star_stat(u, L, keep, far=None):
    """Statistica di R4 PH-1: rapporto pesato per cellula con la maschera vera diviso lo stesso con la maschera spostata.
    L: dict con nuc, exp, nuc_shift, exp_shift (etichette per bin); far: dict con far, far_shift (CP-1) o None."""
    a = cellweighted(per_cell_density(L["nuc"], L["exp"], keep, u, None if far is None else far["far"]))
    b = cellweighted(per_cell_density(L["nuc_shift"], L["exp_shift"], keep, u, None if far is None else far["far_shift"]))
    return a / b

def nc_simulate(S, r, pos_um, kern, d_img, c, U):
    """Conteggi simulati (Poisson per inversione della funzione di ripartizione con uniformi fisse U: n_rep x n_bin)."""
    lam = c * np.maximum(blur_at(nc_field(S, r), pos_um, *kern, d_img), 1e-12)
    return _stats.poisson.ppf(U, lam[None, :])

def invert_r(target, S, pos_um, kern, d_img, c, U, L, keep, far=None, lo=0.05, hi=50.0, iters=22):
    """r per cui la media delle repliche della statistica simulata e' uguale all'osservata (bisezione su log r)."""
    def stat(r):
        Y = nc_simulate(S, r, pos_um, kern, d_img, c, U)
        return float(np.mean([star_stat(y, L, keep, far) for y in Y]))
    a, b = np.log(lo), np.log(hi)
    fa, fb = stat(lo) - target, stat(hi) - target
    if fa * fb > 0:
        return float("nan"), fa, fb
    for _ in range(iters):
        m = 0.5 * (a + b); fm = stat(np.exp(m)) - target
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return float(np.exp(0.5 * (a + b))), fa, fb

def synth_roi(rng, side_um=300.0, upp=0.2737, nuc_r_um=3.0, density=1 / 150.0, d_exp=8.0):
    """ROI sintetica: maschera di etichette nativa con nuclei a disco, posizioni dei bin di 2 µm, etichette dei bin (come R4)."""
    from R6_common import footprint_labels
    H = int(side_um / upp)
    mask = np.zeros((H, H), np.int32)
    n = rng.poisson(density * side_um ** 2)
    cy, cx = rng.uniform(0, side_um, n), rng.uniform(0, side_um, n)
    yy, xx = (np.mgrid[0:H, 0:H] + 0.5) * upp
    R = int(np.ceil(nuc_r_um / upp)) + 1
    for i, (y0, x0) in enumerate(zip(cy, cx), start=1):
        r0, c0 = int(y0 / upp), int(x0 / upp)
        sl = (slice(max(r0 - R, 0), min(r0 + R + 1, H)), slice(max(c0 - R, 0), min(c0 + R + 1, H)))
        d = (yy[sl] - y0) ** 2 + (xx[sl] - x0) ** 2 <= nuc_r_um ** 2
        mask[sl][d] = i
    nb = int(np.floor(H * upp / SUB_UM)) // SUB                  # bin interi sulla griglia di 0.5 µm
    br, bc = np.mgrid[0:nb, 0:nb]
    pos = np.c_[(br.ravel() + 0.5) * 2.0, (bc.ravel() + 0.5) * 2.0]
    rr, cc = np.floor(pos[:, 0] / upp).astype(int), np.floor(pos[:, 1] / upp).astype(int)
    keep = np.unique(mask[mask > 0])
    s = int(round(20.0 / upp)); ms = np.roll(np.roll(mask, s, 0), s, 1)
    L = {"nuc": footprint_labels(mask, rr, cc, upp, d_exp, "nuc"), "exp": footprint_labels(mask, rr, cc, upp, d_exp, "exp"),
         "nuc_shift": footprint_labels(ms, rr, cc, upp, d_exp, "nuc"), "exp_shift": footprint_labels(ms, rr, cc, upp, d_exp, "exp")}
    def far_of(mm, lab_nuc):
        d_in = ndimage.distance_transform_edt(mm > 0) * upp; d_out = ndimage.distance_transform_edt(mm == 0) * upp
        return np.where(lab_nuc > 0, d_in[rr, cc] >= 2.0, d_out[rr, cc] >= 2.0)
    far = {"far": far_of(mask, L["nuc"]), "far_shift": far_of(ms, L["nuc_shift"])}
    return {"mask": mask, "upp": upp, "pos": pos, "rr": rr, "cc": cc, "keep": keep, "L": L, "far": far,
            "rows": br.ravel(), "cols": bc.ravel()}

def observed_cell_density(L, keep, u):
    pc = per_cell_density(L["nuc"], L["exp"], keep, u)
    n = pc["n_nuc"] + pc["n_ring"]
    return {int(k): float((pc["u_nuc"][i] + pc["u_ring"][i]) / n[i]) for i, k in enumerate(keep) if n[i] > 0}

def stat_index(L, keep, far=None):
    """Indici precalcolati (cellula per bin o -1) per la statistica di R4 PH-1, con o senza i bin lontani dal bordo (far)."""
    pos = np.full(int(max(L["exp"].max(), L["exp_shift"].max(), L["nuc"].max(), L["nuc_shift"].max(), np.max(keep))) + 2, -1, np.int64)
    pos[np.asarray(keep, np.int64)] = np.arange(len(keep))
    IX = {"n": len(keep)}
    for mk in ("", "_shift"):
        ln, le = L["nuc" + mk], L["exp" + mk]
        s = np.ones(len(ln), bool) if far is None else far["far" + mk]
        inuc = np.where((pos[ln] >= 0) & s, pos[ln], -1)
        iring = np.where((ln == 0) & (pos[le] >= 0) & s, pos[le], -1)
        IX["nuc" + mk], IX["ring" + mk] = inuc, iring
    return IX

def star_fast(u, IX):
    out = []
    for mk in ("", "_shift"):
        a, b = IX["nuc" + mk], IX["ring" + mk]
        ma, mb = a >= 0, b >= 0
        un = np.bincount(a[ma], weights=u[ma], minlength=IX["n"]); nn = np.bincount(a[ma], minlength=IX["n"])
        ur = np.bincount(b[mb], weights=u[mb], minlength=IX["n"]); nr_ = np.bincount(b[mb], minlength=IX["n"])
        ok = (nn > 0) & (nr_ > 0)
        out.append((un[ok] / nn[ok]).mean() / (ur[ok] / nr_[ok]).mean())
    return float(out[0] / out[1])

def invert_r_fast(target, S, pos_um, kern, d_img, c, U, IX, lo=0.05, hi=50.0, iters=22):
    """Come invert_r, con la statistica a indici precalcolati."""
    def stat(r):
        Y = nc_simulate(S, r, pos_um, kern, d_img, c, U)
        return float(np.mean([star_fast(y, IX) for y in Y]))
    a, b = np.log(lo), np.log(hi)
    fa, fb = stat(lo) - target, stat(hi) - target
    if fa * fb > 0:
        return float("nan"), fa, fb
    for _ in range(iters):
        m = 0.5 * (a + b); fm = stat(np.exp(m)) - target
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return float(np.exp(0.5 * (a + b))), fa, fb

def refine_densities(S, r, pos_um, kern, d_img, c, L, keep, Dobs, acell_obs, acell_mask):
    """Passo autoconsistente: le densita' per cellula (e dell'acellulare) del campo vero si correggono finche' i conteggi attesi
    dopo la sfocatura riproducono quelle osservate (le osservate sono gia' sfocate)."""
    lam = c * np.maximum(blur_at(nc_field(S, r), pos_um, *kern, d_img), 1e-12)
    Dsim = observed_cell_density(L, keep, lam)
    Q = dict(S); D = S["D"].copy()
    for i, l in enumerate(S["labs"]):
        l = int(l)
        if l in Dobs and l in Dsim and Dsim[l] > 0:
            D[i] = D[i] * Dobs[l] / Dsim[l]
    Q["D"] = D
    a_sim = float(lam[acell_mask].mean()) if acell_mask.any() else 0.0
    Q["acell"] = S["acell"] * (acell_obs / a_sim) if a_sim > 0 else S["acell"]
    return Q

def invert_r_refined(target, S, pos_um, kern, d_img, c, U, IX, L, keep, Dobs, acell_obs, acell_mask, n_refine=2):
    r, fa, fb = invert_r_fast(target, S, pos_um, kern, d_img, c, U, IX)
    for _ in range(n_refine):
        if not np.isfinite(r): break
        S = refine_densities(S, r, pos_um, kern, d_img, c, L, keep, Dobs, acell_obs, acell_mask)
        r, fa, fb = invert_r_fast(target, S, pos_um, kern, d_img, c, U, IX)
    return r, S
