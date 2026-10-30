"""tools/R4_common.py — R4: funzioni pure (testate in tools/R4_test.py).
Pre-registrazione: docs/steps/R4.md (cancelli ff04e8f, ok di Luca e0c2341)."""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from R6_common import footprint_labels, aggregate_counts, read_h5_csc  # noqa: F401  (riuso dichiarato: stesse impronte di R6)

SEED = 20261030
B = 2000
D_SHIFT_UM = 20.0        # maschera spostata (come R6 CP-2)
INTERIOR_UM = 20.0       # cellule con centro a >= 20 um dal bordo della ROI
MIN_CELLS = 50           # P-2: tipo con almeno 50 cellule nella ROI
MIN_MEAN = 2.0           # P-2: geni con media >= 2 UMI per cellula
DROPOUT = 0.30           # CP-2
MIN_TYPE_FRAC = 0.05     # P-4: tipi con almeno il 5 % dei nuclei
NON_TYPES = ("non classificabile", "non risolto")
STRATA = {0: "acellular", 1: "nucleus", 2: "ring", 3: "other_nucleus"}

# ---------------------------------------------------------------- 2 -> 8 um
def coarse_index(rows, cols, f=4):
    return np.asarray(rows) // f, np.asarray(cols) // f

def aggregate_totals(rows, cols, umi, f=4):
    """Somma dei totali per bin fine sui bin grossi (riga // f, colonna // f)."""
    r, c = coarse_index(rows, cols, f)
    df = pd.DataFrame({"row": r, "col": c, "umi": np.asarray(umi, dtype=np.int64)})
    return df.groupby(["row", "col"], sort=True, as_index=False)["umi"].sum()

def column_totals(X):
    """X csc (geni x bin): UMI e geni rilevati per bin."""
    X = X.tocsc()
    umi = np.asarray(X.sum(axis=0)).ravel().astype(np.int64)
    ngen = np.diff(X.indptr).astype(np.int64)
    return umi, ngen

# ---------------------------------------------------------------- strati
def shifted_mask(mask, um_per_px, shift_um=D_SHIFT_UM):
    s = int(round(shift_um / um_per_px))
    return np.roll(np.roll(mask, s, axis=0), s, axis=1)

def strata(lab_nuc, lab_exp, keep_labels):
    """1 nucleo (bin in un nucleo tenuto); 2 corona (fuori dai nuclei, attribuito a un nucleo tenuto entro d_exp);
    3 altro nucleo (bin in, o attribuito a, un nucleo non tenuto); 0 acellulare (nessun nucleo entro d_exp)."""
    lab_nuc = np.asarray(lab_nuc); lab_exp = np.asarray(lab_exp)
    kn = np.isin(lab_nuc, keep_labels); ke = np.isin(lab_exp, keep_labels)
    s = np.zeros(len(lab_nuc), np.int8)
    s[kn] = 1
    s[(lab_nuc == 0) & ke] = 2
    s[(lab_nuc > 0) & ~kn] = 3
    s[(lab_nuc == 0) & (lab_exp > 0) & ~ke] = 3
    return s

def density_ratio(umi, s):
    """UMI medi per bin nel nucleo / UMI medi per bin nella corona (bin raggruppati nella ROI)."""
    umi = np.asarray(umi, float); s = np.asarray(s)
    n, r = s == 1, s == 2
    if n.sum() == 0 or r.sum() == 0 or umi[r].mean() == 0:
        return float("nan")
    return float(umi[n].mean() / umi[r].mean())

# ---------------------------------------------------------------- zeri in eccesso
def nb_zero_expected(mean, var):
    """P(0) di una binomiale negativa con media e varianza date (momenti); Poisson se var <= media."""
    mean = np.asarray(mean, float); var = np.asarray(var, float)
    out = np.exp(-mean)
    od = var > mean
    theta = mean[od] ** 2 / (var[od] - mean[od])
    out[od] = np.power(1.0 + mean[od] / theta, -theta)
    return out

def weighted_moments(X, w=None):
    """X csr (cellule x geni). w pesi interi per cellula (bootstrap) o None. -> media, varianza (ddof 0), frazione di zeri."""
    X = X.tocsr()
    n_cells = X.shape[0]
    if w is None:
        w = np.ones(n_cells)
    w = np.asarray(w, float); n = w.sum()
    m1 = np.asarray(X.T @ w).ravel() / n
    X2 = X.multiply(X).tocsr()
    m2 = np.asarray(X2.T @ w).ravel() / n
    nz = X.copy(); nz.data = (nz.data > 0).astype(float)
    pnz = np.asarray(nz.T @ w).ravel() / n
    return m1, np.maximum(m2 - m1 ** 2, 0.0), 1.0 - pnz

def eligible_genes(X, min_mean=MIN_MEAN):
    m1, _, _ = weighted_moments(X)
    return np.flatnonzero(m1 >= min_mean)

def excess_zeros(X, genes, w=None):
    """Frazione osservata di cellule con 0 UMI meno la P(0) binomiale negativa, per i geni dati."""
    Xg = X.tocsc()[:, genes]
    m1, var, z = weighted_moments(Xg, w)
    return z - nb_zero_expected(m1, var)

def nb_sf_excess(X, s):
    """PH-2: per gene, P(0) osservata meno P(0) media di una binomiale negativa con media s_i * mu per cellula
    (s = fattore di taglia, media 1) e theta dai momenti con fattore di taglia; Poisson se il residuo e' <= 0."""
    X = X.tocsc(); s = np.asarray(s, float); n = X.shape[0]
    mu = np.asarray(X.sum(axis=0)).ravel() / s.sum()
    x2 = np.asarray(X.multiply(X).T @ np.ones(n)).ravel(); xs = np.asarray(X.T @ s).ravel()
    resid = x2 - 2 * mu * xs + mu ** 2 * (s ** 2).sum() - mu * s.sum()          # sum_i (x_i - s_i mu)^2 - s_i mu
    den = mu ** 2 * (s ** 2).sum()
    pois = resid <= 0
    theta = np.where(pois, 1.0, den / np.where(pois, 1.0, resid))
    nz = X.copy(); nz.data = (nz.data > 0).astype(float); z_obs = 1 - np.asarray(nz.sum(axis=0)).ravel() / n
    lam = np.outer(s, mu)
    p0 = np.where(pois[None, :], np.exp(-lam), np.power(1 + lam / theta[None, :], -theta[None, :]))
    return z_obs - p0.mean(axis=0)

def apply_dropout(X, frac, rng):
    """Azzera a caso una frazione `frac` dei valori non nulli (CP-2)."""
    Y = X.tocsr().copy()
    Y.data = np.where(rng.random(len(Y.data)) < frac, 0, Y.data).astype(Y.data.dtype)
    Y.eliminate_zeros()
    return Y

# ---------------------------------------------------------------- bootstrap e verdetti
def boot_one_stage(values, stat=np.mean, B=B, rng=None):
    """Bootstrap per blocchi: ricampiona le ROI (valori per ROI). -> (stima, lo, hi, repliche)."""
    v = np.asarray(values, float); k = len(v)
    rng = rng if rng is not None else np.random.default_rng(SEED)
    idx = rng.integers(0, k, size=(B, k))
    reps = np.array([stat(v[i]) for i in idx])
    return float(stat(v)), float(np.nanpercentile(reps, 2.5)), float(np.nanpercentile(reps, 97.5)), reps

def boot_two_stage(groups, stat, B=B, rng=None):
    """groups: lista (una per ROI) di array 1D (cellule). Ricampiona ROI, poi cellule dentro ogni ROI estratta;
    stat applicata alle cellule raggruppate. -> (stima, lo, hi, repliche)."""
    rng = rng if rng is not None else np.random.default_rng(SEED)
    k = len(groups)
    est = stat(np.concatenate(groups))
    reps = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, k, size=k)
        parts = [groups[j][rng.integers(0, len(groups[j]), size=len(groups[j]))] for j in pick if len(groups[j])]
        reps[b] = stat(np.concatenate(parts)) if parts else np.nan
    return float(est), float(np.nanpercentile(reps, 2.5)), float(np.nanpercentile(reps, 97.5)), reps

def verdict(lo, hi, a, b):
    """Esito del protocollo per un intervallo [lo, hi] e un intervallo atteso [a, b] (a o b possono essere infiniti)."""
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return "non valutabile"
    if lo >= a and hi <= b:
        return "confermata"
    if hi < a or lo > b:
        return "smentita"
    return "indeterminata"

# ---------------------------------------------------------------- spettro (D-5, BL-089)
def welch_spectrum(img, valid, tile=64, min_valid=1.0):
    """Periodogramma medio (Welch, finestra di Hann) su tessere tile x tile con frazione di bin validi >= min_valid.
    I bin non validi prendono la media dei validi della tessera. -> (spettro medio con l'origine al centro, n tessere)."""
    H, W = img.shape
    win = np.outer(np.hanning(tile), np.hanning(tile))
    acc = np.zeros((tile, tile)); k = 0
    for r0 in range(0, H - tile + 1, tile):
        for c0 in range(0, W - tile + 1, tile):
            v = valid[r0:r0 + tile, c0:c0 + tile]
            if v.mean() < min_valid:
                continue
            t = img[r0:r0 + tile, c0:c0 + tile].astype(float).copy()
            mu = t[v].mean(); t[~v] = mu; t -= mu
            acc += np.abs(np.fft.fftshift(np.fft.fft2(t * win))) ** 2; k += 1
    return (acc / k if k else acc), k

def spectral_peak(S, r_min=3, thr=10.0, off_axis=False):
    """Rapporto massimo fra potenza e mediana della sua corona di frequenze (raggio intero), raggio >= r_min, fino agli angoli
    (comprese le frequenze di Nyquist, +-0.5 cicli per bin: correzione F01 della revisione). off_axis=True esclude le righe e le colonne
    dello spettro entro 2 campioni da fx = 0, fy = 0 o |f| = 0.5 (assi e Nyquist: strisce per riga/colonna e alternanza pari/dispari), per cercare un
    reticolo diagonale (BL-089, F16). -> dict(ratio, present, fy, fx in cicli per bin)."""
    n = S.shape[0]; cy = n // 2
    yy, xx = np.indices(S.shape); rad = np.rint(np.hypot(yy - cy, xx - cy)).astype(int)
    ok = np.ones(S.shape, bool)
    if off_axis:                                   # banda di +-2 campioni attorno ad assi e linee di Nyquist (dispersione della finestra)
        for k in (-2, -1, 0, 1, 2):
            ok[(cy + k) % n, :] = ok[:, (cy + k) % n] = ok[k % n, :] = ok[:, k % n] = False
    ratio = np.zeros_like(S)
    for r in range(r_min, int(rad.max()) + 1):
        m = (rad == r) & ok
        if m.sum() < 4:
            continue
        med = np.median(S[m])
        if med > 0:
            ratio[m] = S[m] / med
    i = np.unravel_index(np.argmax(ratio), S.shape)
    return {"ratio": float(ratio[i]), "present": bool(ratio[i] > thr), "fy": float((i[0] - cy) / n), "fx": float((i[1] - cy) / n)}
