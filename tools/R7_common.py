"""tools/R7_common.py — R7: funzioni comuni (zone comuni, copertura degli intervalli per area, correzioni, pezzi disgiunti).

Convenzioni: aree di R2 in pixel dell'immagine a piena risoluzione (c0, r0, c1, r1, estremi esclusi a destra);
valori per area in un vettore di lunghezza n (5 aree, 4 finestre). Seme dello step 20261108.
"""
import itertools
import numpy as np
from scipy import stats

SEED = 20261108
NREP = 20000
B = 2000
BRADLEY = (0.925, 0.975)
LEVEL = 0.95


# ---------- geometria ----------
def box_overlap_matrix(boxes):
    """boxes: array (n, 4) di (c0, r0, c1, r1). -> matrice n x n delle quote di area in comune (diagonale 1).
    Per quadrati della stessa area la quota è simmetrica; per aree diverse si usa sqrt(a_i a_j) al denominatore."""
    b = np.asarray(boxes, float); n = len(b)
    F = np.eye(n)
    for i, j in itertools.combinations(range(n), 2):
        w = max(0.0, min(b[i, 2], b[j, 2]) - max(b[i, 0], b[j, 0]))
        h = max(0.0, min(b[i, 3], b[j, 3]) - max(b[i, 1], b[j, 1]))
        ai = (b[i, 2] - b[i, 0]) * (b[i, 3] - b[i, 1]); aj = (b[j, 2] - b[j, 0]) * (b[j, 3] - b[j, 1])
        F[i, j] = F[j, i] = w * h / np.sqrt(ai * aj)
    return F


def unit_overlap_matrix(member):
    """member: array booleano (unità x n) di appartenenza alle aree. -> n_ij / sqrt(n_i n_j)."""
    M = np.asarray(member, float)
    N = M.T @ M
    d = np.sqrt(np.diag(N))
    return N / np.outer(d, d)


def rho_bar(R):
    R = np.asarray(R, float); n = len(R)
    return float(R[~np.eye(n, dtype=bool)].mean())


def k_factor(R):
    """Rapporto fra l'errore standard vero della media di n aree con correlazione R e il valore atteso
    di s/sqrt(n) (E[s^2] = sigma^2 (1 - rho_bar); Var(xbar) = sigma^2 (1 + (n-1) rho_bar) / n)."""
    n = len(R); rb = rho_bar(R)
    return float(np.sqrt((1 + (n - 1) * rb) / (1 - rb)))


# ---------- bootstrap esatto ----------
_W = {}
def boot_weights(n):
    """Pesi delle n^n estrazioni ordinate equiprobabili del bootstrap (righe: estrazioni; colonne: aree)."""
    if n not in _W:
        idx = np.array(list(itertools.product(range(n), repeat=n)))
        _W[n] = np.stack([np.bincount(r, minlength=n) for r in idx]).astype(float) / n
    return _W[n]


def exact_boot_ci(X, level=LEVEL, chunk=2000):
    """X: (nrep, n). Intervallo percentile del bootstrap della media con B infinito (enumerazione esatta),
    stessa interpolazione lineare di numpy.percentile usata in R4_common.boot_one_stage."""
    X = np.atleast_2d(np.asarray(X, float)); W = boot_weights(X.shape[1])
    a = 100 * (1 - level) / 2
    lo = np.empty(len(X)); hi = np.empty(len(X))
    for s in range(0, len(X), chunk):
        M = X[s:s + chunk] @ W.T
        lo[s:s + chunk] = np.percentile(M, a, axis=1); hi[s:s + chunk] = np.percentile(M, 100 - a, axis=1)
    return lo, hi


# ---------- simulazione dei valori per area ----------
def simulate_values(R, nrep, rng, marginal="normal", cv=None, df=3):
    """Valori per area con correlazione gaussiana R e marginale data. -> (X, valore vero del tessuto).
    normal: N(0,1); lognormal: exp(s z) con CV dato, vero = exp(s^2/2); t: t con df gradi (copula gaussiana), vero 0."""
    R = np.asarray(R, float)
    L = np.linalg.cholesky(R)
    Z = rng.standard_normal((nrep, len(R))) @ L.T
    if marginal == "normal":
        return Z, 0.0
    if marginal == "lognormal":
        s = np.sqrt(np.log(1 + cv ** 2))
        return np.exp(s * Z), float(np.exp(s ** 2 / 2))
    if marginal == "t":
        return stats.t.ppf(stats.norm.cdf(Z), df), 0.0
    raise ValueError(marginal)


def needed_kappa(X, truth, lo, hi):
    """Per ogni replica, il fattore minimo che, applicato alle due semi-ampiezze attorno alla stima, copre il vero."""
    est = X.mean(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        left = np.where(truth < est, (est - truth) / (est - lo), 0.0)
        right = np.where(truth > est, (truth - est) / (hi - est), 0.0)
    return np.nan_to_num(np.maximum(left, right), nan=np.inf, posinf=np.inf)


def calibrate_kappa(R, rng, nrep=NREP, level=LEVEL):
    """C-a: kappa tale che la copertura sia esattamente `level` con valori normali di correlazione R."""
    X, truth = simulate_values(R, nrep, rng)
    lo, hi = exact_boot_ci(X)
    return float(np.quantile(needed_kappa(X, truth, lo, hi), level))


def coverage_widened(X, truth, kappa, lo=None, hi=None):
    """Copertura dell'intervallo percentile allargato di kappa attorno alla stima (C-a)."""
    if lo is None: lo, hi = exact_boot_ci(X)
    est = X.mean(1)
    lo2 = est - kappa * (est - lo); hi2 = est + kappa * (hi - est)
    return float(np.mean((lo2 <= truth) & (hi2 >= truth)))


def ci_t_design(X, R, level=LEVEL):
    """C-b: xbar +/- t(n-1) s/sqrt(n) k(R)."""
    n = X.shape[1]; est = X.mean(1); s = X.std(1, ddof=1)
    q = stats.t.ppf(1 - (1 - level) / 2, n - 1)
    h = q * s / np.sqrt(n) * k_factor(R)
    return est - h, est + h


def coverage_interval(lo, hi, truth):
    return float(np.mean((lo <= truth) & (hi >= truth)))


def kappa_range(R, rng, nrep=NREP):
    """Ampiezza attesa max-min di n aree indipendenti diviso quella con correlazione R (valori normali)."""
    n = len(R)
    Xi, _ = simulate_values(np.eye(n), nrep, rng); Xr, _ = simulate_values(R, nrep, rng)
    return float(np.mean(Xi.max(1) - Xi.min(1)) / np.mean(Xr.max(1) - Xr.min(1)))


def mc_se(p, nrep):
    return float(np.sqrt(p * (1 - p) / nrep))


# ---------- pezzi disgiunti ----------
def assign_pieces(x, y, boxes):
    """x, y: coordinate globali (pixel). boxes (n, 4). -> indice dell'area che contiene il punto con il centro
    più vicino (a parità, indice più basso); -1 fuori da tutte. Gli estremi destri sono esclusi."""
    x = np.asarray(x, float); y = np.asarray(y, float); b = np.asarray(boxes, float)
    inside = (x[:, None] >= b[None, :, 0]) & (x[:, None] < b[None, :, 2]) & (y[:, None] >= b[None, :, 1]) & (y[:, None] < b[None, :, 3])
    cx = (b[:, 0] + b[:, 2]) / 2; cy = (b[:, 1] + b[:, 3]) / 2
    d2 = (x[:, None] - cx[None, :]) ** 2 + (y[:, None] - cy[None, :]) ** 2
    d2 = np.where(inside, d2, np.inf)
    out = np.argmin(d2, axis=1)            # argmin restituisce il primo indice in caso di parità
    out[~inside.any(1)] = -1
    return out


def widen(est, lo, hi, kappa):
    """Allarga di kappa le due semi-ampiezze attorno alla stima."""
    return est - kappa * (est - lo), est + kappa * (hi - est)
