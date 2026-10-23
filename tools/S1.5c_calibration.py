"""tools/S1.5c_calibration.py - S1.5c: calibrazione del modello a compartimenti per tipo (definizione in docs/steps/S1.5c.md).
Per ogni tipo principale t: campo gaussiano latente Z_t (covarianza exp(-r^2 / (2 l_t^2)), caratteristiche di Fourier casuali, M = 256),
probabilita' locali p_t(x) = softmax_t(log f_t + a_t + s_t Z_t(x)) (tipi minori s = 0), a_t tali che la media di p_t sui punti = f_t.
Bersaglio: rapporto di connessione dei tipi R_t(b) = P(vicino di tipo t | focale di tipo t, distanza in b) / f_t, bin 0-10-20-40-80-160 um,
dei nuclei veri con un tipo (etichette R6 PH-1). Le etichette simulate si estraggono SULLE POSIZIONI DEI NUCLEI VERI (approssimazione dichiarata:
la calibrazione non include l'accoppiamento fra tipo e densita' locale della semina). Perdita: somma su tipi principali e bin di
w (log R_sim - log R_real)^2, w = min(coppie reali dello stesso tipo nel bin, 1000). Griglia l in {8,16,...,512} um, s in {0,0.5,...,4,5,6};
discesa per coordinate (2 giri, tipi in ordine fisso), 3 repliche dei campi per valutazione (stessi semi per ogni punto della griglia).
Pieghe: lascia-fuori-un-ROI (5) + tutte le 5 ROI (preset). Uso: .venv/bin/python tools/S1.5c_calibration.py A1 [A2 ...]
Uscite: results/S1.5c/S1.5c_compartments_<A>.csv (parametri), results/S1.5c/S1.5c_compartments_curves_<A>.csv (curve reali e adattate)."""
import sys, os, importlib.util, numpy as np, pandas as pd
from scipy.spatial import cKDTree
spec = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(spec); spec.loader.exec_module(cal)
RES = "results/S1.5c"; os.makedirs(RES, exist_ok=True)
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
BINS = np.array([0, 10, 20, 40, 80, 160.]); NB = len(BINS) - 1
ELL = [8, 16, 32, 64, 128, 256, 512]; SS = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6]; M = 256; NREP = 3; NFOC = 2000; SEED = 20261023

def rff(P, ell, rng):
    w = rng.normal(0, 1 / ell, size=(M, 2)); b = rng.uniform(0, 2 * np.pi, M)
    return np.sqrt(2 / M) * np.cos(P @ w.T + b).sum(1)

class ROI:
    def __init__(self, P, lab, types, rng):
        self.P = P; self.types = types; self.K = len(types)
        n = len(P); foc = np.arange(n) if n <= NFOC else np.sort(rng.choice(n, NFOC, replace=False))
        tr = cKDTree(P); lists = tr.query_ball_point(P[foc], BINS[-1])
        fi = np.repeat(foc, [len(l) for l in lists]); j = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists])
        k = fi != j; fi, j = fi[k], j[k]; d = np.hypot(*(P[j] - P[fi]).T); self.fi, self.j, self.b = fi, j, np.digitize(d, BINS) - 1
        self.real = self.counts(lab)
    def counts(self, L):
        idx = L[self.fi] * NB + self.b; same = (L[self.fi] == L[self.j]).astype(float)
        den = np.bincount(idx, minlength=self.K * NB).reshape(self.K, NB); num = np.bincount(idx, weights=same, minlength=self.K * NB).reshape(self.K, NB)
        f = np.bincount(L, minlength=self.K) / len(L); return num, den, f

def labels(roi, f, ell, s, Zcache, rep, rng):
    eta = np.tile(np.log(f), (len(roi.P), 1))
    for t in range(roi.K):
        if s[t] > 0: eta[:, t] += s[t] * Zcache[(rep, t, ell[t])]
    a = np.zeros(roi.K)
    for _ in range(40):
        e = eta + a; e -= e.max(1, keepdims=True); p = np.exp(e); p /= p.sum(1, keepdims=True)
        a += np.log(f) - np.log(p.mean(0))
    u = rng.random(len(roi.P)); return (p.cumsum(1) < u[:, None]).sum(1).clip(0, roi.K - 1)

def ratio(num, den, f): return num / np.maximum(den, 1) / f[:, None]

def loss_of(rois, f, ell, s, main_idx):
    num_r = sum(r.real[0] for r in rois); den_r = sum(r.real[1] for r in rois)
    Rr = ratio(num_r, den_r, f); num_s = 0; den_s = 0
    for k, r in enumerate(rois):
        for rep in range(NREP):
            L = labels(r, f, ell, s, r.Z, rep, np.random.default_rng(SEED + 7919 * rep + 104729 * k))
            n_, d_, _ = r.counts(L); num_s = num_s + n_; den_s = den_s + d_
    Rs = ratio(num_s, den_s, f); w = np.minimum(num_r, 1000)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = np.where((num_r > 0) & (num_s > 0), w * (np.log(Rs) - np.log(Rr)) ** 2, w * 4.0)
    return float(e[main_idx].sum()), Rr, Rs

for A in sys.argv[1:]:
    lab, _ = cal.load(A)
    lab = lab[~lab.type.isin(list(cal.UNTYPED))]
    types = sorted(lab.type.unique()); main_idx = [types.index(t) for t in MAIN[A]]
    f = lab.type.value_counts(normalize=True).reindex(types).values
    rng = np.random.default_rng(SEED); rois = {}
    for r, e in lab.groupby("roi_id"):
        P = e[["x_um", "y_um"]].values; L = np.array([types.index(t) for t in e.type])
        ro = ROI(P, L, types, rng); ro.Z = {}
        for rep in range(NREP):
            for t in main_idx:
                for l in ELL: ro.Z[(rep, t, l)] = rff(P, l, np.random.default_rng([SEED, rep, t, l]))
        rois[r] = ro
    folds = [("loo_" + r, [k for k in rois if k != r]) for r in sorted(rois)] + [("pooled", sorted(rois))]
    prow, crow = [], []
    for fold, keys in folds:
        rs = [rois[k] for k in keys]; ell = [ELL[0]] * len(types); s = [0.0] * len(types)
        best, _, _ = loss_of(rs, f, ell, s, main_idx); trace = [best]
        for it in range(2):
            for t in main_idx:
                for l in ELL:
                    for sv in SS:
                        e2, s2 = list(ell), list(s); e2[t] = l; s2[t] = sv
                        L_, _, _ = loss_of(rs, f, e2, s2, main_idx)
                        if L_ < best - 1e-12: best, ell, s = L_, e2, s2
                trace.append(best)
        L_, Rr, Rs = loss_of(rs, f, ell, s, main_idx)
        L0, _, R0 = loss_of(rs, f, [ELL[0]] * len(types), [0.0] * len(types), main_idx)
        for t in main_idx:
            prow.append(dict(archetype=A, fold=fold, cell_type=types[t], ell_um=ell[t], s=s[t], f_t=f[t], loss=L_, loss_null=L0, n_roi=len(keys)))
            for b in range(NB):
                crow.append(dict(archetype=A, fold=fold, cell_type=types[t], bin=f"{int(BINS[b])}-{int(BINS[b+1])}", real=Rr[t, b], fitted=Rs[t, b], null=R0[t, b]))
        print(A, fold, "loss", round(L_, 2), "null", round(L0, 2), [(types[t], ell[t], s[t]) for t in main_idx], flush=True)
    pd.DataFrame(prow).to_csv(f"{RES}/S1.5c_compartments_{A}.csv", index=False)
    pd.DataFrame(crow).to_csv(f"{RES}/S1.5c_compartments_curves_{A}.csv", index=False)
