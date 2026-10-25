"""tools/S1.5d_calibration_types.py - S1.5d: ricalibrazione dei compartimenti per tipo (BL-159; definizione in docs/steps/S1.5d.md).
Modello di S1.5c (campo gaussiano latente per tipo principale, caratteristiche di Fourier casuali M = 256, p_t = softmax(log f_t + a_t + s_t Z_t),
a_t tali che la media di p_t sui punti = f_t; etichette estratte SULLE POSIZIONI DEI NUCLEI VERI con un tipo). Sei correzioni pre-registrate:
(1) campi indipendenti per ROI (seme [SEED, indice della ROI, replica, tipo, l]); (2) griglia l in {4, 8, ..., 1024} um, s in {0, 0.5, ..., 4, 5, 6, 7, 8};
(3) f = composizione del preset (results/S1.5d/S1.5d_preset_types.csv); (4) bersaglio = metrica di P-6: rapporto di connessione con focali = celle
interne del tipo (al piu' 1 500, campione fisso come tools/S1.5c_metrics.py) e frazione del tipo fra i centri con un tipo della ROI, bin 0-10-20-40-80-160;
(5) perdita sommata per ROI: sum_ROI sum_(tipo, bin) w (log R_sim - log R_vero)^2, w = min(coppie vere dello stesso tipo nel bin della ROI, 1000),
R_sim = media delle repliche dei conteggi; (6) discesa per coordinate in 3 giri: il primo su tutta la griglia con 5 repliche, il secondo e il terzo sui
vicini (un passo in l e in s, compreso il punto attuale) con 20 repliche. Pieghe: lascia-fuori-un-ROI e tutte le ROI.
Uscite: results/S1.5d/S1.5d_compartments_<A>.csv (parametri, perdita, bordo della griglia), S1.5d_compartments_curves_<A>.csv (piega «pooled»,
per ROI, tipo, bin e replica: R simulato con 20 repliche; R vero). Uso: S15D_WORKERS=6 .venv/bin/python tools/S1.5d_calibration_types.py A1 [A2 ...] (una piega per processo)"""
import sys, os, importlib.util, numpy as np, pandas as pd
from scipy.spatial import cKDTree
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py")
RES = "results/S1.5d"; os.makedirs(RES, exist_ok=True)
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
BINS = np.array([0, 10, 20, 40, 80, 160.]); NB = len(BINS) - 1
ELL = [4, 8, 16, 32, 64, 128, 256, 512, 1024]; SS = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8]
M = 256; NREP1 = 5; NREP2 = 20; NFOC = 1500; SEED = 20261025; FOC_SEED = 20261023

def rff(P, ell, seed):
    rng = np.random.default_rng(seed); w = rng.normal(0, 1 / ell, size=(M, 2)); b = rng.uniform(0, 2 * np.pi, M)
    return np.sqrt(2 / M) * np.cos(P @ w.T + b).sum(1)

class ROI:
    def __init__(self, k, P, L, interior, K, main_idx):
        self.k, self.P, self.K, self.main = k, P, K, main_idx; n = len(P); self.I = interior
        rng = np.random.default_rng(FOC_SEED); self.fi, self.j, self.b = [], [], []
        tr = cKDTree(P); self.real = None; self.L = L
        # focali per tipo: celle interne del tipo (al piu' 1500), come la metrica di P-6; le coppie dipendono dall'etichetta del focale,
        # quindi si memorizzano le coppie di TUTTE le celle interne che potrebbero essere focali (campione di 1500 per tipo fatto a ogni valutazione)
        self.cand = np.where(interior)[0]
        lists = tr.query_ball_point(P[self.cand], BINS[-1])
        fi = np.repeat(self.cand, [len(l) for l in lists]); j = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists])
        keep = fi != j; fi, j = fi[keep], j[keep]; d = np.hypot(*(P[j] - P[fi]).T)
        self.pf, self.pj, self.pb = fi.astype(np.int32), j.astype(np.int32), np.minimum(np.digitize(d, BINS) - 1, NB - 1).astype(np.int8)
        del fi, j, d
        self.real = self.ratio(L)[0]; self.wreal = np.minimum(self.ratio(L)[1], 1000)
        self.Z = {}
    def focal_mask(self, L):
        """Per ogni tipo principale: al piu' NFOC focali fra le celle interne del tipo (campione fisso con FOC_SEED, come tools/S1.5c_metrics.py)."""
        m = np.zeros(len(self.P), bool); rng = np.random.default_rng(FOC_SEED)
        for t in self.main:
            foc = np.where(self.I & (L == t))[0]
            if len(foc) < 20: continue
            sub = foc if len(foc) <= NFOC else np.sort(rng.choice(foc, NFOC, replace=False)); m[sub] = True
        return m
    def ratio(self, L):
        fm = self.focal_mask(L); k = fm[self.pf]; fi, j, b = self.pf[k], self.pj[k], self.pb[k]
        idx = L[fi].astype(np.int64) * NB + b; same = (L[fi] == L[j]).astype(float)
        den = np.bincount(idx, minlength=self.K * NB).reshape(self.K, NB); num = np.bincount(idx, weights=same, minlength=self.K * NB).reshape(self.K, NB)
        f = np.bincount(L, minlength=self.K) / len(L)
        with np.errstate(divide="ignore", invalid="ignore"): R = num / den / f[:, None]
        return R, num
    def field(self, rep, t, ell):
        key = (rep, t, ell)
        if key not in self.Z: self.Z[key] = rff(self.P, ell, [SEED, self.k, rep, t, ell])
        return self.Z[key]

def draw(roi, f, ell, s, rep):
    eta = np.tile(np.log(f), (len(roi.P), 1))
    for t in roi.main:
        if s[t] > 0: eta[:, t] += s[t] * roi.field(rep, t, ell[t])
    a = np.zeros(roi.K)
    for _ in range(60):
        e = eta + a; e -= e.max(1, keepdims=True); p = np.exp(e); p /= p.sum(1, keepdims=True)
        a += np.log(f) - np.log(p.mean(0))
    u = np.random.default_rng([SEED + 1, roi.k, rep]).random(len(roi.P)); return (p.cumsum(1) < u[:, None]).sum(1).clip(0, roi.K - 1)

def loss(rois, f, ell, s, nrep, keep_reps=False):
    tot = 0.0; per = []
    for r in rois:
        Rs = np.array([r.ratio(draw(r, f, ell, s, rep))[0] for rep in range(nrep)])
        with np.errstate(divide="ignore", invalid="ignore"):
            Rm = np.nanmean(Rs, axis=0)
            e = np.where((r.real > 0) & (Rm > 0), r.wreal * (np.log(Rm) - np.log(r.real)) ** 2, r.wreal * 4.0)
        e = np.nan_to_num(e[r.main], nan=0.0); tot += float(e.sum())
        if keep_reps: per.append(Rs)
    return (tot, per) if keep_reps else tot

def neighbours(v, grid):
    i = grid.index(v); return [grid[j] for j in (i - 1, i, i + 1) if 0 <= j < len(grid)]

PT = pd.read_csv("results/S1.5d/S1.5d_preset_types.csv")
from multiprocessing import Pool
NW = 6; ROIS_ = None; CTX_ = None
def main(A):
    _, d = cal.load(A); d = d[~d.type.isin(list(cal.UNTYPED))]
    pt = PT[PT.archetype == A]; types = sorted(pt.cell_type); main_idx = [types.index(t) for t in MAIN[A]]
    f = pt.set_index("cell_type").fraction.reindex(types).values; f = f / f.sum()
    rlist = sorted(d.roi_id.unique()); rois = {}
    for k, r in enumerate(rlist):
        e = d[d.roi_id == r]; e = e[e.type.isin(types)]
        rois[r] = ROI(k, e[["x", "y"]].values, np.array([types.index(t) for t in e.type]), e.interior.values.astype(bool), len(types), main_idx)
    RR = pd.read_csv("results/S1.5c/S1.5c_real_refs.csv"); RR = RR[(RR.archetype == A) & (RR.metric == "markconn")]
    worst = 0.0
    for r in rlist:                                                    # controllo: R vero = riferimento di S1.5c (stessa metrica di P-6)
        for t in main_idx:
            for b in range(NB):
                v = RR[(RR.roi_id == r) & (RR.cell_type == types[t]) & (RR.bin == f"{int(BINS[b])}-{int(BINS[b+1])}")].value
                if len(v) and np.isfinite(v.iloc[0]): worst = max(worst, abs(v.iloc[0] - rois[r].real[t, b]))
    if worst > 1e-9: raise SystemExit(f"{A}: rapporto di connessione vero diverso dal riferimento di S1.5c ({worst:.3g})")
    print(A, "R vero = riferimento di S1.5c, scarto massimo", worst, flush=True)
    folds = [("loo_" + r, [q for q in rlist if q != r]) for r in rlist] + [("pooled", rlist)]
    global ROIS_, CTX_; ROIS_ = rois; CTX_ = (A, f, types, main_idx, rlist)
    with Pool(min(NW, len(folds))) as pool: outs = pool.map(run_fold, folds)
    prow = [x for o in outs for x in o[0]]; crow = [x for o in outs for x in o[1]]
    pd.DataFrame(prow).to_csv(f"{RES}/S1.5d_compartments_{A}.csv", index=False)
    pd.DataFrame(crow).to_csv(f"{RES}/S1.5d_compartments_curves_{A}.csv", index=False)

def run_fold(arg):
    fold, keys = arg; A, f, types, main_idx, rlist = CTX_; rois = ROIS_; prow, crow = [], []
    if True:
        rs = [rois[q] for q in keys]; ell = [ELL[0]] * len(types); s = [0.0] * len(types)
        L0 = loss(rs, f, ell, s, NREP1); best = L0
        for t in main_idx:                                             # giro 1: tutta la griglia, 5 repliche
            for l in ELL:
                for sv in SS:
                    e2, s2 = list(ell), list(s); e2[t] = l; s2[t] = sv; L_ = loss(rs, f, e2, s2, NREP1)
                    if L_ < best - 1e-12: best, ell, s = L_, e2, s2
        best = loss(rs, f, ell, s, NREP2)
        for it in range(2):                                            # giri 2 e 3: vicini, 20 repliche
            for t in main_idx:
                for l in neighbours(ell[t], ELL):
                    for sv in neighbours(s[t], SS):
                        e2, s2 = list(ell), list(s); e2[t] = l; s2[t] = sv; L_ = loss(rs, f, e2, s2, NREP2)
                        if L_ < best - 1e-12: best, ell, s = L_, e2, s2
        Lnull = loss(rs, f, [ELL[0]] * len(types), [0.0] * len(types), NREP2)
        for t in main_idx:
            prow.append(dict(archetype=A, fold=fold, cell_type=types[t], ell_um=ell[t], s=s[t], f_t=f[t], loss=best, loss_null=Lnull, n_roi=len(keys),
                             edge_ell=ell[t] in (ELL[0], ELL[-1]) and s[t] > 0, edge_s=s[t] == SS[-1]))
        if fold == "pooled":
            _, per = loss(rs, f, ell, s, NREP2, keep_reps=True)
            for r, Rs in zip(rs, per):
                for t in main_idx:
                    for b in range(NB):
                        base = dict(archetype=A, roi_id=rlist[r.k], cell_type=types[t], bin=f"{int(BINS[b])}-{int(BINS[b+1])}", real=r.real[t, b], n_pairs_real=int(r.ratio(r.L)[1][t, b]))
                        crow += [dict(base, rep=q, sim=Rs[q, t, b]) for q in range(NREP2)]
        print(A, fold, "loss", round(best, 2), "null", round(Lnull, 2), [(types[t], ell[t], s[t]) for t in main_idx], flush=True)
    return prow, crow

if __name__ == "__main__":
    NW = int(os.environ.get("S15D_WORKERS", "6"))
    for A in sys.argv[1:]: main(A)
