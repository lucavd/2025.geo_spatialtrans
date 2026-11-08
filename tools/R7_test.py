"""tools/R7_test.py — R7, cancello G-1: batteria di test con PASS/FAIL per asserzione e mutanti.
Uso: .venv/bin/python tools/R7_test.py   (uscita 0 se tutto passa e tutti i mutanti sono rilevati)"""
import sys, numpy as np
sys.path.insert(0, "tools")
import R7_common as C
import R4_common as R4

RESULTS = []
def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok))); print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True); return bool(ok)

# matrice tipo linfonodo (quote vere arrotondate) per i test
F_A4 = np.array([[1, .002, 0, 0, 0], [.002, 1, .131, .384, .082], [0, .131, 1, 0, 0], [0, .384, 0, 1, .278], [0, .082, 0, .278, 1]])

def t_overlap():
    b = np.array([[0, 0, 10, 10], [5, 0, 15, 10], [5, 5, 15, 15], [20, 20, 30, 30], [0, 0, 10, 10]])
    F = C.box_overlap_matrix(b)
    ok = np.isclose(F[0, 1], .5) and np.isclose(F[0, 2], .25) and F[0, 3] == 0 and np.isclose(F[0, 4], 1) and np.allclose(F, F.T)
    return check("T1 quote di zona comune su quadrati noti", ok, f"{F[0,1]:.3f} {F[0,2]:.3f} {F[0,3]:.3f} {F[0,4]:.3f}")

def t_units():
    M = np.zeros((100, 2), bool); M[:60, 0] = True; M[40:, 1] = True     # 60 e 60 unità, 20 in comune
    U = C.unit_overlap_matrix(M)
    return check("T2 quota di unità condivise", np.isclose(U[0, 1], 20 / 60), f"{U[0,1]:.4f}")

def t_k():
    rng = np.random.default_rng(1); n = 5
    X, _ = C.simulate_values(F_A4, 400000, rng)
    brute = np.sqrt(X.mean(1).var() / np.mean(X.var(1, ddof=1) / n))
    k = C.k_factor(F_A4)
    return check("T3 fattore k contro simulazione bruta", abs(brute / k - 1) < 0.01 and k > 1.05, f"k {k:.4f} bruta {brute:.4f}")

def t_t_interval():
    rng = np.random.default_rng(2)
    X, tr = C.simulate_values(np.eye(5), C.NREP, rng)
    lo, hi = C.ci_t_design(X, np.eye(5)); cov = C.coverage_interval(lo, hi, tr); se = C.mc_se(.95, C.NREP)
    return check("T4 intervallo t con aree normali indipendenti copre 0.95", abs(cov - .95) < 3 * se, f"{cov:.4f} (3 SE {3*se:.4f})")

def t_exact_boot():
    v = np.array([3.1, 4.7, 2.2, 5.9, 4.0])
    lo, hi = C.exact_boot_ci(v[None, :])
    _, blo, bhi, _ = R4.boot_one_stage(v, B=200000, rng=np.random.default_rng(3))
    tol = 0.01 * (v.max() - v.min())
    return check("T5 bootstrap esatto uguale a boot_one_stage di R4 (B 200 000)", abs(lo[0] - blo) < tol and abs(hi[0] - bhi) < tol,
                 f"esatto [{lo[0]:.4f}, {hi[0]:.4f}] R4 [{blo:.4f}, {bhi:.4f}]")

def t_kappa():
    rng = np.random.default_rng(4)
    k = C.calibrate_kappa(F_A4, rng)
    X, tr = C.simulate_values(F_A4, C.NREP, np.random.default_rng(5))
    cov = C.coverage_widened(X, tr, k); se = C.mc_se(.95, C.NREP)
    return check("T6 kappa calibrato copre 0.95 su repliche nuove", abs(cov - .95) < 3 * se and k > 1, f"kappa {k:.4f} copertura {cov:.4f}")

def t_pieces():
    rng = np.random.default_rng(6)
    b = np.array([[0, 0, 100, 100], [60, 0, 160, 100], [40, 50, 140, 150], [300, 300, 400, 400]], float)
    x = rng.uniform(-10, 410, 50000); y = rng.uniform(-10, 410, 50000)
    p = C.assign_pieces(x, y, b)
    inside = (x[:, None] >= b[None, :, 0]) & (x[:, None] < b[None, :, 2]) & (y[:, None] >= b[None, :, 1]) & (y[:, None] < b[None, :, 3])
    in_union = inside.any(1)
    ok1 = np.all((p >= 0) == in_union)                                      # unione dei pezzi = unione delle aree
    ok2 = np.all(inside[np.arange(len(p))[p >= 0], p[p >= 0]])              # ogni pezzo dentro la sua area
    cx = (b[:, 0] + b[:, 2]) / 2; cy = (b[:, 1] + b[:, 3]) / 2
    d2 = np.where(inside, (x[:, None] - cx) ** 2 + (y[:, None] - cy) ** 2, np.inf)
    ok3 = np.all(d2[np.arange(len(p))[p >= 0], p[p >= 0]] <= d2[p >= 0].min(1))   # centro più vicino
    ok4 = np.all(p[inside[:, 3]] == 3)                                      # area disgiunta: pezzo = area
    return check("T7 pezzi disgiunti (partizione, contenimento, centro più vicino, area isolata intatta)", ok1 and ok2 and ok3 and ok4, f"{ok1} {ok2} {ok3} {ok4}")

def t_verdict():
    cases = [(1, 2, 0, 3), (-1, 2, 0, 3), (4, 5, 0, 3), (np.nan, 1, 0, 1), (0, 1, -np.inf, 2), (2.5, 3.5, 0, 3)]
    from R4_common import verdict
    ok = all(verdict(*cs) == R4.verdict(*cs) for cs in cases) and [R4.verdict(*cs) for cs in cases][:3] == [R4.verdict(1, 2, 0, 3), R4.verdict(-1, 2, 0, 3), R4.verdict(4, 5, 0, 3)]
    outs = [R4.verdict(*cs) for cs in cases]
    return check("T8 funzione dei verdetti di R4 (tre esiti e non valutabile)", ok and len(set(outs[:3])) == 3 and outs[3] == "non valutabile", str(outs))

def t_krange():
    rng = np.random.default_rng(7)
    ki = C.kappa_range(np.eye(5), rng); kf = C.kappa_range(F_A4, np.random.default_rng(8))
    return check("T9 kappa_range: 1 senza zone comuni, > 1 con", abs(ki - 1) < 0.01 and kf > 1.01, f"{ki:.4f} {kf:.4f}")

def t_widen():
    lo, hi = C.widen(10.0, 8.0, 13.0, 1.5)
    return check("T10 allargamento attorno alla stima, asimmetria conservata", np.isclose(lo, 7.0) and np.isclose(hi, 14.5), f"{lo} {hi}")

TESTS = [t_overlap, t_units, t_k, t_t_interval, t_exact_boot, t_kappa, t_pieces, t_verdict, t_krange, t_widen]
print("== batteria"); base = [t() for t in TESTS]

# ---------- mutanti ----------
import scipy.stats as st
MUT = []
def mutant(name, patch, tests):
    saved = {k: getattr(C, k) for k in patch}
    for k, f in patch.items(): setattr(C, k, f)
    n0 = len(RESULTS)
    try:
        res = [t() for t in tests]
    except Exception as e:
        res = [False]; print(f"      (eccezione: {e})")
    for k, f in saved.items(): setattr(C, k, f)
    del RESULTS[n0:]
    caught = not all(res); MUT.append((name, caught))
    print(f"{'PASS' if caught else 'FAIL'}  mutante rilevato: {name}", flush=True)

orig_sim, orig_widen = C.simulate_values, C.widen
def m_ci_t(X, R, level=C.LEVEL):
    n = X.shape[1]; est = X.mean(1); s = X.std(1, ddof=1); q = st.t.ppf(1 - (1 - level) / 2, n)
    h = q * s / np.sqrt(n) * C.k_factor(R); return est - h, est + h
def m_pieces(x, y, boxes):
    b = np.asarray(boxes, float); x = np.asarray(x, float); y = np.asarray(y, float)
    inside = (x[:, None] >= b[None, :, 0]) & (x[:, None] < b[None, :, 2]) & (y[:, None] >= b[None, :, 1]) & (y[:, None] < b[None, :, 3])
    cx = (b[:, 0] + b[:, 2]) / 2; cy = (b[:, 1] + b[:, 3]) / 2
    d2 = np.where(inside, (x[:, None] - cx) ** 2 + (y[:, None] - cy) ** 2, -np.inf)
    out = np.argmax(d2, axis=1); out[~inside.any(1)] = -1; return out
def m_cov_one_side(X, truth, kappa, lo=None, hi=None):
    if lo is None: lo, hi = C.exact_boot_ci(X)
    est = X.mean(1); hi2 = est + kappa * (hi - est); return float(np.mean((lo <= truth) & (hi2 >= truth)))
_Wsave = C.boot_weights
def m_weights(n):
    import itertools
    idx = np.array(list(itertools.permutations(range(n))))
    return np.stack([np.bincount(r, minlength=n) for r in idx]).astype(float) / n
def m_krange(R, rng, nrep=C.NREP):
    n = len(R); Xi, _ = orig_sim(np.eye(n), nrep, rng); Xr, _ = orig_sim(R, nrep, rng)
    return float(np.mean(Xr.max(1) - Xr.min(1)) / np.mean(Xi.max(1) - Xi.min(1)))
def m_calib_noF(R, rng, nrep=C.NREP, level=C.LEVEL):
    X, truth = orig_sim(np.eye(len(R)), nrep, rng); lo, hi = C.exact_boot_ci(X)
    return float(np.quantile(C.needed_kappa(X, truth, lo, hi), level))

print("== mutanti")
mutant("M1 correlazione ignorata nel fattore k", {"k_factor": lambda R: 1.0}, [t_k])
mutant("M2 gradi di libertà n invece di n-1", {"ci_t_design": m_ci_t}, [t_t_interval])
mutant("M3 pezzo al centro più lontano", {"assign_pieces": m_pieces}, [t_pieces])
mutant("M4 kappa su un solo lato", {"coverage_widened": m_cov_one_side}, [t_kappa])
mutant("M5 bootstrap senza reinserimento", {"boot_weights": m_weights}, [t_exact_boot])
mutant("M6 kappa_range invertito", {"kappa_range": m_krange}, [t_krange])
mutant("M7 kappa calibrato senza le zone comuni", {"calibrate_kappa": m_calib_noF}, [t_kappa])

nb = sum(ok for _, ok in RESULTS); nm = sum(ok for _, ok in MUT)
print(f"== batteria {nb}/{len(RESULTS)}; mutanti rilevati {nm}/{len(MUT)}")
sys.exit(0 if nb == len(RESULTS) and nm == len(MUT) else 1)
