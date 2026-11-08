"""tools/R7_kappa_types.py — R7: kappa con le matrici specifiche (regola R-2), precisione Monte Carlo di kappa,
analisi post hoc PH-1 (copertura dell'allargamento calibrato in funzione dell'asimmetria dei valori per area).
Uso: .venv/bin/python tools/R7_kappa_types.py -> results/R7/coverage/R7_kappa_types.csv, R7_kappa_mc.csv, R7_PH1_cv.csv"""
import os, sys, itertools
import numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R7_common as C
OUT = "results/R7/coverage"
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
F = pd.read_csv("results/R7/pieces/R7_F_matrices.csv")
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
def mat(A, col):
    ids = list(rois[rois.archetype == A].roi_id); M = np.eye(5)
    for _, r in F[F.archetype == A].iterrows():
        i, j = ids.index(r.i), ids.index(r.j); M[i, j] = M[j, i] = r[col]
    return M
rows, mc = [], []
for ai, A in enumerate(ARCH):
    for ti, col in enumerate(["F_box", "F_bin", "F_nuc"]):
        R = mat(A, col)
        ks = [C.calibrate_kappa(R, np.random.default_rng([C.SEED, 500 + ai, ti, s])) for s in range(5)]
        rows.append(dict(archetype=A, matrix=col, rho_bar=C.rho_bar(R), kappa_mean5=float(np.mean(ks)), kappa_sd5=float(np.std(ks, ddof=1)),
                         kappa_range=C.kappa_range(R, np.random.default_rng([C.SEED, 600 + ai, ti]))))
        print(rows[-1], flush=True)
for nm, R in (("independent_5", np.eye(5)), ("windows_4", np.eye(4))):
    ks = [C.calibrate_kappa(R, np.random.default_rng([C.SEED, 700 + len(R), s])) for s in range(5)]
    mc.append(dict(case=nm, kappa_mean5=float(np.mean(ks)), kappa_sd5=float(np.std(ks, ddof=1))))
pd.DataFrame(rows).to_csv(f"{OUT}/R7_kappa_types.csv", index=False); pd.DataFrame(mc).to_csv(f"{OUT}/R7_kappa_mc.csv", index=False)
# PH-1: copertura di C-a (kappa calibrato su valori normali) con valori lognormali di CV crescente; kappa necessario per lognormale
ph = []
for nm, R in (("independent_5", np.eye(5)), ("A3", mat("A3", "F_box"))):
    k_norm = C.calibrate_kappa(R, np.random.default_rng([C.SEED, 800, len(nm)]))
    for ci, cv in enumerate([0.05, 0.1, 0.2, 0.3, 0.45, 0.62, 1.0]):
        rng = np.random.default_rng([C.SEED, 801, len(nm), ci])
        X, tr = C.simulate_values(R, C.NREP, rng, marginal="lognormal", cv=cv); lo, hi = C.exact_boot_ci(X)
        cov = C.coverage_widened(X, tr, k_norm, lo, hi)
        k_ln = float(np.quantile(C.needed_kappa(X, tr, lo, hi), 0.95))
        ph.append(dict(case=nm, cv=cv, kappa_normal=k_norm, coverage_Ca=cov, se=C.mc_se(cov, C.NREP), kappa_needed_lognormal=k_ln))
        print(ph[-1], flush=True)
pd.DataFrame(ph).to_csv(f"{OUT}/R7_PH1_cv.csv", index=False)
