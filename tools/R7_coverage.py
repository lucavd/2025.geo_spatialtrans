"""tools/R7_coverage.py — R7: simulazione di copertura (M-3), fattori per tessuto (M-a), previsione P-1, regola R-1.
Uso: .venv/bin/python tools/R7_coverage.py  -> results/R7/coverage/*.csv"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R7_common as C
from R4_common import verdict
OUT = "results/R7/coverage"; os.makedirs(OUT, exist_ok=True)
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
R4roi = pd.read_csv("results/R4/R4_roi.csv")

# CV alto per lo scenario lognormale: 90° percentile dei CV fra aree delle quantità di R4_roi.csv (tessuto x quantità)
skip = {"archetype", "roi_id", "n_bins2", "n_bins8", "n_nuclei", "G3c_all_identical"}
qcols = [c for c in R4roi.columns if c not in skip and not c.startswith("spec_") and R4roi[c].dtype.kind == "f"]
cvs = []
for a in ARCH:
    for q in qcols:
        v = R4roi.loc[R4roi.archetype == a, q].to_numpy(float); v = v[np.isfinite(v)]
        if len(v) >= 4 and abs(v.mean()) > 0: cvs.append(dict(archetype=a, quantity=q, cv=float(v.std(ddof=1) / abs(v.mean()))))
cvs = pd.DataFrame(cvs); cvs.to_csv(f"{OUT}/R7_between_area_cv.csv", index=False)
CV_HI = float(np.quantile(cvs.cv, 0.90))
print(f"CV alto (90° percentile su {len(cvs)} quantità x tessuto): {CV_HI:.4f}", flush=True)

cases = []
for a in ARCH:
    g = rois[rois.archetype == a]
    cases.append(dict(case=a, unit="area", n=5, R=C.box_overlap_matrix(g[["c0", "r0", "c1", "r1"]].to_numpy(float))))
cases.append(dict(case="independent_5", unit="area", n=5, R=np.eye(5)))
cases.append(dict(case="windows_4", unit="window", n=4, R=np.eye(4)))

K, P1 = [], []
SCEN = [("normal", {}), ("lognormal", {"cv": CV_HI}), ("t", {"df": 3})]
for ci, cs in enumerate(cases):
    R, n = cs["R"], cs["n"]
    rng = lambda p: np.random.default_rng([C.SEED, ci, p])
    kappa = C.calibrate_kappa(R, rng(1))
    kr = C.kappa_range(R, rng(2))
    X, tr = C.simulate_values(R, C.NREP, rng(3)); lo, hi = C.exact_boot_ci(X)
    cov_now = C.coverage_interval(lo, hi, tr); cov_mm = C.coverage_interval(X.min(1), X.max(1), tr)
    K.append(dict(case=cs["case"], unit=cs["unit"], n=n, rho_bar=C.rho_bar(R), k=C.k_factor(R), kappa=kappa, kappa_range=kr,
                  cover_now=cov_now, cover_now_se=C.mc_se(cov_now, C.NREP), cover_minmax=cov_mm, cover_minmax_se=C.mc_se(cov_mm, C.NREP)))
    for si, (sc, kw) in enumerate(SCEN):
        X, tr = C.simulate_values(R, C.NREP, rng(10 + si), marginal=sc, **kw); lo, hi = C.exact_boot_ci(X)
        for meth, cv_ in (("C-a", C.coverage_widened(X, tr, kappa, lo, hi)), ("C-b", C.coverage_interval(*C.ci_t_design(X, R), tr)), ("now", C.coverage_interval(lo, hi, tr))):
            se = C.mc_se(cv_, C.NREP)
            P1.append(dict(case=cs["case"], scenario=sc, method=meth, coverage=cv_, se=se, lo=cv_ - 2 * se, hi=cv_ + 2 * se,
                           verdict=verdict(cv_ - 2 * se, cv_ + 2 * se, *C.BRADLEY)))
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in K[-1].items()}), flush=True)

K = pd.DataFrame(K); P1 = pd.DataFrame(P1)
K.to_csv(f"{OUT}/R7_kappa.csv", index=False); P1.to_csv(f"{OUT}/R7_P1_cells.csv", index=False)
# R-1: celle tessuto x scenario dentro Bradley (sei tessuti, tre scenari)
tis = P1[P1.case.isin(ARCH)]
cnt = {m: int((tis[tis.method == m].verdict == "confermata").sum()) for m in ("C-a", "C-b")}
choice = "C-a" if cnt["C-a"] >= cnt["C-b"] else "C-b"
json.dump(dict(rule="R-1", cells_within_bradley=cnt, n_cells=int(len(tis) // 3), choice=choice, cv_hi=CV_HI), open(f"{OUT}/R7_R1.json", "w"), indent=1)
print("R-1:", cnt, "->", choice)
