"""tools/R7_p1_v3.py — correzione C-1: P-1 con i kappa di R7_kappa_v3_all.csv (F_box, 10^6 repliche) invece della calibrazione a un seme.
-> results/R7/coverage/R7_P1_cells_v3.csv"""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R7_common as C
from R4_common import verdict
K3 = pd.read_csv("results/R7/coverage/R7_kappa_v3_all.csv"); CV_HI = json.load(open("results/R7/coverage/R7_R1.json"))["cv_hi"]
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
cases = [(a, C.box_overlap_matrix(rois[rois.archetype == a][["c0", "r0", "c1", "r1"]].to_numpy(float)), K3[(K3.archetype == a) & (K3.matrix == "F_box")].kappa_boot.iloc[0]) for a in ["A1", "A2", "A3", "A4", "A5", "A6"]]
cases += [("independent_5", np.eye(5), K3[K3.matrix == "identity_5"].kappa_boot.iloc[0]), ("windows_4", np.eye(4), K3[K3.matrix == "identity_4"].kappa_boot.iloc[0])]
rows = []
for ci, (nm, R, kappa) in enumerate(cases):
    for si, (sc, kw) in enumerate([("normal", {}), ("lognormal", {"cv": CV_HI}), ("t", {"df": 3})]):
        X, tr = C.simulate_values(R, C.NREP, np.random.default_rng([C.SEED, 4000 + ci, si]), marginal=sc, **kw); lo, hi = C.exact_boot_ci(X)
        for meth, cv_ in (("C-a", C.coverage_widened(X, tr, kappa, lo, hi)), ("C-b", C.coverage_interval(*C.ci_t_design(X, R), tr)), ("now", C.coverage_interval(lo, hi, tr))):
            se = C.mc_se(cv_, C.NREP); rows.append(dict(case=nm, scenario=sc, method=meth, kappa=kappa, coverage=cv_, se=se, lo=cv_ - 2 * se, hi=cv_ + 2 * se, verdict=verdict(cv_ - 2 * se, cv_ + 2 * se, *C.BRADLEY)))
P = pd.DataFrame(rows); P.to_csv("results/R7/coverage/R7_P1_cells_v3.csv", index=False)
print(P[P.method == "C-a"].pivot_table(index="case", columns="scenario", values="coverage").round(4))
tis = P[P.case.isin(["A1", "A2", "A3", "A4", "A5", "A6"])]
print("R-1 celle in Bradley:", {m: int((tis[tis.method == m].verdict == "confermata").sum()) for m in ("C-a", "C-b")})
