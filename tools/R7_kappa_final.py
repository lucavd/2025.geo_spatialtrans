"""tools/R7_kappa_final.py — R7: tabella finale dei fattori (media di 5 semi, regola R-2) e fattore per l'asimmetria (PH-1).
-> results/R7/coverage/R7_kappa_final.csv, R7_cv_factor.csv"""
import pandas as pd, numpy as np
K = pd.read_csv("results/R7/coverage/R7_kappa_types.csv"); MC = pd.read_csv("results/R7/coverage/R7_kappa_mc.csv")
rows = []
for A, g in K.groupby("archetype"):
    d = g.set_index("matrix")
    rows.append(dict(archetype=A, unit="area", kappa_box=d.loc["F_box", "kappa_mean5"], kappa_bin=d.loc["F_bin", "kappa_mean5"], kappa_nuc=d.loc["F_nuc", "kappa_mean5"],
                     kappa_se=float(g.kappa_sd5.max() / np.sqrt(5)), kappa_range_box=d.loc["F_box", "kappa_range"], kappa_range_bin=d.loc["F_bin", "kappa_range"],
                     kappa_range_nuc=d.loc["F_nuc", "kappa_range"]))
w = MC.set_index("case").loc["windows_4"]
rows.append(dict(archetype="all", unit="window_4", kappa_box=w.kappa_mean5, kappa_bin=w.kappa_mean5, kappa_nuc=w.kappa_mean5, kappa_se=w.kappa_sd5 / np.sqrt(5),
                 kappa_range_box=1.0, kappa_range_bin=1.0, kappa_range_nuc=1.0))
pd.DataFrame(rows).to_csv("results/R7/coverage/R7_kappa_final.csv", index=False)
P = pd.read_csv("results/R7/coverage/R7_PH1_cv.csv"); P = P[P.case == "independent_5"]
P["cv_factor"] = np.maximum(1.0, P.kappa_needed_lognormal / P.kappa_normal)
P[["cv", "coverage_Ca", "se", "kappa_needed_lognormal", "kappa_normal", "cv_factor"]].to_csv("results/R7/coverage/R7_cv_factor.csv", index=False)
print(pd.DataFrame(rows).round(4).to_string(index=False)); print(P[["cv", "coverage_Ca", "cv_factor"]].round(4).to_string(index=False))
