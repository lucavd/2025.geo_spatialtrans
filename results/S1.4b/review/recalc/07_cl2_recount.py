
import pandas as pd, numpy as np
from scipy.stats import norm
a = pd.read_csv("results/S1.4b/S1.4b_arbiter.csv"); a = a[a.mutant=="none"]
E = 2*norm.cdf(-4)*a.n_cells
thr = E + 3*np.sqrt(E)
ap = a.n_absz_gt4 <= thr
pas = ap & (a.misassigned==0) & (a.uncovered==0) & (a.double==0)
print("n sets", len(a), "PASS", pas.sum(), "sum obs", a.n_absz_gt4.sum(), "sum E %.2f"%E.sum(), "maxdiff E", np.abs(E-a.exp_absz_gt4).max(), "agree pass_", (pas==a.pass_).all())
print("failed:", list(a.set[~pas]))
b = pd.read_csv("results/S1.4b/S1.4b_posthoc_arbiter_PHB.csv")
print("PH-B obs > null_q99:", list(b.set[b.observed>b.null_q99]), " p<0.05:", list(b.set[b.p_value<0.05]))
print("Bonferroni p<0.05/30:", list(b.set[b.p_value<0.05/30]))
