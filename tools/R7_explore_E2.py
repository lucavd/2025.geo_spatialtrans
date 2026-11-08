"""tools/R7_explore_E2.py — R7, esplorazione E-2 (prima della pre-registrazione).
Copertura dell'intervallo bootstrap percentile della media di 5 aree (come R4_common.boot_one_stage)
e dell'intervallo min-max, con valori per area normali di correlazione nota: indipendenti, e con la
correlazione uguale alla quota condivisa (geometria delle aree di R2). Bootstrap esatto (tutte le 5^5
estrazioni ordinate, equivalente a B infinito). Uscita: results/R7/explore/R7_E2_coverage.csv"""
import itertools, numpy as np, pandas as pd, os
OUT = "results/R7/explore"; os.makedirs(OUT, exist_ok=True)
R = pd.read_csv("results/R2/R2_rois_checked.csv")
idx = np.array(list(itertools.product(range(5), repeat=5)))          # 3125 x 5
W = np.stack([np.bincount(r, minlength=5) for r in idx]) / 5.0       # pesi per estrazione
def wpct(v, q):  # percentili esatti della distribuzione bootstrap (3125 medie equiprobabili)
    m = np.sort(W @ v); return np.percentile(m, q)
def Fmat(g):
    n = len(g); F = np.eye(n)
    for i, j in itertools.combinations(range(n), 2):
        a, b = g.iloc[i], g.iloc[j]
        w = max(0, min(a.c1, b.c1) - max(a.c0, b.c0)); h = max(0, min(a.r1, b.r1) - max(a.r0, b.r0))
        F[i, j] = F[j, i] = w * h / ((a.c1 - a.c0) * (a.r1 - a.r0))
    return F
rng = np.random.default_rng(20261108)
NREP = 4000
rows = []
cases = {"independent": np.eye(5)}
for A, g in R.groupby("archetype"): cases[A + "_overlap"] = Fmat(g.reset_index(drop=True))
for name, S in cases.items():
    L = np.linalg.cholesky(S)
    Z = rng.standard_normal((NREP, 5)) @ L.T
    M = Z @ W.T                                   # NREP x 3125
    lo = np.percentile(M, 2.5, axis=1); hi = np.percentile(M, 97.5, axis=1)
    cov_boot = np.mean((lo <= 0) & (hi >= 0))
    cov_mm = np.mean((Z.min(1) <= 0) & (Z.max(1) >= 0))
    n = 5; rb = S[~np.eye(n, dtype=bool)].mean()
    rows.append(dict(case=name, rho_bar=rb, k_theory=np.sqrt((1+(n-1)*rb)/(1-rb)),
                     cover_boot=cov_boot, mc_se=np.sqrt(cov_boot*(1-cov_boot)/NREP), cover_minmax=cov_mm,
                     halfwidth_ratio_needed=None))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/R7_E2_coverage.csv", index=False)
