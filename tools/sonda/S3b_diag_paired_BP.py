"""tools/sonda/S3b_diag_paired_BP.py - S3b, diagnosi (non decide nulla): come S3b_diag_paired.py per i bracci B (media dei 2 semi per area)
e P, appaiati a V per area. Termine del braccio in h^2: indipendente (formula di tools/sonda/S3_D.py per_measure) contro appaiato.
Uscita: results/sonda/S3b/diag_paired_BP.csv"""
import os, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
U = pd.read_csv("results/sonda/S3b/area_values.csv"); H = pd.read_csv("results/sonda/S3b/h.csv"); rows = []
key = lambda s: s.str.extract(r"_(tumour|stroma)_(\d+)")[1].astype(int)
for arm in ("B", "P"):
    for A in ("tumour", "stroma"):
        V = U[(U.arm == "V") & (U.A == A)].copy(); X = U[(U.arm == arm) & (U.A == A)].copy(); V["k"] = key(V.unit); X["k"] = key(X.unit)
        cols = ["M1", "M2", "M3", "M4", "M5", "M6"]; X = X.groupby("k")[cols].mean().sort_index(); V = V.set_index("k").sort_index()
        for m, col in zip(["M-1", "M-2", "M-3", "M-4", "M-5", "M-6"], cols):
            v, x = V[col].to_numpy(float), X[col].to_numpy(float); h = float(H[(H.A == A) & (H.measure == m)].h.iloc[0]); n = len(v); d = x - v
            sp2V = v.var(ddof=1)
            rows.append(dict(arm=arm, A=A, measure=m, r=float(np.corrcoef(x, v)[0, 1]), q_ind=((x.mean() - v.mean()) ** 2 - sp2V / n - x.var(ddof=1) / n) / h ** 2,
                             q_pair=(d.mean() ** 2 - d.var(ddof=1) / n) / h ** 2))
P = pd.DataFrame(rows); P.to_csv("results/sonda/S3b/diag_paired_BP.csv", index=False)
for arm, g in P.groupby("arm"):
    qi, qp = g.q_ind.mean(), g.q_pair.mean()
    print(arm, "median r", round(g.r.median(), 2), "| independent root", round(np.sign(qi) * abs(qi) ** .5, 3), "| paired root", round(np.sign(qp) * abs(qp) ** .5, 3))
