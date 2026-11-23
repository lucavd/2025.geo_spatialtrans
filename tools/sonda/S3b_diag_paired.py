"""tools/sonda/S3b_diag_paired.py - S3b, diagnosi (non decide nulla): T omologo e V sono appaiati per area (k). Per ogni voce:
correlazione fra T_k e V_k sulle 5 aree; termine di rumore dello stimatore (sp^2 / nT + sp^2 / nV, insiemi indipendenti, come tools/riesame_D.py)
contro quello appaiato var(T_k - V_k) / n; termine di T in h^2 nei due modi. Uscita: results/sonda/S3b/diag_paired.csv"""
import os, sys, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); import riesame_D as RD
U = pd.read_csv("results/sonda/S3b/area_values.csv"); H = pd.read_csv("results/sonda/S3b/h.csv"); rows = []
for A in ("tumour", "stroma"):
    V = U[(U.arm == "V") & (U.A == A)].copy(); T = U[(U.arm == "T") & (U.A == A)].copy()
    V["k"] = V.unit.str.rsplit("_", n=1).str[-1].astype(int); T["k"] = T.unit.str.rsplit("_", n=1).str[-1].astype(int)
    V, T = V.sort_values("k"), T.sort_values("k"); assert (V.k.values == T.k.values).all()
    for m, col in zip(["M-1", "M-2", "M-3", "M-4", "M-5", "M-6"], ["M1", "M2", "M3", "M4", "M5", "M6"]):
        v, t = V[col].to_numpy(float), T[col].to_numpy(float); h = float(H[(H.A == A) & (H.measure == m)].h.iloc[0]); n = len(v)
        sp2 = (v.var(ddof=1) + t.var(ddof=1)) / 2; d = t - v
        q_ind = ((t.mean() - v.mean()) ** 2 - sp2 / n - sp2 / n) / h ** 2; q_pair = (d.mean() ** 2 - d.var(ddof=1) / n) / h ** 2
        rows.append(dict(A=A, measure=m, r_TV=float(np.corrcoef(t, v)[0, 1]), noise_ind=2 * sp2 / n, noise_pair=d.var(ddof=1) / n, diff2=float(d.mean() ** 2), q_ind=q_ind, q_pair=q_pair))
P = pd.DataFrame(rows); P.to_csv("results/sonda/S3b/diag_paired.csv", index=False)
print(P.round(3).to_string(index=False))
qi, qp = P.q_ind.mean(), P.q_pair.mean(); print("mean q: independent", round(qi, 3), "signed root", round(np.sign(qi) * abs(qi) ** 0.5, 3), "| paired", round(qp, 3), "signed root", round(np.sign(qp) * abs(qp) ** 0.5, 3))
print("median r_TV", round(P.r_TV.median(), 2), "share noise_pair < noise_ind", (P.noise_pair < P.noise_ind).mean())
