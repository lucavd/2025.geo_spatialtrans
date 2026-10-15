#!/usr/bin/env python
"""review/code/nonoverlap_sensitivity.py — revisore "code": il conteggio PASS di V-a/V-b (regola, passo 2) dipende
dall'algoritmo di non sovrapposizione? 30 fold x 10 semi; algoritmi onepass (= banco), greedy (riduce il piu' grande),
lp_r (LP max sum R s). Metriche e regola come tools/S1.4b_analysis.R (mediana sui semi; P/W/F col metro; archetipo)."""
import glob, os, sys
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from scipy.stats import spearmanr
K32 = np.sqrt(2 * np.pi / (32 * np.sin(2 * np.pi / 32))); G = 0.99
def ks(a, b):
    a = np.sort(a); b = np.sort(b); v = np.unique(np.concatenate([a, b]))
    return np.max(np.abs(np.searchsorted(a, v, "right") / len(a) - np.searchsorted(b, v, "right") / len(b)))
def algs(x, y, R, dreg):
    n = len(x); pr = cKDTree(np.column_stack([x, y])).query_pairs(2 * R.max() / G, output_type="ndarray")
    i, j = pr[:, 0], pr[:, 1]; dij = np.hypot(x[i] - x[j], y[i] - y[j]); bad = R[i] + R[j] > G * dij; i, j, dij = i[bad], j[bad], dij[bad]
    out = {}
    f_ = G * dij / (R[i] + R[j]); sp = np.ones(n); np.minimum.at(sp, i, f_); np.minimum.at(sp, j, f_); out["onepass"] = np.minimum(sp, np.minimum(1, G * dreg / R))
    s = np.minimum(1, G * dreg / R)
    for it in range(2000):
        r = s * R; v = r[i] + r[j] - G * dij; vi = np.where(v > 1e-12)[0]
        if not len(vi): break
        for k in vi[np.argsort(-v[vi])]:
            a, b = i[k], j[k]; ra, rb = s[a] * R[a], s[b] * R[b]
            if ra + rb <= G * dij[k] + 1e-12: continue
            big, sml = (a, b) if ra >= rb else (b, a); rs = s[sml] * R[sml]
            if rs < G * dij[k]: s[big] = (G * dij[k] - rs) / R[big]
            else: f2 = G * dij[k] / (ra + rb); s[a] *= f2; s[b] *= f2
    out["greedy"] = s
    m = len(i); A = coo_matrix((np.concatenate([R[i], R[j]]), (np.r_[np.arange(m), np.arange(m)], np.r_[i, j])), shape=(m, n)).tocsr()
    res = linprog(-R, A_ub=A, b_ub=G * dij, bounds=np.column_stack([np.zeros(n), np.minimum(1, G * dreg / R)]), method="highs")
    out["lp_r"] = np.clip(res.x, 0, 1)
    return out
def fold(f):
    d = pd.read_parquet(f); A, k = os.path.basename(f).replace(".parquet", "").split("_fold"); ev = d.ev.to_numpy()
    rnc = (d.a_n / d.a_t).to_numpy()[ev]; ran = d.a_n.to_numpy()[ev]; at = d.a_t.to_numpy()[ev]; rsp = spearmanr(ran, at)[0]; dg = d.d_gen.to_numpy()[ev]
    rows = []
    for s in range(1, 11):
        tgt = d[f"tgt{s}"].to_numpy(); R = K32 * np.sqrt(tgt / np.pi)
        for nm, sv in algs(d.x.to_numpy(), d.y.to_numpy(), R, d.d_reg.to_numpy()).items():
            af = (sv**2 * tgt)[ev]
            rows.append(dict(A=A, k=int(k), seed=s, alg=nm, ks_nc=ks(af / at, rnc), ks_an=ks(af, ran), d_sp=abs(spearmanr(af, at)[0] - rsp),
                             trunc_all=float(np.mean(sv < 1 - 1e-12)) if nm != "onepass" else float(np.mean(sv < 1)), loss=float(np.mean(1 - sv**2)), b7=float(np.mean(np.sqrt(af / np.pi) > dg))))
    return rows
if __name__ == "__main__":
    fs = sorted(glob.glob("/tmp/code_tgt/*_fold*.parquet"))
    with Pool(30) as p: R = [r for rr in p.map(fold, fs) for r in rr]
    M = pd.DataFrame(R); M.to_csv("results/S1.4b/review/code/nonoverlap_sens_raw.csv", index=False)
    Y = pd.read_csv("/tmp/code_tgt/yards.csv")
    Fd = M.groupby(["alg", "A", "k"]).median(numeric_only=True).reset_index().merge(Y, on=["A", "k"])
    Fd["b7d"] = (M.assign(b7d=lambda z: z.b7).groupby(["alg", "A", "k"]).apply(lambda z: np.median(np.abs(z.b7 - Y.set_index(["A", "k"]).loc[(z.A.iloc[0], z.k.iloc[0]), "b7_real"])), include_groups=False)).reset_index(drop=True)
    cl = lambda v, a, b: np.where(v <= a, "P", np.where(v <= b, "W", "F"))
    Fd["v1"] = cl(Fd.ks_nc, Fd.y_nc_med, Fd.y_nc_max); Fd["v2"] = cl(Fd.ks_an, Fd.y_an_med, Fd.y_an_max); Fd["v3"] = cl(Fd.d_sp, Fd.y_sp_med, Fd.y_sp_max); Fd["v7"] = cl(Fd.b7d, Fd.y_b7_med, Fd.y_b7_max)
    av = lambda o: "PASS" if (o == "P").sum() >= 4 else ("FAIL" if (o == "F").sum() >= 2 else "WARN")
    At = Fd.sort_values("k").groupby(["alg", "A"]).apply(lambda z: pd.Series(dict(B1=av(z.v1), B2=av(z.v2), B3=av(z.v3), B7=av(z.v7), f1="".join(z.v1), f2="".join(z.v2), f3="".join(z.v3), f7="".join(z.v7),
                                                       B4=z.trunc_all.mean(), loss=z.loss.mean())), include_groups=False).reset_index()
    At.to_csv("results/S1.4b/review/code/nonoverlap_sens_archetypes.csv", index=False); Fd.to_csv("results/S1.4b/review/code/nonoverlap_sens_folds.csv", index=False)
    pd.set_option("display.width", 200); print(At.to_string())
    for a in At.alg.unique(): print(a, "score", int((At[At.alg == a][["B1", "B2", "B3", "B7"]] == "PASS").sum().sum()))
