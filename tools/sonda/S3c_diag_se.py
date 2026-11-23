"""tools/sonda/S3c_diag_se.py - S3c, diagnosi (non decide nulla): perche' se di D(B) scende da 0,241 (metro di S3b) a 0,141 (S3c).
1. Reimplementa il bootstrap parametrico di tools/riesame_D.py (RD.draw + RD.q_arm, celle indipendenti, t con df della combinata per V e T,
   t con n-1 per B) e lo verifica sui dati reali contro results/sonda/S3b/D.json (D_se di B e T).
2. Sui dati sintetici del caso «osservato» di S3c (stesso generatore di S3c_D.synthetic, 1000 repliche, seme 20261128): deviazione
   standard vera di D(B), D(T), D(P), D(P)-D(T) fra repliche, contro la media delle se stimate: bootstrap per coppie (i), (ii) = (i) x f,
   parametrico di S3b. Uscita: results/sonda/S3c/diag_se.json"""
import os, sys, json, numpy as np
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools/sonda")
import S3c_D as C
from scipy import stats
TIS = C.TIS; NBP = 1000

def param_se(sd, rng, nb=NBP):
    """se di D(T), D(B), D(P) col bootstrap parametrico di riesame_D (stimatore indipendente, celle indipendenti)."""
    qs = {a: [] for a in ("T", "B", "P")}; obs = {a: [] for a in ("T", "B", "P")}
    for A in TIS:
        X = sd[A]["X"]; h = sd[A]["h"]; V, T = X["V"], X["T"]; n = V.shape[0]
        sp2 = ((V - V.mean(0)) ** 2).sum(0) + ((T - T.mean(0)) ** 2).sum(0); df = 2 * (n - 1); sp = np.sqrt(sp2 / df)
        Vb = V.mean(0) + sp / np.sqrt(n) * rng.standard_t(df, (nb, V.shape[1]))
        Tb = T.mean(0) + sp / np.sqrt(n) * rng.standard_t(df, (nb, V.shape[1]))
        qs["T"].append(((Tb - Vb) ** 2 - 2 * sp ** 2 / n) / h ** 2); obs["T"].append(((T.mean(0) - V.mean(0)) ** 2 - 2 * sp ** 2 / n) / h ** 2)
        for a in ("B", "P"):
            Xa = X[a]; sa = Xa.std(0, ddof=1); Ab = Xa.mean(0) + sa / np.sqrt(n) * rng.standard_t(n - 1, (nb, V.shape[1]))
            qs[a].append(((Ab - Vb) ** 2 - sp ** 2 / n - sa ** 2 / n) / h ** 2); obs[a].append(((Xa.mean(0) - V.mean(0)) ** 2 - sp ** 2 / n - sa ** 2 / n) / h ** 2)
    D = {a: C.ssqrt(np.concatenate(qs[a], axis=-1).mean(-1)) for a in qs}; Do = {a: float(C.ssqrt(np.concatenate(obs[a]).mean())) for a in obs}
    return Do, {a: float(np.std(D[a], ddof=1)) for a in D}, D

data, drop, Vs = C.load("results/sonda/S3b/area_values.csv", "results/sonda/S3b/h.csv")
ref = json.load(open("results/sonda/S3b/D.json"))["by"]["tumour+stroma"]
Do, se_real, _ = param_se(data, np.random.default_rng(20261128), 4000)
out = {"real_reimpl": {"D": Do, "se": se_real}, "real_S3b": {"D": {a: ref["D"][a] for a in ("T", "B", "P")}, "se": {a: ref["D_se"][a] for a in ("T", "B", "P")}}}
print("real reimpl D", {a: round(v, 4) for a, v in Do.items()}, "se", {a: round(v, 4) for a, v in se_real.items()})
print("real S3b    D", {a: round(ref["D"][a], 4) for a in ("T", "B", "P")}, "se", {a: round(ref["D_se"][a], 4) for a in ("T", "B", "P")})
# sintetico, caso osservato
f = float(np.sqrt(5 / 4) * stats.t.ppf(0.975, 4) / stats.norm.ppf(0.975)); base = {}
for A in TIS:
    X = data[A]["X"]; Y = np.hstack([X[a] for a in ("V", "T", "B", "P")]); S = np.cov(Y, rowvar=False)
    base[A] = {"Sig": 0.95 * S + 0.05 * np.diag(np.diag(S)), "mu": Y.mean(0), "m": X["V"].shape[1]}
rng = np.random.default_rng(20261128); rows = []
for r in range(1000):
    sd = {}
    for A in TIS:
        m = base[A]["m"]; Y = rng.multivariate_normal(base[A]["mu"], base[A]["Sig"], size=5)
        sd[A] = {"X": {a: Y[:, j * m:(j + 1) * m] for j, a in enumerate(("V", "T", "B", "P"))}, "h": data[A]["h"]}
    Dp = {a: float(v) for a, v in C.D_of(sd, ("T", "B", "P"), TIS).items()}
    Db = C.D_of(sd, ("T", "B", "P"), TIS, C.boot_idx(np.random.default_rng([20261128, r]), 2000))
    Dind, sep, Dpb = param_se(sd, np.random.default_rng([20261129, r]))
    rows.append({"Dp": Dp, "Dind": Dind, "se_i": {a: float(np.std(Db[a], ddof=1)) for a in Db}, "se_par": sep,
                 "se_i_PT": float(np.std(Db["P"] - Db["T"], ddof=1)), "se_par_PT": float(np.std(Dpb["P"] - Dpb["T"], ddof=1))})
syn = {}
for a in ("T", "B", "P"):
    syn[a] = {"true_sd_paired": float(np.std([x["Dp"][a] for x in rows], ddof=1)), "true_sd_independent": float(np.std([x["Dind"][a] for x in rows], ddof=1)),
              "mean_se_i": float(np.mean([x["se_i"][a] for x in rows])), "mean_se_ii": float(f * np.mean([x["se_i"][a] for x in rows])),
              "mean_se_param_S3b": float(np.mean([x["se_par"][a] for x in rows]))}
syn["P-T"] = {"true_sd_paired": float(np.std([x["Dp"]["P"] - x["Dp"]["T"] for x in rows], ddof=1)), "true_sd_independent": float(np.std([x["Dind"]["P"] - x["Dind"]["T"] for x in rows], ddof=1)),
              "mean_se_i": float(np.mean([x["se_i_PT"] for x in rows])), "mean_se_ii": float(f * np.mean([x["se_i_PT"] for x in rows])), "mean_se_param_S3b": float(np.mean([x["se_par_PT"] for x in rows]))}
out["synthetic_observed"] = syn; out["f"] = f
json.dump(out, open("results/sonda/S3c/diag_se.json", "w"), indent=1)
for k, v in syn.items(): print(k, {kk: round(vv, 3) for kk, vv in v.items()})
