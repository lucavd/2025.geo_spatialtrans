"""tools/sonda/S3c_figure.py - figura di S3c: (a) D per braccio, stimatore di S3b (aree indipendenti) contro appaiato di S3c, con intervalli;
(b) controprova CP-perm: rumore R della permutazione identica contro le 120 permutazioni, per tessuto; (c) taratura dell'errore standard sui dati
sintetici (results/sonda/S3c/diag_se.json). Uscita: results/sonda/S3c/S3c_D.png"""
import os, sys, json, itertools, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools/sonda"); import S3c_D as C
new = json.load(open("results/sonda/S3c/D.json")); old = json.load(open("results/sonda/S3b/D.json")); dg = json.load(open("results/sonda/S3c/diag_se.json"))
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white"})
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8), gridspec_kw={"width_ratios": [1.2, 1, 1]})
arms = ["T", "B", "P"]; lab = {"T": "T (H1, ceiling)", "B": "B (trivial sim, floor)", "P": "P (positive control)"}; x = np.arange(3)
for off, src, col, name in ((-0.15, old, "#999999", "S3b: independent areas"), (0.15, new, "#1f5fa8", "S3c: paired areas")):
    g = src["by"]["tumour+stroma"]; y = [g["D"][a] for a in arms]; ci = [g["D_ci"][a] for a in arms]
    ax[0].errorbar(x + off, y, yerr=[[yy - c[0] for yy, c in zip(y, ci)], [c[1] - yy for yy, c in zip(y, ci)]], fmt="o", color=col, capsize=3, label=name)
fc = new["floor_ceiling"]; ax[0].axhline(fc["tetto"], color="#1f5fa8", lw=0.8, ls=":"); ax[0].axhline(fc["soglia"], color="#c0392b", lw=0.8, ls="--")
ax[0].text(2.45, fc["soglia"], "threshold = ceiling + 1 h", color="#c0392b", va="bottom", ha="right", fontsize=8); ax[0].axhline(0, color="k", lw=0.5)
ax[0].set_xticks(x, [lab[a] for a in arms], fontsize=8); ax[0].set_ylabel("D (h units, lower = closer to H2)"); ax[0].legend(frameon=False, fontsize=8, loc="upper left")
ax[0].set_title("a  D per arm, 95% interval", loc="left")
data, _, _ = C.load("results/sonda/S3b/area_values.csv", "results/sonda/S3b/h.csv")
for i, (A, col) in enumerate((("tumour", "#b5651d"), ("stroma", "#2e8b57"))):
    X = data[A]["X"]; h = data[A]["h"]; R = np.array([((X["T"][list(p)] - X["V"]).var(0, ddof=1) / h ** 2).sum() for p in itertools.permutations(range(5))])
    ax[1].hist(R, bins=30, alpha=0.5, color=col, label=f"{A}: identity rank {int((R < R[0]).sum()) + 1}/120"); ax[1].axvline(R[0], color=col, lw=1.5)
ax[1].set_xlabel("noise R = sum over measures of var(T - V) / h^2"); ax[1].set_ylabel("permutations"); ax[1].legend(frameon=False, fontsize=8)
ax[1].set_title("b  true pairing vs 120 shuffled pairings", loc="left")
s = dg["synthetic_observed"]; keys = ["T", "B", "P", "P-T"]; xx = np.arange(4); w = 0.2
ax[2].bar(xx - 1.5 * w, [s[k]["true_sd_paired"] for k in keys], w, color="k", label="true SD across replicates")
ax[2].bar(xx - 0.5 * w, [s[k]["mean_se_i"] for k in keys], w, color="#9ecae1", label="pair bootstrap (i)")
ax[2].bar(xx + 0.5 * w, [s[k]["mean_se_ii"] for k in keys], w, color="#1f5fa8", label="pair bootstrap x f (ii, used)")
ax[2].bar(xx + 1.5 * w, [s[k]["mean_se_param_S3b"] for k in keys], w, color="#999999", label="S3b parametric")
ax[2].set_xticks(xx, ["D(T)", "D(B)", "D(P)", "D(P) - D(T)"]); ax[2].set_ylabel("standard error (h)"); ax[2].legend(frameon=False, fontsize=7)
ax[2].set_title("c  SE calibration, synthetic data", loc="left")
fig.tight_layout(); fig.savefig("results/sonda/S3c/S3c_D.png", dpi=170); print("saved")
