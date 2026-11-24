"""tools/sonda/S4_figure.py - figura di S4: (a) D per braccio (T, B, Y, G, U, A, AE) con l'intervallo al 95 % (metodo di S3c), tetto, soglia,
pavimento - delta; (b) per voce (tessuto x misura), quota del divario B - T chiusa dalle sostituzioni singole Y, G, U e da A, AE.
Dati: results/sonda/S4/D.json, results/sonda/S4/per_measure.csv. Uscita: results/sonda/S4/S4_D.png"""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
R = json.load(open("results/sonda/S4/D.json")); P = pd.read_csv("results/sonda/S4/per_measure.csv"); g = R["by"]["tumour+stroma"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white"})
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1.15, 1.5]})
arms = ["T", "B", "Y", "G", "U", "A", "AE"]
lab = {"T": "T\nH1\n(ceiling)", "B": "B\ntrivial\n(floor)", "Y": "Y\ntrue\ntype map", "G": "G\ntrue\nnuclei", "U": "U\ntrue\nUMI map",
       "A": "A\nnuclei +\ntypes + UMI", "AE": "AE\nA + true\ncell expr."}
col = {"T": "#1f5fa8", "B": "#555555", "Y": "#d95f02", "G": "#7570b3", "U": "#1b9e77", "A": "#e7298a", "AE": "#a6761d"}
x = np.arange(len(arms)); y = [g["D"][a] for a in arms]; ci = [g["D_ci"][a] for a in arms]
for i, a in enumerate(arms):
    ax[0].errorbar(i, y[i], yerr=[[y[i] - ci[i][0]], [ci[i][1] - y[i]]], fmt="o", color=col[a], capsize=3)
tet, pav = g["D"]["T"], g["D"]["B"]
for v, txt, c, ls in ((tet + 1.0, "threshold = ceiling + 1 h", "#c0392b", "--"), (pav - 1.0, "floor - 1 h", "#777777", ":")):
    ax[0].axhline(v, color=c, lw=0.8, ls=ls); ax[0].text(len(arms) - 0.5, v, txt, color=c, va="bottom", ha="right", fontsize=8)
ax[0].axhline(0, color="k", lw=0.5); ax[0].set_xticks(x, [lab[a] for a in arms], fontsize=7.5)
ax[0].set_ylabel("D (h units, lower = closer to H2)"); ax[0].set_title("a  D per arm, tumour + stroma, 95% interval", loc="left")
P = P[P.B_minus_T_q > 0.5].copy(); P["cell"] = P.A + " " + P.measure; sub = ["Y", "G", "U", "A", "AE"]; w = 0.16; xx = np.arange(len(P))
for j, a in enumerate(sub):
    ax[1].bar(xx + (j - 2) * w, P[f"{a}_closed_of_gap"].clip(-1, 1.2), w, color=col[a], label=lab[a].replace("\n", " "))
ax[1].text(0.99, 0.02, "bars clipped to [-1, 1.2]", transform=ax[1].transAxes, ha="right", fontsize=7, color="#555555")
ax[1].axhline(1, color="k", lw=0.5, ls=":"); ax[1].axhline(0, color="k", lw=0.5)
ax[1].set_xticks(xx, [f"{c}\n({s:.0%} of gap)" for c, s in zip(P.cell, P.B_minus_T_q / P.B_minus_T_q.sum())], fontsize=7.5)
ax[1].set_ylabel("share of the B - T gap closed (q units; 1 = reaches ceiling)"); ax[1].legend(frameon=False, fontsize=7, ncol=2)
ax[1].set_title("b  per measure: which true component closes the gap", loc="left")
fig.tight_layout(); fig.savefig("results/sonda/S4/S4_D.png", dpi=170); print("saved")
