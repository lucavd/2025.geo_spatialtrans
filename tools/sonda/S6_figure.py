"""tools/sonda/S6_figure.py - figura di S6: (a) D per braccio (T, B, Y, Yalt, YH1, Yshuf) con l'intervallo al 95 % (metodo di S3c), tetto + 1 h,
pavimento - 1 h (soglia della regola di RIESAME_C par. 6) e pavimento - 0,5 h; (b) per voce (tessuto x misura) il contributo q a D^2 di B, Y, Yalt,
YH1, Yshuf, T. Dati: results/sonda/S6/D.json, results/sonda/S6/per_measure.csv. Uscita: results/sonda/S6/S6_D.png"""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
R = json.load(open("results/sonda/S6/D.json")); P = pd.read_csv("results/sonda/S6/per_measure.csv"); g = R["by"]["tumour+stroma"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white"})
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1.0, 1.6]})
arms = ["T", "B", "Y", "Yalt", "YH1", "Yshuf"]
lab = {"T": "T\nH1\n(ceiling)", "B": "B\ntrivial\n(floor)", "Y": "Y\ntypes copied\nfrom measured area", "Yalt": "Yalt\ntypes from\nanother area",
       "YH1": "YH1\ntypes from\nhomologous H1 area", "Yshuf": "Yshuf\nmeasured-area types\nshuffled"}
col = {"T": "#1f5fa8", "B": "#555555", "Y": "#d95f02", "Yalt": "#1b9e77", "YH1": "#7570b3", "Yshuf": "#e7298a"}
y = [g["D"][a] for a in arms]; ci = [g["D_ci"][a] for a in arms]
for i, a in enumerate(arms):
    ax[0].errorbar(i, y[i], yerr=[[y[i] - ci[i][0]], [ci[i][1] - y[i]]], fmt="o", color=col[a], capsize=3)
tet, pav = g["D"]["T"], g["D"]["B"]
for v, txt, c, ls in ((tet + 1.0, "threshold = ceiling + 1 h", "#c0392b", "--"), (pav - 1.0, "floor - 1 h (rule, delta = 1)", "#333333", ":"),
                      (pav - 0.5, "floor - 0.5 h (sensitivity)", "#999999", ":")):
    ax[0].axhline(v, color=c, lw=0.8, ls=ls); ax[0].text(len(arms) - 0.5, v, txt, color=c, va="bottom", ha="right", fontsize=7.5)
ax[0].axhline(0, color="k", lw=0.5); ax[0].set_xticks(np.arange(len(arms)), [lab[a] for a in arms], fontsize=7)
ax[0].set_ylabel("D (h units, lower = closer to H2)"); ax[0].set_title("a  D per arm, tumour + stroma, 95% interval", loc="left")
P = P[P.B_minus_T_q > 0.5].copy(); P["cell"] = P.A + " " + P.measure; sub = ["B", "Y", "Yalt", "YH1", "Yshuf", "T"]; w = 0.13; xx = np.arange(len(P))
for j, a in enumerate(sub):
    ax[1].bar(xx + (j - 2.5) * w, P[f"{a}_q"], w, color=col[a], label=lab[a].replace("\n", " "))
ax[1].axhline(0, color="k", lw=0.5); ax[1].set_xticks(xx, list(P.cell), fontsize=8)
ax[1].set_ylabel("contribution q to D^2 (paired, h units squared)"); ax[1].legend(frameon=False, fontsize=7, ncol=2)
ax[1].set_title("b  per measure: the type map acts on M-4; M-1 is what remains", loc="left")
fig.tight_layout(); fig.savefig("results/sonda/S6/S6_D.png", dpi=170); print("saved")
