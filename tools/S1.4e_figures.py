#!/usr/bin/env python
"""tools/S1.4e_figures.py - S1.4e: figure del report (sfondo bianco, un pannello per tessuto dove serve).
Uso (root del repo): .venv/bin/python tools/S1.4e_figures.py   Uscite: results/S1.4e/figures/S1.4e_fig*.png"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
R = Path("results/S1.4e"); F = R / "figures"; F.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white"})
TIS = {"A1": "intestino di topo", "A2": "tumore del colon", "A3": "stroma del tumore", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}
COL = {"K2": "#4d4d4d", "ER": "#3b75af", "ET": "#e08214"}; LBL = {"K2": "disco di oggi", "ER": "ellisse a caso", "ET": "ellisse lungo il territorio"}
V = pd.read_csv(R / "S1.4e_verdicts.csv"); D = pd.read_csv(R / "S1.4e_desc.csv")
def band_plot(pid, tissues, fname, ylab, title):
    fig, ax = plt.subplots(figsize=(6.4, 2.8)); x = np.arange(len(tissues))
    for k, A in enumerate(tissues):
        r = V[(V.id == pid) & (V.archetype == A)].iloc[0]; ax.add_patch(plt.Rectangle((k - 0.35, r.lo), 0.7, r.hi - r.lo, color="0.88", zorder=0))
        for j, arm in enumerate(["K2", "ER", "ET"]):
            q = V[(V.id == pid) & (V.archetype == A) & (V.arm == arm)].iloc[0]
            ax.errorbar(k + (j - 1) * 0.2, q.est, yerr=2 * q.se, fmt="o", ms=4, color=COL[arm], label=LBL[arm] if k == 0 else None, capsize=2)
    ax.set_xticks(x, [TIS[A] for A in tissues], rotation=15); ax.set_ylabel(ylab); ax.set_title(title, loc="left"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(F / fname, dpi=200); plt.close(fig)
band_plot("P-1", list(TIS), "S1.4e_fig1_contacts.png", "nuclei a contatto (frazione)", "Nuclei con almeno un vicino a contatto (banda: riferimento manuale)")
band_plot("P-3", ["A1", "A2", "A3", "A4", "A6"], "S1.4e_fig3_alignment.png", "di traverso al vicino (frazione)",
          "Asse lungo di traverso al vicino a meno di 5 µm (banda: ROI vere)")
c = D[D.desc == "banco A curva"]; LAB = ["<5", "5-6", "6-7", "7-8", "8-10"]
fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.4), sharey=True)
for ax, A in zip(axs, ["A1", "A2", "A3", "A4"]):
    for arm, col, mk in [("manuale", "black", "o"), ("K2", COL["K2"], "s"), ("ER", COL["ER"], "^"), ("ET", COL["ET"], "D")]:
        g = c[(c.archetype == A) & (c.arm == arm)].set_index("bin").reindex(LAB)
        ax.plot(range(5), g.est, marker=mk, ms=3, color=col, lw=1.4 if arm in ("manuale", "ET") else 1, label={"manuale": "giudizio di Luca"}.get(arm, LBL.get(arm)))
    ax.set_xticks(range(5), LAB, fontsize=6.5); ax.set_title(TIS[A], loc="left", fontsize=8); ax.set_xlabel("distanza fra i centri (µm)", fontsize=7)
axs[0].set_ylabel("coppie a contatto (frazione)"); axs[0].legend(frameon=False, fontsize=6.5)
fig.suptitle("Stessi centri dei veri (i punti di Luca): contatto per distanza", x=0.01, ha="left", fontsize=8.5)
fig.tight_layout(); fig.savefig(F / "S1.4e_fig2_same_centres.png", dpi=200); plt.close(fig)
print("figure scritte:", sorted(p.name for p in F.glob("*.png")))
