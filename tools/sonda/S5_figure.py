"""tools/sonda/S5_figure.py - S5 (docs/sonda/REGOLA_S5.md): figura grezza di A3. D_p per braccio (Tp, XB, XV, Bp, B0p, B1p, Pp) con l'intervallo (ii),
fascia del tetto del pannello (intervallo di D_p(T)) e linea tetto + 1 h; un pannello per tessuto e uno combinato; A3 nel titolo del combinato.
Fondo bianco. Uscita: results/sonda/S5/fig_A3.png. Uso: .venv/bin/python -B tools/sonda/S5_figure.py"""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); R = json.load(open("results/sonda/S5/D.json")); LAB = R["labels"]
ORDER = ["T", "XB", "XV", "B", "B0", "B1", "P"]; COL = {"T": "#1b7837", "XB": "#2166ac", "XV": "#67a9cf", "B": "#b2182b", "B0": "#ef8a62", "B1": "#ef8a62", "P": "#555555"}
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True, facecolor="white")
for ax, g in zip(axs, ("tumour", "stroma", "tumour+stroma")):
    r = R["by"][g]; ax.set_facecolor("white"); tl, th = r["D_ci"]["T"]
    ax.axhspan(tl, th, color="#1b7837", alpha=0.12, lw=0, label="tetto del pannello (intervallo di Tp)")
    ax.axhline(r["D"]["T"] + 1.0, color="#1b7837", ls="--", lw=1, label="tetto + 1 h"); ax.axhline(0, color="k", lw=0.6)
    for i, k in enumerate(ORDER):
        lo, hi = r["D_ci"][k]; ax.errorbar(i, r["D"][k], yerr=[[r["D"][k] - lo], [hi - r["D"][k]]], fmt="o", color=COL[k], capsize=3, ms=6)
        ax.text(i, hi + 0.08, f"{r['D'][k]:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([LAB[k] for k in ORDER]); ax.set_title(g if g != "tumour+stroma" else
        f"tumour+stroma   A3 = {r['A3']['value']:.2f} [{r['A3']['lo']:.2f}; {r['A3']['hi']:.2f}], se {r['A3']['se']:.2f}", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("D_p (unita' h_p, intervallo ii)"); axs[0].legend(frameon=False, fontsize=8, loc="upper left")
fig.suptitle("S5: D sul pannello G_p (284 geni), X1 contro Vp; A3 = D_p(XB) - D_p(Bp)", fontsize=11); fig.tight_layout()
fig.savefig("results/sonda/S5/fig_A3.png", dpi=150, facecolor="white"); print("fig ok")
