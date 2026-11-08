"""tools/R7_figures.py — R7: figure del report -> results/R7/figures/*.png"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "results/R7/figures"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
NAMES = {"A1": "intestine (mouse)", "A2": "colon tumour", "A3": "tumour stroma", "A4": "lymph node", "A5": "brain (mouse)", "A6": "heart"}
ARCH = list(NAMES)
# F1: copertura
K = pd.read_csv("results/R7/coverage/R7_kappa.csv").set_index("case"); P1 = pd.read_csv("results/R7/coverage/R7_P1_cells.csv")
P2 = pd.read_csv("results/R7/cp1/R7_P2.csv")
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(6); w = 0.18
ax.axhspan(0.925, 0.975, color="0.9", label="Bradley 0.925-0.975")
ax.bar(x - 1.5 * w, K.loc[ARCH, "cover_now"], w, color="#b2182b", label="bootstrap over 5 areas (current)")
ax.bar(x - 0.5 * w, K.loc[ARCH, "cover_minmax"], w, color="#ef8a62", label="min-max of 5 areas")
ca = P1[(P1.method == "C-a") & (P1.scenario == "normal")].set_index("case").loc[ARCH, "coverage"]
ax.bar(x + 0.5 * w, ca, w, color="#2166ac", label="corrected, geometry only (normal values)")
p2 = P2[P2.field == "umi"].set_index("archetype").reindex(ARCH)
ax.bar(x + 1.5 * w, p2.coverage, w, color="#67a9cf", label="geometry-corrected, under measured tissue structure (UMI)")
ax.errorbar(x + 1.5 * w, p2.coverage, yerr=[p2.coverage - p2.cov_lo, p2.cov_hi - p2.coverage], fmt="none", ecolor="k", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels([NAMES[a] for a in ARCH], rotation=20, ha="right"); ax.set_ylim(0.5, 1.0); ax.set_ylabel("coverage of a nominal 95 % interval")
ax.legend(fontsize=7.5, loc="lower left", frameon=False); fig.tight_layout(); fig.savefig(f"{OUT}/R7_F1_coverage.png", dpi=160)
# F2: correlazione implicita contro quota in comune
pr = pd.read_csv("results/R7/cp1/R7_cp1_pairs.csv")
fig, axs = plt.subplots(1, 5, figsize=(13, 3), sharex=True, sharey=True)
for ax, a in zip(axs, ARCH[:5]):
    for f, mk, col in (("umi", "o", "#2166ac"), ("nuclei", "s", "#b2182b")):
        d = pr[(pr.archetype == a) & (pr.field == f)]
        ax.scatter(d.F_box, d.R_implied, marker=mk, color=col, s=22, label={"umi": "UMI per bin", "nuclei": "nuclei per mm²"}[f])
    ax.plot([0, 1], [0, 1], color="0.5", lw=0.8, ls="--"); ax.set_title(NAMES[a], fontsize=9); ax.set_xlim(-0.02, 0.5); ax.set_ylim(-0.05, 1)
    ax.set_xlabel("shared area fraction")
axs[0].set_ylabel("implied correlation of area means"); axs[0].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/R7_F2_implied_vs_shared.png", dpi=160)
# F3: semivariogrammi (normalizzati alla soglia del modello)
vg = pd.read_csv("results/R7/cp1/R7_cp1_variograms.csv")
fig, axs = plt.subplots(1, 5, figsize=(13, 3), sharey=True)
for ax, a in zip(axs, ARCH[:5]):
    for f, col in (("umi", "#2166ac"), ("nuclei", "#b2182b")):
        d = vg[(vg.archetype == a) & (vg.field == f)]; pp = P2[(P2.archetype == a) & (P2.field == f)].iloc[0]; s = pp.c0 + pp.c1
        ax.plot(d.h_um, d.gamma / s, "o", color=col, ms=3); ax.plot(d.h_um, d.model / s, "-", color=col, lw=1, label=f"{f}: L = {pp.L_um:.0f} µm")
    ax.set_title(NAMES[a], fontsize=9); ax.set_xlabel("distance between 48 µm tiles (µm)"); ax.legend(fontsize=7, frameon=False)
axs[0].set_ylabel("semivariance / model sill")
fig.tight_layout(); fig.savefig(f"{OUT}/R7_F3_variograms.png", dpi=160)
# F4: P-3
P3 = pd.read_csv("results/R7/pieces/R7_P3.csv"); P3 = P3[P3.archetype != "A6"]
fig, ax = plt.subplots(figsize=(8, 3.6)); qs = list(dict.fromkeys(P3.quantity))
for i, q in enumerate(qs):
    d = P3[P3.quantity == q].set_index("archetype").reindex(ARCH[:5])
    ax.scatter(np.arange(5) + (i - 3) * 0.1, np.sign(d["shift"]) * d.shift_over_bound, s=20, label=q)
ax.axhline(1, color="k", lw=0.8, ls="--"); ax.axhline(-1, color="k", lw=0.8, ls="--"); ax.axhline(0, color="0.6", lw=0.6)
ax.set_xticks(range(5)); ax.set_xticklabels([NAMES[a] for a in ARCH[:5]], rotation=15); ax.set_ylim(-1.2, 1.2)
ax.set_ylabel("shift / half of corrected half-width"); ax.legend(fontsize=7, ncol=4, frameon=False, loc="lower center")
fig.tight_layout(); fig.savefig(f"{OUT}/R7_F4_pieces_shift.png", dpi=160)
print("ok")
