#!/usr/bin/env python
"""tools/S1.4f_figures.py - S1.4f: figure del report (sfondo bianco, un pannello per tessuto dove serve).
Uso (root del repo): .venv/bin/python tools/S1.4f_figures.py   Uscite: results/S1.4f/figures/S1.4f_fig*.png"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
R = Path("results/S1.4f"); F = R / "figures"; F.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white"})
TIS = {"A1": "intestino di topo", "A2": "tumore del colon", "A3": "stroma del tumore", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}
COL = {"K2": "#4d4d4d", "ER": "#3b75af", "ET": "#e08214", "EP": "#1b9e77", "V0": "#7570b3", "V1": "#e7298a", "V2": "#1b9e77"}
LBL = {"K2": "disco di oggi", "ER": "ellisse a caso", "ET": "ellisse esattamente lungo il territorio", "EP": "ellisse orientata in parte (compressione)",
       "V0": "dispersione costante", "V1": "legata all'allungamento del territorio", "V2": "legata alla compressione (scelta)"}
V = pd.read_csv(R / "S1.4f_verdicts.csv"); CP = pd.read_csv(R / "S1.4f_cp.csv")
# ---- figura 1: contatti ----
fig, ax = plt.subplots(figsize=(6.4, 2.8)); tissues = list(TIS)
for k, A in enumerate(tissues):
    r = V[(V.id == "P-1") & (V.archetype == A)].iloc[0]; ax.add_patch(plt.Rectangle((k - 0.35, r.lo), 0.7, r.hi - r.lo, color="0.88", zorder=0))
    for j, arm in enumerate(["K2", "EP", "ET"]):
        q = V[(V.id == "P-1") & (V.archetype == A) & (V.arm == arm)].iloc[0]
        ax.errorbar(k + (j - 1) * 0.2, q.est, yerr=2 * q.se, fmt="o", ms=4, color=COL[arm], label=LBL[arm] if k == 0 else None, capsize=2)
ax.set_xticks(range(len(tissues)), [TIS[A] for A in tissues], rotation=15); ax.set_ylabel("nuclei a contatto (frazione)")
ax.set_title("Nuclei con almeno un vicino a contatto (banda: riferimento manuale)", loc="left"); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(F / "S1.4f_fig1_contacts.png", dpi=200); plt.close(fig)
# ---- figura 2: allineamento per classe di distanza ----
SH = pd.read_csv(R / "S1.4f_shape_B.csv"); F2 = pd.read_csv(R / "explore/S1.4f_F2_units.csv"); RF = pd.read_csv(R / "S1.4f_refs.csv")
PH = pd.read_csv(R / "S1.4f_PH1_territory_nn.csv")
CLS = ["nn<5", "nn5-6", "nn6-8"]; XL = ["< 5", "5-6", "6-8"]
def roi_est(d):
    g = d.groupby(["roi_id", "cls"]).p3.mean().groupby("cls"); return g.mean().reindex(CLS), (2 * g.std(ddof=1) / np.sqrt(g.count())).reindex(CLS)
fig, axs = plt.subplots(1, 4, figsize=(7.4, 3.4), sharey=True)
for ax, A in zip(axs, ["A1", "A2", "A3", "A4"]):
    rf = RF[RF.archetype == A].set_index("cls").reindex(CLS)
    ax.fill_between(range(3), rf.lo, rf.hi, color="0.85", lw=0, label="nuclei veri (minimo-massimo dei 5 ROI)")
    ax.plot(range(3), rf.est, color="black", lw=1.5, marker="o", ms=3, label="nuclei veri (media)")
    for arm, src in [("ET", SH[SH.tag == "ET"]), ("V0", SH[SH.tag == "V0"]), ("V1", SH[SH.tag == "V1"]), ("V2", SH[SH.tag == "V2"])]:
        m, e = roi_est(src[src.archetype == A]); ax.errorbar(range(3), m, yerr=e, color=COL[arm], marker="o", ms=2.5, lw=1.2 if arm == "V2" else 0.9, capsize=1.5, label=LBL[arm])
    m, e = roi_est(F2[(F2.arm == "ER") & (F2.archetype == A)]); ax.errorbar(range(3), m, yerr=e, color=COL["ER"], marker="o", ms=2.5, lw=0.9, capsize=1.5, label=LBL["ER"])
    for src, ls, lab in [("real", "--", "territori veri"), ("sim", ":", "territori simulati")]:
        t = PH[(PH.source == src) & (PH.what == "territory") & (PH.archetype == A)].groupby(["roi_id", "cls"]).frac.mean().groupby("cls").mean().reindex(CLS)
        ax.plot(range(3), t, color="0.35", ls=ls, lw=1, label=f"asse del territorio, {lab}")
    ax.set_xticks(range(3), XL); ax.set_xlabel("distanza dal vicino (µm)", fontsize=7); ax.set_title(TIS[A], loc="left", fontsize=8)
axs[0].set_ylabel("di traverso al vicino (frazione)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, fontsize=6.5, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("Allineamento per distanza dal vicino: la taratura usa solo la classe < 5 µm", x=0.01, ha="left", fontsize=8.5)
fig.tight_layout(rect=(0, 0.2, 1, 1)); fig.savefig(F / "S1.4f_fig2_profile.png", dpi=200); plt.close(fig)
# ---- figura 3: stessi centri ----
q = CP[CP.id.isin(["CP-1", "S-1"])]
fig, ax = plt.subplots(figsize=(5.2, 2.6))
for k, A in enumerate(["A1", "A2", "A3", "A4"]):
    for j, (cid, con) in enumerate([("CP-1", "EP-K2"), ("CP-1", "ET-K2"), ("S-1", "EP-K2"), ("S-1", "ET-K2")]):
        r = q[(q.id == cid) & (q.archetype == A) & (q.contrast == con)].iloc[0]; arm = con.split("-")[0]
        ax.errorbar(k + (j - 1.5) * 0.15, r.est, yerr=2 * r.se, fmt="o" if cid == "CP-1" else "s", mfc=COL[arm] if cid == "CP-1" else "white", color=COL[arm], ms=4, capsize=2,
                    label=(f"{LBL[arm]} meno disco" + ("" if cid == "CP-1" else ", senza i centri doppi al bordo")) if k == 0 else None)
ax.axhline(0, color="0.5", lw=0.8); ax.set_xticks(range(4), [TIS[A] for A in ["A1", "A2", "A3", "A4"]], rotation=10)
ax.set_ylabel("differenza di contatto (coppie < 6 µm)"); ax.set_title("Sui punti di Luca: contatti tolti rispetto al disco (± 2 errori standard, 4 finestre)", loc="left", fontsize=7.5)
ax.legend(frameon=False, fontsize=6); fig.tight_layout(); fig.savefig(F / "S1.4f_fig3_same_centres.png", dpi=200); plt.close(fig)
print("figure scritte:", sorted(p.name for p in F.glob("*.png")))
