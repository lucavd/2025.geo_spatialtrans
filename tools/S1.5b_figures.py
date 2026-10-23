"""tools/S1.5b_figures.py - S1.5b: descrittiva D-desc-1 (distanza al vicino piu' vicino per tipo, reale contro simulato) e figure del report.
Uso: .venv/bin/python tools/S1.5b_figures.py <braccio vincente>   Uscite: results/S1.5b/S1.5b_desc_nn.csv, results/S1.5b/figures/*.png"""
import sys, os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import importlib.util
spec = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(spec); spec.loader.exec_module(cal)
WIN = sys.argv[1]; RES = "results/S1.5b"; FIG = f"{RES}/figures"; os.makedirs(FIG, exist_ok=True)
TIS = {"A1": "intestino di topo", "A2": "tumore del colon", "A3": "stroma del tumore", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}
ARMLAB = {"group": "(a) modello del gruppo", "recentred": "(b) ricentrato", "single": "un tipo (S1.5), etichette a caso"}
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 9})
S = pd.read_parquet("/mnt/micron/geo_spatialtrans/S1.5b/cells_seed1.parquet")
V = pd.read_csv(f"{RES}/S1.5b_verdicts.csv"); T = pd.read_csv(f"{RES}/S1.5b_panel_types.csv")
RT = pd.read_csv(f"{RES}/S1.5b_real_refs_roi_type.csv"); RT = RT[RT.level == "type"]
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"], "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"],
        "A4": ["linfociti B", "linfociti T"], "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
# ---- D-desc-1: distanza al vicino piu' vicino (qualsiasi tipo), per tipo; celle interne; reale = centri dei nuclei veri (insieme di valutazione)
rows = []
for A in cal.ARCH:
    lab, d = cal.load(A); ev = d[d.interior]
    for r, e in d.groupby("roi_id"):
        tr = cKDTree(e[["x", "y"]].values); dd, _ = tr.query(e[["x", "y"]].values, k=2); e = e.assign(nn=dd[:, 1])
        for t in MAIN[A]:
            x = e[e.interior & (e.type == t)].nn
            if len(x): rows.append(dict(archetype=A, roi_id=r, source="reale", cell_type=t, nn_median=x.median(), nn_p05=x.quantile(0.05), n=len(x)))
    for arm in [WIN, "single"]:
        for r, e in S[(S.archetype == A) & (S.arm == arm)].groupby("roi_id"):
            tr = cKDTree(e[["x", "y"]].values); dd, _ = tr.query(e[["x", "y"]].values, k=2); e = e.assign(nn=dd[:, 1])
            for t in MAIN[A]:
                x = e[e.interior & (e.cell_type == t)].nn
                if len(x): rows.append(dict(archetype=A, roi_id=r, source=arm, cell_type=t, nn_median=x.median(), nn_p05=x.quantile(0.05), n=len(x)))
NN = pd.DataFrame(rows); NN.to_csv(f"{RES}/S1.5b_desc_nn.csv", index=False)
# ---- figura 1: area mediana del nucleo per tipo, simulato (ROI) contro intervallo reale
fig, axs = plt.subplots(2, 3, figsize=(11, 6.5))
for ax, A in zip(axs.ravel(), cal.ARCH):
    for k, t in enumerate(MAIN[A]):
        rt = RT[(RT.archetype == A) & (RT.unit == t)]
        ax.plot([k, k], [rt.median_area_nuc.min(), rt.median_area_nuc.max()], color="0.7", lw=8, solid_capstyle="butt", zorder=1)
        for j, (arm, col, mk) in enumerate([("group", "#1f77b4", "o"), ("recentred", "#d62728", "s"), ("single", "0.3", "x")]):
            v = V[(V.id == "P-2") & (V.arm == arm) & (V.archetype == A) & (V.unit == t)].iloc[0]
            ax.errorbar(k + (j - 1) * 0.18, v.est, yerr=2 * v.se, fmt=mk, color=col, ms=5, capsize=2, label=ARMLAB[arm] if (A == "A1" and k == 0) else None, zorder=3)
    ax.set_xticks(range(len(MAIN[A]))); ax.set_xticklabels([t.replace("stromale e muscolo liscio", "stromale") for t in MAIN[A]], rotation=20, ha="right")
    ax.set_title(TIS[A]); ax.set_ylabel("area mediana del nucleo (µm²)")
axs[0, 0].legend(fontsize=7, loc="upper right", frameon=False)
fig.suptitle("Area mediana del nucleo per tipo: simulato (media dei 5 ROI ± 2 ES) e intervallo fra ROI reali (barra grigia)")
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5b_fig1_area_by_type.png", dpi=150); plt.close(fig)
# ---- figura 2: cervello, distribuzione dell'area del nucleo (celle interne, seme 1), reale contro simulato
lab, d = cal.load("A5"); ev = d[d.interior & d.matched & (d.a_n > 0)]
thr = float(pd.read_csv(f"{RES}/S1.5b_real_refs_roi.csv").small_threshold_um2.dropna().iloc[0])
fig, ax = plt.subplots(figsize=(6.5, 4))
bins = np.logspace(np.log10(3), np.log10(400), 50)
ax.hist(ev.a_n, bins=bins, density=True, histtype="stepfilled", color="0.8", label="reale (5 ROI)")
for arm, col in [(WIN, "#d62728"), ("single", "0.2")]:
    s = S[(S.archetype == "A5") & (S.arm == arm) & S.interior]
    ax.hist(s.a_n, bins=bins, density=True, histtype="step", lw=1.5, color=col, label=ARMLAB[arm] + " (seme 1)")
ax.axvspan(thr / np.sqrt(2), thr * np.sqrt(2), color="#ffcc66", alpha=0.3, label="banda fra le due taglie")
ax.set_xscale("log"); ax.set_xlabel("area del nucleo (µm², scala logaritmica)"); ax.set_ylabel("densità"); ax.legend(fontsize=7, frameon=False)
ax.set_title("Cervello di topo: area del nucleo, celle interne")
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5b_fig2_brain_area.png", dpi=150); plt.close(fig)
# ---- figura 3: rapporto nucleo/territorio per tipo
fig, axs = plt.subplots(2, 3, figsize=(11, 6.5))
for ax, A in zip(axs.ravel(), cal.ARCH):
    for k, t in enumerate(MAIN[A]):
        rt = RT[(RT.archetype == A) & (RT.unit == t)]
        ax.plot([k, k], [rt.median_nt.min(), rt.median_nt.max()], color="0.7", lw=8, solid_capstyle="butt", zorder=1)
        v = V[(V.id == "P-7") & (V.arm == WIN) & (V.archetype == A) & (V.unit == t)].iloc[0]
        ax.errorbar(k, v.est, yerr=2 * v.se, fmt="s", color="#d62728" if WIN == "recentred" else "#1f77b4", ms=5, capsize=2, zorder=3)
    ax.set_xticks(range(len(MAIN[A]))); ax.set_xticklabels([t.replace("stromale e muscolo liscio", "stromale") for t in MAIN[A]], rotation=20, ha="right")
    ax.set_title(TIS[A]); ax.set_ylabel("rapporto nucleo/territorio")
fig.suptitle(f"Rapporto mediano nucleo/territorio per tipo: simulato {ARMLAB[WIN]} contro intervallo fra ROI reali")
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5b_fig3_nt_by_type.png", dpi=150); plt.close(fig)
print(NN.groupby(["archetype", "cell_type", "source"]).nn_median.mean().unstack("source").round(2).to_string())

# ---- D-desc-2: tempo di place_cells() per unita' (s), bracci per tipo contro il braccio a un tipo
U = pd.read_csv(f"{RES}/S1.5b_panel_units.csv")
tm = U.groupby(["archetype", "arm"]).elapsed_s.median().unstack("arm").round(1); tm.to_csv(f"{RES}/S1.5b_desc_time.csv"); print(tm.to_string())
# ---- D-desc-3: un ROI per tessuto, H&E accanto ai centri simulati colorati per tipo (braccio vincente, seme 1)
import tifffile
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
fig, axs = plt.subplots(6, 2, figsize=(8, 23))
for k, A in enumerate(cal.ARCH):
    r = rois[rois.archetype == A].iloc[0]; img = tifffile.imread(f"/mnt/micron/geo_spatialtrans/R2/rois/{A}_{r.roi_id}.tif")
    if img.ndim == 3 and img.shape[0] in (3, 4): img = np.moveaxis(img, 0, -1)
    step = max(1, img.shape[0] // 900); im = img[::step, ::step, :3]; side = r.side_px * r.um_per_px
    axs[k, 0].imshow(im, extent=[0, side, side, 0]); axs[k, 0].set_title(f"{TIS[A]} ({r.roi_id}): H&E"); axs[k, 0].set_xlabel("µm")
    s = S[(S.archetype == A) & (S.roi_id == r.roi_id) & (S.arm == WIN)]
    for t, g in s.groupby("cell_type"): axs[k, 1].scatter(g.x, g.y, s=1.5 if A in ("A1", "A4") else 4, label=f"{t} ({len(g)})")
    axs[k, 1].set_xlim(0, side); axs[k, 1].set_ylim(side, 0); axs[k, 1].set_aspect("equal"); axs[k, 1].set_title("centri simulati per tipo (seme 1)")
    axs[k, 1].legend(fontsize=6, markerscale=3, loc="upper right", frameon=True)
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5b_fig4_types_vs_he.png", dpi=110); plt.close(fig)
