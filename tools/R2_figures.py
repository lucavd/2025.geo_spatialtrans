#!/usr/bin/env python
"""tools/R2_figures.py — figure del report R2 (sfondo bianco, un pannello per archetipo dove pertinente). Output: results/R2/figures/."""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from skimage.segmentation import find_boundaries
sys.path.insert(0, str(Path(__file__).resolve().parent)); import R2_common as C
R = C.ROOT / "results/R2"; F = R / "figures"; F.mkdir(exist_ok=True)
ARCS = ["A1", "A2", "A3", "A4", "A5", "A6"]
NAMES = {"A1": "A1 epitelio (intestino tenue, mouse)", "A2": "A2 tumore (CRC)", "A3": "A3 stroma (CRC)", "A4": "A4 linfoide (linfonodo)", "A5": "A5 corteccia (mouse)", "A6": "A6 miocardio (cuore)"}
MCOL = {"cellpose_rgb": "#1b7837", "stardist_he": "#762a83", "spaceranger": "#e08214", "cellpose_hed": "#999999"}
MLAB = {"cellpose_rgb": "Cellpose-SAM (RGB)", "stardist_he": "StarDist HE", "spaceranger": "Space Ranger 4 (StarDist custom)", "cellpose_hed": "Cellpose-SAM ematossilina (esplorativa)"}
summ = pd.read_csv(R / "R2_roi_summary.csv"); nat = summ[summ.scale == 1]
nuc = pd.read_parquet(R / "R2_nuclei_all.parquet"); nuc = nuc[(nuc.scale == 1) & nuc.keep]
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.size": 10})

# Fig 1: densita' per archetipo e metodo (punti = ROI, barra = mediana)
fig, ax = plt.subplots(figsize=(11, 5))
methods = ["cellpose_rgb", "stardist_he", "spaceranger"]
for i, A in enumerate(ARCS):
    for j, m in enumerate(methods):
        d = nat[(nat.archetype == A) & (nat.method == m)].density_per_mm2
        if len(d) == 0: continue
        x = i + (j - 1) * 0.25
        ax.scatter(np.full(len(d), x) + np.random.uniform(-0.05, 0.05, len(d)), d, color=MCOL[m], s=22, alpha=0.85, label=MLAB[m] if i == 0 else None)
        ax.hlines(d.median(), x - 0.1, x + 0.1, color=MCOL[m], lw=2.5)
ax.set_xticks(range(6)); ax.set_xticklabels(ARCS); ax.set_ylabel("nuclei / mm² di tessuto (ROI 1 mm²)"); ax.set_yscale("log")
ax.set_title("Densità nucleare per archetipo — punti = ROI, barra = mediana"); ax.legend(frameon=False, loc="upper right"); ax.grid(axis="y", alpha=0.3)
fig.tight_layout(); fig.savefig(F / "fig1_density.png", dpi=150); plt.close(fig)

# Fig 2: distribuzione area nucleare per archetipo (ECDF per metodo), un pannello per archetipo
fig, axs = plt.subplots(2, 3, figsize=(14, 8)); axs = axs.ravel()
for ax, A in zip(axs, ARCS):
    for m in methods:
        a = nuc[(nuc.archetype == A) & (nuc.method == m)].area_um2.values
        if len(a) == 0: continue
        a = np.sort(a); ax.plot(a, np.arange(1, len(a) + 1) / len(a), color=MCOL[m], lw=1.6, label=f"{MLAB[m]} (n={len(a):,}, med {np.median(a):.0f})")
    ax.set_xlim(0, 120); ax.set_title(NAMES[A]); ax.set_xlabel("area nucleare (µm²)"); ax.set_ylabel("ECDF"); ax.legend(fontsize=7.5, frameon=False, loc="lower right"); ax.grid(alpha=0.3)
fig.suptitle("Area nucleare — ECDF per archetipo e metodo (scala nativa, nuclei validi)"); fig.tight_layout(); fig.savefig(F / "fig2_area_ecdf.png", dpi=150); plt.close(fig)

# Fig 3: istogrammi area (log) per archetipo, Cellpose-RGB — per la bimodalita' (B-R2.4)
fig, axs = plt.subplots(2, 3, figsize=(14, 7)); axs = axs.ravel()
for ax, A in zip(axs, ARCS):
    a = nuc[(nuc.archetype == A) & (nuc.method == "cellpose_rgb")].area_um2.values
    ax.hist(np.log10(a), bins=80, color=MCOL["cellpose_rgb"], alpha=0.8); ax.set_title(NAMES[A]); ax.set_xlabel("log10 area nucleare (µm²)")
    for v in [10, 20, 50, 100]: ax.axvline(np.log10(v), color="k", lw=0.5, ls=":")
fig.suptitle("Distribuzione dell'area nucleare (Cellpose-SAM RGB) — linee: 10, 20, 50, 100 µm²"); fig.tight_layout(); fig.savefig(F / "fig3_area_hist.png", dpi=150); plt.close(fig)

# Fig 4: eccentricita' e NN distance per archetipo (Cellpose-RGB), punti = ROI
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, col, lab in zip(axs, ["ecc_median", "nn_median_um", "nuclear_area_fraction"], ["eccentricità nucleare mediana", "distanza al primo vicino mediana (µm)", "frazione di area di tessuto occupata dai nuclei"]):
    for i, A in enumerate(ARCS):
        d = nat[(nat.archetype == A) & (nat.method == "cellpose_rgb")][col]
        ax.scatter(np.full(len(d), i) + np.random.uniform(-0.08, 0.08, len(d)), d, color=MCOL["cellpose_rgb"], s=25); ax.hlines(d.median(), i - 0.2, i + 0.2, color="k", lw=2)
        d2 = nat[(nat.archetype == A) & (nat.method == "stardist_he")][col]
        ax.scatter(np.full(len(d2), i + 0.3), d2, color=MCOL["stardist_he"], s=14, alpha=0.7)
    ax.set_xticks(range(6)); ax.set_xticklabels(ARCS); ax.set_title(lab, fontsize=10); ax.grid(axis="y", alpha=0.3)
fig.suptitle("Forma, spaziatura e occupazione — verde Cellpose-RGB (barra = mediana), viola StarDist"); fig.tight_layout(); fig.savefig(F / "fig4_shape_nn.png", dpi=150); plt.close(fig)

# Fig 5: g(r) per archetipo con envelope CSR
pcf = pd.read_csv(R / "spatial/R2_pcf.csv")
fig, axs = plt.subplots(2, 3, figsize=(14, 8)); axs = axs.ravel()
for ax, A in zip(axs, ARCS):
    d = pcf[pcf.archetype == A]
    for roi, g in d.groupby("roi_id"):
        ax.fill_between(g.r_um, g.g_lo, g.g_hi, color="grey", alpha=0.15, lw=0)
        ax.plot(g.r_um, g.g_obs, lw=1.2, label=roi)
    ax.axhline(1, color="k", lw=0.8, ls="--"); ax.set_xlim(0, 60); ax.set_ylim(0, 2.5); ax.set_title(NAMES[A]); ax.set_xlabel("r (µm)"); ax.set_ylabel("g(r)"); ax.legend(fontsize=7, frameon=False, ncol=2); ax.grid(alpha=0.3)
fig.suptitle("Funzione di correlazione di coppia g(r) dei centroidi nucleari (Cellpose-RGB) — grigio: envelope CSR 95 % (39 sim.)"); fig.tight_layout(); fig.savefig(F / "fig5_pcf.png", dpi=150); plt.close(fig)

# Fig 6: concordanza fra metodi (F1@IoU0.5 e Δ densita') per archetipo
mt = pd.read_csv(R / "R2_matching.csv")
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
pairs = [("cellpose_rgb", "stardist_he", "Cellpose vs StarDist", "#1f78b4"), ("cellpose_rgb", "spaceranger", "Cellpose vs Space Ranger", "#e08214"), ("stardist_he", "spaceranger", "StarDist vs Space Ranger", "#762a83")]
for j, (a, b, lab, col) in enumerate(pairs):
    d = mt[(mt.method_a == a) & (mt.method_b == b)]
    for i, A in enumerate(ARCS):
        e = d[d.archetype == A]
        if len(e) == 0: continue
        x = i + (j - 1) * 0.25
        axs[0].scatter(np.full(len(e), x), e.f1_iou05, color=col, s=22, label=lab if i == 0 else None); axs[0].hlines(e.f1_iou05.median(), x - 0.1, x + 0.1, color=col, lw=2.5)
        axs[1].scatter(np.full(len(e), x), e.delta_density_rel * 100, color=col, s=22); axs[1].hlines(e.delta_density_rel.median() * 100, x - 0.1, x + 0.1, color=col, lw=2.5)
axs[0].axhline(0.7, color="green", ls=":", lw=1); axs[0].axhline(0.5, color="red", ls=":", lw=1); axs[0].set_ylim(0, 1); axs[0].set_ylabel("F1 a IoU ≥ 0.5"); axs[0].set_title("Concordanza oggetto per oggetto"); axs[0].legend(frameon=False, fontsize=8)
axs[1].axhline(0, color="k", lw=0.8); axs[1].axhspan(-10, 10, color="green", alpha=0.08); axs[1].axhspan(-25, 25, color="orange", alpha=0.06); axs[1].set_ylabel("Δ densità del secondo metodo vs primo (%)"); axs[1].set_title("Differenza di densità")
for ax in axs: ax.set_xticks(range(6)); ax.set_xticklabels(ARCS); ax.grid(axis="y", alpha=0.3)
fig.tight_layout(); fig.savefig(F / "fig6_matching.png", dpi=150); plt.close(fig)

# Fig 7: effetto scala (CP-R2.3)
sc = pd.read_csv(R / "R2_scale_effect.csv")
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
for ax, m in zip(axs, ["cellpose_rgb", "stardist_he"]):
    for s, col in [(2, "#4393c3"), (4, "#d6604d")]:
        d = sc[(sc.method == m) & (sc.scale == s)]
        for i, A in enumerate(ARCS):
            e = d[d.archetype == A].density_rel * 100 - 100
            ax.scatter(np.full(len(e), i + (s - 3) * 0.12), e, color=col, s=22, label=f"{s}× ({0.274*s:.2f} µm/px)" if i == 0 else None); ax.hlines(e.median(), i + (s - 3) * 0.12 - 0.1, i + (s - 3) * 0.12 + 0.1, color=col, lw=2.5)
    ax.axhline(0, color="k", lw=0.8); ax.axhspan(-10, 10, color="green", alpha=0.08); ax.set_xticks(range(6)); ax.set_xticklabels(ARCS); ax.set_ylabel("Δ densità vs nativa (%)"); ax.set_title(MLAB[m]); ax.legend(frameon=False, fontsize=8); ax.grid(axis="y", alpha=0.3)
fig.suptitle("Effetto della risoluzione (BL-015 / CP-R2.3): densità a 2× e 4× rispetto alla nativa 0.274 µm/px"); fig.tight_layout(); fig.savefig(F / "fig7_scale.png", dpi=150); plt.close(fig)

# Fig 8: galleria overlay 200 um per archetipo (ROI 1), Cellpose-RGB verde / StarDist viola / SR arancio
fig, axs = plt.subplots(2, 3, figsize=(15, 10)); axs = axs.ravel()
for ax, A in zip(axs, ARCS):
    roi = C.rois_table().query("archetype==@A").roi_id.iloc[0]; rgb = C.load_roi(A, roi); upp = C.rois_table().query("archetype==@A and roi_id==@roi").um_per_px.iloc[0]
    w = int(200 / upp); r0 = rgb.shape[0] // 2 - w // 2; c0 = rgb.shape[1] // 2 - w // 2; im = rgb[r0:r0 + w, c0:c0 + w].copy()
    for m, col in [("stardist_he", [118, 42, 131]), ("spaceranger", [224, 130, 20]), ("cellpose_rgb", [0, 200, 0])]:
        f = C.MASK_DIR / f"{A}_{roi}_{m}_native.npz"
        if f.exists():
            L = C.load_labels(A, roi, m, 1)[r0:r0 + w, c0:c0 + w]; im[find_boundaries(L, mode="inner")] = col
    ax.imshow(im); ax.set_title(f"{NAMES[A]} — {roi}, 200 µm"); ax.set_axis_off()
fig.suptitle("Overlay delle tre segmentazioni: verde Cellpose-RGB, viola StarDist, arancio Space Ranger (dove disponibile)"); fig.tight_layout(); fig.savefig(F / "fig8_overlays.png", dpi=130); plt.close(fig)

# Fig 9: maschere tessuto/bolle per A6 (BL-017) e un ROI A1 (lume)
fig, axs = plt.subplots(1, 4, figsize=(18, 4.8))
for k, (A, roi) in enumerate([("A6", "r1"), ("A1", "r2")]):
    rgb = C.load_roi(A, roi); upp = C.rois_table().query("archetype==@A and roi_id==@roi").um_per_px.iloc[0]
    tm, bm = C.masks_for(rgb, upp, 4, A, roi); rgb4 = C.downsample(rgb, 4)
    axs[2 * k].imshow(rgb4); axs[2 * k].set_title(f"{A}_{roi} H&E (1/4)"); over = rgb4.copy(); over[~tm] = (over[~tm] * 0.3).astype(np.uint8); over[bm] = [255, 0, 0]
    axs[2 * k + 1].imshow(over); axs[2 * k + 1].set_title(f"maschera: scuro = non tessuto ({100*(1-tm.mean()):.0f} %), rosso = bolle ({100*bm.mean():.1f} %)")
for ax in axs: ax.set_axis_off()
fig.tight_layout(); fig.savefig(F / "fig9_masks.png", dpi=120); plt.close(fig)
print("figures ok")
