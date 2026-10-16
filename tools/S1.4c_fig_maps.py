#!/usr/bin/env python
"""tools/S1.4c_fig_maps.py — S1.4c, figura 1: riquadro di 40 x 40 µm al centro del ROI di test del fold 1 (seme 1):
maschera reale (segmentatore di forma, solo etichette-generatore) contro K-c, K-1, K-2 rasterizzati con l'operatore di
S1.4c; pixel di bordo a contatto (Chebyshev <= 2 px) in nero. Uso: .venv/bin/python tools/S1.4c_fig_maps.py"""
import importlib.util, numpy as np, pandas as pd, shapely, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
spec = importlib.util.spec_from_file_location("ras", "tools/S1.4c_raster.py"); ras = importlib.util.module_from_spec(spec); spec.loader.exec_module(ras)
OUT = "/mnt/micron/geo_spatialtrans/S1.4c"; SIDE = 40.0
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
fold1 = {"A1": "r1", "A2": "r1", "A3": "r1", "A4": "f1", "A5": "r1", "A6": "r1"}
import glob, os
def roi_of(A):
    f = sorted(glob.glob(f"{OUT}/gen_in/{A}_*_fold1_s01_Kc.parquet"))[0]; return os.path.basename(f).split("_")[1]
fig, ax = plt.subplots(6, 4, figsize=(11, 16.5), facecolor="white")
rng = np.random.default_rng(7)
for i, A in enumerate(["A1", "A2", "A3", "A4", "A5", "A6"]):
    roi = roi_of(A); rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]; upp = float(rw.um_per_px); side = int(rw.side_px)
    c0 = r0 = int(side / 2 - SIDE / upp / 2); n = int(SIDE / upp)
    seg = ras.SHAPE_SEG[A]
    gen = pd.read_parquet(f"{ras.R3}/real/{A}_{roi}_{seg}_gen.parquet")
    lab = np.load(f"{ras.MASKS}/{A}_{roi}_{seg}_native.npz")["labels"].astype(np.int64)
    keep = np.zeros(int(lab.max()) + 1, bool); keep[gen["label"].to_numpy().astype(np.int64)] = True; lab[~keep[lab]] = 0
    imgs = [("reale (" + seg.split("_")[0] + ")", lab[r0:r0 + n, c0:c0 + n])]
    for V in ["Kc", "K1", "K2"]:
        t = pd.read_parquet(f"{OUT}/gen_in/{A}_{roi}_fold1_s01_{V}.parquet"); g = shapely.from_wkb(t["wkb"].to_numpy())
        box = shapely.box(c0 * upp - 20, r0 * upp - 20, (c0 + n) * upp + 20, (r0 + n) * upp + 20); sel = shapely.intersects(g, box)
        L, _ = ras.rasterize(g[sel], t["idx"].to_numpy()[sel], (side, side), upp); imgs.append(({"Kc": "K-c (scala)", "K1": "K-1 (ritaglio)", "K2": "K-2 (ritaglio + pari area)"}[V], L[r0:r0 + n, c0:c0 + n]))
    for j, (ttl, L) in enumerate(imgs):
        u, inv = np.unique(L, return_inverse=True); cols = rng.uniform(0.35, 0.95, (len(u), 3)); cols[u == 0] = 1.0
        rgb = cols[inv.reshape(L.shape)]
        cont = ras.boundary(L) & ras.other_within(L, 2); rgb[cont] = 0.0
        ax[i, j].imshow(rgb, extent=(0, SIDE, SIDE, 0), interpolation="nearest"); ax[i, j].set_xticks([]); ax[i, j].set_yticks([])
        ax[i, j].set_title(f"{A} {roi} — {ttl}" if j == 0 else ttl, fontsize=9)
fig.suptitle("Nuclei in un riquadro di 40 × 40 µm (fold 1, seme 1); nero = bordo a contatto (≤ 2 px)", fontsize=11)
fig.tight_layout(); fig.savefig("results/S1.4c/fig/S1.4c_fig1_maps.png", dpi=130, facecolor="white"); print("ok")
