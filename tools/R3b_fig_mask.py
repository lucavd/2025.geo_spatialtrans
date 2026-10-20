#!/usr/bin/env python
"""tools/R3b_fig_mask.py — R3b: figura di controllo visivo della maschera CO su una finestra di linfonodo e una di cuore:
H&E, contorni dei nuclei di forma StarDist (verde) e dei poligoni Space Ranger rimasti (magenta), punti di Luca (giallo).
Uscita: reports/figures/R3b_mask_windows.png.  Uso (root del repo): .venv/bin/python tools/R3b_fig_mask.py"""
import sys
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from skimage.segmentation import find_boundaries
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "tools")); import R2_common as C
M = Path("/mnt/micron/geo_spatialtrans/R3b"); OUT = ROOT / "reports/figures"; OUT.mkdir(parents=True, exist_ok=True)
wins = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv").set_index("win_id")
pts = pd.read_csv(ROOT / "results/R2b/R2b_points.csv"); pts = pts[pts.rater == "Luca"]
fig, axs = plt.subplots(1, 2, figsize=(12, 6.2), facecolor="white")
for ax, wid, title in zip(axs, ["A4_w2", "A6_w2"], ["Linfonodo (A4_w2, 60 µm)", "Cuore (A6_w2, 250 µm)"]):
    w = wins.loc[wid]; A, roi = w.archetype, w.roi_id; c0, r0, s = int(w.c0_roi), int(w.r0_roi), int(w.side_px)
    rgb = C.load_roi(A, roi)[r0:r0 + s, c0:c0 + s]
    lab = np.load(M / "masks" / f"{A}_{roi}_consensus_native.npz")["labels"][r0:r0 + s, c0:c0 + s].astype(np.int64)
    mp = pd.read_parquet(M / "masks" / f"{A}_{roi}_consensus_map.parquet")
    is_sr = np.zeros(int(max(lab.max(), mp.label.max())) + 1, bool); is_sr[mp.label[mp.src == "sr"].to_numpy()] = True
    b = find_boundaries(lab, mode="inner")
    img = rgb.astype(float) / 255
    img[b & ~is_sr[lab]] = [0.0, 0.75, 0.2]; img[b & is_sr[lab]] = [0.9, 0.0, 0.8]
    ax.imshow(img); p = pts[pts.win_id == wid]
    ax.scatter(p.x_win_px, p.y_win_px, s=9, c="yellow", edgecolors="black", linewidths=0.3)
    n_sr = len(np.intersect1d(np.unique(lab[lab > 0]), mp.label[mp.src == "sr"].to_numpy())); n_all = len(np.unique(lab[lab > 0]))
    ax.set_title(f"{title}\nforma StarDist (verde) {n_all - n_sr}, poligono Space Ranger (magenta) {n_sr}; punti di Luca {len(p)}", fontsize=10)
    ax.set_axis_off()
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(OUT / "R3b_mask_windows.png", dpi=150, facecolor="white"); print("ok")
