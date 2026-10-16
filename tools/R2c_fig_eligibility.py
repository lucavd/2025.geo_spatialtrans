#!/usr/bin/env python
"""tools/R2c_fig_eligibility.py — R2c, C-5: una finestra per archetipo (w1) con punti di Luca, focali, coppie da giudicare,
coppie di fascia campionate e cerchio d_max attorno a un focale. Output results/R2c/fig/R2c_eligibility.png (+ .pdf)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/R2c"; R2B = ROOT / "results/R2b"
(RES / "fig").mkdir(exist_ok=True)
wins = pd.read_csv(R2B / "R2b_windows.csv").set_index("win_id")
pts = pd.read_csv(R2B / "R2b_points.csv"); pts = pts[pts.rater == "Luca"]
pr = pd.read_csv(RES / "R2c_pairs.csv"); foc = pd.read_csv(RES / "R2c_focal.csv"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype")
fig, axs = plt.subplots(2, 3, figsize=(15, 10.4), facecolor="white")
for ax, A in zip(axs.ravel(), ["A1", "A2", "A3", "A4", "A5", "A6"]):
    wid = f"{A}_w1"; w = wins.loc[wid]; upp = float(w.um_per_px)
    im = np.asarray(Image.open(w.png)); ax.imshow(im, extent=[0, im.shape[1] * upp, im.shape[0] * upp, 0])
    g = pts[pts.win_id == wid].set_index("idx"); X = g.x_um - w.x0_um; Y = g.y_um - w.y0_um
    f = foc[(foc.win_id == wid)]; fi = set(f.idx[f.focal]); ei = set(f.idx[f.eligible])
    for q, col, lw in [(pr[(pr.win_id == wid) & pr.judged], "#00bcd4", 1.0), (pr[(pr.win_id == wid) & pr.band_sampled], "#e91e63", 1.6)]:
        for r in q.itertuples():
            ax.plot([X[r.i], X[r.j]], [Y[r.i], Y[r.j]], color=col, lw=lw, alpha=0.9)
    ax.scatter(X, Y, s=6, c="#424242", zorder=3)
    ax.scatter(X[list(ei - fi)], Y[list(ei - fi)], s=10, facecolors="none", edgecolors="#9e9e9e", lw=0.6, zorder=3)
    ax.scatter(X[list(fi)], Y[list(fi)], s=28, c="#ffc107", edgecolors="k", lw=0.5, zorder=4)
    d = float(dm.d_max_um[A]); s = int(w.side_px) * upp
    ax.add_patch(plt.Rectangle((d, d), s - 2 * d, s - 2 * d, fill=False, ls="--", ec="#ffc107", lw=1))
    k0 = sorted(fi)[0]; ax.add_patch(plt.Circle((X[k0], Y[k0]), d, fill=False, ec="#ffc107", lw=0.8, ls=":"))
    nj = int(((pr.win_id == wid) & pr.judged).sum()); nbs = int(((pr.win_id == wid) & pr.band_sampled).sum())
    ax.set_title(f"{wid}  d_max {d:.1f} µm  focali {len(fi)}  coppie {nj}  fascia {nbs}", fontsize=10)
    ax.set_xlim(0, s); ax.set_ylim(s, 0); ax.set_xlabel("µm"); ax.set_ylabel("µm")
fig.suptitle("R2c C-5: punti di Luca (grigio), eleggibili (anello), focali (giallo), coppie da giudicare (ciano), fascia campionata (magenta); "
             "tratteggio = zona eleggibile, punteggiato = d_max", fontsize=10)
fig.tight_layout(); fig.savefig(RES / "fig/R2c_eligibility.png", dpi=130, facecolor="white"); fig.savefig(RES / "fig/R2c_eligibility.pdf", facecolor="white")
print("ok")
