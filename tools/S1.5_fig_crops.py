#!/usr/bin/env python
"""tools/S1.5_fig_crops.py - S1.5: ritaglio H&E della prima finestra di R2b di ogni tessuto (per la figura di esempio del report).
Uscita: /mnt/micron/geo_spatialtrans/S1.5/fig/<A>_<roi>_<win>.png e crops.csv (finestra in µm nel sistema del ROI)."""
import os, numpy as np, pandas as pd, tifffile
from PIL import Image
OUT = "/mnt/micron/geo_spatialtrans/S1.5/fig"; os.makedirs(OUT, exist_ok=True)
w = pd.read_csv("results/R2b/R2b_windows.csv"); rows = []
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    r = w[w.archetype == A].iloc[0]
    img = tifffile.imread(f"/mnt/micron/geo_spatialtrans/R2/rois/{A}_{r.roi_id}.tif")[..., :3]
    upp = float(r.um_per_px); c0 = int(round(r.x0_um / upp)); r0 = int(round(r.y0_um / upp)); n = int(r.side_px)
    crop = img[r0:r0 + n, c0:c0 + n]
    f = f"{OUT}/{A}_{r.roi_id}_{r.win_id}.png"; Image.fromarray(crop).save(f)
    rows.append(dict(archetype=A, roi_id=r.roi_id, win_id=r.win_id, x0_um=r.x0_um, y0_um=r.y0_um, side_um=n * upp, png=f))
pd.DataFrame(rows).to_csv(f"{OUT}/crops.csv", index=False); print(pd.DataFrame(rows))
