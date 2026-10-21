#!/usr/bin/env python
"""tools/S1.5_cp1_resample.py - S1.5, controprova CP-1: le 30 immagini H&E dei ROI di R2 ricampionate per media d'area
(PIL, filtro BOX) a 1, 2, 5, 10 µm/px, partendo dalla risoluzione nativa (um_per_px di results/R2/R2_rois_checked.csv).
Per ogni ROI e scala scrive /mnt/micron/geo_spatialtrans/S1.5/cp1/img/<A>_<roi>_<s>um.parquet con x, y (colonna e riga,
da 1), value (grigio = media dei tre canali / 255, come prepare_image()) e valid (maschera del tessuto valido di R2,
ricampionata con lo stesso filtro e soglia 0.5). Il .npz vuoto accanto segna il file come fatto.
Uso (root del repo): .venv/bin/python tools/S1.5_cp1_resample.py
"""
import os, numpy as np, pandas as pd, tifffile
from PIL import Image
OUT = "/mnt/micron/geo_spatialtrans/S1.5/cp1/img"; os.makedirs(OUT, exist_ok=True)
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
for r in rois.itertuples():
    img = tifffile.imread(f"/mnt/micron/geo_spatialtrans/R2/rois/{r.archetype}_{r.roi_id}.tif")
    if img.ndim == 3 and img.shape[0] in (3, 4) and img.shape[-1] not in (3, 4): img = np.moveaxis(img, 0, -1)
    img = img[..., :3]
    assert img.dtype == np.uint8 and img.shape[0] == img.shape[1] == r.side_px, (r.archetype, r.roi_id, img.shape, img.dtype)
    m = np.array(Image.open(f"results/R2/tissue_masks/{r.archetype}_{r.roi_id}_valid_ds4.png").convert("L")) > 127
    side_um = r.side_px * r.um_per_px
    for s in (1, 2, 5, 10):
        b = f"{r.archetype}_{r.roi_id}_{s}um"
        if os.path.exists(os.path.join(OUT, b + ".npz")): continue
        n = int(round(side_um / s))
        rgb = np.asarray(Image.fromarray(img).resize((n, n), Image.BOX), dtype=np.float64)
        mv = np.asarray(Image.fromarray((m * 255).astype(np.uint8)).resize((n, n), Image.BOX), dtype=np.float64) / 255 > 0.5
        gray = rgb.mean(axis=2) / 255.0
        yy, xx = np.mgrid[1:n + 1, 1:n + 1]
        pd.DataFrame(dict(x=xx.ravel().astype(np.int32), y=yy.ravel().astype(np.int32), value=gray.ravel(), valid=mv.ravel())).to_parquet(os.path.join(OUT, b + ".parquet"), index=False)
        np.savez(os.path.join(OUT, b + ".npz"), side_um=side_um, n=n)
    print(r.archetype, r.roi_id, img.shape, flush=True)
