# tools/R3_export_tissue.py — esporta tissue_ds4 (maschera di tessuto R2, senza sottrazione delle bolle) per P4.
import numpy as np, glob, os
from skimage import io
out = "/mnt/micron/geo_spatialtrans/R3/posthoc/masks_tissue"; os.makedirs(out, exist_ok=True)
for f in sorted(glob.glob("/mnt/micron/geo_spatialtrans/R2/masks/*_masks.npz")):
    z = np.load(f); b = os.path.basename(f).replace("_masks.npz", "")
    io.imsave(f"{out}/{b}_tissue_ds4.png", (z["tissue_ds4"] * 255).astype(np.uint8), check_contrast=False)
print(len(glob.glob(out + "/*.png")))
