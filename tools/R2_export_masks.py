
# esporta maschere valide (tessuto & ~bolle) a 1/4 come PNG per R
import sys, numpy as np; sys.path.insert(0, "tools"); import R2_common as C
from skimage.io import imsave
out = C.ROOT / "results/R2/tissue_masks"; out.mkdir(exist_ok=True)
for _, r in C.rois_table().iterrows():
    rgb = C.load_roi(r.archetype, r.roi_id); tm, bm = C.masks_for(rgb, r.um_per_px, 4, r.archetype, r.roi_id)
    imsave(out / f"{r.archetype}_{r.roi_id}_valid_ds4.png", ((tm & ~bm) * 255).astype(np.uint8), check_contrast=False)
print("masks ok")
