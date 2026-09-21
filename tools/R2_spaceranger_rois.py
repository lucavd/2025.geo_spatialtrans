#!/usr/bin/env python
"""tools/R2_spaceranger_rois.py — terza segmentazione: poligoni nucleari di Space Ranger (nucleus_segmentations.geojson,
coordinate pixel full-res x=col,y=row) rasterizzati nei ROI a risoluzione nativa. Metodo 'spaceranger', scala 1.
C-R2.5: per ogni ROI, frazione di poligoni SR che si sovrappongono (>0 px) a un nucleo Cellpose-RGB, e OD ematossilina
media dentro i poligoni vs anello esterno (3 um): se l'orientamento fosse sbagliato (BL-020) entrambe crollerebbero.
"""
import sys, json, numpy as np, pandas as pd
from pathlib import Path
from skimage.draw import polygon as draw_polygon
from skimage import color, morphology
sys.path.insert(0, str(Path(__file__).resolve().parent)); import R2_common as C
DATA = Path("/mnt/micron/geo_spatialtrans/data_real")
rois = C.rois_table(); rois = rois[rois.inside_image]
rows = []
for ds, g in rois.groupby("dataset", sort=False):
    p = DATA / ds / "segmented_outputs_x/segmented_outputs/nucleus_segmentations.geojson"
    if not p.exists():
        print(ds, "nessuna segmentazione Space Ranger (BL-019)", flush=True); continue
    feats = json.load(open(p))["features"]
    # bbox per feature per filtrare velocemente
    polys = [np.asarray(f["geometry"]["coordinates"][0], dtype=np.float64) for f in feats]
    bb = np.array([[q[:, 0].min(), q[:, 0].max(), q[:, 1].min(), q[:, 1].max()] for q in polys])
    print(ds, "poligoni SR:", len(polys), flush=True)
    for _, r in g.iterrows():
        sel = np.where((bb[:, 0] >= r.c0) & (bb[:, 1] < r.c1) & (bb[:, 2] >= r.r0) & (bb[:, 3] < r.r1))[0]   # interamente nel ROI
        lab = np.zeros((r.side_px, r.side_px), np.uint32)
        for k, i in enumerate(sel, start=1):
            q = polys[i]; rr, cc = draw_polygon(q[:, 1] - r.r0, q[:, 0] - r.c0, shape=lab.shape)
            lab[rr, cc] = k
        rgb = C.load_roi(r.archetype, r.roi_id); upp = r.um_per_px
        tm, bm = C.masks_for(rgb, upp, 1, r.archetype, r.roi_id)
        tab = C.nuclei_table(lab, upp, tm, bm, r.archetype, r.roi_id, "spaceranger", 1)
        C.save_labels(lab, r.archetype, r.roi_id, "spaceranger", 1)
        tab.to_parquet(C.NUC_DIR / f"{r.archetype}_{r.roi_id}_spaceranger_native.parquet", index=False)
        # C-R2.5
        cp = C.load_labels(r.archetype, r.roi_id, "cellpose_rgb", 1)
        n = lab.max()
        overlap = np.bincount(lab[(lab > 0) & (cp > 0)], minlength=n + 1)[1:] > 0
        hod = color.rgb2hed(rgb)[..., 0]
        inside = np.bincount(lab.ravel(), weights=hod.ravel(), minlength=n + 1)[1:] / np.maximum(np.bincount(lab.ravel(), minlength=n + 1)[1:], 1)
        ring = morphology.dilation(lab, morphology.disk(int(round(3 / upp)))); ring[lab > 0] = 0
        rin = np.bincount(ring.ravel(), weights=hod.ravel(), minlength=n + 1)[1:] / np.maximum(np.bincount(ring.ravel(), minlength=n + 1)[1:], 1)
        rows.append(dict(archetype=r.archetype, roi_id=r.roi_id, n_sr=int(n), n_keep=int(tab.keep.sum()) if len(tab) else 0,
                         frac_overlap_cellpose=float(overlap.mean()) if n else np.nan, frac_od_inside_gt_ring=float((inside > rin).mean()) if n else np.nan,
                         od_inside_median=float(np.median(inside)) if n else np.nan, od_ring_median=float(np.median(rin)) if n else np.nan))
        print(f"  {r.archetype}_{r.roi_id}: SR n={n} overlap_cp={overlap.mean():.3f} OD_in>ring={(inside>rin).mean():.3f}", flush=True)
pd.DataFrame(rows).to_csv(C.ROOT / "results/R2/R2_spaceranger_check.csv", index=False)
