#!/usr/bin/env python
"""tools/R2_exclusive_check.py — arbitro per il disaccordo fra segmentatori (scala nativa).
Per ogni coppia (A,B) e ROI: oggetti di A appaiati a B (IoU>=0.5) vs esclusivi di A. Per ciascun oggetto: OD ematossilina media
dentro e in un anello esterno di 2 um; 'evidenza nucleare' = OD_in > OD_ring * 1.3 e OD_in > OD mediana del tessuto.
Se gli esclusivi hanno evidenza nucleare come gli appaiati -> nuclei veri persi da B; se no -> falsi positivi di A.
Output: results/R2/R2_exclusive_check.csv
"""
import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
from skimage import color, morphology
sys.path.insert(0, str(Path(__file__).resolve().parent)); import R2_common as C
rois = C.rois_table()
if len(sys.argv) > 1: rois = rois[(rois.archetype + '_' + rois.roi_id).isin(sys.argv[1].split(','))]
import time
pairs = [("cellpose_rgb", "stardist_he"), ("stardist_he", "cellpose_rgb"), ("cellpose_rgb", "spaceranger"), ("spaceranger", "cellpose_rgb"), ("stardist_he", "spaceranger"), ("spaceranger", "stardist_he")]
rows = []
def matched_ids(la, lb):
    a = la.ravel().astype(np.int64); b = lb.ravel().astype(np.int64); m = (a > 0) & (b > 0)
    inter = sparse.coo_matrix((np.ones(m.sum(), np.int64), (a[m], b[m])), shape=(la.max() + 1, lb.max() + 1)).tocsr().tocoo()   # tocsr somma i duplicati
    area_a = np.bincount(la.ravel(), minlength=la.max() + 1); area_b = np.bincount(lb.ravel(), minlength=lb.max() + 1)
    iou = inter.data / (area_a[inter.row] + area_b[inter.col] - inter.data)
    return set(inter.row[iou >= 0.5].tolist())
for _, r in rois.iterrows():
    A, roi, upp = r.archetype, r.roi_id, r.um_per_px
    rgb = C.load_roi(A, roi); hod = color.rgb2hed(rgb)[..., 0]
    tm, bm = C.masks_for(rgb, upp, 1, A, roi); valid = tm & ~bm; od_tissue = np.median(hod[valid & (hod > 0)])   # mediana dei pixel con OD > 0 (rgb2hed tronca a 0: in A5/A6 meta' dei pixel e' 0; rilievo revisore R2)
    labs = {}
    for m in ["cellpose_rgb", "stardist_he", "spaceranger"]:
        f = C.MASK_DIR / f"{A}_{roi}_{m}_native.npz"
        if f.exists(): labs[m] = np.where(valid, C.load_labels(A, roi, m, 1), 0)
    ring_cache = {}
    for ma, mb in pairs:
        if ma not in labs or mb not in labs: continue
        t0 = time.time()
        la, lb = labs[ma], labs[mb]; n = la.max()
        if n == 0: continue
        if ma not in ring_cache:
            from scipy import ndimage as ndi
            dist, (iy, ix) = ndi.distance_transform_edt(la == 0, return_indices=True)     # distanza dal nucleo piu' vicino e suo indice
            ring = np.where((la == 0) & (dist <= 2 / upp), la[iy, ix], 0)                  # anello di 2 um attribuito al nucleo piu' vicino
            cnt = np.maximum(np.bincount(la.ravel(), minlength=n + 1)[1:], 1); rc = np.maximum(np.bincount(ring.ravel(), minlength=n + 1)[1:], 1)
            od_in = np.bincount(la.ravel(), weights=hod.ravel(), minlength=n + 1)[1:] / cnt
            od_rg = np.bincount(ring.ravel(), weights=hod.ravel(), minlength=n + 1)[1:] / rc
            area = cnt * upp ** 2
            ring_cache[ma] = (od_in, od_rg, area)
        od_in, od_rg, area = ring_cache[ma]
        ids = np.arange(1, n + 1); mset = matched_ids(la, lb); mt = np.isin(ids, np.fromiter(mset, dtype=np.int64, count=len(mset)))
        present = np.bincount(la.ravel(), minlength=n + 1)[1:] > 0          # label con area > 0 nell'area valida (niente oggetti fantasma)
        ev = (od_in > od_rg * 1.3) & (od_in > od_tissue)
        for lab, sel in [("matched", mt & present), ("exclusive", (~mt) & present)]:
            if sel.sum() == 0: continue
            print(f'   {ma} vs {mb} {lab} {time.time()-t0:.1f}s', flush=True) if lab=='exclusive' else None
            rows.append(dict(archetype=A, roi_id=roi, method=ma, other=mb, cls=lab, n=int(sel.sum()), frac_of_method=float(sel.mean()),
                             od_in_median=float(np.median(od_in[sel])), od_ring_median=float(np.median(od_rg[sel])), od_tissue_median=float(od_tissue),
                             frac_nuclear_evidence=float(ev[sel].mean()), area_median=float(np.median(area[sel]))))
    print(A, roi, "ok", flush=True)
d = pd.DataFrame(rows); d.to_csv(C.ROOT / "results/R2/R2_exclusive_check.csv", index=False)
s = d.groupby(["method", "other", "cls", "archetype"]).agg(n=("n", "sum"), frac_of_method=("frac_of_method", "median"), od_in=("od_in_median", "median"), od_ring=("od_ring_median", "median"), nuclear_evidence=("frac_nuclear_evidence", "median")).reset_index()
s.to_csv(C.ROOT / "results/R2/R2_exclusive_summary.csv", index=False); print(s.round(3).to_string())
