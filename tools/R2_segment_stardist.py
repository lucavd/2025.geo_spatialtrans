#!/usr/bin/env python
"""tools/R2_segment_stardist.py — segmentazione nucleare con StarDist 2D_versatile_he sui ROI R2 (eseguire con tools/py_stardist.sh).
Uso: tools/py_stardist.sh tools/R2_segment_stardist.py [--scales 1,2,4] [--rois A4_f2] [--tag X] [--seed 0]
"""
import sys, time, argparse, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import R2_common as C
import tensorflow as tf
from stardist.models import StarDist2D
from csbdeep.utils import normalize

ap = argparse.ArgumentParser()
ap.add_argument("--scales", default="1,2,4"); ap.add_argument("--rois", default=None); ap.add_argument("--tag", default=None)
ap.add_argument("--seed", type=int, default=0); ap.add_argument("--force", action="store_true"); ap.add_argument("--prob", type=float, default=None); ap.add_argument("--nms", type=float, default=None)
a = ap.parse_args()
method = a.tag or "stardist_he"
scales = [int(s) for s in a.scales.split(",")]
tf.random.set_seed(a.seed); np.random.seed(a.seed)
model = StarDist2D.from_pretrained("2D_versatile_he")
rois = C.rois_table(); rois = rois[rois.inside_image]
if a.rois: rois = rois[(rois.archetype + "_" + rois.roi_id).isin(a.rois.split(","))]
runs = []
for _, r in rois.iterrows():
    rgb = C.load_roi(r.archetype, r.roi_id); upp0 = r.um_per_px
    for s in scales:
        outp = C.NUC_DIR / f"{r.archetype}_{r.roi_id}_{method}_{C.SCALES[s]}.parquet"
        if outp.exists() and not a.force:
            print(f"{r.archetype}_{r.roi_id} {method} s={s} gia' presente, salto", flush=True); continue
        img = C.downsample(rgb, s); upp = upp0 * s
        x = normalize(img, 1, 99.8, axis=(0, 1))
        t0 = time.time()
        kw = {}
        if a.prob is not None: kw["prob_thresh"] = a.prob
        if a.nms is not None: kw["nms_thresh"] = a.nms
        if max(img.shape[:2]) > 2048:
            labels, _ = model.predict_instances_big(x, axes="YXC", block_size=2048, min_overlap=128, context=128, n_tiles=(2, 2, 1), **kw)
        else:
            labels, _ = model.predict_instances(x, **kw)
        dt = time.time() - t0
        tm, bm = C.masks_for(rgb, upp0, s, r.archetype, r.roi_id)
        tab = C.nuclei_table(labels, upp, tm, bm, r.archetype, r.roi_id, method, s)
        C.save_labels(labels, r.archetype, r.roi_id, method, s)
        tab.to_parquet(C.NUC_DIR / f"{r.archetype}_{r.roi_id}_{method}_{C.SCALES[s]}.parquet", index=False)
        n_keep = int(tab.keep.sum()) if len(tab) else 0
        runs.append(dict(archetype=r.archetype, roi_id=r.roi_id, method=method, scale=s, um_per_px=upp, n_labels=int(labels.max()),
                         n_keep=n_keep, label_px_sum=int((labels > 0).sum()), tissue_frac=float(tm.mean()), bubble_frac=float(bm.mean()),
                         seconds=round(dt, 1), seed=a.seed, prob_thresh=model.thresholds.prob if a.prob is None else a.prob,
                         nms_thresh=model.thresholds.nms if a.nms is None else a.nms))
        print(f"{r.archetype}_{r.roi_id} {method} s={s} n={labels.max()} keep={n_keep} tissue={tm.mean():.2f} bubble={bm.mean():.3f} {dt:.0f}s", flush=True)
out = C.ROOT / "results/R2/R2_seg_runs.csv"
pd.DataFrame(runs).to_csv(out, mode="a", header=not out.exists(), index=False)
