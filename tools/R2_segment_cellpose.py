#!/usr/bin/env python
"""tools/R2_segment_cellpose.py — segmentazione nucleare con Cellpose-SAM (cpsam_v2) sui ROI R2.
Uso: .venv/bin/python tools/R2_segment_cellpose.py [--variant hed|rgb] [--scales 1,2,4] [--rois A4_f2,A1_r1] [--tag X] [--seed 0]
Input al modello: variante 'hed' = canale ematossilina (nuclei chiari, 1 canale); 'rgb' = immagine H&E cosi' com'e'.
Output: maschere npz su micron, tabella nuclei parquet in results/R2/nuclei/, log tempi in results/R2/R2_seg_runs.csv.
"""
import sys, time, argparse, numpy as np, pandas as pd, torch
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import R2_common as C
from cellpose import models

ap = argparse.ArgumentParser()
ap.add_argument("--variant", default="hed"); ap.add_argument("--scales", default="1,2,4")
ap.add_argument("--rois", default=None); ap.add_argument("--tag", default=None); ap.add_argument("--seed", type=int, default=0); ap.add_argument("--force", action="store_true")
ap.add_argument("--diameter", type=float, default=None); ap.add_argument("--flow_threshold", type=float, default=0.4); ap.add_argument("--cellprob_threshold", type=float, default=0.0)
a = ap.parse_args()
method = a.tag or f"cellpose_{a.variant}"
scales = [int(s) for s in a.scales.split(",")]
torch.manual_seed(a.seed); np.random.seed(a.seed)
model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")
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
        x = C.hematoxylin(img) if a.variant == "hed" else img
        t0 = time.time()
        labels, flows, styles = model.eval(x, diameter=a.diameter, flow_threshold=a.flow_threshold, cellprob_threshold=a.cellprob_threshold,
                                           batch_size=32, normalize=True, channel_axis=(None if x.ndim == 2 else 2))
        dt = time.time() - t0
        tm, bm = C.masks_for(rgb, upp0, s, r.archetype, r.roi_id)
        tab = C.nuclei_table(labels, upp, tm, bm, r.archetype, r.roi_id, method, s)
        C.save_labels(labels, r.archetype, r.roi_id, method, s)
        tab.to_parquet(C.NUC_DIR / f"{r.archetype}_{r.roi_id}_{method}_{C.SCALES[s]}.parquet", index=False)
        n_keep = int(tab.keep.sum()) if len(tab) else 0
        runs.append(dict(archetype=r.archetype, roi_id=r.roi_id, method=method, scale=s, um_per_px=upp, n_labels=int(labels.max()),
                         n_keep=n_keep, label_px_sum=int((labels > 0).sum()), tissue_frac=float(tm.mean()), bubble_frac=float(bm.mean()),
                         seconds=round(dt, 1), seed=a.seed, flow_threshold=a.flow_threshold, cellprob_threshold=a.cellprob_threshold))
        print(f"{r.archetype}_{r.roi_id} {method} s={s} n={labels.max()} keep={n_keep} tissue={tm.mean():.2f} bubble={bm.mean():.3f} {dt:.0f}s", flush=True)
out = C.ROOT / "results/R2/R2_seg_runs.csv"
pd.DataFrame(runs).to_csv(out, mode="a", header=not out.exists(), index=False)
