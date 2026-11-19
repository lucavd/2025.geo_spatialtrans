#!/usr/bin/env python
"""tools/F3.3_explore_seg.py - F3.3, esplorazione E-1 (solo dati di taratura): la catena delle aree nuove riproduce i ritagli e le
maschere di R2? Ritaglio dall'immagine a piena risoluzione con il riquadro di results/R2/R2_rois_checked.csv (R1_verify.read_region,
come tools/R2_rois.py) contro il TIFF salvato; Cellpose-SAM RGB (cpsam_v2, parametri di R2) contro la maschera nativa salvata.
Uscita: /mnt/micron/geo_spatialtrans/F3.3/explore/E1_seg.csv. Nessun dato di verifica letto.
Uso: .venv/bin/python tools/F3.3_explore_seg.py A1_r1,A2_r1,...
"""
import sys, time, numpy as np, pandas as pd, torch, tifffile
from pathlib import Path
sys.path.insert(0, "tools")
import R2_common as C
from R1_verify import read_region
from R2_rois import img_path, DATA
from cellpose import models
OUT = Path("/mnt/micron/geo_spatialtrans/F3.3/explore"); OUT.mkdir(parents=True, exist_ok=True)
tags = sys.argv[1].split(",")
rois = C.rois_table(); rois["tag"] = rois.archetype + "_" + rois.roi_id
torch.manual_seed(0); np.random.seed(0)
model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")
def match_frac(a, b):
    """quota dei nuclei di a con un nucleo di b a IoU >= 0.5 (appaiamento per sovrapposizione massima)"""
    m = (a > 0) & (b > 0)
    if a.max() == 0: return np.nan
    pa, pb = a[m].astype(np.int64), b[m].astype(np.int64)
    key = pa * (int(b.max()) + 1) + pb
    u, cnt = np.unique(key, return_counts=True); ia, ib = u // (int(b.max()) + 1), u % (int(b.max()) + 1)
    na = np.bincount(a.ravel(), minlength=a.max() + 1); nb = np.bincount(b.ravel(), minlength=b.max() + 1)
    iou = cnt / (na[ia] + nb[ib] - cnt)
    best = pd.DataFrame({"a": ia, "iou": iou}).groupby("a").iou.max()
    return float((best >= 0.5).sum() / (np.unique(a[a > 0]).size))
rows = []
for tag in tags:
    r = rois.set_index("tag").loc[tag]
    saved = C.load_roi(r.archetype, r.roi_id)
    crop = read_region(img_path(r.dataset), int(r.r0), int(r.r1), int(r.c0), int(r.c1))
    same_crop = bool(crop.shape == saved.shape and np.array_equal(crop, saved))
    t0 = time.time()
    lab, _, _ = model.eval(saved, diameter=None, flow_threshold=0.4, cellprob_threshold=0.0, batch_size=32, normalize=True, channel_axis=2)
    dt = time.time() - t0
    old = C.load_labels(r.archetype, r.roi_id, "cellpose_rgb", 1)
    same_mask = bool(lab.shape == old.shape and np.array_equal(lab.astype(np.uint32), old))
    rows.append(dict(tag=tag, same_crop=same_crop, same_mask=same_mask, n_new=int(lab.max()), n_old=int(old.max()),
                     match_old_in_new=match_frac(old.astype(np.int64), lab.astype(np.int64)),
                     match_new_in_old=match_frac(lab.astype(np.int64), old.astype(np.int64)), sec=round(dt, 1)))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "E1_seg.csv", index=False)
