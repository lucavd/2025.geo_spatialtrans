"""tools/F3.3_segment.py - F3.3: per ogni area della tabella results/F3.3/F3.3_areas.csv, ritaglio dell'immagine a piena risoluzione
(R1_verify.read_region, come tools/R2_rois.py), maschere di tessuto e bolle a 1/4 (tools/R2_common.py, come masks_for) e maschera valida
a 1/4 in PNG (come tools/R2_export_masks.py); nuclei: Cellpose-SAM RGB con i parametri di R2 (aree 'cellpose_rgb') oppure i poligoni di
Space Ranger interamente nel riquadro rasterizzati come tools/R2_spaceranger_rois.py (aree 'consensus'; StarDist e consenso in
tools/F3.3_stardist.py e tools/F3.3_consensus.py). Tabella dei nuclei e flag keep con R2_common.nuclei_table.
Uscite in /mnt/micron/geo_spatialtrans/F3.3/seg/<area>/. Uso: .venv/bin/python tools/F3.3_segment.py <area,area,...|IN|holdout>"""
import sys, os, json, time, numpy as np, pandas as pd, tifffile
from pathlib import Path
sys.path.insert(0, "tools")
import R2_common as C
from R1_verify import read_region
from skimage.io import imsave
from skimage.draw import polygon as draw_polygon
OUT = Path("/mnt/micron/geo_spatialtrans/F3.3/seg")
A = pd.read_csv("results/F3.3/F3.3_areas.csv")
sel = sys.argv[1]
A = A[A.level == "IN"] if sel == "IN" else A[A.level != "IN"] if sel == "holdout" else A[A.area.isin(sel.split(","))]
model = None; _sr = {}
def cellpose(rgb):
    global model
    import torch
    from cellpose import models
    if model is None:
        torch.manual_seed(0); np.random.seed(0); model = models.CellposeModel(gpu=True, pretrained_model="cpsam_v2")
    lab, _, _ = model.eval(rgb, diameter=None, flow_threshold=0.4, cellprob_threshold=0.0, batch_size=32, normalize=True, channel_axis=2)
    return lab
def sr_polys(d):
    if d not in _sr:
        _sr.clear(); feats = json.load(open(f"{d}/segmented_outputs_x/segmented_outputs/nucleus_segmentations.geojson"))["features"]
        polys = [np.asarray(f["geometry"]["coordinates"][0], dtype=np.float64) for f in feats]
        bb = np.array([[q[:, 0].min(), q[:, 0].max(), q[:, 1].min(), q[:, 1].max()] for q in polys]); _sr[d] = (polys, bb)
    return _sr[d]
for _, r in A.iterrows():
    t0 = time.time(); o = OUT / r.area; o.mkdir(parents=True, exist_ok=True)
    if (o / "seg.json").exists(): print(r.area, "presente, salto", flush=True); continue
    rgb = read_region(r.image, int(r.r0), int(r.r1), int(r.c0), int(r.c1))
    tifffile.imwrite(o / "roi.tif", rgb, photometric="rgb", compression="zlib")
    rgb4 = C.downsample(rgb, 4); upp = float(r.upp); t4 = C.tissue_mask(rgb4, upp * 4); b4 = C.bubble_mask(rgb4, upp * 4)
    np.savez_compressed(o / "masks.npz", tissue_ds4=t4, bubble_ds4=b4)
    tm4, bm4 = C.masks_for(rgb, upp, 4); assert np.array_equal(tm4, t4[:tm4.shape[0], :tm4.shape[1]])
    imsave(o / "valid_ds4.png", ((tm4 & ~bm4) * 255).astype(np.uint8), check_contrast=False)
    tm, bm = C.masks_for(rgb, upp, 1)
    info = {"area": r.area, "seg": r.seg, "side_px": int(r.side_px)}
    if r.seg == "cellpose_rgb":
        lab = cellpose(rgb).astype(np.uint32); meth = "cellpose_rgb"
    else:
        polys, bb = sr_polys(r.ds_dir)
        s = np.where((bb[:, 0] >= r.c0) & (bb[:, 1] < r.c1) & (bb[:, 2] >= r.r0) & (bb[:, 3] < r.r1))[0]
        lab = np.zeros((int(r.side_px), int(r.side_px)), np.uint32)
        for k, i in enumerate(s, start=1):
            q = polys[i]; rr, cc = draw_polygon(q[:, 1] - r.r0, q[:, 0] - r.c0, shape=lab.shape); lab[rr, cc] = k
        meth = "spaceranger"
    np.savez_compressed(o / f"{meth}_native.npz", labels=lab)
    tab = C.nuclei_table(lab, upp, tm, bm, r.A, r.area, meth, 1); tab.to_parquet(o / f"{meth}_native.parquet", index=False)
    info.update(method=meth, n_labels=int(lab.max()), n_keep=int(tab.keep.sum()) if len(tab) else 0, sec=round(time.time() - t0, 1))
    json.dump(info, open(o / "seg.json", "w")); print(info, flush=True)
