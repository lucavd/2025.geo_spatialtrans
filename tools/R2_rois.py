#!/usr/bin/env python
"""tools/R2_rois.py — ROI 1x1 mm per archetipo (sessione R2).
Input : results/R2/R2_rois.csv  (archetype,dataset,roi_id,label,x_mm,y_mm)   centri in mm dell'immagine full-res (x=col,y=row)
Output: /mnt/micron/geo_spatialtrans/R2/rois/<archetype>_<roi_id>.tif   ritaglio nativo RGB uint8 (fuori git)
        results/R2/rois_thumbs/<archetype>_<roi_id>.png                   miniatura 1/4
        results/R2/rois_thumbs/<archetype>_sheet.png                      tavola per archetipo
        results/R2/R2_rois_checked.csv                                     bbox px, um/px, copertura in_tissue (C-R2.4)
Copertura in_tissue: frazione dei bin 8 um attesi (125x125) dentro il bbox che sono in_tissue==1.
"""
import sys, json, numpy as np, pandas as pd, tifffile
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R1_verify import read_region
ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/mnt/micron/geo_spatialtrans/data_real")
ROI_DIR = Path("/mnt/micron/geo_spatialtrans/R2/rois"); ROI_DIR.mkdir(parents=True, exist_ok=True)
TH = ROOT / "results/R2/rois_thumbs"; TH.mkdir(parents=True, exist_ok=True)
SIDE_UM = 1000.0

def img_path(ds):
    c = [p for p in (DATA / ds).glob(f"{ds}_tissue_image.*") if p.suffix in (".btf", ".tif", ".tiff")]
    assert len(c) == 1, c; return c[0]

def main(only=None):
    rois = pd.read_csv(ROOT / "results/R2/R2_rois.csv")
    if only: rois = rois[rois.archetype.isin(only)]
    rows = []
    for ds, g in rois.groupby("dataset", sort=False):
        sq = DATA / ds / "binned_outputs_x/binned_outputs/square_008um"
        upp = json.load(open(sq / "spatial/scalefactors_json.json"))["microns_per_pixel"]
        pos = pd.read_parquet(sq / "spatial/tissue_positions.parquet")
        pos = pos.dropna(subset=["pxl_col_in_fullres"])
        p = img_path(ds)
        with tifffile.TiffFile(p) as tf:
            sh = tf.series[0].shape; H, W = (sh[1], sh[2]) if tf.series[0].axes[0] in "SC" else (sh[0], sh[1])
        side = int(round(SIDE_UM / upp))
        for _, r in g.iterrows():
            cx, cy = r.x_mm * 1000 / upp, r.y_mm * 1000 / upp
            c0, r0 = int(round(cx - side / 2)), int(round(cy - side / 2)); c1, r1 = c0 + side, r0 + side
            inside = (c0 >= 0) and (r0 >= 0) and (c1 <= W) and (r1 <= H)
            m = pos[pos.pxl_col_in_fullres.between(c0, c1) & pos.pxl_row_in_fullres.between(r0, r1)]
            exp_bins = (SIDE_UM / 8.0) ** 2
            cov_bins = len(m) / exp_bins; cov_tissue = m.in_tissue.sum() / exp_bins
            crop = read_region(p, r0, r1, c0, c1) if inside else None
            tag = f"{r.archetype}_{r.roi_id}"
            if inside:
                tifffile.imwrite(ROI_DIR / f"{tag}.tif", crop, photometric="rgb", compression="zlib")
                th = crop[:side // 4 * 4, :side // 4 * 4].reshape(side // 4, 4, side // 4, 4, 3).mean(axis=(1, 3)).astype(np.uint8)
                fig, ax = plt.subplots(figsize=(6, 6), facecolor="white"); ax.imshow(th); ax.set_axis_off()
                ax.set_title(f"{tag} — {r.label}\n({r.x_mm:.2f},{r.y_mm:.2f}) mm · tissue {cov_tissue:.2f}", fontsize=9)
                fig.tight_layout(); fig.savefig(TH / f"{tag}.png", dpi=100); plt.close(fig)
            rows.append(dict(r, um_per_px=upp, side_px=side, c0=c0, r0=r0, c1=c1, r1=r1, img_H=H, img_W=W,
                             inside_image=inside, cov_bins=round(cov_bins, 3), cov_in_tissue=round(cov_tissue, 3)))
            print(tag, r.label, "inside", inside, "cov_in_tissue", round(cov_tissue, 3), flush=True)
    out = pd.DataFrame(rows); out.to_csv(ROOT / "results/R2/R2_rois_checked.csv", index=False)
    for A, g in out.groupby("archetype"):
        n = len(g); fig, axs = plt.subplots(1, n, figsize=(4.2 * n, 4.6), facecolor="white")
        for ax, (_, r) in zip(np.atleast_1d(axs), g.iterrows()):
            f = TH / f"{r.archetype}_{r.roi_id}.png"
            if f.exists(): ax.imshow(plt.imread(f))
            ax.set_axis_off()
        fig.tight_layout(); fig.savefig(TH / f"{A}_sheet.png", dpi=100); plt.close(fig)

if __name__ == "__main__":
    main(sys.argv[1:] or None)
