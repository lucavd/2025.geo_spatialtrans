#!/usr/bin/env python
"""tools/R2_overview.py — panoramiche per la scelta dei ROI (sessione R2).
Per ogni dataset: H&E sottocampionato (fattore F, media a blocchi) e mappa dei cluster graphclust 8 um
proiettata nelle stesse coordinate pixel (pxl_row/pxl_col in fullres / F). Griglia ogni 1 mm.
Output: results/R2/overview/<A>_<dataset>_overview.png e <..>_lowres.npy (H&E ridotto, per riuso).
"""
import sys, json, numpy as np, pandas as pd, tifffile, zarr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/mnt/micron/geo_spatialtrans/data_real")
OUT = ROOT / "results/R2/overview"; OUT.mkdir(parents=True, exist_ok=True)
F = 16
MAP = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer",
       "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}

def img_path(ds):
    c = list((DATA / ds).glob(f"{ds}_tissue_image.*"))
    c = [p for p in c if p.suffix in (".btf", ".tif", ".tiff")]
    assert len(c) == 1, c
    return c[0]

def lowres(path, F=F, chunk=4096):
    store = tifffile.imread(path, aszarr=True)
    try:
        z = zarr.open(store, mode="r")
        if isinstance(z, zarr.Group): z = z["0"]
        with tifffile.TiffFile(path) as tf: axes = tf.series[0].axes
        cyx = not axes.startswith("YX")
        H, W = (z.shape[1], z.shape[2]) if cyx else (z.shape[0], z.shape[1])
        h, w = H // F, W // F
        out = np.zeros((h, w, 3), np.float32)
        for r0 in range(0, h * F, chunk):
            r1 = min(r0 + chunk, h * F)
            blk = z[:, r0:r1, :w * F] if cyx else z[r0:r1, :w * F]
            if cyx: blk = np.moveaxis(blk, 0, -1)
            blk = np.asarray(blk)[..., :3].astype(np.float32)
            nb = (r1 - r0) // F
            blk = blk[:nb * F].reshape(nb, F, w, F, 3).mean(axis=(1, 3))
            out[r0 // F:r0 // F + nb] = blk
            print(f"  rows {r1}/{H}", flush=True)
    finally:
        store.close()
    return out.astype(np.uint8), (H, W)

def main(arcs):
    for A in arcs:
        ds = MAP[A]; d = DATA / ds
        p = img_path(ds); print(A, ds, p.name, flush=True)
        lr, (H, W) = lowres(p)
        np.save(OUT / f"{A}_{ds}_lowres_F{F}.npy", lr)
        sq = d / "binned_outputs_x/binned_outputs/square_008um"
        pos = pd.read_parquet(sq / "spatial/tissue_positions.parquet")
        cl = pd.read_csv(sq / "analysis/clustering/gene_expression_graphclust/clusters.csv")
        pos = pos.merge(cl, left_on="barcode", right_on="Barcode", how="left")
        sf = json.load(open(sq / "spatial/scalefactors_json.json")); upp = sf["microns_per_pixel"]
        pos = pos[pos.in_tissue == 1]
        n_cl = int(pos.Cluster.max())
        fig, ax = plt.subplots(1, 2, figsize=(22, 11 * H / W + 1), facecolor="white")
        ext = [0, W, H, 0]
        ax[0].imshow(lr, extent=ext); ax[0].set_title(f"{A} {ds} — H&E (1/{F})")
        ax[1].imshow(np.full_like(lr, 255), extent=ext)
        cmap = plt.get_cmap("tab20", n_cl)
        sc = ax[1].scatter(pos.pxl_col_in_fullres, pos.pxl_row_in_fullres, c=pos.Cluster, cmap=cmap, s=0.15, vmin=0.5, vmax=n_cl + 0.5, rasterized=True)
        ax[1].set_title(f"graphclust 8 um ({n_cl} cluster, n={len(pos)})")
        cb = fig.colorbar(sc, ax=ax[1], ticks=range(1, n_cl + 1), fraction=0.03); cb.set_label("cluster")
        mm_px = 1000 / upp
        for a in ax:
            a.set_xticks(np.arange(0, W, mm_px)); a.set_xticklabels([f"{i}" for i in range(len(a.get_xticks()))])
            a.set_yticks(np.arange(0, H, mm_px)); a.set_yticklabels([f"{i}" for i in range(len(a.get_yticks()))])
            a.grid(color="k", alpha=0.25, lw=0.5); a.set_xlabel("mm (x = col)"); a.set_ylabel("mm (y = row)")
            a.set_xlim(0, W); a.set_ylim(H, 0)
        fig.tight_layout(); fig.savefig(OUT / f"{A}_{ds}_overview.png", dpi=110); plt.close(fig)
        pos[["barcode", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres", "Cluster"]].to_parquet(OUT / f"{A}_{ds}_pos8_clusters.parquet")
        print(f"  saved; image {H}x{W}, um/px {upp:.5f}, clusters {n_cl}", flush=True)

if __name__ == "__main__":
    main(sys.argv[1:] or list(MAP))
