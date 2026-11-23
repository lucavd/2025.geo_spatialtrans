"""tools/sonda/S3b_explore_masks.py - S3b, esplorazione prima della regola: maschere in_tissue (8 um) e H&E hires di H1 e H2,
metà di sviluppo (x < m - 50 um) e aree di S3. Nessuna espressione letta. Uscita: /mnt/micron/geo_spatialtrans/S3b_sonda/explore/masks.png"""
import os, sys, json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
from PIL import Image
OUT = "/mnt/micron/geo_spatialtrans/S3b_sonda/explore"; os.makedirs(OUT, exist_ok=True)
fig, ax = plt.subplots(2, 2, figsize=(12, 12)); info = {}
for j, sec in enumerate(("H1", "H2")):
    tp = C.positions(sec, "008"); m = C.split_m(sec)
    R, Cc = int(tp.array_row.max()) + 1, int(tp.array_col.max()) + 1
    M = np.zeros((R, Cc)); M[tp.array_row, tp.array_col] = tp.in_tissue
    sf = json.load(open(f"{C.bdir(sec, '008')}/spatial/scalefactors_json.json"))
    img = np.asarray(Image.open(f"{C.bdir(sec, '008')}/spatial/tissue_hires_image.png"))
    A = C.areas(sec)
    a = ax[0, j]; a.imshow(M, cmap="Greys", extent=(0, Cc * 8, R * 8, 0)); a.axvline(m - 50, color="r"); a.set_title(f"{sec} in_tissue 8um, m={m:.0f}")
    for r_ in A.itertuples():
        a.add_patch(plt.Rectangle((r_.x0, r_.y0), r_.side, r_.side, fill=False, color="tab:orange" if r_.tissue == "tumour" else "tab:blue")); a.text(r_.x0, r_.y0, f"{r_.tissue[0]}{r_.k}", color="k", fontsize=8)
    ax[1, j].imshow(img); ax[1, j].set_title(f"{sec} hires H&E {img.shape}")
    # pixel di piena risoluzione dei quadratini, per la relazione array -> immagine
    cols = [c for c in tp.columns if "pxl" in c]
    info[sec] = dict(R=R, C=Cc, m=m, n_tissue=int(tp.in_tissue.sum()), sf=sf, img_shape=list(img.shape), pxl_cols=cols,
                     corner_px=tp.loc[[tp.index[0], tp.index[-1]], ["array_row", "array_col"] + cols].values.tolist())
fig.savefig(f"{OUT}/masks.png", dpi=80); json.dump(info, open(f"{OUT}/masks.json", "w"), indent=1, default=float); print(json.dumps(info, default=float))
