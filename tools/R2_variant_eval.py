#!/usr/bin/env python
"""tools/R2_variant_eval.py — rivalutazione delle varianti di input di Cellpose (richiesta di Luca, 2026-09-20).
Per ROI di prova e per variante: area degli oggetti, copertura del tessuto, contenuto nucleare di ogni oggetto
(OD ematossilina media/massima; numero di centroidi StarDist e Cellpose-RGB contenuti), nuclei RGB/StarDist non coperti.
Figure: overlay 300 um con titoli; oggetti hed colorati per contenuto nucleare.
"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from skimage import color, measure
from skimage.segmentation import find_boundaries
sys.path.insert(0, str(Path(__file__).resolve().parent)); import R2_common as C
OUT = C.ROOT / "results/R2/variant_eval"; OUT.mkdir(exist_ok=True)
METHODS = {"test_cp_hed": "Cellpose-SAM ematossilina", "test_cp_hed_d25": "Cellpose-SAM ematossilina, diam 25 px", "test_cp_rgb": "Cellpose-SAM RGB", "test_sd": "StarDist HE"}
rows = []; cover_rows = []
for A, roi in [("A4", "f2"), ("A1", "r1")]:
    rgb = C.load_roi(A, roi); rr = C.rois_table().query("archetype==@A and roi_id==@roi").iloc[0]; upp = rr.um_per_px
    hod = color.rgb2hed(rgb)[..., 0]                      # OD ematossilina (continua)
    tm, bm = C.masks_for(rgb, upp, 1); tissue_px = tm.sum()
    labs = {m: C.load_labels(A, roi, m, 1) for m in METHODS if (C.MASK_DIR / f"{A}_{roi}_{m}_native.npz").exists()}
    cents = {}
    for m, L in labs.items():
        p = measure.regionprops_table(L, properties=("label", "centroid", "area")); d = pd.DataFrame(p)
        cents[m] = (d["centroid-0"].round().astype(int).values, d["centroid-1"].round().astype(int).values)
    # OD di riferimento: mediana dell'OD nei nuclei StarDist (nuclei "certi") e nel tessuto fuori da ogni nucleo StarDist
    sd = labs["test_sd"]; od_nuc = np.median(hod[sd > 0]); od_bg = np.median(hod[(sd == 0) & tm])
    thr = (od_nuc + od_bg) / 2
    for m, L in labs.items():
        p = measure.regionprops_table(L, intensity_image=hod, properties=("label", "area", "intensity_mean", "intensity_max"))
        d = pd.DataFrame(p); d["area_um2"] = d.area * upp ** 2
        # frazione di pixel "nucleari" (OD > soglia) in ogni oggetto
        nucpx = np.bincount(L.ravel(), weights=(hod > thr).ravel(), minlength=L.max() + 1)[1:]
        d["frac_nuclear_px"] = nucpx / d.area
        for other in [k for k in labs if k != m]:
            ci = np.bincount(L[cents[other]], minlength=L.max() + 1)[1:]     # quanti centroidi di 'other' cadono in ogni oggetto
            d[f"n_{other}_centroids"] = ci
        d["archetype"] = A; d["roi_id"] = roi; d["method"] = m
        rows.append(d)
        # nuclei di 'other' non coperti da alcun oggetto di m
        for other in [k for k in labs if k != m]:
            covered = L[cents[other]] > 0
            cover_rows.append(dict(archetype=A, roi_id=roi, method=m, other=other, n_other=len(covered), frac_other_covered=covered.mean()))
        cover_rows.append(dict(archetype=A, roi_id=roi, method=m, other="tissue", n_other=int(L.max()), frac_other_covered=(L > 0)[tm].mean()))
    # figure: overlay 300 um, 2 finestre (centro e quarto in alto a sinistra)
    w = int(300 / upp)
    wins = {"centro": (rgb.shape[0] // 2 - w // 2, rgb.shape[1] // 2 - w // 2), "NW": (rgb.shape[0] // 4 - w // 2, rgb.shape[1] // 4 - w // 2)}
    for wname, (r0, c0) in wins.items():
        sub = rgb[r0:r0 + w, c0:c0 + w]
        fig, axs = plt.subplots(2, 3, figsize=(21, 14), facecolor="white"); axs = axs.ravel()
        axs[0].imshow(sub); axs[0].set_title(f"{A}_{roi} H&E — finestra {wname} 300 µm", fontsize=13)
        hsub = C.hematoxylin(sub); axs[1].imshow(hsub, cmap="gray"); axs[1].set_title("canale ematossilina (input variante hed)", fontsize=13)
        for ax, m in zip(axs[2:], [k for k in METHODS if k in labs]):
            L = labs[m][r0:r0 + w, c0:c0 + w]; b = find_boundaries(L, mode="outer"); im = sub.copy(); im[b] = [0, 255, 0]
            ax.imshow(im); ax.set_title(f"{METHODS[m]}\nn oggetti nella finestra = {len(np.unique(L)) - 1}", fontsize=13)
        for ax in axs: ax.set_axis_off()
        fig.tight_layout(); fig.savefig(OUT / f"{A}_{roi}_{wname}_overlay300.png", dpi=90); plt.close(fig)
    # figura: oggetti hed colorati per contenuto nucleare (verde: >=1 centroide StarDist; giallo: solo RGB; rosso: nessuno)
    r0, c0 = wins["centro"]; sub = rgb[r0:r0 + w, c0:c0 + w]
    for m in ["test_cp_hed", "test_cp_rgb"]:
        if m not in labs: continue
        L = labs[m]; d = pd.concat([x for x in rows if x.method.iloc[0] == m and x.archetype.iloc[0] == A])
        cls = np.zeros(L.max() + 1, np.uint8)
        has_sd = d["n_test_sd_centroids"].values > 0
        has_other = (d["n_test_cp_rgb_centroids"].values > 0) if m == "test_cp_hed" else (d["n_test_cp_hed_centroids"].values > 0)
        cls[d.label.values] = np.where(has_sd, 1, np.where(has_other, 2, 3))
        cw = cls[L[r0:r0 + w, c0:c0 + w]]
        pal = np.array([[0, 0, 0], [0, 200, 0], [230, 200, 0], [230, 0, 0]], np.uint8)
        over = sub.copy().astype(float); col = pal[cw].astype(float); mk = cw > 0
        over[mk] = 0.45 * over[mk] + 0.55 * col[mk]
        fig, ax = plt.subplots(1, 2, figsize=(16, 8.4), facecolor="white")
        ax[0].imshow(sub); ax[0].set_title(f"{A}_{roi} H&E — centro 300 µm", fontsize=13)
        ax[1].imshow(over.astype(np.uint8)); ax[1].set_title(f"{METHODS[m]}: verde = contiene un nucleo StarDist; giallo = contiene solo un nucleo dell'altra variante Cellpose; rosso = nessun nucleo", fontsize=11)
        for a in ax: a.set_axis_off()
        fig.tight_layout(); fig.savefig(OUT / f"{A}_{roi}_{m}_nucleus_content.png", dpi=90); plt.close(fig)
    print(A, roi, "OD nucleo SD", round(od_nuc, 3), "OD fondo", round(od_bg, 3), "soglia", round(thr, 3), flush=True)
obj = pd.concat(rows); obj.to_parquet(OUT / "objects.parquet", index=False)
cov = pd.DataFrame(cover_rows); cov.to_csv(OUT / "coverage.csv", index=False)
# tabella di sintesi
def q(s): return f"{s.median():.0f} [{s.quantile(.25):.0f}–{s.quantile(.75):.0f}]"
summ = []
for (A, roi, m), g in obj.groupby(["archetype", "roi_id", "method"]):
    c_tis = cov.query("archetype==@A and roi_id==@roi and method==@m and other=='tissue'").frac_other_covered.iloc[0]
    row = dict(archetype=A, roi_id=roi, method=m, n=len(g), area_um2_median_IQR=q(g.area_um2), frac_tissue_covered=round(c_tis, 3),
               frac_obj_nuclear_px_median=round(g.frac_nuclear_px.median(), 2), frac_obj_with_sd_nucleus=round((g.get("n_test_sd_centroids", pd.Series(np.nan)) > 0).mean(), 3),
               frac_obj_with_ge2_sd_nuclei=round((g.get("n_test_sd_centroids", pd.Series(np.nan)) >= 2).mean(), 3))
    if "n_test_cp_rgb_centroids" in g: row["frac_obj_with_rgb_nucleus"] = round((g.n_test_cp_rgb_centroids > 0).mean(), 3)
    if "n_test_cp_hed_centroids" in g: row["frac_obj_with_hed_centroid"] = round((g.n_test_cp_hed_centroids > 0).mean(), 3)
    summ.append(row)
pd.DataFrame(summ).to_csv(OUT / "summary.csv", index=False)
print(pd.DataFrame(summ).to_string()); print(cov.to_string())
