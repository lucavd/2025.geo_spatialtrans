#!/usr/bin/env python
"""tools/R2_metrics.py — metriche R2 dalle tabelle nuclei e dalle maschere.
Output (results/R2/):
  R2_roi_summary.csv      una riga per ROI x metodo x scala: n nuclei validi, area tessuto (mm2), densita' (nuclei/mm2 di tessuto),
                          area nucleare mediana/IQR, eq_diam, eccentricita', NN distance mediana, frazioni flag
  R2_archetype_summary.csv  per archetipo x metodo x scala: mediana fra ROI e min-max (n ROI)
  R2_matching.csv         concordanza fra coppie di metodi a scala nativa per ROI: F1@IoU0.5, precision, recall, Dice binario, Δ densita'
  R2_scale_effect.csv     densita' e area a 2x/4x relative alla nativa (CP-R2.3)
  R2_repro.csv            C-R2.2: confronto run ripetuti
  R2_nuclei_all.parquet   tabella nuclei concatenata (fuori git se > 50 MB)
"""
import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy.spatial import cKDTree
from scipy import sparse
sys.path.insert(0, str(Path(__file__).resolve().parent)); import R2_common as C
R = C.ROOT / "results/R2"
rois = C.rois_table().set_index(["archetype", "roi_id"])

# --- tabella nuclei
parts = [pd.read_parquet(p) for p in sorted(C.NUC_DIR.glob("*.parquet")) if "test_" not in p.name]
nuc = pd.concat(parts, ignore_index=True)
nuc.to_parquet(R / "R2_nuclei_all.parquet", index=False)

def nn_median(d):
    if len(d) < 2: return np.nan
    t = cKDTree(d[["x_um", "y_um"]].values); dist, _ = t.query(d[["x_um", "y_um"]].values, k=2); return float(np.median(dist[:, 1]))

runs = pd.read_csv(R / "R2_seg_runs.csv").drop_duplicates(["archetype", "roi_id", "method", "scale"], keep="last")
rows = []; MASKFRAC = {}
for (A, roi, m, s), d in nuc.groupby(["archetype", "roi_id", "method", "scale"]):
    upp = rois.loc[(A, roi), "um_per_px"]; side_mm2 = (rois.loc[(A, roi), "side_px"] * upp / 1000) ** 2
    if (A, roi) not in MASKFRAC:
        tm_, bm_ = C.masks_for(C.load_roi(A, roi), upp, 4, A, roi); MASKFRAC[(A, roi)] = (float(tm_.mean()), float((tm_ & bm_).mean()))
    tissue_frac, bubble_frac = MASKFRAC[(A, roi)]      # frazione di tessuto e frazione (del ROI) di bolle dentro il tessuto
    tissue_mm2 = side_mm2 * (tissue_frac - bubble_frac)     # area valida = tessuto meno bolle nel tessuto
    k = d[d.keep]
    rows.append(dict(archetype=A, roi_id=roi, method=m, scale=s, um_per_px=k.um_per_px.iloc[0] if len(k) else upp * s,
                     n_all=len(d), n_keep=len(k), frac_small=d.flag_small.mean(), frac_large=d.flag_large.mean(), frac_out_tissue=(~d.in_tissue).mean(), frac_bubble=d.in_bubble.mean(),
                     tissue_frac=tissue_frac, bubble_frac=bubble_frac, tissue_mm2=tissue_mm2, density_per_mm2=len(k) / tissue_mm2 if tissue_mm2 > 0 else np.nan,
                     area_median=k.area_um2.median(), area_q25=k.area_um2.quantile(.25), area_q75=k.area_um2.quantile(.75), area_mean=k.area_um2.mean(), area_cv=k.area_um2.std() / k.area_um2.mean(),
                     eq_diam_median=k.eq_diam_um.median(), ecc_median=k.eccentricity.median(), ecc_q75=k.eccentricity.quantile(.75), major_minor_median=(k.major_um / k.minor_um.clip(lower=1e-3)).median(),
                     nn_median_um=nn_median(k), nuclear_area_fraction=k.area_um2.sum() / (tissue_mm2 * 1e6) if tissue_mm2 > 0 else np.nan))
summ = pd.DataFrame(rows).sort_values(["method", "scale", "archetype", "roi_id"]); summ.to_csv(R / "R2_roi_summary.csv", index=False)

agg = summ.groupby(["archetype", "method", "scale"]).agg(n_roi=("roi_id", "count"), n_nuclei=("n_keep", "sum"),
      density_median=("density_per_mm2", "median"), density_min=("density_per_mm2", "min"), density_max=("density_per_mm2", "max"),
      area_median=("area_median", "median"), area_min=("area_median", "min"), area_max=("area_median", "max"),
      eq_diam_median=("eq_diam_median", "median"), ecc_median=("ecc_median", "median"), nn_median_um=("nn_median_um", "median"),
      nuclear_area_fraction=("nuclear_area_fraction", "median"), frac_small=("frac_small", "mean"), frac_large=("frac_large", "mean")).reset_index()
agg.to_csv(R / "R2_archetype_summary.csv", index=False)

# --- matching fra metodi (scala nativa)
def match(la, lb):
    """F1 a IoU>=0.5 fra due label image (stessa forma). Coppie IoU>=0.5 sono automaticamente uno-a-uno."""
    a = la.ravel().astype(np.int64); b = lb.ravel().astype(np.int64); m = (a > 0) | (b > 0)
    a = a[m]; b = b[m]
    inter = sparse.coo_matrix((np.ones(len(a), np.int64), (a, b)), shape=(la.max() + 1, lb.max() + 1)).tocsr()
    area_a = np.bincount(la.ravel(), minlength=la.max() + 1); area_b = np.bincount(lb.ravel(), minlength=lb.max() + 1)
    inter = inter.tocoo(); keep = (inter.row > 0) & (inter.col > 0)
    r_, c_, v = inter.row[keep], inter.col[keep], inter.data[keep]
    iou = v / (area_a[r_] + area_b[c_] - v)
    tp = int((iou >= 0.5).sum()); na, nb = int((area_a[1:] > 0).sum()), int((area_b[1:] > 0).sum())   # oggetti con area > 0 (le label azzerate fuori area valida non contano; rilievo revisore R2)
    prec = tp / nb if nb else np.nan; rec = tp / na if na else np.nan
    f1 = 2 * tp / (na + nb) if (na + nb) else np.nan
    dice = 2 * float(((la > 0) & (lb > 0)).sum()) / float((la > 0).sum() + (lb > 0).sum())
    return dict(tp=tp, n_a=na, n_b=nb, precision_b_vs_a=prec, recall_b_vs_a=rec, f1_iou05=f1, dice_binary=dice, median_iou_matched=float(np.median(iou[iou >= 0.5])) if tp else np.nan)

pairs = [("cellpose_rgb", "stardist_he"), ("cellpose_rgb", "spaceranger"), ("stardist_he", "spaceranger"), ("cellpose_rgb", "cellpose_hed")]
mrows = []
for (A, roi), r in rois.iterrows():
    for ma, mb in pairs:
        fa, fb = C.MASK_DIR / f"{A}_{roi}_{ma}_native.npz", C.MASK_DIR / f"{A}_{roi}_{mb}_native.npz"
        if not (fa.exists() and fb.exists()): continue
        la, lb = C.load_labels(A, roi, ma, 1), C.load_labels(A, roi, mb, 1)
        tm, bm = C.masks_for(C.load_roi(A, roi), r.um_per_px, 1, A, roi); valid = tm & ~bm
        la = np.where(valid, la, 0); lb = np.where(valid, lb, 0)    # confronto solo nell'area valida
        d = match(la, lb)
        da = summ.query("archetype==@A and roi_id==@roi and method==@ma and scale==1").density_per_mm2; db = summ.query("archetype==@A and roi_id==@roi and method==@mb and scale==1").density_per_mm2
        d.update(archetype=A, roi_id=roi, method_a=ma, method_b=mb, density_a=float(da.iloc[0]) if len(da) else np.nan, density_b=float(db.iloc[0]) if len(db) else np.nan)
        d["delta_density_rel"] = (d["density_b"] - d["density_a"]) / d["density_a"] if d["density_a"] else np.nan
        mrows.append(d); print(f"match {A}_{roi} {ma} vs {mb}: F1={d['f1_iou05']:.3f} dice={d['dice_binary']:.3f} dDens={d['delta_density_rel']:+.3f}", flush=True)
pd.DataFrame(mrows).to_csv(R / "R2_matching.csv", index=False)

# --- effetto scala (CP-R2.3)
base = summ[summ.scale == 1].set_index(["archetype", "roi_id", "method"])
sc = summ[summ.scale > 1].join(base[["density_per_mm2", "area_median", "n_keep"]], on=["archetype", "roi_id", "method"], rsuffix="_native")
sc["density_rel"] = sc.density_per_mm2 / sc.density_per_mm2_native; sc["area_rel"] = sc.area_median / sc.area_median_native
sc[["archetype", "roi_id", "method", "scale", "um_per_px", "n_keep", "n_keep_native", "density_rel", "area_rel"]].to_csv(R / "R2_scale_effect.csv", index=False)

# --- C-R2.2 riproducibilita'
rep = []
for m in ["cellpose_rgb", "stardist_he"]:
    for (A, roi) in [("A4", "f2"), ("A1", "r1")]:
        f1_, f2_ = C.MASK_DIR / f"{A}_{roi}_{m}_native.npz", C.MASK_DIR / f"{A}_{roi}_{m}_repro_native.npz"
        if f1_.exists() and f2_.exists():
            a, b = C.load_labels(A, roi, m, 1), C.load_labels(A, roi, m + "_repro", 1)
            rep.append(dict(archetype=A, roi_id=roi, method=m, n_run1=int(a.max()), n_run2=int(b.max()), px_run1=int((a > 0).sum()), px_run2=int((b > 0).sum()),
                            identical_labels=bool(np.array_equal(a, b)), identical_binary=bool(np.array_equal(a > 0, b > 0)), f1_iou05=match(a, b)["f1_iou05"]))
pd.DataFrame(rep).to_csv(R / "R2_repro.csv", index=False); print(pd.DataFrame(rep).to_string())
print(agg[agg.scale == 1].to_string())
