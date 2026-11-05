"""S2.2, esplorazione E-1 (lato vero del banco della gara): per ogni ROI, i quadratini di 2 um in tessuto che cadono
nell'impronta (nucleo + corona fino a d_exp, come R6) di un nucleo vero con un tipo del preset (criterio corretto di R6, PH-1).
Uso (root del repo): .venv/bin/python tools/S2.2_real.py <KEY>   KEY = A1 | A23 | A4 | A5 | A6
Scrive /mnt/micron/geo_spatialtrans/S2.2/real/<A>_<roi>_bins.npz (+ _cells.parquet) e results/S2.2/E1/<A>_<roi>_real.json;
salta le ROI già scritte. Pre-registrazione: docs/steps/S2.2.md."""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow.parquet as pq
sys.path.insert(0, "tools")
from R6_common import footprint_labels, read_h5_csc
import S22_common as S

M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"; OUT = f"{M}/S2.2/real"; RES = "results/S2.2/E1"
os.makedirs(OUT, exist_ok=True); os.makedirs(RES, exist_ok=True)
KEY = sys.argv[1]
DS = {"A1": ("Visium_HD_Mouse_Small_Intestine", ["A1"]), "A23": ("Visium_HD_Human_Colon_Cancer", ["A2", "A3"]),
      "A4": ("Visium_HD_Human_Lymph_Node_FFPE", ["A4"]), "A5": ("Visium_HD_6p5mm_Mouse_Brain", ["A5"]), "A6": ("Visium_HD_6p5mm_Human_Heart", ["A6"])}
SRC = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
# tipi dei preset dello Step 1 (R/04b0_tissue_presets.R, righe 42-43); controllati contro R da R/testing/test_S2.2.R
PRESET_TYPES = {"A1": ["epiteliale", "immunitario", "stromale e muscolo liscio"], "A2": ["epiteliale", "fibroblasti"],
                "A3": ["endoteliale", "epiteliale", "fibroblasti", "immunitario"],
                "A4": ["endoteliali", "linfociti B", "linfociti T", "mieloidi", "plasmacellule", "stromali"],
                "A5": ["astrociti", "microglia", "oligodendrociti", "vascolari", "neuroni"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
MERGE = {"epiteliale normale": "epiteliale", "epiteliale tumorale": "epiteliale"}   # A2, A3: un solo epitelio nel preset (DEC-014)
SEED_M4 = 20261106
ds, archs = DS[KEY]
d_exp = float(pd.read_csv("results/R6/R6_dexp.csv").query("dataset == 'pooled'")["d_exp_um"].iloc[0])
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype.isin(archs)]
t0 = time.time()
X, bc, gnames, gids = read_h5_csc(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
tp = pq.read_table(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet").to_pandas()
col_of = pd.Series(np.arange(len(bc)), index=bc)
print(f"[{KEY}] h5 {X.shape} in {time.time()-t0:.0f}s; d_exp {d_exp}", flush=True)
for _, R in rois.iterrows():
    A, roi, um = R.archetype, R.roi_id, R.um_per_px
    tag = f"{A}_{roi}"
    if os.path.exists(f"{RES}/{tag}_real.json"):
        print(f"[{tag}] già fatto", flush=True); continue
    mask = np.load(f"{M}/{'R3b' if SRC[A]=='consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
    lab6 = pd.read_parquet(f"{M}/R6/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet")
    lab6 = lab6[lab6.roi_id == roi].copy()
    lab6["ptype"] = lab6.cell_type.replace(MERGE)
    ev = lab6[lab6.ptype.isin(PRESET_TYPES[A])]
    box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)].copy()
    box["rr"] = np.floor(box.pxl_row_in_fullres - R.r0).astype(int); box["cc"] = np.floor(box.pxl_col_in_fullres - R.c0).astype(int)
    filt = box[box.barcode.isin(col_of.index)]
    lab = footprint_labels(mask, filt.rr.to_numpy(), filt.cc.to_numpy(), um, d_exp, mode="exp")
    labn = footprint_labels(mask, filt.rr.to_numpy(), filt.cc.to_numpy(), um, d_exp, mode="nuc")
    sel = np.isin(lab, ev.label.to_numpy())
    cols = col_of.reindex(filt.barcode[sel]).to_numpy()
    Xb = X[:, cols].T.tocsr().astype(np.int32)                    # quadratini di valutazione x geni
    rows_, cols_ = filt.array_row.to_numpy()[sel], filt.array_col.to_numpy()[sel]
    lb = lab[sel]; innuc = labn[sel] == lb
    umi_c = pd.Series(np.asarray(Xb.sum(axis=1)).ravel()).groupby(lb).sum()
    nb_c = pd.Series(lb).value_counts()
    cells = ev[["label", "x_um", "y_um", "cluster", "cell_type", "ptype"]].copy()
    cells["umi_eval"] = umi_c.reindex(cells.label).fillna(0).to_numpy(); cells["n_bins_eval"] = nb_c.reindex(cells.label).fillna(0).astype(int).to_numpy()
    cells = cells[cells.n_bins_eval > 0]
    sp.save_npz(f"{OUT}/{tag}_X.npz", Xb)
    np.savez_compressed(f"{OUT}/{tag}_bins.npz", rows=rows_, cols=cols_, label=lb, in_nucleus=innuc)
    cells.to_parquet(f"{OUT}/{tag}_cells.parquet", index=False)
    ctype = dict(zip(cells.label, cells.ptype)); cxy = {l: (x, y) for l, x, y in zip(cells.label, cells.x_um, cells.y_um)}
    met = S.bin_metrics(Xb, rows_, cols_, lb, ctype, cxy, SEED_M4)
    met.update(archetype=A, roi_id=roi, n_cells=int(len(cells)), n_bins_tissue=int(len(filt)), frac_bins_eval=float(sel.mean()),
               n_types=int(cells.ptype.nunique()), types=";".join(f"{t}:{n}" for t, n in cells.ptype.value_counts().items()))
    json.dump(met, open(f"{RES}/{tag}_real.json", "w"), indent=1)
    print(f"[{tag}] bins {Xb.shape[0]} cells {len(cells)} m1 {met['m1']:.2f} m2 {met['m2']:.1f} m3 {met['m3']:.3f} m4r {met['m4_ratio']:.3f} ({time.time()-t0:.0f}s)", flush=True)
