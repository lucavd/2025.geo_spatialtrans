"""R6: aggregate 2 um bins onto the real nuclei of the 30 ROIs (footprints exp / nuc / shift) and match nuclei to
Space Ranger 4 cells where available. Usage: python tools/R6_aggregate.py <dataset_key>  (A1|A23|A4|A5|A6)
Outputs: /mnt/micron/geo_spatialtrans/R6/counts/<A>_<roi>_<fp>.npz (+ _genes.tsv once per dataset) and <A>_<roi>_nuclei.parquet"""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow.parquet as pq
from scipy.spatial import cKDTree
sys.path.insert(0, "tools"); from R6_common import footprint_labels, aggregate_counts, read_h5_csc
M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"; OUT = f"{M}/R6/counts"; os.makedirs(OUT, exist_ok=True)
KEY = sys.argv[1]
DS = {"A1": ("Visium_HD_Mouse_Small_Intestine", ["A1"]), "A23": ("Visium_HD_Human_Colon_Cancer", ["A2", "A3"]),
      "A4": ("Visium_HD_Human_Lymph_Node_FFPE", ["A4"]), "A5": ("Visium_HD_6p5mm_Mouse_Brain", ["A5"]), "A6": ("Visium_HD_6p5mm_Human_Heart", ["A6"])}
ds, archs = DS[KEY]
SRC = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
d_exp = float(pd.read_csv("results/R6/R6_dexp.csv").query("dataset == 'pooled'")["d_exp_um"].iloc[0])
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype.isin(archs)]
t0 = time.time()
X, bc, gnames, gids = read_h5_csc(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
pd.DataFrame({"gene_id": gids, "gene_name": gnames}).to_csv(f"{OUT}/{KEY}_genes.tsv", sep="\t", index=False)
tp = pq.read_table(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet").to_pandas()
col_of = pd.Series(np.arange(len(bc)), index=bc)
bmf = f"{D}/{ds}/{ds}_barcode_mappings.parquet"
bm = pq.read_table(bmf, columns=["square_002um", "cell_id", "in_nucleus"]).to_pandas().set_index("square_002um") if os.path.exists(bmf) else None
print(f"[{KEY}] loaded h5 {X.shape} nnz {X.nnz} in {time.time()-t0:.0f}s; d_exp = {d_exp} um", flush=True)
for _, R in rois.iterrows():
    A, roi, um = R.archetype, R.roi_id, R.um_per_px
    tag = f"{A}_{roi}"
    mask = np.load(f"{M}/{'R3b' if SRC[A]=='consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
    nuc = pd.read_parquet(f"{M}/R3b/nuclei/{tag}_consensus_native.parquet" if SRC[A] == "consensus" else f"results/R2/nuclei/{tag}_{SRC[A]}_native.parquet")
    keep_labels = nuc.loc[nuc.keep, "label"].to_numpy()
    box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)].copy()
    box["rr"] = np.floor(box.pxl_row_in_fullres - R.r0).astype(int); box["cc"] = np.floor(box.pxl_col_in_fullres - R.c0).astype(int)
    filt = box[box.barcode.isin(col_of.index)]          # bins present in the filtered matrix (in tissue)
    cols = col_of.reindex(filt.barcode).to_numpy()
    Xb = X[:, cols]
    tab = nuc[nuc.keep][["label", "x_um", "y_um", "area_um2", "eccentricity", "solidity", "in_tissue"]].copy()
    tab["archetype"] = A; tab["roi_id"] = roi
    info = {"archetype": A, "roi_id": roi, "n_bins_box": len(box), "n_bins_filtered": len(filt), "umi_roi": float(Xb.sum()), "n_nuclei_keep": len(keep_labels), "d_exp_um": d_exp}
    for fp in ["exp", "nuc", "shift"]:
        lab = footprint_labels(mask, filt.rr.to_numpy(), filt.cc.to_numpy(), um, d_exp, mode=fp)
        C, nb, umi = aggregate_counts(Xb, lab, keep_labels)
        sp.save_npz(f"{OUT}/{tag}_{fp}.npz", C.astype(np.int32))
        tab[f"n_bins_{fp}"] = nb; tab[f"umi_{fp}"] = umi
        info[f"bins_assigned_{fp}"] = int((lab > 0).sum()); info[f"umi_assigned_{fp}"] = float(umi.sum())
    # match to Space Ranger cell: the 2 um bin containing the nucleus centroid, if Space Ranger puts it in a nucleus
    if bm is not None:
        tree = cKDTree(np.c_[box.pxl_row_in_fullres - R.r0, box.pxl_col_in_fullres - R.c0])
        d, j = tree.query(np.c_[tab.y_um / um, tab.x_um / um])
        bcs = box.barcode.to_numpy()[j]
        m = bm.reindex(bcs)
        tab["sr_bin_dist_px"] = d; tab["sr_in_nucleus"] = m.in_nucleus.fillna(False).to_numpy(); tab["sr_cell_id"] = m.cell_id.to_numpy()
        tab.loc[~tab.sr_in_nucleus, "sr_cell_id"] = np.nan
        info["n_matched_sr"] = int(tab.sr_in_nucleus.sum())
    tab.to_parquet(f"{OUT}/{tag}_nuclei.parquet", index=False)
    json.dump(info, open(f"{OUT}/{tag}_info.json", "w"))
    print(f"[{KEY}] {tag}: bins box {len(box)} filtered {len(filt)} assigned exp {info['bins_assigned_exp']} nuc {info['bins_assigned_nuc']}; umi exp {info['umi_assigned_exp']:.0f}/{info['umi_roi']:.0f}; nuclei {len(keep_labels)}; {time.time()-t0:.0f}s", flush=True)
print(f"[{KEY}] done in {time.time()-t0:.0f}s")
