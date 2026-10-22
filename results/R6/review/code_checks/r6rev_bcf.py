
import sys, os, json, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow.parquet as pq, h5py
from scipy import ndimage
from scipy.spatial import cKDTree
M="/mnt/micron/geo_spatialtrans"; D=f"{M}/data_real"
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
SRC = {"A1": "cellpose_rgb", "A5": "cellpose_rgb", "A6": "consensus"}
d_exp = 8.0
def read_h5(path):
    with h5py.File(path, "r") as h:
        m = h["matrix"]; shape = tuple(int(v) for v in m["shape"][:])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape); bc = m["barcodes"][:].astype(str)
    return X, bc
for A, roi in [("A1","r1"), ("A5","r4"), ("A6","r1")]:
    ds = DS[A]; R = rois[(rois.archetype==A)&(rois.roi_id==roi)].iloc[0]; um = R.um_per_px; tag=f"{A}_{roi}"
    print(f"\n===== {tag} ({ds}) um_per_px={um:.6f} d_exp={d_exp}")
    X, bc = read_h5(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
    tp = pq.read_table(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet").to_pandas()
    print("tissue_positions dtypes:", tp.pxl_row_in_fullres.dtype, "example", tp.pxl_row_in_fullres.iloc[0], tp.pxl_col_in_fullres.iloc[0])
    mask = np.load(f"{M}/{'R3b' if SRC[A]=='consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
    nuc = pd.read_parquet(f"{M}/R3b/nuclei/{tag}_consensus_native.parquet" if SRC[A]=="consensus" else f"results/R2/nuclei/{tag}_{SRC[A]}_native.parquet")
    tab = pd.read_parquet(f"{M}/R6/counts/{tag}_nuclei.parquet"); info = json.load(open(f"{M}/R6/counts/{tag}_info.json"))
    print("mask shape", mask.shape, "dtype", mask.dtype, "box rows", R.r1-R.r0, "cols", R.c1-R.c0)
    labs_mask = set(np.unique(mask)) - {0}; print("labels in mask:", len(labs_mask), "nuclei table:", len(nuc), "keep:", int(nuc.keep.sum()), "mask labels == table labels:", labs_mask == set(nuc.label))
    keep = nuc.loc[nuc.keep, "label"].to_numpy(); nonkeep = nuc.loc[~nuc.keep, "label"].to_numpy()
    # --- (b) centroid convention: center_of_mass(mask) (row, col) vs (y_um/um, x_um/um)
    com = ndimage.center_of_mass(mask > 0, mask, nuc.label.to_numpy()); com = np.array(com)
    d_std = np.hypot(com[:,0]*um - nuc.y_um.values, com[:,1]*um - nuc.x_um.values)
    d_swp = np.hypot(com[:,0]*um - nuc.x_um.values, com[:,1]*um - nuc.y_um.values)
    print(f"(b) centroid: median |COM(row,col)*um - (y_um,x_um)| = {np.median(d_std):.4f} um (max {d_std.max():.4f}); swapped axes: median {np.median(d_swp):.2f} um")
    # bins in box
    box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)].copy()
    box["rr"] = np.floor(box.pxl_row_in_fullres - R.r0).astype(int); box["cc"] = np.floor(box.pxl_col_in_fullres - R.c0).astype(int)
    colof = pd.Series(np.arange(len(bc)), index=bc); filt = box[box.barcode.isin(colof.index)].copy()
    print("bins box", len(box), "(info", info["n_bins_box"], ") filtered", len(filt), "(info", info["n_bins_filtered"], ")")
    rr = filt.rr.to_numpy(); cc = filt.cc.to_numpy()
    lab_nuc = mask[rr, cc].astype(np.int64)
    # mean bin position per nucleus (nuc mode) vs centroid
    df = pd.DataFrame({"lab": lab_nuc, "x": cc*um, "y": rr*um}); g = df[df.lab>0].groupby("lab")[["x","y"]].mean()
    g = g.join(nuc.set_index("label")[["x_um","y_um"]], how="inner")
    dd = np.hypot(g.x - g.x_um, g.y - g.y_um); dd_sw = np.hypot(g.x - g.y_um, g.y - g.x_um)
    print(f"(b) nuc-mode bins: nuclei with >=1 bin {len(g)}; median dist mean-bin-position vs centroid = {np.median(dd):.3f} um, q90 {np.quantile(dd,.9):.3f}, max {dd.max():.2f}; swapped axes median {np.median(dd_sw):.1f} um")
    # --- (c) independent exp footprint: KD-tree on nucleus pixels
    pr, pc = np.nonzero(mask); tree = cKDTree(np.c_[pr, pc]); plab = mask[pr, pc].astype(np.int64)
    amb = lab_nuc == 0
    dist, j = tree.query(np.c_[rr[amb], cc[amb]])
    near = plab[j]; near[dist*um > d_exp] = 0
    lab_exp = lab_nuc.copy(); lab_exp[amb] = near
    print(f"(c) exp (independent): bins assigned {int((lab_exp>0).sum())} (info {info['bins_assigned_exp']}); nuc: {int((lab_nuc>0).sum())} (info {info['bins_assigned_nuc']})")
    print(f"(c) expanded ambient bins: max dist*um = {np.max(dist[near>0]*um) if (near>0).any() else np.nan:.3f} um (<= d_exp)")
    nk = np.isin(lab_exp, nonkeep); print(f"(c) bins assigned to keep=False nuclei (exp): {int(nk.sum())} = info bins_assigned_exp - sum(n_bins_exp parquet) = {info['bins_assigned_exp'] - int(tab.n_bins_exp.sum())}")
    # contrast: remove non-keep nuclei from the mask -> expansion no longer blocked
    mask2 = mask.copy(); mask2[np.isin(mask2, nonkeep)] = 0
    pr2, pc2 = np.nonzero(mask2); tree2 = cKDTree(np.c_[pr2, pc2]); plab2 = mask2[pr2, pc2].astype(np.int64)
    lab_nuc2 = mask2[rr, cc].astype(np.int64); amb2 = lab_nuc2 == 0
    dist2, j2 = tree2.query(np.c_[rr[amb2], cc[amb2]]); near2 = plab2[j2]; near2[dist2*um > d_exp] = 0
    lab_exp2 = lab_nuc2.copy(); lab_exp2[amb2] = near2
    print(f"(c) contrast without keep=False in mask: bins to keep nuclei {int(np.isin(lab_exp2, keep).sum())} vs with blocking {int(np.isin(lab_exp, keep).sum())} (sum n_bins_exp parquet {int(tab.n_bins_exp.sum())})")
    # --- (a) independent aggregation vs npz
    cols = colof.reindex(filt.barcode).to_numpy(); Xb = X[:, cols].tocsc()
    pos = pd.Series(np.arange(len(tab)), index=tab.label.values)
    for name, L in [("nuc", lab_nuc), ("exp", lab_exp)]:
        sel = np.isin(L, tab.label.values); idx_n = pos.reindex(L[sel]).to_numpy(); idx_b = np.flatnonzero(sel)
        I = sp.csr_matrix((np.ones(len(idx_b)), (idx_n, idx_b)), shape=(len(tab), Xb.shape[1])); C = (I @ Xb.T.tocsr()).tocsr()
        Cs = sp.load_npz(f"{M}/R6/counts/{tag}_{name}.npz")
        diff = (C - Cs); nz = diff.nnz; rowsdiff = int((np.asarray(abs(diff).sum(axis=1)).ravel() > 0).sum())
        print(f"(a) {name}: independent counts vs npz: shape {C.shape} vs {Cs.shape}; nonzero differences {nz}; nuclei with any difference {rowsdiff}; umi parquet == independent row sums: {np.array_equal(np.asarray(C.sum(axis=1)).ravel(), tab[f'umi_{name}'].values)}")
    # --- (f) Space Ranger matching
    bmf = f"{D}/{ds}/{ds}_barcode_mappings.parquet"
    if os.path.exists(bmf):
        bm = pq.read_table(bmf, columns=["square_002um","cell_id","in_nucleus","in_cell"]).to_pandas()
        print("(f) barcode_mappings rows", len(bm), "unique square_002um", bm.square_002um.nunique(), "cell_id example", bm.cell_id.dropna().iloc[0], "in_nucleus dtype", bm.in_nucleus.dtype)
        bm = bm.set_index("square_002um")
        treeb = cKDTree(np.c_[box.pxl_row_in_fullres - R.r0, box.pxl_col_in_fullres - R.c0])
        d, jj = treeb.query(np.c_[tab.y_um/um, tab.x_um/um]); bcs = box.barcode.to_numpy()[jj]; m = bm.reindex(bcs)
        in_n = m.in_nucleus.fillna(False).to_numpy().astype(bool); cid = m.cell_id.to_numpy()
        print(f"(f) sr_in_nucleus agree: {(in_n == tab.sr_in_nucleus.values).mean():.4f}; sr_cell_id agree (where in_nucleus): {(pd.Series(cid[in_n]) == pd.Series(tab.sr_cell_id.values[in_n])).mean():.4f}; n in_nucleus {in_n.sum()} (parquet {int(tab.sr_in_nucleus.sum())})")
        print(f"(f) sr_bin_dist_px: median {np.median(tab.sr_bin_dist_px):.2f}, max {tab.sr_bin_dist_px.max():.2f} px; bin pitch px = {2/um:.2f}; half-diagonal = {np.sqrt(2)/um:.2f}")
        # bin containing the centroid (floor to grid) vs nearest centre: compare via array_row/col
        bx = box.set_index("barcode"); near_bins = bx.loc[bcs]
        # the centroid in pixels of the box; the bin whose [row, row+pitch) contains the centroid
        print(f"(f) fraction of nuclei whose centroid is within half-pitch of the matched bin centre: {np.mean(d <= 1.0/um + 1e-6):.4f}")
        # in_cell but not in_nucleus
        print(f"(f) nearest bin in_cell but not in_nucleus: {int((m.in_cell.fillna(False).to_numpy().astype(bool) & ~in_n).sum())}; in_nucleus with cell_id NaN: {int((in_n & pd.isna(cid)).sum())}")
        ct = pd.read_csv(f"{D}/{ds}/segmented_outputs_x/segmented_outputs/cell_types/10x_Cloud/cell_types.csv")
        print("(f) cell_types.csv rows", len(ct), "unique barcode", ct.barcode.nunique(), "example", ct.barcode.iloc[0], "| matched sr_cell_id found in cell_types:", pd.Series(tab.sr_cell_id.dropna().unique()).isin(ct.barcode).mean().round(4))
        print("(f) coarse types:", ct.coarse_cell_type.value_counts().head(12).to_dict())
