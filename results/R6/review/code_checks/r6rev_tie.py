
import numpy as np, pandas as pd, pyarrow.parquet as pq, os, sys, time
from scipy.spatial import cKDTree
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0,"tools"); from R6_common import footprint_labels
M="/mnt/micron/geo_spatialtrans"; D=f"{M}/data_real"
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
for A, roi, ds in [("A1","r1","Visium_HD_Mouse_Small_Intestine"), ("A6","r1","Visium_HD_6p5mm_Human_Heart")]:
    R = rois[(rois.archetype==A)&(rois.roi_id==roi)].iloc[0]; um=R.um_per_px
    mask = np.load(f"{M}/{'R3b' if A=='A6' else 'R2'}/masks/{A}_{roi}_{'consensus' if A=='A6' else 'cellpose_rgb'}_native.npz")["labels"]
    tp = pq.read_table(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet", columns=["pxl_row_in_fullres","pxl_col_in_fullres","in_tissue"]).to_pandas()
    box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)]
    rr = np.floor(box.pxl_row_in_fullres.values - R.r0).astype(int); cc = np.floor(box.pxl_col_in_fullres.values - R.c0).astype(int)
    L = footprint_labels(mask, rr, cc, um, 8.0, mode="exp")
    pr, pc = np.nonzero(mask); tree = cKDTree(np.c_[pr,pc]); plab = mask[pr,pc].astype(np.int64)
    amb = mask[rr,cc]==0; dist, j = tree.query(np.c_[rr[amb],cc[amb]]); near = plab[j]; near[dist*um>8.0]=0
    L2 = mask[rr,cc].astype(np.int64); L2[amb]=near
    diff = np.flatnonzero(L != L2)
    print(f"{A}_{roi}: box bins {len(rr)}; exp labels differ on {len(diff)} bins ({len(diff)/len(rr):.4%}); assigned both {(L>0).sum()} vs {(L2>0).sum()}", flush=True)
    ties = 0; dd=[]
    for i in diff[:1500]:
        a, b = L[i], L2[i]
        if a==0 or b==0: dd.append(np.nan); continue
        pa = np.c_[np.nonzero(mask==a)]; pb = np.c_[np.nonzero(mask==b)]
        da = np.min(np.hypot(pa[:,0]-rr[i], pa[:,1]-cc[i])); db = np.min(np.hypot(pb[:,0]-rr[i], pb[:,1]-cc[i]))
        dd.append(da-db); ties += abs(da-db) < 1e-9
    dd = np.array(dd); print(f"   among first {min(1500,len(diff))} differing bins: exact distance ties {ties}; one side 0 (d_exp boundary) {int(np.isnan(dd).sum())}; max |d_EDT - d_KD| = {np.nanmax(np.abs(dd)) if len(dd) else 0:.3g} px", flush=True)
