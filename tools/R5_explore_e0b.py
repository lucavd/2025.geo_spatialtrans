
"""R5 exploration E0b (geometry only): do the H&E images cover the out-of-tissue bins? Raw matrix size."""
import numpy as np, pandas as pd, h5py, tifffile, sys
DS = ["Visium_HD_Mouse_Small_Intestine","Visium_HD_Human_Colon_Cancer","Visium_HD_Human_Lymph_Node_FFPE","Visium_HD_6p5mm_Mouse_Brain","Visium_HD_6p5mm_Human_Heart"]
IMG = {"Visium_HD_Human_Lymph_Node_FFPE": "tissue_image.tif"}
rows=[]
for ds in DS:
    base=f"data/real/datasets/{ds}"
    t=pd.read_parquet(f"{base}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet")
    fn=f"{base}/{ds}_"+IMG.get(ds,"tissue_image.btf")
    with tifffile.TiffFile(fn) as tf: H,W=tf.pages[0].shape[:2]
    inside=(t.pxl_row_in_fullres>=0)&(t.pxl_row_in_fullres<H)&(t.pxl_col_in_fullres>=0)&(t.pxl_col_in_fullres<W)
    with h5py.File(f"{base}/binned_outputs_x/binned_outputs/square_002um/raw_feature_bc_matrix.h5","r") as h: nraw=h["matrix/shape"][1]
    out=t.in_tissue==0
    rows.append(dict(dataset=ds, img=f"{H}x{W}", n_bins=len(t), n_raw=int(nraw), frac_all_in_img=round(inside.mean(),4),
                     frac_out_tissue_in_img=round(inside[out].mean(),4), cols=",".join(t.columns)))
df=pd.DataFrame(rows); df.to_csv(sys.argv[1],index=False); print(df.drop(columns="cols").to_string()); print(rows[0]["cols"])
