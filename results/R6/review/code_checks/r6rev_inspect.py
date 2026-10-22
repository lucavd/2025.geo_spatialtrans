
import pandas as pd, numpy as np, scipy.sparse as sp, glob, os, pyarrow.parquet as pq
M="/mnt/micron/geo_spatialtrans"
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
print(rois.columns.tolist())
print(rois[["archetype","roi_id","r0","r1","c0","c1","um_per_px"]].to_string())
f = f"{M}/R6/counts/A1_"+rois.query("archetype=='A1'").roi_id.iloc[0]+"_nuclei.parquet"
t = pd.read_parquet(f); print(f, t.shape); print(t.dtypes); print(t.head(3).to_string())
lab = pd.read_parquet(f"{M}/R6/labels/A1_exp_umi20_seed20261022_nuclei.parquet"); print(lab.dtypes); print(lab.head(3).to_string())
lab1 = pd.read_parquet(f"{M}/R6/labels/A1_exp_umi20_seed20261022_nuclei_ph1.parquet"); print(lab1.columns.tolist())
print(os.listdir(f"{M}/R6/terr")[:5]); tt = pd.read_parquet(glob.glob(f"{M}/R6/terr/A1_*.parquet")[0]); print(tt.dtypes); print(tt.head(2).to_string())
print(pq.read_schema(f"{M}/data_real/Visium_HD_6p5mm_Mouse_Brain/Visium_HD_6p5mm_Mouse_Brain_barcode_mappings.parquet"))
print(os.path.getsize(f"{M}/data_real/Visium_HD_6p5mm_Mouse_Brain/Visium_HD_6p5mm_Mouse_Brain_barcode_mappings.parquet")/1e6,"MB")
nuc = pd.read_parquet("results/R2/nuclei/A1_"+rois.query("archetype=='A1'").roi_id.iloc[0]+"_cellpose_rgb_native.parquet"); print(nuc.dtypes); print(nuc.keep.value_counts()); print(nuc.head(2).to_string())
