
import json, os, pandas as pd, numpy as np, h5py
R=os.path.expanduser('~/2025.geo_spatialtrans/results/S1.3c')
ds='Visium_HD_Mouse_Small_Intestine'
me=json.load(open(f'{R}/review/recalc/recalc_bg_{ds}.json')); rep=json.load(open(f'{R}/S1.3c_bg_{ds}_bg.json'))
print('mine', me['n_umi_missing'], me['n_bins'], me['bg']['com']); print('rep', rep['bg']['com'], rep['log'])
m=pd.read_parquet(f'/tmp/recalc_s13c/{ds}_bins.parquet')
na=m.umi_com.isna(); print('NaN', na.sum(), m[na].head(3).to_string(), m[na].in_tissue.value_counts().to_dict())
p=f'/mnt/micron/geo_spatialtrans/data_real/{ds}/binned_outputs_x/binned_outputs/square_002um/raw_feature_bc_matrix.h5'
with h5py.File(p) as f: bc=f['matrix/barcodes'][:]; print('raw n bc', len(bc), bc[:2])
print('tp n', len(m), m.barcode.head(2).tolist(), 'grid', m.array_row.max()+1, m.array_col.max()+1)
