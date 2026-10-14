
import pandas as pd, numpy as np, glob, os
R='/home/'+os.environ.get('USER','')+'/2025.geo_spatialtrans'
R=os.path.expanduser('~/2025.geo_spatialtrans')
rois=pd.read_csv(R+'/results/R2/R2_rois_checked.csv')
sel=[('A1','r1'),('A2','r1'),('A3','r1'),('A4','f1'),('A5','r1'),('A6','r1')]
print(rois[[tuple(x) in sel for x in rois[['archetype','roi_id']].values]].to_string())
nuc=pd.read_parquet(R+'/results/R2/R2_nuclei_all.parquet')
print(nuc.columns.tolist(), nuc.shape); print(nuc.head(3).to_string())
print(nuc.groupby(['archetype','method']).size() if 'method' in nuc else '')
z=np.load('/mnt/micron/geo_spatialtrans/R2/masks/A3_r1_cellpose_rgb_native.npz'); print({k:(z[k].shape,z[k].dtype) for k in z.files})
z=np.load('/mnt/micron/geo_spatialtrans/R2/masks/A3_r1_masks.npz'); print({k:(z[k].shape,z[k].dtype) for k in z.files})
g=pd.read_parquet('/mnt/micron/geo_spatialtrans/R3/real/A3_r1_cellpose_rgb_gen.parquet'); print(g.columns.tolist(), g.shape); print(g.head(3).to_string())
print(sorted(os.path.basename(p) for p in glob.glob('/mnt/micron/geo_spatialtrans/R3/real/A4_f1*')+glob.glob('/mnt/micron/geo_spatialtrans/R3/real/A6_r1*')))
print(os.listdir('/mnt/micron/geo_spatialtrans/data_real/'))
print(os.listdir(R+'/results/R2/tissue_masks')[:5])
from skimage import io
v=io.imread(R+'/results/R2/tissue_masks/A3_r1_valid_ds4.png'); print(v.shape,v.dtype,np.unique(v))
