
import pandas as pd, numpy as np, os
V=os.path.expanduser('~/2025.geo_spatialtrans/results/S1.3c/review/recalc')
T=pd.read_csv(V+'/recalc_vs_reported.csv')
n=T[T.quantity.str.contains('phi_null \\(|^E', regex=True)].copy()
n['z_px4']=n.diff_px4/n.seed_se
pd.set_option('display.width',250)
print(n[['unit','quantity','reported','recalc_px4','diff_px4','seed_se','z_px4']].to_string(float_format=lambda v:f'{v:.4g}'))
print('max |z|', n.z_px4.abs().max())
# n righe deterministiche, quante identiche in px4
d=T[~T.quantity.str.contains('phi_null|^E|E \\(', regex=True) & ~T.unit.str.contains('rank')]
print('deterministic rows', len(d), 'px4 identical', int((d.diff_px4.abs()<1e-9).sum() + d.reported.isna().sum()), 'NaN reported', int(d.reported.isna().sum()))
print('floor_clip identical', int((d.diff_floor.abs()<1e-9).sum()), 'rint identical', int((d.diff_rint.abs()<1e-9).sum()))
# max abs diff (floor_clip, rint) per quantity type vs reported
for q in ['n_tissue_bins (PH1)','n_lumen_bins CP-2 (PH1)','n_tissue_windows','phi_obs','R [4,8)/[0,2)','BA','CP-2 n_windows','CP-2 n_empty']:
    x=d[d.quantity==q]; print(q, 'max|diff| floor', x.diff_floor.abs().max(), 'rint', x.diff_rint.abs().max())
