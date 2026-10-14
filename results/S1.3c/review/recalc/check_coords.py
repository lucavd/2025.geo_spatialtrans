
import numpy as np, pandas as pd, os
from scipy import ndimage as ndi
R=os.path.expanduser('~/2025.geo_spatialtrans')
nuc=pd.read_parquet(R+'/results/R2/R2_nuclei_all.parquet')
for A,roi,meth,gm in [('A3','r1','cellpose_rgb','cellpose_rgb'),('A4','f1','stardist_he','spaceranger')]:
    lab=np.load(f'/mnt/micron/geo_spatialtrans/R2/masks/{A}_{roi}_{meth}_native.npz')['labels']
    n=nuc[(nuc.archetype==A)&(nuc.roi_id==roi)&(nuc.method==meth)]
    print(A, meth, n.scale.value_counts().to_dict(), n.keep.sum(), len(n), lab.max())
    n1=n[n.scale==1].head(200)
    cm=np.array(ndi.center_of_mass(np.ones_like(lab), lab, n1.label.values))
    u=n1.um_per_px.iloc[0]
    print(' max |x_um - col*u|', np.abs(n1.x_um.values-(cm[:,1]+0.5)*u).max(), np.abs(n1.x_um.values-cm[:,1]*u).max(), ' y', np.abs(n1.y_um.values-cm[:,0]*u).max())
    g=pd.read_parquet(f'/mnt/micron/geo_spatialtrans/R3/real/{A}_{roi}_{gm}_gen.parquet')
    print(' gen', len(g), g[['x','y']].describe().loc[['min','max']].to_dict())
    if gm==meth:
        mm=g.merge(n[n.scale==1][['label','x_um','y_um','keep','in_tissue','in_bubble']],on='label')
        print(' gen vs nuc', len(mm), np.abs(mm.x-mm.x_um).max(), np.abs(mm.y-mm.y_um).max(), mm.keep.mean())
