
import numpy as np, pandas as pd
from skimage import io, measure
for A,r in [("A1","r1"),("A5","r1"),("A5","r3"),("A6","r1"),("A6","r3")]:
    v=io.imread(f"results/R2/tissue_masks/{A}_{r}_valid_ds4.png"); v=(v[...,0] if v.ndim==3 else v)>127
    t=io.imread(f"/mnt/micron/geo_spatialtrans/R3/posthoc/masks_tissue/{A}_{r}_tissue_ds4.png"); t=(t[...,0] if t.ndim==3 else t)>127
    hv=measure.label(~v,connectivity=1); ht=measure.label(~t,connectivity=1)
    print(A,r,"valid px",v.sum(),"tissue px",t.sum(),"valid&~tissue",(v&~t).sum(),"holes(valid)",hv.max(),"holes(tissue)",ht.max())
c=pd.read_parquet("/mnt/micron/geo_spatialtrans/R3/cells_all.parquet",columns=["archetype","method","theta_N","ecc_N","area_nuc"])
for A in ["A4","A6"]:
    z=c[(c.archetype==A)&(c.method=="spaceranger")&(c.ecc_N>=.8)]
    print(A, (np.round(np.degrees(z.theta_N),1)).value_counts().head(6).to_dict(), "area_nuc median", z.area_nuc.median())
    zz=c[(c.archetype==A)&(c.method=="spaceranger")]; print(A,"SR frac ecc>=0.8", (zz.ecc_N>=.8).mean().round(3), "SD frac", (c[(c.archetype==A)&(c.method=="stardist_he")].ecc_N>=.8).mean().round(3))

