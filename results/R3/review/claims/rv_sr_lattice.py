
import numpy as np, pandas as pd
c=pd.read_parquet("/mnt/micron/geo_spatialtrans/R3/cells_all.parquet",columns=["archetype","roi_id","method","x","y","theta_T","ecc_T","interior"])
for A,me in [("A4","spaceranger"),("A4","stardist_he"),("A6","spaceranger")]:
    z=c[(c.archetype==A)&(c.method==me)&c.interior]
    fx=np.mod(z.x,0.5); fy=np.mod(z.y,0.5)
    th=np.degrees(np.abs(z.theta_T)); ax=np.minimum(th, np.abs(90-th))  # distance from 0/90
    print(A,me,"x mod 0.5 unique(0.01):",np.unique(np.round(fx,2)).size,"| theta_T within 10deg of 0/90:",(ax<10).mean().round(3),"(uniform 0.444)")

