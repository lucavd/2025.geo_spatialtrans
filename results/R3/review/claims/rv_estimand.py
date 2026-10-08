
import numpy as np, pandas as pd
RV="results/R3/review/claims"
print(pd.read_csv(f"{RV}/rv_B34_spearman_roi.csv").round(3).to_string())
v=pd.read_csv(f"{RV}/rv_verdicts_recomputed.csv"); print(v[v.id.isin(["CP-1a","CP-2b","CP-3","CP-4"])][["id","archetype","note"]].to_string())
PRIM=dict(A1="cellpose_rgb",A2="cellpose_rgb",A3="cellpose_rgb",A4="spaceranger",A5="cellpose_rgb",A6="spaceranger")
c=pd.read_parquet("/mnt/micron/geo_spatialtrans/R3/cells_all.parquet",columns=["archetype","roi_id","method","area","area_nuc","interior"])
rr=pd.read_csv("results/R3/R3_roi_summary.csv")
out=[]
for A in PRIM:
    z=c[(c.archetype==A)&(c.method==PRIM[A])]; r=rr[(rr.archetype==A)&(rr.method==PRIM[A])]
    lam=r.n.sum()/r.area_poly_um2.sum()
    zi=z[z.interior]
    out.append(dict(archetype=A, sqrt_meanAnuc_x_lambda=np.sqrt(z.area_nuc.mean()*lam), median_ratio_interior=np.median(np.sqrt(zi.area_nuc/zi.area)),
                    sqrt_ratio_of_medians=np.sqrt(zi.area_nuc.median()/zi.area.median()), nc_of_means=zi.area_nuc.mean()/zi.area.mean(), nc_median=np.median(zi.area_nuc/zi.area)))
o=pd.DataFrame(out); print(o.round(3).to_string()); o.to_csv(f"{RV}/rv_B32_estimand.csv",index=False)

