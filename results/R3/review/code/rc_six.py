import pandas as pd, numpy as np
C=pd.read_csv("results/R3/review/code/rc_compare_roi.csv")
sel=[("A1","r1","cellpose_rgb"),("A2","r4","cellpose_rgb"),("A3","r5","cellpose_rgb"),("A4","p1","spaceranger"),("A5","r2","cellpose_rgb"),("A6","r3","spaceranger")]
mets=["n","sum_area","frac_interior","median_area","median_eq_r","cv","median_ecc_T","median_nc","n_eligible","median_dtheta","cv_loc"]
J=pd.read_csv("results/R3/review/code/rc_roi_all.csv"); J=J[J.jitter_um==0]
ca=pd.read_parquet("/mnt/micron/geo_spatialtrans/R3/cells_all.parquet"); ci=ca[ca.interior]
rows=[]
for A,r,m in sel:
    d=C[(C.archetype==A)&(C.roi_id==r)&(C.method==m)].set_index("metric")
    row={"ROI":f"{A} {r} {m}"}
    for k in mets: row[k]=f"{d.loc[k,'declared']:.6g} / {d.loc[k,'recomputed']:.6g}"
    t=ci[(ci.archetype==A)&(ci.roi_id==r)&(ci.method==m)]; j=J[(J.archetype==A)&(J.roi_id==r)&(J.method==m)].iloc[0]
    row["frac_cut"]=f"{np.mean(t.frac_out>0.05):.6g} / {j.frac_cut:.6g}"
    rows.append(row)
o=pd.DataFrame(rows); o.to_csv("results/R3/review/code/rc_six_roi_table.csv",index=False)
print(o.T.to_string())
