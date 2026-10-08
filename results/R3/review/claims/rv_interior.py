
import pandas as pd, numpy as np
RV="results/R3/review/claims"; OUT="/mnt/micron/geo_spatialtrans/R3"
cp1=pd.read_csv(f"{RV}/rv_cp1_roi.csv"); p4=pd.read_csv("results/R3/R3_posthoc_P4_tissue_only.csv")
m=cp1.merge(p4[["archetype","roi_id","n","frac_interior","median_eq_r"]].rename(columns={"n":"n_p4","frac_interior":"fi_p4","median_eq_r":"eqr_p4"}))
ns=pd.read_csv("results/R3/R3_null_summary.csv"); rr=pd.read_csv("results/R3/R3_roi_summary.csv")
lam=rr.set_index(["archetype","roi_id","method"])
cs=ns[ns.model=="CSR"].groupby(["archetype","roi_id"]).agg(csr_meanA=("mean_area","median"),csr_fi=("frac_interior","median"),csr_eqr=("median_eq_r","median")).reset_index()
m=m.merge(cs)
m["csr_meanA_x_lam"]=[r.csr_meanA*lam.loc[(r.archetype,r.roi_id,"cellpose_rgb" if r.archetype not in("A4","A6") else "spaceranger")].lambda_mm2/1e6 for r in m.itertuples()]
print(m[["archetype","roi_id","n_real","n_p4","frac_int_real","fi_p4","frac_int_csr","frac_int_rsa","csr_meanA_x_lam"]].round(3).to_string())
e=pd.read_csv(f"{RV}/rv_P2_by_ecc.csv"); print(e.round(2).to_string())
c=pd.read_parquet(f"{OUT}/cells_all.parquet",columns=["archetype","roi_id","method","x","y","area","interior","theta_N","ecc_N"])
for A,me in [("A4","spaceranger"),("A4","stardist_he"),("A1","cellpose_rgb"),("A6","spaceranger")]:
    z=c[(c.archetype==A)&(c.method==me)&(c.ecc_N>=.8)]
    print(A,me,"unique theta_N (deg, 0.01 rounding):",np.unique(np.round(np.degrees(z.theta_N),2)).size,"of",len(z), "top5 share:", (np.round(np.degrees(z.theta_N),2).value_counts().head(5).sum()/len(z)).round(3))
# eq_r of territories away from the ROI frame (tissue-mask clipping kept as real geometry)
rois=pd.read_csv("results/R2/R2_rois_checked.csv"); rois["side"]=rois.side_px*rois.um_per_px
PRIM=dict(A1="cellpose_rgb",A2="cellpose_rgb",A3="cellpose_rgb",A4="spaceranger",A5="cellpose_rgb",A6="spaceranger")
out=[]
for A in PRIM:
    z=c[(c.archetype==A)&(c.method==PRIM[A])].merge(rois[["archetype","roi_id","side"]])
    mg=60 if A in("A5","A6") else 25
    far=(z.x>mg)&(z.x<z.side-mg)&(z.y>mg)&(z.y<z.side-mg)
    out.append(dict(archetype=A,margin=mg,eqr_med_interior=np.median(np.sqrt(z.area[z.interior]/np.pi)),eqr_med_far_all=np.median(np.sqrt(z.area[far]/np.pi)),
               frac_far_interior=z.interior[far].mean(),n_far=far.sum()))
o=pd.DataFrame(out); print(o.round(3).to_string()); o.to_csv(f"{RV}/rv_eqr_far_from_frame.csv",index=False)
m.to_csv(f"{RV}/rv_interior_fraction.csv",index=False)

