import pandas as pd, numpy as np
D="results/R3/review/code"; R="/mnt/micron/geo_spatialtrans/R3/real"
P={'A1':'cellpose_rgb','A2':'cellpose_rgb','A3':'cellpose_rgb','A4':'spaceranger','A5':'cellpose_rgb','A6':'spaceranger'}
def axd(a,b):
    d=np.abs(a-b)%np.pi; return np.minimum(d,np.pi-d)
out=[]
for A,roi in [("A1","r1"),("A2","r4"),("A3","r5"),("A4","p1"),("A5","r2"),("A6","r3"),("A6","r1")]:
    me=P[A]; mine=pd.read_parquet(f"{D}/cells/{A}_{roi}_{me}.parquet")
    tv=pd.read_parquet(f"{R}/{A}_{roi}_{me}_cells.parquet"); nu=pd.read_parquet(f"{R}/{A}_{roi}_{me}_nuc.parquet")
    th=tv.merge(nu[["idx","label","ecc_N","theta_N","frac_out"]],on=["idx","label"])
    mm=mine.merge(th,on="label",suffixes=("","_R3"))
    I=mm.interior&mm.interior_R3
    fc_m=np.mean(mine.frac_out[mine.interior]>0.05); fc_t=np.mean(th.frac_out[th.interior]>0.05)
    out.append(dict(A=A,roi=roi,n_mine=len(mine),n_R3=len(th),n_matched=len(mm),interior_disagree=int((mm.interior!=mm.interior_R3).sum()),
        max_d_area=float(np.abs(mm.area-mm.area_R3).max()),max_d_eccT=float(np.abs(mm.ecc_T-mm.ecc_T_R3)[I].max()),
        max_d_thetaT_deg=float(np.degrees(axd(mm.theta_T,mm.theta_T_R3)[I&(mm.ecc_T>0.05)]).max()),
        max_d_eccN=float(np.nanmax(np.abs(mm.ecc_N-mm.ecc_N_R3))),
        max_d_thetaN_deg=float(np.nanmax(np.degrees(axd(mm.theta_N,mm.theta_N_R3))[mm.ecc_N>0.2])),
        max_d_fracout=float(np.nanmax(np.abs(mm.frac_out-mm.frac_out_R3))), n_fracout_diff=int((np.abs(mm.frac_out-mm.frac_out_R3)>1e-12).sum()),
        frac_cut_mine=fc_m, frac_cut_R3=fc_t, nsides_disagree=int((mm.nsides!=mm.nsides_R3)[I].sum())))
    if A=="A6" and roi=="r3":
        extra=mine[~mine.label.isin(th.label)]
        print("A6 r3 extra in mine:\n", extra[["label","x","y"]].to_string())
        px=999.9/913
pd.set_option("display.width",250)
o=pd.DataFrame(out); print(o.to_string()); o.to_csv(f"{D}/rc_compare_cells.csv",index=False)
