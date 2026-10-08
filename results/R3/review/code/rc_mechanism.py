import pandas as pd, numpy as np
from scipy.spatial import cKDTree
pd.set_option("display.width",250)
D="results/R3/review/code/cells"
def axd(a,b):
    d=np.abs(a-b)%np.pi; return np.minimum(d,np.pi-d)
rows=[]
for tag in ["A1_r1_cellpose_rgb","A2_r4_cellpose_rgb","A3_r5_cellpose_rgb","A4_p1_spaceranger","A4_p1_stardist_he","A5_r2_cellpose_rgb","A6_r3_spaceranger","A6_r3_stardist_he"]:
    d=pd.read_parquet(f"{D}/{tag}.parquet")
    T=cKDTree(d[["x","y"]].values); dist,ix=T.query(d[["x","y"]].values,k=7)
    vx=d.x.values[ix[:,1]]-d.x.values; vy=d.y.values[ix[:,1]]-d.y.values
    phi_nn=np.degrees(axd(np.arctan2(vy,vx),d.theta_N.values))   # angolo fra direzione del 1° vicino e asse maggiore del nucleo
    # tutti i 6 vicini: angolo medio pesato
    el=(d.interior&(d.ecc_N>=0.8)&(d.ecc_T>=0.5)).values
    dth=np.degrees(axd(d.theta_T.values,d.theta_N.values))
    # lunghezza dell'asse maggiore del nucleo ~ sqrt(area*4/pi / sqrt(1-e^2)) (ellisse equivalente)
    b_over_a=np.sqrt(1-np.clip(d.ecc_N.values,0,0.9999)**2); a_len=2*np.sqrt(d.area_nuc.values/(np.pi*b_over_a))
    gap=dist[:,1]/a_len
    t=np.nanquantile(gap[el],[1/3,2/3]); g1=el&(gap<=t[0]); g3=el&(gap>t[1])
    # contrasto: nuclei con ecc_N bassa (<0.5) -> direzione del 1° vicino rispetto all'asse (dovrebbe essere ~45)
    lo=d.interior.values&(d.ecc_N.values<0.5)
    rows.append(dict(tag=tag,n_el=int(el.sum()),med_dth=np.median(dth[el]),med_phi_nn_el=np.median(phi_nn[el]),frac_phi_nn_gt60=np.mean(phi_nn[el]>60),
        med_phi_nn_round=np.median(phi_nn[lo]),med_dth_closest_tertile=np.median(dth[g1]),med_dth_farthest_tertile=np.median(dth[g3]),
        gap_tertiles=f"{t[0]:.2f}/{t[1]:.2f}",nc_med=np.median((d.area_nuc/d.area)[d.interior]),frac_cut=np.mean(d.frac_out[d.interior]>0.05),
        frac_nn_lt_sum_eqr=np.mean(dist[:,1]<2*np.sqrt(d.area_nuc.values/np.pi)),
        theta_T_near_axes=np.mean(np.minimum(np.degrees(d.theta_T[el])%90,90-np.degrees(d.theta_T[el])%90)<5),
        theta_N_near_axes=np.mean(np.minimum(np.degrees(d.theta_N[el])%90,90-np.degrees(d.theta_N[el])%90)<5)))
o=pd.DataFrame(rows); print(o.round(3).to_string()); o.to_csv("results/R3/review/code/rc_alignment_mechanism.csv",index=False)
