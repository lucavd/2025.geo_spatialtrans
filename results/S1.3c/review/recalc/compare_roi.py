
import pandas as pd, numpy as np, json, os, re
R=os.path.expanduser('~/2025.geo_spatialtrans/results/S1.3c'); V=R+'/review/recalc'
roi=pd.read_csv(R+'/S1.3c_roi.csv'); roi=roi[roi.variant=='primary']
ck=pd.read_csv(R+'/S1.3c_checkC.csv'); ph1=pd.read_csv(R+'/S1.3c_PH1_levels.csv'); cp=pd.read_csv(R+'/S1.3c_cp_roi.csv')
c5=pd.read_csv(R+'/S1.3c_C5.csv').rename(columns={'Unnamed: 0':'archetype'}); ph3=pd.read_csv(R+'/S1.3c_PH3_calibration.csv')
rows=[]
def add(unit,q,rep,a,b,se=None,c=None):
    rows.append(dict(unit=unit,quantity=q,reported=rep,recalc_rint=a,recalc_floor=b,recalc_px4=c,diff_px4=(c-rep) if c is not None and rep is not None and not pd.isna(rep) else None,diff_rint=(a-rep) if a is not None and rep is not None and not pd.isna(rep) else None,
                     diff_floor=(b-rep) if b is not None and rep is not None and not pd.isna(rep) else None, seed_se=se))
M={m:pd.concat([pd.read_csv(f'{V}/recalc_roi_{m}_{A}.csv') for A in ['A1','A2','A3','A4','A5','A6']]) for m in ('rint_noclip','floor_clip','px4')}
for _,x in M['floor_clip'].iterrows():
    A,ro,k=x.archetype,x.roi_id,int(x.k); y=M['rint_noclip'].query('archetype==@A and roi_id==@ro and k==@k').iloc[0]; z=M['px4'].query('archetype==@A and roi_id==@ro and k==@k').iloc[0]
    u=f'{A}_{ro}_k{k}'; rr=roi.query('archetype==@A and roi_id==@ro and k==@k').iloc[0]
    if k==3:
        m=re.match(r'(\d+) / (\d+)', ck.query('check=="C-1" and unit==@A+"_"+@ro').value.iloc[0])
        add(f'{A}_{ro}','n_bins_roi (C-1)',int(m.group(1)),y.n_bins_roi,x.n_bins_roi,c=z.n_bins_roi)
        p=ph1.query('archetype==@A and roi_id==@ro').iloc[0]
        add(f'{A}_{ro}','n_tissue_bins (PH1)',p.n_tissue_bins,y.n_tissue_bins,x.n_tissue_bins,c=z.n_tissue_bins)
        add(f'{A}_{ro}','umi_tissue_mean (PH1)',p.umi_tissue_mean,y.umi_tissue_mean,x.umi_tissue_mean,c=z.umi_tissue_mean)
        add(f'{A}_{ro}','n_lumen_bins CP-2 (PH1)',p.n_lumen_bins,y.n_lumen_bins,x.n_lumen_bins,c=z.n_lumen_bins)
        add(f'{A}_{ro}','umi_lumen_mean (PH1)',p.umi_lumen_mean,y.umi_lumen_mean,x.umi_lumen_mean,c=z.umi_lumen_mean)
        add(f'{A}_{ro}','n_gen',rr.n_gen,y.n_gen,x.n_gen,c=z.n_gen)
        add(f'{A}_{ro}','R [4,8)/[0,2)',rr.R_profile,y.R,x.R,c=z.R)
    add(u,'thr',rr.thr,y.thr,x.thr,c=z.thr)
    add(u,'n_tissue_windows',rr.n_tissue_windows,y.n_tissue_windows,x.n_tissue_windows,c=z.n_tissue_windows)
    add(u,'n_empty',rr.n_empty,y.n_empty,x.n_empty,c=z.n_empty)
    add(u,'phi_obs',rr.phi_obs,y.phi_obs,x.phi_obs,c=z.phi_obs)
    se=x.phi_null_sd*np.sqrt(1/20+1/30)
    add(u,'phi_null (media)',rr.phi_null,y.phi_null_mean,x.phi_null_mean,se,c=z.phi_null_mean)
    add(u,'phi_null_min',rr.phi_null_min,y.phi_null_min,x.phi_null_min,c=z.phi_null_min)
    add(u,'phi_null_max',rr.phi_null_max,y.phi_null_max,x.phi_null_max,c=z.phi_null_max)
    add(u,'E',rr.E,y.E,x.E,se,c=z.E)
    add(u,'r_star',rr.r_star,y.r_star,x.r_star,c=z.r_star)
    add(u,'BA',rr.BA,y.BA if not pd.isna(y.BA) else None,x.BA if not pd.isna(x.BA) else None,c=z.BA)
    c=cp.query('cp=="CP-2" and archetype==@A and roi_id==@ro and k==@k').iloc[0]
    add(u,'CP-2 n_windows',c.n_windows,y.cp2_n_windows,x.cp2_n_windows,c=z.cp2_n_windows); add(u,'CP-2 n_empty',c.n_empty,y.cp2_n_empty,x.cp2_n_empty,c=z.cp2_n_empty)
C={m:pd.read_csv(f'{V}/recalc_ctrlA_A3r1_{m}.csv') for m in ('rint_noclip','floor_clip','px4')}
c5r=c5.query('archetype=="A3"').iloc[0]; p3=ph3.query('archetype=="A3" and roi_id=="r1" and type=="A" and f_target==0.1').iloc[0]
for mode in ('thr','rank'):
    y=C['rint_noclip'][C['rint_noclip'].zero.str.contains(mode)].iloc[0]; x=C['floor_clip'][C['floor_clip'].zero.str.contains(mode)].iloc[0]
    z=C['px4'][C['px4'].zero.str.contains(mode)].iloc[0]
    u=f'A3_r1_ctrlA_{mode}'
    add(u,'frac_zeroed (C5 A_frac_zeroed)',c5r.A_frac_zeroed,json.loads(y.zero)['frac_zeroed'],json.loads(x.zero)['frac_zeroed'],c=json.loads(z.zero)['frac_zeroed'])
    add(u,'phi_obs (PH3 A f=0.1)',p3.phi_obs,y.phi_obs,x.phi_obs,c=z.phi_obs)
    se=x.phi_null_sd*np.sqrt(1/20+1/30)
    add(u,'phi_null (PH3 A f=0.1)',p3.phi_null,y.phi_null_mean,x.phi_null_mean,se,c=z.phi_null_mean)
    add(u,'E (C5 A_E)',c5r.A_E,y.E,x.E,se,c=z.E); add(u,'E (PH3 A f=0.1)',p3.E,y.E,x.E,se,c=z.E)
    add(u,'BA (C5 A_BA)',c5r.A_BA,y.BA,x.BA,c=z.BA)
T=pd.DataFrame(rows); T.to_csv(V+'/recalc_vs_reported.csv',index=False)
pd.set_option('display.width',250)
print(T.to_string(float_format=lambda v:f'{v:.6g}'))
det=T[~T.quantity.str.contains('phi_null|^E|E \\(', regex=True)]
print('DETERMINISTIC px4 nonzero diffs:'); print(det[(det.diff_px4.abs()>1e-9)].to_string())
for A in ('A2','A6'):
    n=pd.read_csv(f'{V}/recalc_roi_floor_clip_n200_{A}.csv')
    for _,x in n.iterrows():
        rr=roi.query('archetype==@x.archetype and roi_id==@x.roi_id and k==@x.k').iloc[0]
        zz=(rr.phi_null-x.phi_null_mean)/(x.phi_null_sd/np.sqrt(20))
        print('N200', A, int(x.k), 'rep', rr.phi_null, 'mine200', x.phi_null_mean, 'sd', x.phi_null_sd, 'z(rep 20-mean)', round(zz,2), 'range200', x.phi_null_min, x.phi_null_max, 'rep range', rr.phi_null_min, rr.phi_null_max, 'n_win', x.n_tissue_windows)
