
import sys, numpy as np, pandas as pd, importlib.util
sys.argv=['x','edge']
spec=importlib.util.spec_from_file_location('rv','results/S1.3c/review/method/rev_method.py'); RV=importlib.util.module_from_spec(spec); spec.loader.exec_module(RV)
E=RV.E; rois=E.rois_table(); out=[]
for a,roi in (('A1','r1'),('A6','r1'),('A5','r1')):
    r=rois[(rois.archetype==a)&(rois.roi_id==roi)].iloc[0]
    s,valid4,px4,vb,nm,dmap,dnuc,terr,ng=RV.setup(r); u=s.umi.values.astype(float); idx=np.where(vb)[0]
    Rs=[]
    for p in range(1,6):
        u2=u.copy(); u2[idx]=E.permute_within(terr[idx],u[idx],np.random.default_rng(p)); Rs.append(RV.ratio(dnuc[vb],u2[vb]))
    out.append((a,roi,round(RV.ratio(dnuc[vb],u[vb]),3),np.round(Rs,3).tolist()))
for o in out: print(o)
