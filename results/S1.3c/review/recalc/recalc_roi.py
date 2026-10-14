
"""Revisione avversariale S1.3c (recalc): metriche per ROI, codice indipendente (non importa tools/S1.3c_ecm.py)."""
import numpy as np, pandas as pd, json, os, sys, time
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import io
R=os.path.expanduser('~/2025.geo_spatialtrans'); REV=R+'/results/S1.3c/review/recalc'; TMP='/tmp/recalc_s13c'
MNT='/mnt/micron/geo_spatialtrans'
ROIS=[('A1','r1','cellpose_rgb','cellpose_rgb'),('A2','r1','cellpose_rgb','cellpose_rgb'),('A3','r1','cellpose_rgb','cellpose_rgb'),
      ('A4','f1','spaceranger','stardist_he'),('A5','r1','cellpose_rgb','cellpose_rgb'),('A6','r1','spaceranger','stardist_he')]
SEEDS=list(range(9001,9001+int(os.environ.get('RECALC_NPERM','30'))))   # seed del revisore (30), diversi da 1..20 della pipeline
rois=pd.read_csv(R+'/results/R2/R2_rois_checked.csv'); nuc=pd.read_parquet(R+'/results/R2/R2_nuclei_all.parquet')
MAP=os.environ.get('RECALC_MAP','rint_noclip')   # rint_noclip (revisore) | floor_clip (convenzione della pipeline, per isolare la causa)
def lookup(mask, iy, ix, f=4):
    a=iy//f; b=ix//f; ok=(a<mask.shape[0])&(b<mask.shape[1])
    if MAP=='px4':   # riproduce la pipeline: px4 = side_um / shape[1] (fattore di scala, non 4 px nativi)
        return mask[np.clip((iy_f/(SIDE_UM/mask.shape[1])).astype(int),0,mask.shape[0]-1), np.clip((ix_f/(SIDE_UM/mask.shape[1])).astype(int),0,mask.shape[1]-1)], int((~ok).sum())
    if MAP=='floor_clip':
        return mask[np.clip(a,0,mask.shape[0]-1), np.clip(b,0,mask.shape[1]-1)], int((~ok).sum())
    out=np.zeros(len(iy),bool); out[ok]=mask[a[ok],b[ok]]; return out, int((~ok).sum())
def windows(wid, sel_full, umi, k):
    """wid: id finestra per bin del ROI; sel_full: bin che qualificano; ritorna id finestre complete, totale UMI, mappa bin->finestra"""
    u,inv=np.unique(wid, return_inverse=True)
    n=np.bincount(inv); nq=np.bincount(inv, weights=sel_full)
    full=(n==k*k)&(nq==k*k)
    return inv, full
def phi(inv, full, umi, thr):
    tot=np.bincount(inv, weights=umi, minlength=len(full))
    return float((tot[full]<=thr).mean()), tot
def permute_within(umi, terr, mask, rng):
    out=umi.copy(); ii=np.where(mask)[0]; t=terr[ii]
    a=np.argsort(t, kind='stable'); b=np.lexsort((rng.random(len(ii)), t))
    new=np.empty(len(ii), umi.dtype); new[a]=umi[ii][b]; out[ii]=new; return out
def analyse(A, roi, terr_m, mask_m, k_list=(3,5), zero=None):
    row=rois[(rois.archetype==A)&(rois.roi_id==roi)].iloc[0]; ds=row.dataset; u=row.um_per_px; S=int(row.side_px)
    bins=pd.read_parquet(f'{TMP}/{ds}_bins.parquet')
    b=bins[(bins.pxl_col_in_fullres>=row.c0)&(bins.pxl_col_in_fullres<row.c0+S)&(bins.pxl_row_in_fullres>=row.r0)&(bins.pxl_row_in_fullres<row.r0+S)].reset_index(drop=True)
    px=b.pxl_col_in_fullres.values-row.c0; py=b.pxl_row_in_fullres.values-row.r0
    rnd=np.floor if MAP in ('floor_clip','px4') else np.rint
    global iy_f, ix_f, SIDE_UM
    SIDE_UM=S*u
    iy_f=py*u; ix_f=px*u
    ix=np.clip(rnd(px).astype(int),0,S-1); iy=np.clip(rnd(py).astype(int),0,S-1)
    valid=io.imread(f'{R}/results/R2/tissue_masks/{A}_{roi}_valid_ds4.png')>0
    tis, n_oob=lookup(valid, iy, ix)
    mk=np.load(f'{MNT}/R2/masks/{A}_{roi}_masks.npz'); tds=mk['tissue_ds4']; bub=mk['bubble_ds4']
    in_t,_=lookup(tds,iy,ix); in_b,_=lookup(bub,iy,ix)
    lab=np.load(f'{MNT}/R2/masks/{A}_{roi}_{mask_m}_native.npz')['labels']
    kl=nuc[(nuc.archetype==A)&(nuc.roi_id==roi)&(nuc.method==mask_m)&(nuc.scale==1)&(nuc.keep)].label.values
    nm=np.isin(lab, kl)
    dmap=ndi.distance_transform_edt(~nm)*u
    d=dmap[iy,ix]
    g=pd.read_parquet(f'{MNT}/R3/real/{A}_{roi}_{terr_m}_gen.parquet')
    _,terr=cKDTree(g[['x','y']].values).query(np.c_[px*u, py*u])
    umi=b.umi_com.values.astype(np.int64)
    out=dict(archetype=A, roi_id=roi, dataset=ds, n_bins_roi=len(b), n_tissue_bins=int(tis.sum()), n_bins_outside_ds4=n_oob, n_gen=len(g), n_keep_nuclei=len(kl),
             n_terr_used=int(len(np.unique(terr[tis]))), nuc_frac_px=float(nm.mean()))
    out['R']=float(umi[tis&(d>=4)&(d<8)].mean()/umi[tis&(d<2)].mean())
    out['umi_tissue_mean']=float(umi[tis].mean())
    lum=(b.in_tissue.values==1)&(~in_t)&(~in_b)
    out['n_lumen_bins']=int(lum.sum()); out['umi_lumen_mean']=float(umi[lum].mean()) if lum.any() else None
    # alternativa: fuori dal tessuto VALIDO (non dalla tissue_ds4) e fuori bolle
    out['n_lumen_bins_altvalid']=int(((b.in_tissue.values==1)&(~tis)&(~in_b)).sum())
    zero_info=None
    if zero is not None:
        dt=d[tis]; r90=np.percentile(dt,90)
        if zero=='thr': z=tis&(d>r90)
        else:  # rank: esattamente il 10 % dei bin di tessuto con d_nuc piu' alto (pareggi rotti a caso, seed 777)
            ii=np.where(tis)[0]; o=np.lexsort((np.random.default_rng(777).random(len(ii)), -d[ii])); z=np.zeros(len(b),bool); z[ii[o[:int(round(0.1*len(ii)))]]]=True
        umi=umi.copy(); umi[z]=0; zero_info=dict(mode=zero, r90=float(r90), frac_zeroed=float(z.sum()/tis.sum()), n_zeroed=int(z.sum()))
    thrs=json.load(open(f'{REV}/recalc_bg_{ds}.json'))['bg']['com']
    res=[]
    for k in k_list:
        thr=thrs[str(k)]['thr']
        wid=(b.array_row.values//k).astype(np.int64)*100000+(b.array_col.values//k)
        inv,full=windows(wid, tis, umi, k)
        po,tot=phi(inv,full,umi,thr)
        nulls=[]
        for sd in SEEDS:
            up=permute_within(umi, terr, tis, np.random.default_rng(sd))
            if sd==SEEDS[0]:
                t0=np.bincount(terr[tis], weights=umi[tis]); t1=np.bincount(terr[tis], weights=up[tis]); assert np.array_equal(t0,t1)
                assert np.array_equal(np.sort(up[tis]), np.sort(umi[tis]))
            nulls.append(phi(inv,full,up,thr)[0])
        nulls=np.array(nulls)
        # modello A
        dmin=np.full(len(full), np.inf); np.minimum.at(dmin, inv, d)
        dm=dmin[full]; emp=tot[full]<=thr
        grid=np.arange(0, dm.max()+0.25, 0.25); pf=np.array([(dm>r).mean() for r in grid])
        rs=grid[np.argmin(np.abs(pf-po))]; pred=dm>rs
        BA=float(((pred&emp).sum()/emp.sum() + ((~pred)&(~emp)).sum()/(~emp).sum())/2) if emp.any() else None
        # CP-2 finestre
        invl,fulll=windows(wid, lum, umi, k); totl=np.bincount(invl, weights=umi, minlength=len(fulll))
        res.append(dict(out, k=k, thr=thr, n_tissue_windows=int(full.sum()), n_empty=int(emp.sum()), phi_obs=po, phi_null_mean=float(nulls.mean()),
                        phi_null_min=float(nulls.min()), phi_null_max=float(nulls.max()), phi_null_sd=float(nulls.std(ddof=1)),
                        phi_null_mean_first20=float(nulls[:20].mean()), E=po-float(nulls.mean()), r_star=float(rs), BA=BA,
                        cp2_n_windows=int(fulll.sum()), cp2_n_empty=int((totl[fulll]<=thr).sum()), zero=json.dumps(zero_info)))
    return res
if __name__=='__main__':
    mode=sys.argv[1]
    allr=[]
    if mode=='roi':
        only=sys.argv[2:]
        for A,roi,tm,mm in ROIS:
            if only and A not in only: continue
            t0=time.time(); rr=analyse(A,roi,tm,mm); allr+=rr; print(A,roi,round(time.time()-t0),flush=True)
        pd.DataFrame(allr).to_csv(f'{REV}/recalc_roi_'+MAP+('_n'+str(len(SEEDS)) if len(SEEDS)!=30 else '')+('_'+'_'.join(only) if only else '')+'.csv', index=False)
    else:
        for zm in ('thr','rank'):
            allr+=analyse('A3','r1','cellpose_rgb','cellpose_rgb',k_list=(3,),zero=zm)
        pd.DataFrame(allr).to_csv(f'{REV}/recalc_ctrlA_A3r1_'+MAP+'.csv', index=False)
    print('done')
