
"""Revisione avversariale S1.3c (recalc): UMI per bin e soglie di sfondo, codice indipendente."""
import h5py, numpy as np, pandas as pd, json, os, sys, time
from scipy import ndimage as ndi
BASE='/mnt/micron/geo_spatialtrans/data_real'
OUT='/tmp/recalc_s13c'; os.makedirs(OUT, exist_ok=True)
REV=os.path.expanduser('~/2025.geo_spatialtrans/results/S1.3c/review/recalc')
def per_bc_sum(path, keep_ids=None):
    with h5py.File(path,'r') as f:
        m=f['matrix']; ids=m['features/id'][:]; ft=m['features/feature_type'][:]
        bc=m['barcodes'][:]; indptr=m['indptr'][:]; idx=m['indices'][:]; dat=m['data'][:]
    ge = ft==b'Gene Expression'
    col=np.repeat(np.arange(len(bc)), np.diff(indptr))
    out={'all': np.bincount(col, weights=dat*ge[idx], minlength=len(bc)).astype(np.int64)}
    if keep_ids is not None:
        keep=np.isin(ids, keep_ids)&ge
        out['com']=np.bincount(col, weights=dat*keep[idx], minlength=len(bc)).astype(np.int64)
    return bc, ids, ft, out
def bg_thr(tp, umi, k):
    nr=tp.array_row.max()+1; nc=tp.array_col.max()+1
    g=np.zeros((nr,nc),bool); g[tp.array_row.values[tp.in_tissue.values==1], tp.array_col.values[tp.in_tissue.values==1]]=True
    d=ndi.distance_transform_edt(~g)*2.0
    dist=d[tp.array_row.values, tp.array_col.values]
    isbg=(tp.in_tissue.values==0)&(dist>=50)
    wr=tp.array_row.values//k; wc=tp.array_col.values//k
    wid=wr*(nc//k+2)+wc
    u,inv=np.unique(wid, return_inverse=True)
    n=np.bincount(inv); nbg=np.bincount(inv, weights=isbg); tot=np.bincount(inv, weights=umi)
    ok=(n==k*k)&(nbg==k*k)
    t=tot[ok]
    return dict(thr=float(np.percentile(t,95)), n_windows=int(ok.sum()), mean=float(t.mean()), frac_zero=float((t==0).mean()))
res={}
for ds in sys.argv[1:]:
    t0=time.time(); p=f'{BASE}/{ds}/binned_outputs_x/binned_outputs/square_002um/'
    with h5py.File(p+'filtered_feature_bc_matrix.h5','r') as f: fids=f['matrix/features/id'][:]
    bc,ids,ft,u=per_bc_sum(p+'raw_feature_bc_matrix.h5', keep_ids=fids)
    fbc,_,_,fu=per_bc_sum(p+'filtered_feature_bc_matrix.h5')
    tp=pd.read_parquet(p+'spatial/tissue_positions.parquet')
    tp['barcode']=tp.barcode.astype(str)
    raw=pd.DataFrame({'barcode':bc.astype(str),'umi_com':u['com'],'umi_all':u['all']})
    m=tp.merge(raw,on='barcode',how='left',validate='1:1')
    n_missing=int(m.umi_com.isna().sum())
    n_missing_in_tissue=int((m.umi_com.isna()&(m.in_tissue==1)).sum())
    m[['umi_com','umi_all']]=m[['umi_com','umi_all']].fillna(0).astype(np.int64)  # barcode assenti dalla grezza (SR 3.0.0 omette i bin a 0 UMI) = 0 UMI
    filt=pd.DataFrame({'barcode':fbc.astype(str),'umi_filt':fu['all']})
    m=m.merge(filt,on='barcode',how='left',validate='1:1')
    tis=m.in_tissue.values==1
    c3=dict(n_in_tissue=int(tis.sum()), n_filtered=len(fbc), n_filtered_not_in_tissue=int((~tis & m.umi_filt.notna()).sum()),
            n_tissue_not_in_filtered=int((tis & m.umi_filt.isna()).sum()),
            n_diff_com=int((m.umi_com[tis]!=m.umi_filt[tis]).sum()), n_diff_all=int((m.umi_all[tis]!=m.umi_filt[tis]).sum()),
            frac_umi_excl_tissue=float(1-m.umi_com[tis].sum()/m.umi_all[tis].sum()))
    bg={gs:{str(k):bg_thr(m, m['umi_'+gs].values, k) for k in (3,5)} for gs in ('com','all')}
    res[ds]=dict(C3=c3, bg=bg, n_bins=len(m), n_umi_missing=n_missing, n_umi_missing_in_tissue=n_missing_in_tissue, n_genes_all=int((ft==b'Gene Expression').sum()), n_genes_com=int(np.isin(ids,fids).sum()),
                 feature_types={k.decode():int(v) for k,v in zip(*np.unique(ft,return_counts=True))},
                 umi_total=dict(com=int(m.umi_com.sum()), all=int(m.umi_all.sum())), umi_tissue=dict(com=int(m.umi_com[tis].sum()), all=int(m.umi_all[tis].sum())))
    m[['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres','umi_com','umi_all']].to_parquet(f'{OUT}/{ds}_bins.parquet')
    print(ds, round(time.time()-t0), json.dumps(res[ds]), flush=True)
    json.dump(res[ds], open(f'{REV}/recalc_bg_{ds}.json','w'), indent=1)
