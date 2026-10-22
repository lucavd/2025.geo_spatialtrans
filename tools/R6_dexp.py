"""R6: measure Space Ranger 4 nucleus->cell expansion distance d_exp from barcode_mappings.
For a random sample of cells per dataset: distance (um) from each cell-only 2um bin centre to the
nearest nucleus bin centre of the same cell; report quantiles per dataset and pooled."""
import pyarrow.parquet as pq, pandas as pd, numpy as np, sys
from scipy.spatial import cKDTree
D='/mnt/micron/geo_spatialtrans/data_real'
DS={'A2A3':'Visium_HD_Human_Colon_Cancer','A4':'Visium_HD_Human_Lymph_Node_FFPE','A5':'Visium_HD_6p5mm_Mouse_Brain','A6':'Visium_HD_6p5mm_Human_Heart'}
rng=np.random.default_rng(20261022); rows=[]; alld=[]
for arch,ds in DS.items():
    bm=pq.read_table(f'{D}/{ds}/{ds}_barcode_mappings.parquet',columns=['square_002um','cell_id','in_nucleus','in_cell']).to_pandas()
    bm=bm[bm.in_cell]
    tp=pq.read_table(f'{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet',columns=['barcode','pxl_row_in_fullres','pxl_col_in_fullres']).to_pandas()
    import json; sf=json.load(open(f'{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/scalefactors_json.json'))
    um=sf['microns_per_pixel']
    bm=bm.merge(tp,left_on='square_002um',right_on='barcode')
    cells=bm.cell_id.unique(); samp=rng.choice(cells,size=min(20000,len(cells)),replace=False)
    sub=bm[bm.cell_id.isin(samp)]
    d_all=[]; n_cells=0; n_nonuc=0
    for cid,g in sub.groupby('cell_id'):
        nuc=g[g.in_nucleus]; cyt=g[~g.in_nucleus]
        if len(nuc)==0: n_nonuc+=1; continue
        n_cells+=1
        if len(cyt)==0: continue
        t=cKDTree(nuc[['pxl_row_in_fullres','pxl_col_in_fullres']].values)
        d,_=t.query(cyt[['pxl_row_in_fullres','pxl_col_in_fullres']].values); d_all.append(d*um)
    d_all=np.concatenate(d_all); alld.append(d_all)
    q=np.quantile(d_all,[.5,.9,.95,.99])
    rows.append(dict(dataset=ds,n_cells_sampled=n_cells,n_cells_without_nucleus_bins=n_nonuc,n_cell_only_bins=len(d_all),frac_bins_in_nucleus=sub.in_nucleus.mean(),
                     d_q50=q[0],d_q90=q[1],d_q95=q[2],d_q99=q[3],d_max=d_all.max(),um_per_px=um))
    print(rows[-1],flush=True)
alld=np.concatenate(alld); q95=np.quantile(alld,.95)
rows.append(dict(dataset='pooled',n_cells_sampled=sum(r['n_cells_sampled'] for r in rows),n_cell_only_bins=len(alld),d_q50=np.quantile(alld,.5),d_q90=np.quantile(alld,.9),d_q95=q95,d_q99=np.quantile(alld,.99),d_max=alld.max(),d_exp_um=round(q95*2)/2))
pd.DataFrame(rows).to_csv('results/R6/R6_dexp.csv',index=False); print(rows[-1])
