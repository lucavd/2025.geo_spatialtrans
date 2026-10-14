
import json, glob, os
R=os.path.expanduser('~/2025.geo_spatialtrans/results/S1.3c')
for p in sorted(glob.glob(R+'/review/recalc/recalc_bg_*.json')):
    ds=os.path.basename(p)[len('recalc_bg_'):-5]
    me=json.load(open(p)); rep=json.load(open(f'{R}/S1.3c_bg_{ds}_bg.json'))
    diffs=[]
    for gs in ('com','all'):
        for k in ('3','5'):
            for f in ('thr','n_windows','mean','frac_zero'):
                a=me['bg'][gs][k][f]; b=rep['bg'][gs][k][f]
                if abs(a-b)>1e-12: diffs.append((gs,k,f,a,b))
    c3m=me['C3']; c3r=rep['C3']
    print(ds, 'bg diffs:', diffs, '| thr com3/5', me['bg']['com']['3']['thr'], me['bg']['com']['5']['thr'], 'nwin', me['bg']['com']['3']['n_windows'],
          '| C3 mine diff_com', c3m['n_diff_com'], 'diff_all', c3m['n_diff_all'], 'rep', c3r['common_genes']['n_diff'], c3r['prereg_all_genes']['n_diff'],
          'ntis', c3m['n_in_tissue'], c3r['n_in_tissue'], 'nfilt', c3m['n_filtered'], c3r['n_filtered'], 'excl', round(c3m['frac_umi_excl_tissue'],5), round(c3r.get('frac_umi_excluded_genes_in_tissue',float('nan')),5),
          'umi', me['umi_tissue'], rep['umi_in_tissue'], me['feature_types'])
