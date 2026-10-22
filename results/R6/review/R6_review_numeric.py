
import pandas as pd, numpy as np, os, json
LAB = '/mnt/micron/geo_spatialtrans/R6/labels'
RES = os.path.expanduser('~/2025.geo_spatialtrans/results/R6')
OUT = os.path.join(RES, 'review'); os.makedirs(OUT, exist_ok=True)
SR = {'A5': '/mnt/micron/geo_spatialtrans/data_real/Visium_HD_6p5mm_Mouse_Brain/segmented_outputs_x/segmented_outputs/cell_types/10x_Cloud/cell_types.csv',
      'A6': '/mnt/micron/geo_spatialtrans/data_real/Visium_HD_6p5mm_Human_Heart/segmented_outputs_x/segmented_outputs/cell_types/10x_Cloud/cell_types.csv'}
EXCL = {'non risolto', 'non classificabile'}
TOL = 1e-6

def map_a6(coarse, fine):
    s = (str(coarse) + ' | ' + str(fine)).lower(); c = str(coarse).lower()
    if any(k in s for k in ['cardiac muscle', 'cardiocyte', 'cardiac myocyte', 'myoblast']) or c == 'muscle cell': return 'cardiomiociti'
    if any(k in s for k in ['smooth muscle', 'pericyte', 'mural']): return 'muscolo liscio e periciti'
    if 'endothelial' in s: return 'endoteliali'
    if any(k in s for k in ['fibroblast', 'stromal', 'connective']): return 'fibroblasti'
    if any(k in s for k in ['macrophage', 'b cell', 't cell', 'lymphocyte', 'mast', 'monocyte', 'granulocyte', 'plasma cell', 'lymphoid', 'hematopoietic', 'leukocyte', 'dendritic']): return 'immunitari'
    return 'altro'

def map_a5(coarse, fine):
    s = (str(coarse) + ' | ' + str(fine)).lower()
    if 'astrocyte' in s: return 'astrociti'
    if 'oligodendrocyte' in s: return 'oligodendrociti'
    if 'microglia' in s: return 'microglia'
    if 'neuron associated' in s: return 'altro'
    if any(k in s for k in ['neuron', 'interneuron', 'pyramidal']): return 'neuroni'
    if any(k in s for k in ['endothelial', 'mural', 'pericyte', 'smooth muscle']): return 'vascolari'
    return 'altro'

rows = []
def add(variant, arch, item, mine, theirs):
    d = None if (mine is None or theirs is None or pd.isna(mine) or pd.isna(theirs)) else float(mine) - float(theirs)
    ok = (d is not None) and abs(d) <= TOL
    rows.append(dict(variant=variant, archetype=arch, item=item, mine=mine, theirs=theirs, diff=d, ok=ok))

extra = {}
for variant, suf in [('prereg', ''), ('ph1', '_ph1')]:
    refs = pd.read_csv(f'{RES}/R6_refs_type{suf}.csv')
    pred = pd.read_csv(f'{RES}/R6_predictions{suf}.csv')
    for A in ['A1','A2','A3','A4','A5','A6']:
        d = pd.read_parquet(f'{LAB}/{A}_exp_umi20_seed20261022_nuclei{suf}.parquet')
        assert d['roi_id'].nunique() == 5, (A, d['roi_id'].nunique())
        assert d['cell_type'].notna().all()
        n = len(d)
        g = d.groupby('cell_type')
        frac = g.size() / n
        med_area = g['area_um2'].median()
        med_ecc = g['eccentricity'].median()
        rA = refs[refs.archetype == A].set_index('cell_type')
        for ct in sorted(set(frac.index) | set(rA.index)):
            th = rA.loc[ct] if ct in rA.index else None
            add(variant, A, f'n[{ct}]', int(g.size().get(ct, 0)), None if th is None else th['n'])
            add(variant, A, f'fraction[{ct}]', frac.get(ct, np.nan), None if th is None else th['fraction'])
            add(variant, A, f'median_area_nuc[{ct}]', med_area.get(ct, np.nan), None if th is None else th['median_area_nuc'])
            add(variant, A, f'median_ecc_nuc[{ct}]', med_ecc.get(ct, np.nan), None if th is None else th['median_ecc_nuc'])
        # P-2
        keep = [ct for ct in frac.index if frac[ct] >= 0.05 and ct not in EXCL]
        p2 = med_area[keep].max() / med_area[keep].min()
        th = pred[(pred.id == 'P-2') & (pred.archetype == A)]['estimate']
        add(variant, A, 'P-2 max/min median area (types>=5%)', p2, th.iloc[0] if len(th) else None)
        extra[f'{variant}/{A}/P-2 types kept'] = {ct: round(float(med_area[ct]), 4) for ct in keep}
        if A == 'A5':
            glia = d[d.cell_type.isin(['astrociti', 'oligodendrociti', 'microglia'])]['area_um2'].median()
            neu = d[d.cell_type == 'neuroni']['area_um2'].median()
            th = pred[(pred.id == 'P-1') & (pred.archetype == A)]['estimate']
            add(variant, A, 'P-1 median area neuroni / glia', neu / glia, th.iloc[0] if len(th) else None)
            extra[f'{variant}/A5/P-1'] = dict(neu=float(neu), glia=float(glia))
        if A == 'A6':
            cm = d[d.cell_type == 'cardiomiociti']['eccentricity'].median()
            oth = d[(d.cell_type != 'cardiomiociti') & ~d.cell_type.isin(EXCL)]['eccentricity'].median()
            th = pred[(pred.id == 'P-3') & (pred.archetype == A)]['estimate']
            add(variant, A, 'P-3 median ecc cardiomiociti - others', cm - oth, th.iloc[0] if len(th) else None)
            extra[f'{variant}/A6/P-3'] = dict(cm=float(cm), others=float(oth))
        if A in SR:
            sr = pd.read_csv(SR[A]).rename(columns={'barcode': 'sr_cell_id'})
            assert sr['sr_cell_id'].is_unique
            m = d[d.sr_cell_id.notna()].merge(sr[['sr_cell_id', 'coarse_cell_type', 'fine_cell_type']], on='sr_cell_id', how='left')
            n_unmatched = int(m['coarse_cell_type'].isna().sum())
            m = m[m['coarse_cell_type'].notna()]
            m = m[~m['coarse_cell_type'].isin(['Low UMI Barcode', 'cell'])]
            f = map_a6 if A == 'A6' else map_a5
            m['sr_type'] = [f(a, b) for a, b in zip(m['coarse_cell_type'], m['fine_cell_type'])]
            all_matched = float((m['sr_type'] == m['cell_type']).mean())
            inf = m[(m['sr_type'] != 'altro') & ~m['cell_type'].isin(EXCL)]
            agree = float((inf['sr_type'] == inf['cell_type']).mean())
            th = pred[(pred.id == 'P-6') & (pred.archetype == A)]['estimate']
            add(variant, A, 'P-6 agreement with SR (matched, informative)', agree, th.iloc[0] if len(th) else None)
            th2 = pred[(pred.id == 'P-6 (all matched)') & (pred.archetype == A)]['estimate']
            extra[f'{variant}/{A}/P-6'] = dict(n_sr_nonnull=int(d.sr_cell_id.notna().sum()), n_unmatched_in_csv=n_unmatched,
                 n_after_excl_lowumi_cell=int(len(m)), n_informative=int(len(inf)), agree=agree,
                 all_matched_mine_guess=all_matched, all_matched_theirs=float(th2.iloc[0]) if len(th2) else None,
                 sr_type_counts=m['sr_type'].value_counts().to_dict())

out = pd.DataFrame(rows)
out.to_csv(f'{OUT}/R6_review_numeric.csv', index=False)
print('n_checks', len(out), 'n_discrepancies', int((~out.ok).sum()))
print(out[~out.ok].to_string())
print(json.dumps(extra, indent=0, default=str))
