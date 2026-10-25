"""tools/S1.5d_verify_prereg.py - S1.5d: ricalcolo dai file di ogni numero citato nella pre-registrazione (aggiunta sul modello di taglia ed esplorazione).
Uscita: results/S1.5d/S1.5d_verification_log_prereg.csv"""
import pandas as pd, numpy as np
rows = []
def add(claim, file, sel, val, cited):
    rows.append(dict(claim=claim, source_file=file, selector=sel, recomputed=val, cited=cited, match=abs(float(val) - float(cited)) <= 0.5 * 10 ** (-(len(str(cited).split('.')[-1]) if '.' in str(cited) else 0)) + 1e-12))
C = pd.read_csv('results/S1.5d/S1.5d_size_check.csv'); C3 = pd.read_csv('results/S1.5d/S1.5d_size_check_m3.csv'); F = pd.read_csv('results/S1.5d/S1.5d_explore_sizeform.csv')
c = C[(C.archetype=='A2')&(C.cell_type=='epiteliale')].mean(numeric_only=True); c3 = C3[(C3.archetype=='A2')&(C3.cell_type=='epiteliale')].mean(numeric_only=True)
f = F[(F.archetype=='A2')&(F.cell_type=='epiteliale')].mean(numeric_only=True)
for k, v in [('real_nc095','0.0264'),('M0_nc095','0.0687'),('M1_nc095','0.0581'),('real_med','32.65'),('M0_med','32.57'),('M1_med','31.48'),('real_q90','68.09'),('M0_q90','67.80'),('M1_q90','69.58')]:
    add(f'A2 epithelium {k} (LOO, real territories, group model)', 'S1.5d_size_check.csv', 'A2 epiteliale, mean of 5 folds', c[k], v)
for k, v in [('M3_nc095','0.0293'),('M3_med','32.65'),('M3_q90','68.42')]:
    add(f'A2 epithelium {k}', 'S1.5d_size_check_m3.csv', 'A2 epiteliale, mean of 5 folds', c3[k], v)
add('A2 epithelium logit model nc>=0.95', 'S1.5d_explore_sizeform.csv', 'A2 epiteliale M2_nc095', f['M2_nc095'], '0.0150')
add('A2 epithelium logit model q90', 'S1.5d_explore_sizeform.csv', 'A2 epiteliale M2_q90_an', f['M2_q90_an'], '79.9')
add('A2 epithelium per-type heteroscedastic nc>=0.95', 'S1.5d_explore_sizeform.csv', 'A2 epiteliale M1_nc095', f['M1_nc095'], '0.047')
R = pd.read_csv('results/S1.5d/S1.5d_size_resid.csv'); r = R[(R.archetype=='A2')&(R.size_group==1)].iloc[0]
add('A2 group 1 residual skewness', 'S1.5d_size_resid.csv', 'A2 g1 skew', r['skew'], '-0.45'); add('A2 group 1 q99 of standardized residuals', 'S1.5d_size_resid.csv', 'A2 g1 q99_z', r.q99_z, '2.10')
add('groups with negative skew (of 10)', 'S1.5d_size_resid.csv', 'count skew<0', (R['skew'] < 0).sum(), '8')
G = pd.read_csv('results/S1.5d/S1.5d_size_groups.csv'); add('A2 group 1 slope of log sigma (pooled)', 'S1.5d_size_groups.csv', 'A2 g1 pooled c1', G[(G.archetype=='A2')&(G.size_group==1)&(G.fold=='pooled')].c1.iloc[0], '0.267')
add('share of M0 excess removed by M1', 'S1.5d_size_check.csv', '(M0-M1)/(M0-real)', (c.M0_nc095-c.M1_nc095)/(c.M0_nc095-c.real_nc095), '0.25')
add('M0 / real tail ratio', 'S1.5d_size_check.csv', 'M0/real', c.M0_nc095/c.real_nc095, '2.6')
X = pd.read_csv('results/S1.5d/S1.5d_explore_size_real.csv'); x = X[(X.archetype=='A2')&(X.cell_type=='epiteliale')].iloc[0]
add('A2 epithelium sigma, smallest territory quartile', 'S1.5d_explore_size_real.csv', 'sigma_q1', x.sigma_q1, '0.39'); add('A2 epithelium sigma, largest quartile', 'S1.5d_explore_size_real.csv', 'sigma_q4', x.sigma_q4, '0.60')
add('A2 epithelium homoscedastic lognormal tail, per type', 'S1.5d_explore_size_real.csv', 'frac_nc095_lognormal', x.frac_nc095_lognormal, '0.070'); add('A2 epithelium real tail', 'S1.5d_explore_size_real.csv', 'frac_nc095', x.frac_nc095, '0.027')
pos = X[(X.cell_type!='_all')&(X.cell_type!='non risolto')&(X.cell_type!='non classificabile')]
add('types where sigma rises from Q1 to Q4 (of 20)', 'S1.5d_explore_size_real.csv', 'count sigma_q4>sigma_q1, typed', (pos.sigma_q4 > pos.sigma_q1).sum(), '18')  # in chat era scritto 17: errore di lettura a occhio, corretto (2026-10-25)
V = pd.read_csv('results/S1.5c/S1.5c_verdicts.csv'); v = V[(V.id=='P-8c')]
add('S1.5c P-8c uniform', 'results/S1.5c/S1.5c_verdicts.csv', 'P-8c U est', v[v.arm=='U'].est.iloc[0], '0.022'); add('S1.5c P-8c compartments N-a', 'results/S1.5c/S1.5c_verdicts.csv', 'P-8c N-a est', v[v.arm=='N-a'].est.iloc[0], '0.034')
Q = pd.read_csv('results/S1.5d/S1.5d_real_q90.csv')
def qi(A, t): q = Q[(Q.archetype==A)&(Q.cell_type==t)].iloc[0]; return q.q90_min, q.q90_max
for A, t, lo, hi in [('A5','astrociti','58.96','65.06'),('A3','endoteliale','54.5','63.5'),('A6','immunitari','27.7','33.9'),('A1','immunitario','17.2','18.5'),('A3','fibroblasti','62.8','72.2')]:
    a, b = qi(A, t); add(f'real q90 interval {A} {t} min', 'S1.5d_real_q90.csv', 'q90_min', a, lo); add(f'real q90 interval {A} {t} max', 'S1.5d_real_q90.csv', 'q90_max', b, hi)
inn = 0
for A, t in zip(Q.archetype, Q.cell_type):
    m = C3[(C3.archetype==A)&(C3.cell_type==t)].M3_q90.mean(); a, b = qi(A, t); inn += (a <= m <= b)
add('M3 q90 inside real interval (of 17), real territories', 'S1.5d_size_check_m3.csv + S1.5d_real_q90.csv', 'count', inn, '12')
add('S1.5c bench minutes-worker per unit', 'docs/steps/S1.5c.md (248.7 min, 3 workers, 1800 units)', '248.7*3/1800', 248.7*3/1800, '0.41')
S = pd.read_csv('results/S1.5d/S1.5d_explore_space.csv'); Sr = S.groupby(['archetype','roi_id','source']).mean(numeric_only=True).reset_index().groupby(['archetype','source']).mean(numeric_only=True)
for A, src, k, v in [('A1','real','H_rel_3.0','0.079'),('A1','N-a','H_rel_3.0','0.0011'),('A2','real','xcv2_50','0.220'),('A6','real','H_rel_3.0','0.0095'),('A6','N-a','H_rel_3.0','0.0067'),('A4','real','nn_cv','0.20') if False else ('A3','real','xcv2_50','0.493')]:
    add(f'exploration {A} {src} {k}', 'S1.5d_explore_space.csv', 'mean of ROI means', Sr.loc[(A, src), k], v)
I = pd.read_csv('results/S1.5d/S1.5d_real_intervals.csv')
for A, mt, lo, hi in [('A1','frac_gt3','0.0463','0.0604'),('A1','H_rel_3.0','0.0638','0.0943')]:
    r_ = I[(I.archetype==A)&(I.metric==mt)].iloc[0]; add(f'real interval {A} {mt} min', 'S1.5d_real_intervals.csv', 'min', r_['min'], lo); add(f'real interval {A} {mt} max', 'S1.5d_real_intervals.csv', 'max', r_['max'], hi)
D = pd.DataFrame(rows); D.to_csv('results/S1.5d/S1.5d_verification_log_prereg.csv', index=False)
print(D[['claim','recomputed','cited','match']].to_string(index=False)); print('mismatches:', (~D.match).sum(), 'of', len(D))
