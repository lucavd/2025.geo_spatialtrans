
import re, json
import numpy as np, pandas as pd
OUT = '/tmp/r7_S15_S16'
T = pd.read_pickle(OUT + '/rows.pkl')
T['arm'] = T.notes.str.extract(r'arm=([^ ;]+)')[0]
T['pid'] = T.prediction_id
T['_va5'] = np.where(T.verdict_geom_A5.astype(str) != '', T.verdict_geom_A5, T.verdict_new)
COLS = {'old': '_vold', 'rec': 'verdict_recorded', 'new': 'verdict_new', 'stress': 'verdict_new_stress', 'a5': '_va5'}
def cnt(step, arm, ids, vc, arch=None, what='confermata'):
    s = T[(T.step == step) & (T.arm == arm) & (T.pid.isin(ids))]
    if arch is not None: s = s[s.archetype.isin(arch)]
    return int((s[vc] == what).sum())
SP = ['P-1','P-2','P-3','P-4','P-5','P-6','P-7']
agg = []
def rule(step, rid, desc, fn, touch):
    res = {k: fn(v) for k, v in COLS.items()}
    agg.append(dict(step=step, prediction_id=rid, archetype='panel', quantity=desc, verdict_recorded=res['rec'][0], _vold=res['old'][0],
                    verdict_new=res['new'][0], verdict_new_stress=res['stress'][0], verdict_geom_A5=res['a5'][0], reproduces=res['old'][0] == res['rec'][0],
                    notes=f"counts recorded: {res['rec'][1]} | new: {res['new'][1]} | stress: {res['stress'][1]} | geometric A5 factors: {res['a5'][1]}", touches_conclusion=touch))
# S1.5 choice 2
def r15(vc):
    a = {arm: (cnt('S1.5', arm, ['P-1'], vc), cnt('S1.5', arm, ['P-2'], vc)) for arm in ['preserve_area', 'none']}
    w = 'preserve_area' if a['preserve_area'] >= a['none'] else 'none'
    return w, f"P-1/P-2 confermate preserve_area {a['preserve_area']}, none {a['none']}"
rule('S1.5', 'RULE choice 2', 'size-model conversion arm (P-1 tissues, tie P-2, tie preserve_area)', r15, "rho' preserve_area conversion in the 6 presets")
def r15cp(vc):
    s = T[(T.step == 'S1.5') & (T.pid == 'CP-1') & T.notes.str.contains('variant=pipeline, scale=1 ')]
    n = int((s[vc] == 'confermata').sum())
    return ('method useful' if n >= 3 else 'method removed from recommended use'), f'1 um/px pipeline confermata in {n}/6'
rule('S1.5', 'AGG CP-1 1 um/px', 'image-based scale estimate: ratio inside 0.8-1.2 in >= 3 of 6 tissues', r15cp, 'estimate_pixel_size() out of recommended use (allow_estimate=TRUE)')
# S1.5b
def r15b(vc):
    a = {arm: (cnt('S1.5b', arm, ['P-2'], vc), cnt('S1.5b', arm, ['P-1'], vc)) for arm in ['group', 'recentred']}
    w = 'recentred' if a['recentred'] > a['group'] else 'group'
    return w, f"P-2/P-1 confermate group {a['group']}, recentred {a['recentred']}"
rule('S1.5b', 'RULE choice 2', 'size-model arm (P-2 types of 17, tie P-1 tissues, tie group)', r15b, 'size_arm="recentred" default (DEC-019)')
def r15bcp(vc):
    w = r15b(vc)[0]
    nnull = cnt('S1.5b', 'single', ['P-2'], vc); nw = cnt('S1.5b', w, ['P-2'], vc)
    sm = {a: cnt('S1.5b', 'single', ['P-2'], vc, [a], 'smentita') for a in ['A1', 'A2', 'A5']}
    p3 = T[(T.step == 'S1.5b') & (T.arm == 'single') & (T.pid == 'P-3')].iloc[0]
    hi = p3.hi if vc in ('_vold', 'verdict_recorded') else p3.hi_new
    ok = nnull < nw and all(v >= 1 for v in sm.values()) and hi < 2.1
    return ('confermata' if ok else 'smentita'), f'null P-2 {nnull} vs winner({w}) {nw}; null smentite A1/A2/A5 {sm}; null P-3 hi {hi:.3f} (<2.1)'
rule('S1.5b', 'AGG CP-1', 'controprova nessun tipo (recorded "superata" mapped to confermata)', r15bcp, 'types="multi" default (types add information)')
for x in agg[-1:]:
    x['verdict_recorded'] = 'confermata'; x['reproduces'] = x['_vold'] == 'confermata'
# S1.5c
def c_rules(vc):
    r1 = {a: cnt('S1.5c', a, ['P-1', 'P-5'], vc) for a in ['C-a', 'N-a']}
    d = 'N' if r1['N-a'] > r1['C-a'] else 'C'
    a_, b_ = d + '-a', d + '-b'
    r2 = {a: (cnt('S1.5c', a, ['P-8b'], vc), cnt('S1.5c', a, ['P-8a'], vc)) for a in [a_, b_]}
    m = b_ if r2[b_] > r2[a_] else a_
    sp = {a: cnt('S1.5c', a, SP, vc) for a in [m, 'U']}; sz = {a: cnt('S1.5c', a, ['P-8a', 'P-8b', 'P-8c'], vc) for a in [m, 'U']}
    r3 = 'adottata' if (sp[m] > sp['U'] and sz[m] >= sz['U']) else 'non adottata'
    cp = {p: (cnt('S1.5c', 'G', [p], vc), cnt('S1.5c', m, [p], vc), cnt('S1.5c', 'C-b', [p], vc)) for p in ['P-4', 'P-5', 'P-6']}
    cpok = all(g < w for g, w, _ in cp.values())
    return dict(r1=(d, f'C-a {r1["C-a"]}, N-a {r1["N-a"]}'), r2=(m, f'{a_} {r2[a_]}, {b_} {r2[b_]} (P-8b, P-8a)'),
                r3=(r3, f'winner {m}: spatial {sp[m]} vs U {sp["U"]}; size {sz[m]} vs U {sz["U"]}'),
                cp=('confermata' if cpok else 'smentita', f'(G, winner {m}, C-b same positions) per P-4/P-5/P-6: {cp}'))
rule('S1.5c', 'RULE 1', 'minimum distance: C-a vs N-a on P-1 + P-5 confermate (tie tissue)', lambda vc: c_rules(vc)['r1'], 'min_dist_mode="nucleus" for compartments (DEC-020)')
rule('S1.5c', 'RULE 2', 'size model on winning distance: P-8b (tie P-8a, tie group)', lambda vc: c_rules(vc)['r2'], 'size model "group" with compartments (DEC-020)')
rule('S1.5c', 'RULE 3', 'default: winner vs U, spatial P-1..P-7 more AND size P-8 not fewer', lambda vc: c_rules(vc)['r3'], 'compartments option, not default (DEC-020)')
rule('S1.5c', 'AGG CP-1', 'G worse than winner on P-4, P-5, P-6 (recorded "superata" mapped to confermata)', lambda vc: c_rules(vc)['cp'], 'CP-1 supports type compartments')
agg[-4]['verdict_recorded'] = 'N (per nucleo)'.split(' ')[0]; agg[-3]['verdict_recorded'] = 'N-a'; agg[-2]['verdict_recorded'] = 'non adottata'; agg[-1]['verdict_recorded'] = 'confermata'
for x in agg[-4:]: x['reproduces'] = x['_vold'] == x['verdict_recorded']
# S1.5d
SZC = ['P-8a', 'P-8b', 'P-8c', 'P-10', 'P-14']; SZB = ['P-8a', 'P-8b', 'P-8c']
def d_rules(vc):
    A = {a: cnt('S1.5d', a, ['P-1', 'P-2', 'P-7', 'P-13'], vc) for a in ['RV', 'RD', 'RVD']}
    best = max(A.values()); wA = [a for a in ['RV', 'RD', 'RVD'] if A[a] == best][0]
    out = dict(A=(wA, f'RV {A["RV"]}, RD {A["RD"]}, RVD {A["RVD"]}'))
    def C(base, emp):
        sb, se_ = cnt('S1.5d', base, SZC, vc), cnt('S1.5d', emp, SZC, vc)
        pb, pe = cnt('S1.5d', base, SP, vc), cnt('S1.5d', emp, SP, vc)
        return ('empirical' if (se_ > sb and pe >= pb) else 'gaussian'), f'{base} size {sb} sp {pb}; {emp} size {se_} sp {pe}'
    cC = C('RV', 'W3'); cU = C('U', 'U3'); cUc = C('U', 'U3c')
    out['C'] = (cC[0] + ('' if wA == 'RV' else ' (W3 built on RV; rule A winner changed to ' + wA + ')'), cC[1])
    out['CU'] = cU; out['CUc'] = cUc
    fin = 'W3' if cC[0] == 'empirical' else wA
    ref = 'U3' if cU[0] == 'empirical' else 'U'
    sp = {a: cnt('S1.5d', a, SP, vc) for a in [fin, ref]}; sz = {a: cnt('S1.5d', a, SZB, vc) for a in [fin, ref]}
    out['B'] = ('adottata' if (sp[fin] > sp[ref] and sz[fin] >= sz[ref]) else 'non adottata', f'{fin} sp {sp[fin]} size {sz[fin]} vs {ref} sp {sp[ref]} size {sz[ref]}')
    v11 = (cnt('S1.5d', 'VF', ['P-11'], vc), cnt('S1.5d', 'RV', ['P-11'], vc)); v13 = (cnt('S1.5d', 'VF', ['P-13'], vc), cnt('S1.5d', 'RV', ['P-13'], vc))
    sA1 = cnt('S1.5d', 'VF', ['P-11'], vc, ['A1'], 'smentita')
    out['CP2'] = ('confermata' if (v11[0] < v11[1] and v13[0] < v13[1] and sA1 == 1) else 'smentita', f'P-11 VF/RV {v11}; P-13 VF/RV {v13}; VF P-11 A1 smentita {sA1}')
    ty = {a: cnt('S1.5d', a, ['P-4', 'P-5', 'P-6'], vc) for a in ['N', 'R']}
    out['RN'] = ('nessun dubbio' if ty['R'] >= ty['N'] else 'dubbio', f'P-4+P-5+P-6 N {ty["N"]}, R {ty["R"]}')
    cp3 = T[(T.step == 'S1.5d') & (T.pid == 'CP-3') & T.quantity.str.contains('W3p - W3')]
    p8 = T[(T.step == 'S1.5d') & (T.arm == 'W3p') & (T.pid == 'P-8c')].iloc[0]
    n0 = int((cp3[vc] == 'confermata').sum())
    out['CP3'] = ('confermata' if (p8[vc] == 'confermata' and n0 == len(cp3)) else 'smentita', f"W3p P-8c' {p8[vc]}; W3p-W3 differences containing 0: {n0}/{len(cp3)}")
    cp1 = T[(T.step == 'S1.5d') & (T.pid == 'CP-1')]
    out['CP1'] = ('confermata' if (cp1[vc] == 'confermata').all() else 'smentita', f'{int((cp1[vc]=="confermata").sum())}/{len(cp1)} items pass')
    return out
for k, rid, desc, tc, rec in [('A', 'RULE A', 'heterogeneity family RV/RD/RVD on P-1,P-2,P-7,P-13 (tie fewer parameters, RV before RD)', 'preset heterogeneity table = RV (voids only)', 'RV'),
                              ('C', 'RULE C', 'size model on rule-A winner: W3 (empirical residuals) wins if more size (P-8a,b,c,P-10,P-14) and not fewer spatial', 'gaussian size model kept (DEC-023)', 'gaussian'),
                              ('CU', 'RULE C uniform', 'size model on uniform seeding: U vs U3', 'gaussian size model kept in default uniform seeding (DEC-023)', 'gaussian'),
                              ('CUc', 'RULE C uniform sens. U3c', 'sensitivity D-14: U vs U3c (not the rule)', 'none (Luca chose model of today, DEC-023)', 'empirical'),
                              ('B', 'RULE B', 'default: final vs reference, spatial more AND size P-8 not fewer', 'compartments option, not default', 'non adottata'),
                              ('CP1', 'AGG CP-1', 'heart negative control (parameters + paired diffs)', 'CP-1 calibration not fitting noise', 'confermata'),
                              ('CP2', 'AGG CP-2', 'fine voids VF worse than RV in P-11 and P-13, P-11 A1 smentita', 'CP-2 void size matters', 'confermata'),
                              ('CP3', 'AGG CP-3', "permuted residuals W3p: P-8c' confermata and P-8a,b equal to W3", 'explanation of size tail (territory dependence)', 'smentita'),
                              ('RN', 'AGG R vs N', 'recalibration doubt rule: R fewer type confirmations than N', 'recalibrated compartments adopted (B-232..B-248)', 'nessun dubbio')]:
    rule('S1.5d', rid, desc, (lambda k: (lambda vc: d_rules(vc)[k]))(k), tc)
    agg[-1]['verdict_recorded'] = rec; agg[-1]['reproduces'] = agg[-1]['_vold'] == rec
# S1.5e
def e_rules(vc):
    ids = ['P-1', 'P-2', 'P-7', 'P-13', 'P-15', 'P-16', 'P-17']
    A = {a: cnt('S1.5e', a, ids, vc) for a in ['RV', 'X', 'XA', 'XADs', 'XAD']}
    order = ['RV', 'X', 'XA', 'XADs']; best = max(A[a] for a in order); w = [a for a in order if A[a] == best][0]
    out = dict(A=(w, f'RV {A["RV"]}, X {A["X"]}, XA {A["XA"]}, XADs {A["XADs"]} (XAD sens. {A["XAD"]})'))
    sp = {a: cnt('S1.5e', a, SP, vc) for a in ['Q', w, 'U']}; ty = {a: cnt('S1.5e', a, ['P-4', 'P-5', 'P-6'], vc) for a in ['Q', w]}
    if sp['Q'] == sp[w] and ty['Q'] == ty[w]: q = 'con quote'
    elif sp['Q'] >= sp[w] and ty['Q'] >= ty[w]: q = 'senza quote'
    else: q = 'con quote'
    out['Q'] = (q + ('' if w == 'RV' else ' (Q built on RV; rule A winner changed to ' + w + ')'), f'spatial Q {sp["Q"]} vs {w} {sp[w]}; types Q {ty["Q"]} vs {w} {ty[w]}')
    sz = {a: cnt('S1.5e', a, SZB, vc) for a in [w, 'U']}
    out['B'] = ('adottata' if (sp[w] > sp['U'] and sz[w] >= sz['U']) else 'non adottata', f'{w} sp {sp[w]} size {sz[w]} vs U sp {sp["U"]} size {sz["U"]}')
    return out
for k, rid, desc, tc, rec in [('A', 'RULE A', 'void model RV/X/XA/XADs on P-1,P-2,P-7,P-13,P-15,P-16,P-17 (tie fewer params, order RV,X,XA,XADs)', 'presets keep RV voids (moved cells)', 'RV'),
                              ('Q', 'RULE Q', 'quotas: Q not fewer than winner in spatial and in types (P-4..P-6); tie in both keeps quotas', 'quotas=FALSE default (BL-170, applied in S1.6; DEC-025)', 'senza quote'),
                              ('B', 'RULE B', 'default: winner vs U', 'compartments option, not default', 'non adottata')]:
    rule('S1.5e', rid, desc, (lambda k: (lambda vc: e_rules(vc)[k]))(k), tc)
    agg[-1]['verdict_recorded'] = rec; agg[-1]['reproduces'] = agg[-1]['_vold'] == rec
# S1.6
def s16(vc, which):
    if which == 'P1':
        s = T[(T.step == 'S1.6') & T.pid.str.startswith('P-1[')]
        nc, ns = int((s[vc] == 'confermata').sum()), int((s[vc] == 'smentita').sum())
        return ('smentita' if ns > 0 else ('confermata' if nc >= 167 else 'indeterminata')), f'{nc} confermate, {int((s[vc]=="indeterminata").sum())} indeterminate, {ns} smentite of {len(s)}'
    s = T[(T.step == 'S1.6') & T.pid.str.startswith('CP-1[')]
    j = s[s.notes.str.contains('judged row')]
    nj = int((j[vc] == 'smentita').sum())
    return ('confermata' if nj == len(j) else 'smentita'), f'judged rows smentite {nj}/{len(j)}; all rows smentite {int((s[vc]=="smentita").sum())}/{len(s)}'
rule('S1.6', 'AGG P-1', 'no smentita and >= 167 confermate of 179', lambda vc: s16(vc, 'P1'), 'S1.6 closure: new compartment defaults (quotas=FALSE, exact bound) same law as S1.5e Q')
rule('S1.6', 'AGG CP-1', 'XA - RV smentita in the 5 judged rows', lambda vc: s16(vc, 'CP1'), 'DEC-031: CP-1 failure judged non-blocking')
agg[-1]['verdict_recorded'] = 'smentita'; agg[-1]['reproduces'] = agg[-1]['_vold'] == 'smentita'
A = pd.DataFrame(agg)
for ccol in ['est','lo','hi','a','b','lo_new','hi_new','a_new','b_new','kappa_est','kappa_margin']: A[ccol] = np.nan
A['est_unit'] = 'count'; A['margin_unit'] = 'rule'; A['qtype'] = 'rule'; A['interval_construction'] = 'other'; A['side'] = ''; A['interval_construction_margin'] = 'other'; A['side_margin'] = ''
A['transition'] = np.where(A.verdict_new == A.verdict_recorded, 'unchanged', A.verdict_recorded + '->' + A.verdict_new)
A['notes'] = '[RULE/AGG: count of per-row verdicts, each corrected with its own tissue kappa; no pooled interval, so no panel-wide kappa needed] ' + A.notes + np.where(A.verdict_new_stress != A.verdict_new, '; stress outcome differs: ' + A.verdict_new_stress, '')
A['notes'] = A.notes + np.where(A.verdict_geom_A5 != A.verdict_new, '; geometric-A5 outcome differs: ' + A.verdict_geom_A5, '')
A['notes'] = A.notes + np.where(~A.reproduces, '; NOT REPRODUCED: recomputed old = ' + A._vold, '')
cols = ['step','prediction_id','archetype','quantity','est','lo','hi','a','b','est_unit','margin_unit','qtype','verdict_recorded','reproduces','lo_new','hi_new','a_new','b_new','kappa_est','kappa_margin','verdict_new','verdict_new_stress','verdict_geom_A5','transition','interval_construction','side','interval_construction_margin','side_margin','touches_conclusion','notes']
A['reproduces'] = np.where(A.reproduces.astype(bool), 'TRUE', 'FALSE')
F = pd.concat([T[cols], A[cols]], ignore_index=True)
F.to_csv(OUT + '/R7_M5_S15_S16.csv', index=False)
print(F.shape)
print(A[['step','prediction_id','verdict_recorded','reproduces','verdict_new','verdict_new_stress','notes']].to_string())
