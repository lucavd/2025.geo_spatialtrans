
import sys, json, re, math
import numpy as np, pandas as pd
sys.path.insert(0, 'tools')
from R4_common import verdict as V
OUT = '/tmp/r7_S15_S16'
KT = pd.read_csv(OUT + '/R7_kappa_final_v4.csv').set_index('archetype')
TIS = KT[KT.unit == 'area']; SIM = KT[KT.unit == 'sim_5'].iloc[0]; WIN = KT[KT.unit == 'window_4'].iloc[0]
GEO = {'boot': {'nuc': 1.751, 'bin': 1.760, 'box': 1.743}, 'se': {'nuc': 1.496, 'bin': 1.503, 'box': 1.491}, 'range': {'nuc': 1.020, 'bin': 1.019, 'box': 1.018}}
STRESS = 1.343
M = 2.0   # all estimate intervals in S1.5-S1.6 are est +- 2 SE
S1 = pd.read_csv('results/R2c/R2c_s1.csv'); S1 = S1[S1.set == 'Luca'].set_index('archetype')
S1 = {a: (float(S1.loc[a, 's1_primary']), float(S1.loc[a, 'lo_primary']), float(S1.loc[a, 'hi_primary'])) for a in S1.index}
LOAO_NOTE = 'simulated with leave-one-ROI-out fold parameters (held-out area shares up to 60% of nuclei with training areas in A3/A4); LOAO dependence not corrected (rule e), flagged; across-ROI SE of simulated areas corrected with sim_5'
REAL_EST = {('S1.5', 'P-6'), ('S1.5', 'P-7'), ('S1.5', 'CP-1')}
def side_of(step, pid, eu):
    if eu not in ('area',): return ''
    if (step, pid) in REAL_EST: return 'real'
    if step in ('S1.5c', 'S1.5d', 'S1.5e', 'S1.6') and pid in ('P-2', 'P-9'): return 'real'
    return 'simulated'

def kest(eu, side, q, a, mode):
    if eu != 'area': return 1.0
    if side == 'simulated': k = SIM['kappa_se_' + q] * 2 / M
    elif mode == 'g' and a == 'A5': k = GEO['se'][q] * 2 / M
    else: k = TIS.loc[a, 'kappa_se_' + q] * 2 / M
    return k * STRESS if mode == 's' else k
def krange(q, a, mode):
    if mode == 'g' and a == 'A5': return GEO['range'][q]
    return TIS.loc[a, 'kappa_range_' + q]
def kwin(mode):
    k = WIN['kappa_boot_nuc']; return k * STRESS if mode == 's' else k

def scale_iv(lo, hi, centre, k):
    if k == 1.0: return lo, hi
    lo2 = centre - (centre - lo)*k if np.isfinite(lo) else lo
    hi2 = centre + (hi - centre)*k if np.isfinite(hi) else hi
    return lo2, hi2
WPT = {a: v[0] for a, v in S1.items()}  # centre unused for window (one-stage replacement)
CONSTR_M = {'range': 'range', 'window': 'window_one_stage', 'fixed': 'fixed', 'other': 'other'}

rows = []
def add(step, pid, arch, quantity, est, lo, hi, a, b, eu, mu, q, vrec, notes, touch, clip=False, a_floor=None, rule='verdict', centre_m=None, role='', quarter=False):
    base_pid = re.sub(r'.*\[|\]', '', pid) if quarter else pid
    side = 'simulated' if (quarter and eu == 'area') else side_of(step, pid, eu)
    def vfun(l, h, aa, bb):
        if rule == 'verdict': return V(l, h, aa, bb)
        if rule == 'contains0': return ('confermata' if (l <= 0 <= h) else 'smentita') if np.isfinite(l) and np.isfinite(h) else 'non valutabile'
        return vrec
    vold = vfun(lo, hi, a, b)
    out = dict(step=step, prediction_id=pid, archetype=arch, quantity=quantity, est=est, lo=lo, hi=hi, a=a, b=b,
               est_unit=eu, margin_unit=mu, qtype=q, verdict_recorded=vrec,
               reproduces=('NA' if rule == 'copy' else ('TRUE' if vold == vrec else 'FALSE')),
               interval_construction={'area': 'se_m', 'seed': 'seed', 'fixed': 'fixed', 'other': 'other'}.get(eu, eu),
               side=side, interval_construction_margin=CONSTR_M.get(mu, mu), side_margin=('real' if mu in ('range', 'window') else ''))
    res = {}
    arch_ok = arch in TIS.index
    for mode in ['n', 's', 'g']:
        ke = kest(eu, side, q, arch, mode) if arch_ok else 1.0
        c_e = est if np.isfinite(est) else (lo + hi) / 2
        ln, hn = scale_iv(lo, hi, c_e, ke)
        if clip:
            ln = max(ln, 0.0) if np.isfinite(ln) else ln; hn = min(hn, 1.0) if np.isfinite(hn) else hn
        km = 1.0
        if quarter:
            if mu == 'window' and arch_ok:
                km = kwin(mode); s1, l1, h1 = S1[arch]; m = (h1 - l1) * km / 4
            elif mu == 'range' and arch_ok:
                km = krange(q, arch, mode); m = b * km
            else: m = b
            an, bn = -m, m
        elif mu == 'window' and arch_ok:
            km = kwin(mode); s1, l1, h1 = S1[arch]
            an, bn = s1 - (s1 - l1) * km, s1 + (h1 - s1) * km
        elif mu == 'range' and arch_ok:
            km = krange(q, arch, mode)
            cm = (a + b) / 2 if (np.isfinite(a) and np.isfinite(b)) else (a if np.isfinite(a) else b)
            an, bn = scale_iv(a, b, cm, km)
        else:
            an, bn = a, b
        if a_floor is not None: an = max(an, a_floor)
        res[mode] = (ln, hn, an, bn, ke, km, vfun(ln, hn, an, bn))
    ln, hn, an, bn, ke, km, vn = res['n']
    out.update(lo_new=ln, hi_new=hn, a_new=an, b_new=bn, kappa_est=ke, kappa_margin=km, verdict_new=vn, verdict_new_stress=res['s'][6],
               verdict_geom_A5=(res['g'][6] if arch == 'A5' else ''))
    affected = (eu == 'area') or (mu in ('range', 'window'))
    out['transition'] = 'unaffected' if not affected else ('unchanged' if vn == vrec else f'{vrec}->{vn}')
    nt = list(filter(None, [role, notes]))
    if mu == 'window':
        nt.append('margin: recorded two-stage window interval (B-127..B-132) NOT widened; replaced by the one-stage cluster percentile interval over the 4 windows (R2c_s1.csv, Luca, primary) x kappa_boot window_4' + (' (then /4 as S1.6 margin)' if quarter else ''))
    if side == 'real' and (step, base_pid) not in {('S1.5', 'CP-1')}:
        nt.append('paired simulated - real per ROI: real-area side (tissue kappa_se)')
    if res['s'][6] != vn: nt.append(f'stress verdict differs: {res["s"][6]}')
    if vold != vrec and rule != 'copy': nt.append(f'NOT REPRODUCED: recomputed old = {vold}')
    out['notes'] = '; '.join(nt); out['touches_conclusion'] = touch
    out['_vold'] = vold
    rows.append(out)
def est_iv(est, se):
    return est - 2*se, est + 2*se

def pick_clip(est, se, a, b, vrec, rule='verdict'):
    lo, hi = est_iv(est, se)
    v0 = V(lo, hi, a, b)
    if v0 == vrec: return lo, hi, False
    lc, hc = max(lo, 0.0), min(hi, 1.0)
    if V(lc, hc, a, b) == vrec and (lc != lo or hc != hi): return lc, hc, True
    return lo, hi, False

# ---------------------------------------------------------------- S1.5
d = pd.read_csv('results/S1.5/S1.5_verdicts.csv')
def t15(pid, arm, a):
    win = arm == 'preserve_area'
    if pid == 'P-1': return 'choice 2 (size-model conversion rho\' preserve_area in the 6 presets): winner by number of tissues with P-1 confermata'
    if pid == 'P-2': return 'choice 2 tie-break' + ('; P-2 smentita deferred to BL-143 (step closure)' if win and a in ('A1','A2') else '')
    if pid == 'P-4' and a == 'A2' and win: return 'deferred to BL-138 (DEC-005, no nucleus cap in colon tumour; step closure)'
    if pid == 'P-5' and win and a in ('A2','A3','A4','A6'): return 'deferred to BL-143 (homogeneous seeding; step closure)'
    if pid == 'P-6' and win: return 'deferred to BL-143' + (' (F-1 fragility of the heart)' if a == 'A6' else '') + '; step closure'
    if pid == 'P-8' and a == 'A4': return 'deferred to BL-144 (optimisation before Step 2; step closure)'
    return 'none'
for _, r in d.iterrows():
    src = r.interval_source; pid = r.id; a = r.archetype
    if pid in ('P-1','P-2','P-4'): mu = 'range'
    elif pid == 'P-5': mu = 'window'
    else: mu = 'fixed'
    if pid == 'P-5': eu = 'seed'
    elif pid == 'P-8': eu = 'seed' if r.n_units > 1 else 'fixed'
    else: eu = 'area'
    q = 'box' if pid == 'P-8' else 'nuc'
    clip = (abs(r.lo2 - (r.estimate - 2*r.se)) > 1e-12) or (abs(r.hi2 - (r.estimate + 2*r.se)) > 1e-12)
    role = f'arm={r.arm}' + (' (winner)' if r.winner_arm else ' (losing arm of choice 2)')
    add('S1.5', pid, a, r.quantity, r.estimate, r.lo2, r.hi2, r.int_lo, r.int_hi, eu, mu, q, r.outcome,
        ('margin one-sided [0, max of real ROIs], corrected literally around midpoint' if pid == 'P-4' else ''),
        t15(pid, r.arm, a), clip=clip, centre_m=WPT.get(a) if mu == 'window' else None, role=role)
c1 = pd.read_csv('results/S1.5/S1.5_CP1_verdicts.csv')
for _, r in c1.iterrows():
    vrec = r.outcome; vmap = 'non valutabile' if str(vrec).startswith('rifiuto') else vrec
    est = r.ratio_mean; lo, hi = (r.lo2, r.hi2) if np.isfinite(est) else (np.nan, np.nan)
    role = f'CP-1 variant={r.variant}, scale={r.scale_true} um/px' + (' (descriptive D-7, no verdict of record)' if r.variant == 'full_mask' else '')
    tc = 'estimate_pixel_size() removed from recommended use (D-11; >=3 of 6 tissues needed at 1 um/px)' if (r.variant == 'pipeline' and r.scale_true == 1) else 'none'
    add('S1.5', 'CP-1', r.archetype, f'scale ratio estimated/true at {r.scale_true} um/px ({r.variant})', est, lo, hi, 0.8, 1.2,
        'area' if np.isfinite(est) else 'fixed', 'fixed', 'box', vmap,
        ('recorded "' + vrec + '" mapped to non valutabile (no number returned)') if vmap != vrec else '', tc, role=role)

# ---------------------------------------------------------------- S1.5b
d = pd.read_csv('results/S1.5b/S1.5b_verdicts.csv')
def t15b(pid, arm, a, unit):
    if arm == 'single': return 'CP-1 (types add information: types="multi" default)' if pid in ('P-2','P-3') else 'none'
    if pid == 'P-2': return 'arm rule (choice 2): size_arm="recentred" default (DEC-019), P-2 count'
    if pid == 'P-1': return 'arm rule (choice 2) tie-break on P-1: recentred default (DEC-019)'
    if pid == 'P-5' and arm == 'recentred': return 'closes BL-138 (cytoplasm-free tail corrected by types)'
    if pid in ('P-7','P-8') and arm == 'recentred': return 'deferred to BL-153 (type-dependent seeding; step closure)'
    if pid == 'P-4' and arm == 'recentred': return 'deferred to BL-154 (glia dispersion; step closure)'
    return 'none'
for _, r in d.iterrows():
    pid = r.id; a = r.archetype
    mu = 'fixed' if pid == 'P-6' else 'range'
    lo, hi, clip = pick_clip(r.est, r.se, r.lo, r.hi, r.verdict)
    role = f'arm={r.arm}' + {'recentred':' (winner)', 'group':' (losing arm)', 'single':' (CP-1 null arm)'}[r.arm]
    nt = ''
    if pid == 'P-3': nt = 'margin = range of real ROIs AND >= 2.1 (P-1 of R6, fixed): a_new floored at 2.1'
    if pid == 'P-5': nt = 'margin one-sided [0, max], corrected literally around midpoint'
    add('S1.5b', pid, a, f'{pid} {r.unit}', r.est, lo, hi, r.lo, r.hi, 'area', mu, 'nuc', r.verdict, nt, t15b(pid, r.arm, a, r.unit),
        clip=clip, a_floor=2.1 if pid == 'P-3' else None, role=role)
s15b_rows = [x for x in rows if x['step'] == 'S1.5b']

# ---------------------------------------------------------------- S1.5c, S1.5d, S1.5e, S1.6 verdict tables
RANGE_SRC = ('fra ROI', 'P-7 di S1.5b', 'P-2 di S1.5b', 'P-2 di S1.5', 'P-1 di S1.5b', 'P-5 di S1.5b', 'real_q90')
def mu_of(src):
    s = str(src)
    if s.startswith('B-127'): return 'window'
    if s.startswith('B-2') or s.startswith('B-4') or s == 'S1.4': return 'fixed'
    if any(x in s for x in RANGE_SRC): return 'range'
    return 'other'

def touch_c(pid, arm):
    t = []
    if arm in ('C-a','N-a') and pid in ('P-1','P-5'): t.append('rule 1: min_dist_mode="nucleus" (DEC-020)')
    if arm in ('N-a','N-b') and pid in ('P-8a','P-8b'): t.append('rule 2: size model "group" with compartments (DEC-020)')
    if arm in ('N-a','U') and pid in ('P-1','P-2','P-3','P-4','P-5','P-6','P-7','P-8a','P-8b','P-8c'): t.append('rule 3: compartments not the default (DEC-020)')
    if arm in ('G','N-a') and pid in ('P-4','P-5','P-6'): t.append('CP-1 (types matter)')
    return '; '.join(t) or 'none'
def touch_d(pid, arm):
    t = []
    if arm in ('RV','RD','RVD') and pid in ('P-1','P-2','P-7','P-13'): t.append('rule A: heterogeneity family RV (preset void table)')
    if arm in ('RV','W3') and pid in ('P-1','P-2','P-3','P-4','P-5','P-6','P-7','P-8a','P-8b','P-8c','P-10','P-14'): t.append('rule C: gaussian size model kept (DEC-023)')
    if arm in ('U','U3') and pid in ('P-1','P-2','P-3','P-4','P-5','P-6','P-7','P-8a','P-8b','P-8c','P-10','P-14'): t.append('rule C on uniform seeding: gaussian kept (DEC-023)')
    if arm in ('RV','U') and pid in ('P-1','P-2','P-3','P-4','P-5','P-6','P-7','P-8a','P-8b','P-8c'): t.append('rule B: compartments not the default')
    if arm in ('VF','RV') and pid in ('P-11','P-13'): t.append('CP-2 (void size matters)')
    if arm in ('N','R') and pid in ('P-4','P-5','P-6'): t.append('recalibration doubt rule (R vs N, D-desc-6)')
    if arm in ('W3p','W3') and pid == 'P-8c': t.append("CP-3 (P-8c')")
    return '; '.join(t) or 'none'
def touch_e(pid, arm):
    t = []
    if arm in ('RV','X','XA','XADs') and pid in ('P-1','P-2','P-7','P-13','P-15','P-16','P-17'): t.append('rule A: void model RV kept (presets unchanged)')
    if arm in ('Q','RV') and pid in ('P-1','P-2','P-3','P-4','P-5','P-6','P-7'): t.append('rule Q: quotas=FALSE (BL-170, default since S1.6, DEC-025)')
    if arm in ('RV','U') and pid in ('P-1','P-2','P-3','P-4','P-5','P-6','P-7','P-8a','P-8b','P-8c'): t.append('rule B: compartments not the default')
    return '; '.join(t) or 'none'
ROLE = {'S1.5c': {'U':'reference (=S1.5b)','C-a':'candidate','C-b':'candidate','N-a':'winner','N-b':'candidate','G':'CP-1 arm'},
        'S1.5d': {'U':'default reference','N':'reference (=N-a S1.5c)','R':'recalibration','RV':'winner (rule A)','RD':'candidate','RVD':'candidate','VF':'CP-2 arm','Q':'sensitivity BL-160','W3':'rule C candidate','U3':'rule C uniform candidate','W3p':'CP-3 arm','U3c':'sensitivity D-14'},
        'S1.5e': {'U':'default reference','RV':'winner (rule A)','X':'candidate','XA':'candidate','XAD':'sensitivity (sigma_d=0 allowed)','XADs':'candidate (as pre-registered)','Q':'rule Q arm'},
        'S1.6': {}}
for st, tf in [('S1.5c', touch_c), ('S1.5d', touch_d), ('S1.5e', touch_e), ('S1.6', lambda p, a: 'none (descriptive D-2: state of the cell layer at end of Step 1)')]:
    d = pd.read_csv(f'results/{st}/{st}_verdicts.csv')
    for _, r in d.iterrows():
        pid = r.id; a = r.archetype; mu = mu_of(r.source)
        eu = 'seed' if pid == 'P-1' else 'area'
        lo, hi, clip = pick_clip(r.est, r.se, r.lo, r.hi, r.verdict)
        nt = []
        if not str(r.arm).startswith('U'): nt.append(LOAO_NOTE)
        if pid == 'P-8c': nt.append('margin one-sided [0, max], corrected literally around midpoint')
        if clip: nt.append('estimate interval clipped to [0,1] as in the step')
        if mu == 'other': nt.append(f'margin source not classified: {r.source}')
        role = f"arm={r.arm}" + (f" ({ROLE[st].get(r.arm,'')})" if ROLE[st].get(r.arm) else '')
        if st == 'S1.6': role += ' (descriptive D-2)'
        add(st, pid, a, f'{pid} {r.unit}', r.est, lo, hi, r.lo, r.hi, eu, mu, 'nuc', r.verdict, '; '.join(nt), tf(pid, r.arm),
            clip=clip, centre_m=WPT.get(a) if mu == 'window' else None, role=role)

# S1.5d / S1.5e controprove tables
def parse_pm(s):
    m = re.match(r'\s*([-0-9.e]+)\s*\+-\s*([-0-9.e]+)', str(s));
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)
cp = pd.read_csv('results/S1.5d/S1.5d_cp.csv')
for _, r in cp.iterrows():
    rec = 'confermata' if bool(r['pass']) else 'smentita'
    est, h = parse_pm(r.value)
    cpid = r.cp.split(' ')[0]
    if np.isfinite(est):
        m = re.search(r'(A[1-6])', r['item']); a = m.group(1) if m else 'A6'
        add('S1.5d', cpid, a, f"{r.cp}: {r['item']}", est, est - h, est + h, np.nan, np.nan, 'area', 'fixed', 'nuc', rec,
            'rule: paired difference +- 2 SE contains 0 (D-8); pass mapped to confermata, fail to smentita; values rounded as stored in S1.5d_cp.csv' + ('' if r.judged else '; not judged'),
            'CP-1 (calibration does not fit noise in the heart)' if cpid == 'CP-1' else 'CP-3 (explanation of the size tail: territory dependence)', rule='contains0', role='controprova')
    else:
        add('S1.5d', cpid, 'A6' if cpid == 'CP-1' else 'panel', f"{r.cp}: {r['item']} = {r.value}", np.nan, np.nan, np.nan, np.nan, np.nan, 'fixed' if cpid == 'CP-1' else 'other', 'fixed', 'nuc', rec,
            'calibration parameter check (fixed)' if cpid == 'CP-1' else 'count-based: recomputed in AGG rows', 'CP-1' if cpid == 'CP-1' else cpid, rule='copy', role='controprova' + ('' if r.judged else ' (not judged)'))
cp = pd.read_csv('results/S1.5e/S1.5e_cp.csv')
for _, r in cp.iterrows():
    rec = 'confermata' if bool(r['pass']) else 'smentita'
    cpid = r.cp.split(' ')[0]
    m = re.search(r'(A[1-6])', str(r['item'])); a = m.group(1) if m else 'A6'
    est, h = parse_pm(r.value)
    if np.isfinite(est):
        eu = 'seed' if cpid == 'CP-3' else 'area'
        add('S1.5e', cpid, a, f"{r.cp}: {r['item']}", est, est - h, est + h, np.nan, np.nan, eu, 'fixed', 'nuc', rec,
            'rule: difference +- 2 SE contains 0; pass mapped to confermata, fail to smentita' + ('' if r.judged else '; not judged (arm not winner, reported)'),
            'none (DEC-025: CP failures do not block)', rule='contains0', role='controprova')
    else:
        add('S1.5e', cpid, a, f"{r.cp}: {r['item']} = {r.value}", np.nan, np.nan, np.nan, np.nan, np.nan, 'fixed', 'fixed', 'nuc', rec,
            'calibration parameter check (fixed); pass mapped to confermata, fail to smentita', 'none (DEC-025: CP failures do not block)', rule='copy', role='controprova')

# S1.6 P-1 and CP-1 (paired differences, margin = width/4)
src16 = pd.read_csv('results/S1.5e/S1.5e_verdicts.csv')[['id','source']].drop_duplicates().set_index('id')['source'].to_dict()
CP1_JUDGED = {('P-15','A1'),('P-15','A4'),('P-15','A5'),('P-17','A1'),('P-3','A4')}
for fname, pidp in [('S1.6_P1.csv','P-1'), ('S1.6_CP1.csv','CP-1')]:
    d = pd.read_csv(f'results/S1.6/{fname}')
    for _, r in d.iterrows():
        mu = mu_of(src16.get(r.id, ''))
        eu = 'seed' if r.id == 'P-1' else 'area'
        lo, hi = r['diff'] - 2*r.se, r['diff'] + 2*r.se
        judged = (r.id, r.archetype) in CP1_JUDGED
        if pidp == 'P-1':
            tc = 'S1.6 closure: quotas=FALSE and thinning_bound="exact" defaults reproduce S1.5e Q (P-1)'
            role = 'P-1 Q6 - Q(S1.5e)' + ('; precision insufficient by construction (E3)' if r.by_construction else '')
        else:
            tc = 'CP-1 (DEC-031): comparison not blind' if judged else 'none'
            role = 'CP-1 XA - RV (S1.5e)' + ('; judged row (expected smentita)' if judged else '; reported, not judged')
        add('S1.6', f'{pidp}[{r.id}]', r.archetype, f'{r.id} {r.unit} paired difference', r['diff'], lo, hi, -r.margin, r.margin, eu, mu, 'nuc', r.verdict,
            f'margin = +-1/4 of width of the source interval ({src16.get(r.id,"")}); unit inherited from that source', tc, role=role, quarter=True)
c2 = pd.read_csv('results/S1.6/S1.6_CP2.csv')
for _, r in c2.iterrows():
    rec = 'confermata' if r.ok else 'smentita'
    add('S1.6', 'CP-2', r.archetype, f'accepted-proposal ratio Q6/Q (expected {r.expected})', r.ratio, r.ratio, r.ratio, 1.05 if '>' in r.expected else -np.inf, np.inf if '>' in r.expected else 0.95,
        'fixed', 'fixed', 'box', rec, 'ratio of sums over all units, no uncertainty; pass mapped to confermata', 'S1.6 closure: exact bound active', role='controprova')
add('S1.6b', 'controprova', 'panel', 'doc check rejects the 2025 README/CLAUDE.md (39/83 refs missing)', np.nan, np.nan, np.nan, np.nan, np.nan, 'fixed', 'fixed', 'box', 'confermata',
    'S1.6b has no prediction with a relevance interval (step file, Previsioni); gates are binary; PASS mapped to confermata', 'none', rule='copy', role='controprova (binary)')

T = pd.DataFrame(rows)
T.to_pickle(OUT + '/rows.pkl')
print(T.shape, T.reproduces.value_counts().to_dict())
print(T[T.reproduces == 'FALSE'][['step','prediction_id','archetype','quantity','verdict_recorded','_vold','lo','hi','a','b']].to_string())
