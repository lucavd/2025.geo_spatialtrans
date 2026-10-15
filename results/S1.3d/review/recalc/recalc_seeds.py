# results/S1.3d/review/recalc/recalc_seeds.py — revisione «ricalcolo» S1.3d: i seed di H (rep 61..64) contro tutti i seed di R3, S1.3, S1.3b.
# Formule copiate dal codice (grep di set.seed / random_seed in tools/R3_*.R, tools/S1.3*_*.R, results/S1.3b/review/*): offset rispetto a BASE = 20261008.
import itertools, collections
BASE = 20261008; SEL = 20261010 - BASE   # SEL_SEED = BASE + 2
NROI = 30
def job13b(i, mi, rp): return (i - 1) * 60 + mi * 20 + rp       # tools/S1.3b_run.R job_index (rp 1..20)
H = {}
for i in range(1, NROI + 1):
    for rp in (61, 62, 63, 64): H[(i, 'CSR', rp)] = 1000 * i + 100 + rp
    for rp in (61, 62, 63): H[(i, 'RSA', rp)] = 1000 * i + 200 + rp; H[(i, 'RSArule', rp)] = 1000 * i + 300 + rp
assert len(H) == 300 and len(set(H.values())) == 300
gen = collections.defaultdict(set)   # seed -> usi (generatori di punti dei nulli)
oth = collections.defaultdict(set)   # seed -> altri usi (campionamenti, perturbazioni, permutazioni)
for i in range(1, NROI + 1):
    for rp in range(1, 21):
        gen[1000 * i + 100 + rp].add('R3/S1.3 CSR'); gen[1000 * i + 200 + rp].add('R3/S1.3 RSA'); gen[1000 * i + 300 + rp].add('S1.3 RSArule')
    oth[1000 * i + 300].add('R3_analysis CP-2b permutazioni'); oth[1000 * i + 400].add('R3_posthoc')
for ki in range(1, 6):
    for rp in range(1, 21): oth[50000 + 100 * ki + rp].add('S1.3_cp'); oth[51000 + 100 * ki + rp].add('S1.3_cp')
for rp in range(1, 21): oth[60000 + rp].add('S1.3_cp D-1')
oth[SEL].add('S1.3b sel_grid'); oth[SEL + 1].add('S1.3b sel_cp2'); oth[1].add('S1.3_review_checks (20261009)')
for sd in range(1, 21): oth[SEL + 100 + sd].add('S1.3b placebo / regen_audit')
for i in range(1, NROI + 1):
    for mi in range(3):
        for rp in range(1, 21):
            j = job13b(i, mi, rp); oth[SEL + j].add('S1.3b C-3 st_sample')
            for k in (8, 10, 12): oth[SEL + 7 * j + k].add('S1.3b CP-2 s3b_perturb')
coll_gen = {k: v for k, v in H.items() if v in gen}
coll_oth = {k: (v, sorted(oth[v])) for k, v in H.items() if v in oth}
print('H: 300 seed distinti; range offset', min(H.values()), '..', max(H.values()))
print('collisioni con seed di GENERAZIONE dei nulli (R3, S1.3, S1.3b):', len(coll_gen))
print('collisioni numeriche con seed di ALTRI usi (campionamento/perturbazione/permutazione):', len(coll_oth))
c = collections.Counter(u for _, (_, us) in coll_oth.items() for u in us); print('  per uso:', dict(c))
for k, (v, us) in list(sorted(coll_oth.items()))[:8]: print('  es.', k, 'seed BASE+%d' % v, us)
