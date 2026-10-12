"""R/testing/test_S1.3c.py — asserzioni unitarie su tools/S1.3c_ecm.py (casi a risultato noto). Uso: .venv/bin/python R/testing/test_S1.3c.py"""
import importlib.util, sys, numpy as np, pandas as pd
from pathlib import Path
sys.argv = sys.argv[:1]
spec = importlib.util.spec_from_file_location("ecm", Path(__file__).resolve().parents[2] / "tools/S1.3c_ecm.py")
E = importlib.util.module_from_spec(spec); spec.loader.exec_module(E)
NP = NF = 0
def check(name, cond):
    global NP, NF
    print(("PASS " if cond else "FAIL ") + name); NP += bool(cond); NF += (not cond)

# T1 totali CSC, colonne vuote, filtro feature
indptr = np.array([0, 2, 2, 5, 6]); data = np.array([1, 2, 3, 4, 5, 6]); ind = np.array([0, 1, 0, 1, 2, 2])
check("T1a umi_totals colonne (vuota inclusa)", np.array_equal(E.umi_totals(indptr, data), [3, 0, 12, 6]))
keep = np.array([True, True, False])
check("T1b umi_totals con filtro feature", np.array_equal(E.umi_totals(indptr, data, ind, keep), [3, 0, 7, 0]))
# T2 finestre complete
ar, ac = np.meshgrid(np.arange(6), np.arange(6), indexing="ij"); ar, ac = ar.ravel(), ac.ravel()
ok = np.ones(36, bool); ok[0] = False
wr, wc = E.window_ids(ar, ac, 3); g = E.window_table(wr, wc, np.ones(36), ok, 3)
check("T2a 4 finestre", len(g) == 4); check("T2b 3 complete", int(g.complete.sum()) == 3)
check("T2c somme", np.allclose(g.total, 9))
# T3 permutazione nei gruppi
rng = np.random.default_rng(0); grp = rng.integers(0, 50, 5000); val = rng.poisson(3, 5000)
p = E.permute_within(grp, val, np.random.default_rng(1))
check("T3a somme per gruppo conservate", np.array_equal(np.bincount(grp, p), np.bincount(grp, val)))
check("T3b multinsieme conservato", np.array_equal(np.sort(p), np.sort(val)))
check("T3c ordine cambiato", (p != val).mean() > 0.5)
pm = np.random.default_rng(1).permutation(val)
check("T3d mutante globale viola le somme", not np.array_equal(np.bincount(grp, pm), np.bincount(grp, val)))
check("T3e seed riproducibile", np.array_equal(p, E.permute_within(grp, val, np.random.default_rng(1))))
# T4 modello A
d = np.linspace(0, 10, 4001); emp = d > 6
r, ba = E.fit_model_a(d, emp)
check("T4a r* = 6 (+-0.25)", abs(r - 6) <= 0.25); check("T4b BA ~ 1", ba > 0.99)
r2, ba2 = E.fit_model_a(d, np.random.default_rng(3).random(len(d)) < 0.4)
check("T4c BA ~ 0.5 se vuoto indipendente da d", abs(ba2 - 0.5) < 0.05)
# T5 BA
check("T5 BA nota", abs(E.balanced_accuracy([1, 1, 0, 0], [1, 0, 0, 0]) - (1 + 2 / 3) / 2) < 1e-12)
# T6 componenti
wr = np.repeat(np.arange(5), 5); wc = np.tile(np.arange(5), 5); emp = np.zeros(25, bool)
emp[[0, 1, 5]] = True; emp[24] = True
a, al, n = E.components_large_frac(wr, wc, emp, 6.0, 72.0)
check("T6a 2 componenti", n == 2); check("T6b area vuota 4*36", a == 144)
check("T6c area grande = 108 (3 finestre > 72)", al == 108)
emp2 = np.zeros(25, bool); emp2[[0, 6]] = True
check("T6d diagonale = 2 componenti (4-conn)", E.components_large_frac(wr, wc, emp2, 6.0, 72.0)[2] == 2)
# T7 regola
check("T7a nessuno", E.verdict(0.04, 0.9, 0.9) == "nessuno"); check("T7b A", E.verdict(0.1, 0.8, 0.9) == "A")
check("T7c B", E.verdict(0.1, 0.6, 0.6) == "B"); check("T7d indet. (BA bassa, piccole)", E.verdict(0.1, 0.6, 0.4) == "indeterminato")
check("T7e indet. (BA intermedia)", E.verdict(0.1, 0.7, 0.9) == "indeterminato")
# T8 profilo
check("T8 rapporto profilo", abs(E.profile_ratio(np.array([1, 1, 5, 5]), np.array([4., 4, 2, 2])) - 0.5) < 1e-12)
# T9 frazione nucleare con spostamento
nm = np.zeros((40, 40), bool); nm[10:20, 10:20] = True
f = E.nuc_fraction_shift(nm, np.array([15, 15, 30, 1]), np.array([15, 30, 30, 1]), 0, 0, 2)
check("T9a dentro = 1", f[0] == 1); check("T9b fuori = 0", f[1] == 0 and f[2] == 0); check("T9c bordo = NaN", np.isnan(f[3]))
f2 = E.nuc_fraction_shift(nm, np.array([15]), np.array([30]), -15, 0, 2)
check("T9d spostamento riporta dentro", f2[0] == 1)
# T10 soglia di sfondo su dataset sintetico
R = C = 120; ar, ac = np.meshgrid(np.arange(R), np.arange(C), indexing="ij")
intis = ((ar >= 40) & (ar < 80) & (ac >= 40) & (ac < 80)).ravel()
umi = np.where(intis, 50, np.random.default_rng(5).poisson(0.2, R * C)); umi[(ar.ravel() == 0) & (ac.ravel() == 0)] = 0
b = pd.DataFrame({"array_row": ar.ravel(), "array_col": ac.ravel(), "in_tissue": intis.astype(int), "umi": umi})
bg = E.background_thresholds(b)
dist = np.maximum(np.maximum(40 - ar, ar - 79), np.maximum(40 - ac, ac - 79)).clip(min=0)  # Chebyshev, limite inferiore
check("T10a finestre di sfondo trovate", bg[3]["n_windows"] > 100)
w = E.window_table(*E.window_ids(ar.ravel(), ac.ravel(), 3), umi, np.ones(R * C, bool), 3)
check("T10b soglia = q95 di finestre poisson (<5)", bg[3]["thr"] < 5)
check("T10c nessuna finestra di tessuto nello sfondo", bg[3]["mean"] < 5)
print(f"\n{NP} PASS, {NF} FAIL"); sys.exit(1 if NF else 0)
