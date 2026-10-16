#!/usr/bin/env python
"""tools/R2c_analysis_test.py — test delle funzioni di tools/R2c_analysis.py su dati costruiti (C-4: i mutanti MU1-MU3 devono
far fallire almeno un test). Uso: [R2C_MUTANT=MU1] .venv/bin/python tools/R2c_analysis_test.py ; exit 0 = tutti PASS."""
import sys, importlib.util
import numpy as np, pandas as pd
from pathlib import Path
spec = importlib.util.spec_from_file_location("ana", Path(__file__).resolve().parent / "R2c_analysis.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
res = []
def ck(name, ok, val=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + " " + str(val))

# kappa: a=[t,t,t,s], b=[t,s,s,s] -> po 0.5, pe 0.375, kappa 0.2 (a mano)
k = A.kappa(["t", "t", "t", "s"], ["t", "s", "s", "s"]); ck("K1 kappa con marginali diverse = 0.2", abs(k - 0.2) < 1e-12, k)
ck("K2 kappa accordo perfetto = 1", A.kappa(["t", "s", "t"], ["t", "s", "t"]) == 1.0)
k3 = A.kappa(["t", "s", "t", "s"], ["s", "t", "s", "t"]); ck("K3 kappa disaccordo totale = -1", abs(k3 + 1) < 1e-12, k3)
# bootstrap a grappolo: 4 finestre omogenee 0,0,1,1 (15 nuclei) -> IC [0, 1]; senza reinserimento -> ampiezza 0; iid -> stretto
w = np.repeat(["a", "b", "c", "d"], 15); y = np.repeat([0.0, 0.0, 1.0, 1.0], 15)
est, lo, hi = A.cluster_boot(w, y, B=2000, seed=1)
ck("B1 stima = media", abs(est - 0.5) < 1e-12, est)
ck("B2 IC a grappolo ampio (>= 0.9)", hi - lo >= 0.9, (lo, hi))
y2 = np.r_[np.repeat(0.2, 15), np.repeat(0.3, 15), np.repeat(0.6, 15), np.repeat(0.7, 15)]
e2, l2, h2 = A.cluster_boot(w, y2, B=2000, seed=1)
ck("B3 IC non degenere con medie diverse", h2 - l2 > 0.1, (l2, h2))
ck("B4 NaN esclusi dalla media", abs(A.cluster_boot(["a", "a", "b"], [1, np.nan, 0], B=50, seed=1)[0] - 0.5) < 1e-12)
# S-1 per focale
pa = pd.DataFrame(dict(win_id=["w", "w", "w", "w"], i=[0, 0, 2, 4], j=[1, 2, 3, 5], answer=["separate", "touch", "unsure", "unsure"]))
fo = pd.DataFrame(dict(win_id=["w"] * 5, idx=[0, 1, 3, 4, 9]))
ck("S1 primaria", np.array_equal(A.focal_s1(pa, fo, "primary"), np.array([1, 0, np.nan, np.nan, 0]), equal_nan=True),
   A.focal_s1(pa, fo, "primary"))
ck("S2 non giudicabile = contatto", np.array_equal(A.focal_s1(pa, fo, "unsure_touch"), np.array([1, 0, 1, 1, 0.])))
ck("S3 non giudicabile = separati", np.array_equal(A.focal_s1(pa, fo, "unsure_sep"), np.array([1, 0, 0, 0, 0.])))
pb = pd.DataFrame(dict(win_id=["w", "w"], i=[0, 0], j=[1, 2], answer=["unsure", "separate"]))
ck("S4 unsure + separate -> 0 (non escluso)", A.focal_s1(pb, pd.DataFrame(dict(win_id=["w"], idx=[0])))[0] == 0.0)
# Clopper-Pearson bilaterale: k = 0 -> 1 - 0.025^(1/n)
ck("C1 CP superiore k=0 n=90", abs(A.cp_upper(0, 90) - (1 - 0.025 ** (1 / 90))) < 1e-10, round(A.cp_upper(0, 90), 4))
ck("C2 CP superiore k=1 n=90 > 0.05", A.cp_upper(1, 90) > 0.05, round(A.cp_upper(1, 90), 4))
ck("C3 fascia 0/90 PASS, 1/90 WARN, 2/90 FAIL", (A.band_verdict(0, 90), A.band_verdict(1, 90), A.band_verdict(2, 90)) == ("PASS", "WARN", "FAIL"))
ck("C4 fascia 0/30: limite > 5 % senza contatti = non valutabile (D-6)", A.band_verdict(0, 30) == "non valutabile", round(A.cp_upper(0, 30), 4))
ck("C5 soglia di potenza: 0/72 PASS, 0/71 non valutabile", (A.band_verdict(0, 72), A.band_verdict(0, 71)) == ("PASS", "non valutabile"))
# AUC
ck("U1 AUC perfetta", A.auc([3, 2, 1, 0], [1, 1, 0, 0]) == 1.0)
ck("U2 AUC con pareggio", A.auc([1, 1], [1, 0]) == 0.5)
# B-dec
ck("D1 entrambi dentro, cellpose piu' vicino", A.b_dec(0.2, 0.5, 0.35, 0.36, 0.25) == "cellpose")
ck("D2 entrambi dentro, |d| simili -> indeciso", A.b_dec(0.2, 0.5, 0.35, 0.40, 0.31) == "indeciso")
ck("D3 nessuno", A.b_dec(0.2, 0.3, 0.25, 0.5, 0.1) == "nessuno")
ck("D4 solo StarDist", A.b_dec(0.2, 0.3, 0.25, 0.5, 0.22) == "StarDist")
ck("D5 K-2 PASS/WARN/FAIL", (A.k2_verdict(0.3, 0.2, 0.4, 0.05), A.k2_verdict(0.43, 0.2, 0.4, 0.05), A.k2_verdict(0.5, 0.2, 0.4, 0.05)) == ("PASS", "WARN", "FAIL"))
ck("D6 soglie kappa", (A.kappa_verdict(0.6), A.kappa_verdict(0.45), A.kappa_verdict(0.39)) == ("PASS", "WARN", "FAIL"))
n_fail = len(res) - sum(res)
print("ALL PASS" if n_fail == 0 else f"{n_fail} FAIL", "| mutante", A.MUT)
sys.exit(0 if n_fail == 0 else 1)
