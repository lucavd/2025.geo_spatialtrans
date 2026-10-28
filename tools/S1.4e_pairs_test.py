#!/usr/bin/env python
"""tools/S1.4e_pairs_test.py - S1.4e: test dell'angolo di allineamento e del contatto di coppia (cancello G-1, T5 e T6; mutanti G-2).
Uso (root del repo): [S14E_PYSRC=<cartella con tools/>] .venv/bin/python tools/S1.4e_pairs_test.py"""
import importlib.util, os, sys
import numpy as np, shapely
from shapely import affinity
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; SRC = Path(os.environ.get("S14E_PYSRC", ROOT))
sp = importlib.util.spec_from_file_location("s14e_measure", SRC / "tools/S1.4e_measure.py"); MS = importlib.util.module_from_spec(sp); sp.loader.exec_module(MS)
res = []
def ck(name, ok, val=""):
    res.append(bool(ok)); print(f"{'PASS' if ok else 'FAIL'} {name} {val}")
def ell(cx, cy, a, b, rot=0):
    return affinity.translate(affinity.rotate(affinity.scale(shapely.Point(0, 0).buffer(1, 64), a, b), rot, origin=(0, 0)), cx, cy)
# ---- T5: psi = 90 di fianco, 0 in fila (ellissi con l'asse maggiore lungo y) ----
for name, other, want in [("di fianco", (5, 0), 90.0), ("in fila", (0, 7), 0.0)]:
    G = [ell(0, 0, 2, 3.2), ell(other[0], other[1], 2, 3.2)]; S = np.array([MS.poly_shape(g) for g in G])
    d, j, u = MS.nn_dir(S[:, :2]); psi = MS.psi_deg(S[:, 2], u[:, 0], u[:, 1])
    ck(f"T5 psi {name} = {want:.0f} gradi", np.all(np.abs(psi - want) < 1e-6), f"psi {np.round(psi, 6)}")
G = [ell(10, 10, 2, 3.2, 37)]; S = MS.poly_shape(G[0])
ck("T5 allungamento ed eccentricita' dai momenti (ellisse 2 x 3.2 ruotata di 37 gradi)", abs(S[3] - 1.6) < 1e-3 and abs(S[4] - np.sqrt(1 - 1 / 1.6 ** 2)) < 1e-3
   and abs(((np.degrees(S[2]) - (90 + 37)) + 90) % 180 - 90) < 0.1, f"allungamento {S[3]:.4f}, ecc {S[4]:.4f}, angolo {np.degrees(S[2]):.2f}")
# ---- T6: contatto di coppia: Chebyshev 2 px si', 3 px no; per nucleo coincide con n_c2 > 0 ----
lab = np.zeros((60, 80), np.int64); lab[10:20, 10:20] = 1; lab[10:20, 21:30] = 2; lab[40:50, 10:20] = 3; lab[40:50, 22:30] = 4
pr = MS.contacts(lab); got = set(zip(pr.label_i, pr.label_j))
ck("T6 rettangoli a 2 px si toccano, a 3 px no", got == {(1, 2)}, f"coppie {sorted(got)}")
rng = np.random.default_rng(5); polys, ids = [], []
while len(polys) < 120:
    c = rng.uniform(10, 290, 2); r = rng.uniform(2.5, 5)
    g = shapely.Point(*c).buffer(r, 32)
    if all(g.distance(q) > 0.05 for q in polys): polys.append(g); ids.append(len(polys))
labr, ncf = MS.RAS.rasterize(polys, np.array(ids), (300, 300), 1.0); m = MS.RAS.measure(labr, 1.0, solidity=False, frac_out=False)
pr = MS.contacts(labr); anyc = np.isin(m.label.to_numpy(), np.r_[pr.label_i.to_numpy(), pr.label_j.to_numpy()])
ck("T6 contatto di coppia = n_c2 > 0 per ogni nucleo (120 dischi casuali)", np.array_equal(anyc, m.n_c2.to_numpy() > 0) and anyc.sum() > 5,
   f"{int(anyc.sum())} nuclei a contatto su {len(m)}")
print(f"RESULT {sum(res)}/{len(res)} PASS")
