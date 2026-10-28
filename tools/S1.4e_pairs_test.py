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
# ---- T5b: poligono in piu' parti (revisione G-5): centroide e allungamento con un'origine unica ----
mp2 = shapely.MultiPolygon([shapely.box(0, 0, 1, 1), shapely.box(10, 0, 11, 1)]); S2 = MS.poly_shape(mp2)
ck("T5b due quadrati separati: centroide (5.5, 0.5) e asse lungo x, molto allungato", abs(S2[0] - 5.5) < 1e-9 and abs(S2[1] - 0.5) < 1e-9 and S2[3] > 5
   and abs(np.sin(S2[2])) < 1e-9, f"centroide ({S2[0]:.3f}, {S2[1]:.3f}), allungamento {S2[3]:.2f}")
# ---- T9: regola dei focali e coppie con un estremo focale ----
import pandas as pd
wr = pd.Series(dict(x0_um=0.0, y0_um=0.0, side_px=100, um_per_px=1.0)); P = np.array([[50, 50], [54, 50], [5, 50], [9, 50], [50, 90], [53, 92], [80, 80], [84, 80]], float)
ex = shapely.box(70, 70, 95, 95)
foc = MS.focal_mask(P, wr, 10.0, ex)
ck("T9 focali: a >= d_max dai lati e dalle esclusioni", foc.tolist() == [True, True, False, False, True, False, False, False], f"{foc.astype(int).tolist()}")
ids = np.arange(1, 9); M_ = 10; cset = np.array([1 * M_ + 2, 3 * M_ + 4, 7 * M_ + 8])
pc = MS.pair_counts(P, ids, foc, cset, 10.0, M_)
ck("T9 coppie con almeno un estremo focale (< 6 um), contatti dalle coppie del raster", pc["lt6"] == (2, 1) and pc["6to10"] == (0, 0), f"{pc}")
# ---- T10: verdetti e regola N su una tabella costruita a mano ----
sp2 = importlib.util.spec_from_file_location("s14e_analysis", SRC / "tools/S1.4e_analysis.py"); AN = importlib.util.module_from_spec(sp2); sp2.loader.exec_module(AN)
ck("T10 verdetti: dentro, fuori, a cavallo", [AN.verdict(0.5, 0.01, 0.4, 0.6), AN.verdict(0.8, 0.01, 0.4, 0.6), AN.verdict(0.59, 0.01, 0.4, 0.6)] == ["confermata", "smentita", "indeterminata"])
def tab(conf):
    r = []
    for arm, n in conf.items():
        for k in range(10): r.append(dict(id="P-1" if k < 6 else "P-2", archetype=f"A{k % 6 + 1}", arm=arm, verdict="confermata" if k < n else "smentita"))
    return pd.DataFrame(r)
w1 = AN.rule_n(tab({"K2": 3, "ER": 3, "ET": 4}))[0]; w2 = AN.rule_n(tab({"K2": 3, "ER": 3, "ET": 5}))[0]; w3 = AN.rule_n(tab({"K2": 3, "ER": 6, "ET": 5}))[0]
ck("T10 regola N: +1 resta K2; +2 vince ET; il punteggio piu' alto vince", (w1, w2, w3) == ("K2", "ET", "ER"), f"{w1} {w2} {w3}")
# ---- T11: mappa del banco A (finestra del cuore con esclusioni): il nucleo di ogni punto di Luca contiene il punto ----
mp = pd.read_parquet(MS.M / "A" / "map" / "A6_w4.parquet"); mm = MS.luca_to_cell(mp); pts = pd.read_csv(ROOT / "results/R2b/R2b_points.csv")
pts = pts[(pts.rater == "Luca") & (pts.win_id == "A6_w4")].set_index("idx")
nuc = pd.read_parquet(sorted((MS.M / "A" / "nuc" / "K2").glob("A6_*_A6_w4_s01.parquet"))[0]); G = shapely.from_wkb(nuc.wkb.to_numpy())
inside = [bool(shapely.contains_xy(G[int(cid) - 1], pts.loc[i].x_um, pts.loc[i].y_um)) for i, cid in mm.items()]
ck("T11 mappa punto di Luca -> nucleo (cuore, finestra 4, 11 esclusioni)", len(inside) > 50 and all(inside), f"{sum(inside)}/{len(inside)}")
print(f"RESULT {sum(res)}/{len(res)} PASS")
