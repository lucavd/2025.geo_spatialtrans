"""review/exact/sliver.py — punto (d): geometria esatta dei pezzi «scheggia» (estensione x/y del pezzo esatto, larghezza media = area/(perimetro/2)),
e riproduzione del valore del candidato F2 arrotondando la cella alla griglia g come GEOS (round(v*1e6)/1e6) e ricalcolando il pezzo in esatto."""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import pandas as pd, numpy as np, rx
D = "/mnt/micron/geo_spatialtrans/S1.3d"; g = 1e-6
def extent_x(rings, x0, x1, y0, y1):
    E = rx.edges(rings); x0, x1, y0, y1 = F(x0), F(x1), F(y0), F(y1); xs = {x0, x1}
    for a, b in E:
        if x0 < a[0] < x1: xs.add(a[0])
    xs.update(rx._cross_x_hline(E, y0, x0, x1)); xs.update(rx._cross_x_hline(E, y1, x0, x1)); xs = sorted(xs); on = []
    for q in range(len(xs) - 1):
        if rx._len_clip(rx._ys_at(E, (xs[q] + xs[q + 1]) / 2), y0, y1) > 0: on.append((xs[q], xs[q + 1]))
    return (on[0][0], on[-1][1]) if on else (None, None)
def swap(rings): return [[(y, x) for x, y in r] for r in rings]
cases = [("A4_p2_CSR_62", 8601, 380, 469, 2.081084e-14, "F2/F3 (H)"), ("A4_f2_RSArule_16", 9871, 447, 434, None, "F5 (N)"),
         ("A4_p1_RSArule_06", 27810, 496, 442, None, "F5 (N)"), ("A4_p1_RSArule_18", 9781, 419, 439, None, "F5 (N)")]
rows = []
for key, ci, ix, iy, cand, lab in cases:
    gs = json.load(open(f"{D}/geom/{key}_meta.json"))["gs"]; s = 2
    cell = rx.read_polys(f"{D}/geom/{key}_cells.parquet", want={ci})[ci]
    bx = (F(gs["X0"]) + s * ix, F(gs["X0"]) + s * (ix + 1), F(gs["Y0"]) + s * iy, F(gs["Y0"]) + s * (iy + 1))
    a = rx.area_box(cell, *bx); ex = extent_x(cell, *bx); ey = extent_x(swap(cell), bx[2], bx[3], bx[0], bx[1])
    T = pd.read_parquet(f"{D}/truth/{key}_o2_s2.parquet"); tr = T[(T.i == ci) & (T.ix == ix) & (T.iy == iy)]
    per = float(tr.perim.iloc[0]); w = float(a) / (per / 2)
    # celle vicine che condividono il pezzo di bin: lato/vertice della cella più vicino all'angolo del bin
    nv_in = sum(1 for r_ in cell for p in r_ if bx[0] <= p[0] <= bx[1] and bx[2] <= p[1] <= bx[3])
    rnd = [[(F(round(float(x) * 1e6) / 1e6), F(round(float(y) * 1e6) / 1e6)) for x, y in r_] for r_ in cell]
    a_r = rx.area_box(rnd, *bx)
    rows.append(dict(key=key, cell=ci, ix=ix, iy=iy, label=lab, exact=float(a), exact_eq_truth=(float(a) == float(tr.area.iloc[0])), perim_exact=per,
                     ext_x=float(ex[1] - ex[0]), ext_y=float(ey[1] - ey[0]), ext_x_over_g=float(ex[1] - ex[0]) / g, ext_y_over_g=float(ey[1] - ey[0]) / g,
                     mean_width=w, mean_width_over_g=w / g, area_over_g2=float(a) / g**2, verts_in_bin=nv_in, cand=cand,
                     rounded_cell_piece=float(a_r), rounded_matches_cand=(None if cand is None else abs(float(a_r) - cand) <= 1e-6 * cand), tau=g * per))
R = pd.DataFrame(rows); R.to_csv(f"{D}/review_exact/sliver.csv", index=False); pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(R.T.to_string())
