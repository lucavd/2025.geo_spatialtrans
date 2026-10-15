# review/exact/test_rx.py — test del modulo rx (casi analitici + confronto con shapely in virgola mobile su celle reali)
import sys, os, random; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import rx, shapely, numpy as np
from shapely.geometry import Polygon, box
ok = 0; bad = 0
def chk(name, cond):
    global ok, bad
    print(("PASS " if cond else "FAIL ") + name); ok += cond; bad += (not cond)
sq = lambda x0, y0, x1, y1: [[(F(x0), F(y0)), (F(x1), F(y0)), (F(x1), F(y1)), (F(x0), F(y1))]]
chk("area quadrato 1", rx.area(sq(0, 0, 1, 1)) == 1)
chk("quadrato ∩ quadrato traslato = 1/4", rx.area_inter(sq(0, 0, 1, 1), sq(F(1, 2), F(1, 2), 2, 2)) == F(1, 4))
chk("quadrati adiacenti = 0", rx.area_inter(sq(0, 0, 1, 1), sq(1, 0, 2, 1)) == 0)
chk("identici = 1", rx.area_inter(sq(0, 0, 1, 1), sq(0, 0, 1, 1)) == 1)
tri = [[(F(0), F(0)), (F(4), F(0)), (F(0), F(4))]]
chk("triangolo ∩ bin [0,2]^2 = 4", rx.area_box(tri, 0, 2, 0, 2) == 4)
chk("triangolo ∩ bin [2,4]x[0,2] = 2", rx.area_box(tri, 2, 4, 0, 2) == 2)
chk("triangolo ∩ bin [2,4]x[2,4] = 0", rx.area_box(tri, 2, 4, 2, 4) == 0)
hole = sq(0, 0, 4, 4) + [[(F(1), F(1)), (F(1), F(3)), (F(3), F(3)), (F(3), F(1))]]
chk("bucato area 12", rx.area(hole) == 12 and rx.shoelace(hole) == 12)
chk("bucato ∩ bin [0,2]^2 = 3", rx.area_box(hole, 0, 2, 0, 2) == 3)
chk("bucato ∩ quadrato nel buco = 0", rx.area_inter(hole, sq(F(3, 2), F(3, 2), 2, 2)) == 0)
rot = [[(F(0), F(-1)), (F(1), F(0)), (F(0), F(1)), (F(-1), F(0))]]
chk("rombo ∩ rombo ruotato (ottagono)", rx.area_inter(rot, sq(F(-1, 2), F(-1, 2), F(1, 2), F(1, 2))) == 1)
star = [[(F(0), F(0)), (F(2), F(0)), (F(2), F(2)), (F(1), F(1, 2)), (F(0), F(2))]]          # non convesso (incavo)
chk("incavo: shoelace = slab", rx.shoelace(star) == rx.area(star))
chk("incavo ∩ bin [0,1]x[0,1] = 11/12 (analitico)", rx.area_box(star, 0, 1, 0, 1) == F(11, 12))
v = rx.area_box(star, 0, 1, 0, 1)
P = Polygon([(float(a), float(b)) for a, b in star[0]]); chk("incavo ∩ bin vs shapely", abs(float(v) - P.intersection(box(0, 0, 1, 1)).area) < 1e-15)
# celle reali: slab vs shoelace, e intersezioni/pezzi vs shapely (tolleranza 1e-9 relativa: shapely è in virgola mobile)
D = "/mnt/micron/geo_spatialtrans/S1.3d/geom"; key = "A3_r3_CSR_03"
cells = rx.read_polys(f"{D}/{key}_cells.parquet", want=set(range(1, 60)))
nsl = sum(rx.area(r) == abs(rx.shoelace(r)) for r in cells.values()); chk(f"slab = |shoelace| su {len(cells)} celle reali", nsl == len(cells))
random.seed(1); worst = 0
for ci in list(cells)[:20]:
    pc = rx.pieces(cells[ci], 0, 0, 2, 10**6, 10**6); Pg = Polygon([(float(a), float(b)) for a, b in cells[ci][0]])
    tot = sum(pc.values()); worst = max(worst, abs(float(tot - rx.area(cells[ci]))))
    for (ix, iy), a in pc.items():
        sa = Pg.intersection(box(2 * ix, 2 * iy, 2 * ix + 2, 2 * iy + 2)).area; worst = max(worst, abs(float(a) - sa))
chk(f"pezzi vs shapely e somma pezzi = area (max |d| = {worst:.2e})", worst < 1e-11)
import pyarrow.parquet as pq
pt = pq.read_table(f"{D}/{key}_pairs.parquet").to_pandas(); pt = pt[(pt.i < 60) & (pt.j < 60)]
w2 = 0; nsym = 0
for i, j in zip(pt.i, pt.j):
    a = rx.area_inter(cells[i], cells[j]); b = rx.area_inter(cells[j], cells[i]); nsym += (a == b)
    sa = Polygon([(float(p), float(q)) for p, q in cells[i][0]]).intersection(Polygon([(float(p), float(q)) for p, q in cells[j][0]])).area; w2 = max(w2, abs(float(a) - sa))
chk(f"coppie: simmetria esatta {nsym}/{len(pt)}, vs shapely max |d| = {w2:.2e}", nsym == len(pt) and w2 < 1e-11)
print(f"TOTALE PASS {ok} FAIL {bad}")
