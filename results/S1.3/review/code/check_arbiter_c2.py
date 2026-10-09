#!/usr/bin/env python3
"""results/S1.3/review/code/check_arbiter_c2.py — revisione del codice S1.3 (RA-code), compiti 4 e 5.
(4) ricalcolo con Qhull (scipy.spatial.Voronoi, insieme completo dei punti + 8 fittizi lontani) delle 26 righe di
    results/S1.3/S1.3_c10a_arbitration.csv: area, numero di vertici, lato minimo della cella vera (non ritagliata).
(5) A3 r3 CSR replica 3, celle 766 e 1987: geometria, predicati e intersezione con shapely (GEOS 3.14), campionamento.
"""
import numpy as np, pandas as pd, shapely, pyarrow.parquet as pq
from shapely.geometry import Polygon
from scipy.spatial import Voronoi
REPO = "/home/user/2025.geo_spatialtrans"; D = "/mnt/micron/geo_spatialtrans/S1.3/review_code_data"
OUT = REPO + "/results/S1.3/review/code"

def full_cells(x, y, cells):
    L = max(np.ptp(x), np.ptp(y)); cx, cy = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2; R = 3 * L + 10
    ang = np.arange(8) * np.pi / 4
    pts = np.vstack([np.c_[x, y], np.c_[cx + R * np.cos(ang), cy + R * np.sin(ang)]])
    vor = Voronoi(pts); res = {}
    for c in cells:
        i = c - 1; reg = vor.regions[vor.point_region[i]]; V = vor.vertices[reg]; m = V.mean(0)
        V = V[np.argsort(np.arctan2(V[:, 1] - m[1], V[:, 0] - m[0]))]
        e = np.hypot(*(np.roll(V, -1, 0) - V).T)
        a = 0.5 * abs(np.dot(V[:, 0], np.roll(V[:, 1], -1)) - np.dot(np.roll(V[:, 0], -1), V[:, 1]))
        res[c] = (a, len(V), e.min(), int((e > 1e-9).sum()))
    return res

ar = pd.read_csv(REPO + "/results/S1.3/S1.3_c10a_arbitration.csv", keep_default_na=False); rows = []   # "null" non e NA; rows = []
for (st, A, roi, mod), g in ar.groupby(["set", "archetype", "roi_id", "model"]):
    if st == "real":
        t = pq.read_table(f"/mnt/micron/geo_spatialtrans/R3/real/{A}_{roi}_{mod}_gen.parquet").to_pandas(); x, y = t.x.values, t.y.values
    else:
        t = pd.read_csv(f"{D}/arb_null_{A}_{roi}_{mod}.csv"); x, y = t.x.values, t.y.values
    res = full_cells(x, y, g.cell.values)
    for _, r in g.iterrows():
        a, nv, me, ne = res[r.cell]
        rows.append(dict(set=st, archetype=A, roi_id=roi, model=mod, cell=r.cell, a_qhull=a, nvert_qhull=nv, n_edges_gt1e9=ne,
                         min_edge_qhull=me, a_geos=r.a_geos, a_deldir=r.a_deldir, a_hp=r.a_hp,
                         rel_geos=abs(r.a_geos - a) / a, rel_deldir=abs(r.a_deldir - a) / a, rel_hp=abs(r.a_hp - a) / a,
                         ns_geos=r.ns_geos, ns_deldir=r.ns_deldir, ns_hp=r.ns_hp, min_edge_hp=r.min_edge,
                         geos_ok_report=r.geos_ok, deldir_ok_report=r.deldir_ok))
out = pd.DataFrame(rows); out.to_csv(OUT + "/arbiter_recompute_qhull.csv", index=False)
print("righe ricalcolate:", len(out))
print("max rel |a_geos - a_qhull|:", out.rel_geos.max(), " max rel |a_hp - a_qhull|:", out.rel_hp.max(), " max rel |a_deldir - a_qhull|:", out.rel_deldir.max())
print("righe con rel_deldir > 1e-9:", int((out.rel_deldir > 1e-9).sum()), " rel_geos > 1e-9:", int((out.rel_geos > 1e-9).sum()))
print(out[["archetype", "roi_id", "model", "cell", "nvert_qhull", "n_edges_gt1e9", "ns_geos", "ns_deldir", "ns_hp", "min_edge_qhull", "min_edge_hp"]].to_string())

# ---- (5) A3 r3 CSR 3 ---------------------------------------------------------------------------------------
tv = pd.read_csv(f"{D}/A3_r3_CSR03_tv.csv"); reg = pd.read_csv(f"{D}/A3_r3_CSR03_reg.csv"); cen = pd.read_csv(f"{D}/A3_r3_CSR03_cen.csv")
T = shapely.from_wkb(tv.terr_wkb.values); G = shapely.from_wkb(reg.wkb.values[0])
a, b = T[765], T[1986]
print("\n(5) cella 766: valida", a.is_valid, "anelli interni", len(a.interiors), "vertici", len(a.exterior.coords), "area", repr(a.area))
print("    cella 1987: valida", b.is_valid, "anelli interni", len(b.interiors), "vertici", len(b.exterior.coords), "area", repr(b.area))
print("    shapely relate(766,1987):", a.relate(b), " intersection area:", repr(a.intersection(b).area), " within:", a.within(b))
print("    766 dentro un buco di 1987?", [i for i, h in enumerate(b.interiors) if Polygon(h).buffer(1e-9).contains(a)],
      " area(766 Δ buco):", [repr(Polygon(h).symmetric_difference(a).area) for h in b.interiors if Polygon(h).intersects(a)])
rng = np.random.default_rng(1); bx = a.bounds; P = []
while len(P) < 5000:
    q = np.c_[rng.uniform(bx[0], bx[2], 20000), rng.uniform(bx[1], bx[3], 20000)]
    q = q[shapely.contains_xy(a, q[:, 0], q[:, 1])]; P.extend(q.tolist())
P = np.array(P[:5000]); inb = shapely.contains_xy(b, P[:, 0], P[:, 1])
print("    5000 punti in 766: dentro 1987 =", int(inb.sum()))
# copertura: quanti territori coprono 2e5 punti uniformi nella regione (predicato, non operazioni d'insieme)
bx = G.bounds; P = []
while len(P) < 200000:
    q = np.c_[rng.uniform(bx[0], bx[2], 400000), rng.uniform(bx[1], bx[3], 400000)]
    q = q[shapely.contains_xy(G, q[:, 0], q[:, 1])]; P.extend(q.tolist())
P = np.array(P[:200000]); tree = shapely.STRtree(T); qi, ti = tree.query(shapely.points(P), predicate="within")
cnt = np.bincount(qi, minlength=len(P))
print("    copertura 2e5 punti: 0 ->", int((cnt == 0).sum()), " 1 ->", int((cnt == 1).sum()), " >=2 ->", int((cnt >= 2).sum()))
s = shapely.union_all(T); print("    shapely union_all: area unione", repr(s.area), " area regione", repr(G.area), " Σ aree", repr(float(shapely.area(T).sum())),
      " rel (Σ - unione)/regione", (float(shapely.area(T).sum()) - s.area) / G.area, " rel area(unione Δ regione)", s.symmetric_difference(G).area / G.area)
pairs = tree.query(T, predicate="overlaps"); pr = [(i, j) for i, j in zip(*pairs) if i < j]
print("    coppie con predicato overlaps (shapely):", len(pr), [(i + 1, j + 1, repr(T[i].intersection(T[j]).area)) for i, j in pr][:10])
