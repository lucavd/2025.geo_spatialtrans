"""review/exact/f4_check.py — punto (b): celle F4 invalide. Per ogni caso con c5$info$n_invalid > 0 (c5_info.csv da dump_viol.R) trova la cella
invalida in geom/<key>_f4cells.parquet (validità GEOS via shapely; il revisore non usa l'elenco del run perché non è salvato), la ragione,
il punto; verifica in esatto (Fraction) che: (1) l'anello della cella passa due volte per quel punto; (2) il punto è un vertice della regione
in cui il bordo della regione tocca sé stesso (stesso punto due volte nell'anello della regione, o due anelli che si toccano); (3) area esatta
della cella F4 vs area esatta della cella V0 (geom/<key>_cells.parquet) e tau_i = g*P_i. Uscita: review_exact/f4/f4_invalid.csv"""
import sys, os, json, collections; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import pandas as pd, numpy as np, pyarrow.parquet as pq, shapely, rx
from shapely.geometry import Polygon, MultiPolygon
D = "/mnt/micron/geo_spatialtrans/S1.3d"; O = f"{D}/review_exact/f4"; os.makedirs(O, exist_ok=True)
C5 = pd.read_csv(f"{D}/review_exact/c5_info.csv"); C5 = C5[C5.n_invalid > 0]
def shp(t):
    out = {}
    for (ci, pa), d in t.groupby(["cell", "part"]):
        rings = [g_[["x", "y"]].values for _, g_ in d.sort_values(["ring", "k"]).groupby("ring")]
        out.setdefault(int(ci), []).append(Polygon(rings[0], rings[1:]))
    return {c: (p[0] if len(p) == 1 else MultiPolygon(p)) for c, p in out.items()}
rows = []
for key in C5.key:
    t4 = pq.read_table(f"{D}/geom/{key}_f4cells.parquet").to_pandas(); P4 = shp(t4)
    inv = [c for c, p in P4.items() if not shapely.is_valid(p)]
    meta = json.load(open(f"{D}/geom/{key}_meta.json"))
    reg = pq.read_table(f"{D}/geom/{key}_region.parquet").to_pandas()
    regpts = collections.Counter(zip(reg.x.values, reg.y.values))
    for ci in inv:
        reason = shapely.is_valid_reason(P4[ci]); xy = reason[reason.find("[") + 1: reason.find("]")].split()
        px, py = float(xy[0]), float(xy[1])
        d = t4[t4.cell == ci]; cellpts = collections.Counter(zip(d.x.values, d.y.values))
        # vertici ripetuti nell'anello della cella F4 (anelli non chiusi nel parquet: un punto ripetuto = passaggio doppio)
        rep = [(p, n) for p, n in cellpts.items() if n > 1]
        rp = rep[0][0] if len(rep) == 1 else (px, py)
        on_reg = regpts.get(rp, 0); rr = reg[(reg.x == rp[0]) & (reg.y == rp[1])]
        reg_rings_at = str(sorted(set(zip(rr.part.tolist(), rr.ring.tolist())))); near = abs(rp[0] - px) + abs(rp[1] - py)
        R4 = rx.read_polys(f"{D}/geom/{key}_f4cells.parquet", want={ci})[ci]; R0 = rx.read_polys(f"{D}/geom/{key}_cells.parquet", want={ci})[ci]
        a4 = abs(rx.shoelace(R4)); a0 = rx.area(R0); tau = 1e-6 * meta["P"][ci - 1]
        rows.append(dict(key=key, cell=ci, n_invalid_found=len(inv), n_invalid_run=int(C5[C5.key == key].n_invalid.iloc[0]), reason=reason,
                         n_rings=int(d.groupby(["part", "ring"]).ngroups), n_parts=int(d.part.nunique()), n_vert=len(d),
                         repeated_vertex=str(rep), n_repeated_vertices=len(rep), rep_x=rp[0], rep_y=rp[1], dist_reason_point=near, rep_times_in_cell=cellpts.get(rp, 0), rep_times_in_region=on_reg, region_part_ring_at_rep=reg_rings_at,
                         region_rings=int(reg.groupby(["part", "ring"]).ngroups), region_has_repeated_pts=sum(1 for n in regpts.values() if n > 1),
                         area_f4_exact=float(a4), area_v0_exact=float(a0), area_equal_exact=(a4 == a0), d_area=float(a4 - a0), tau_i=tau, d_over_tau=float(abs(a4 - a0)) / tau))
R = pd.DataFrame(rows); R.to_csv(f"{O}/f4_invalid.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 60)
print(R.drop(columns=["repeated_vertex"]).to_string()); print(R.repeated_vertex.tolist())
