"""review/exact/f4_makevalid_exact.py — area esatta (rx) delle celle F4 dopo st_make_valid vs area esatta V0 e tau_i; stessi vertici?"""
import sys, os, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, pyarrow.parquet as pq, rx, shapely
from shapely.geometry import Polygon
D = "/mnt/micron/geo_spatialtrans/S1.3d"; O = f"{D}/review_exact/f4"
R = pd.read_csv(f"{O}/f4_invalid.csv"); M = pq.read_table(f"{O}/f4_makevalid_coords.parquet").to_pandas(); S = pd.read_csv(f"{O}/f4_makevalid_summary.csv")
rows = []
for _, r in R.iterrows():
    m = M[(M.key == r.key) & (M.cell == r.cell)].rename(columns={})
    m = m.assign(cell=1); tmp = f"{O}/_tmp_{r.key}.parquet"; pq.write_table(__import__("pyarrow").Table.from_pandas(m[["cell", "part", "ring", "k", "x", "y"]].astype({"cell": "int32", "part": "int32", "ring": "int32", "k": "int32"})), tmp)
    Rm = rx.read_polys(tmp)[1]; os.remove(tmp)
    a_mv = rx.area(Rm); R0 = rx.read_polys(f"{D}/geom/{r.key}_cells.parquet", want={int(r.cell)})[int(r.cell)]; a0 = rx.area(R0)
    t4 = pq.read_table(f"{D}/geom/{r.key}_f4cells.parquet").to_pandas(); t4 = t4[t4.cell == r.cell]
    same_pts = set(zip(t4.x, t4.y)) == set(zip(m.x, m.y))
    rings = [g[["x", "y"]].values for _, g in m.groupby("ring")]; P = Polygon(rings[0], rings[1:])
    rows.append(dict(key=r.key, cell=r.cell, valid_after_shapely=bool(shapely.is_valid(P)), n_rings_after=len(rings), same_vertex_set=same_pts,
                     area_exact_after=float(a_mv), area_exact_v0=float(a0), equal_exact=(a_mv == a0), d_over_tau=float(abs(a_mv - a0)) / r.tau_i, tau_i=r.tau_i))
X = pd.DataFrame(rows); X.to_csv(f"{O}/f4_makevalid_exact.csv", index=False); print(X.to_string())
