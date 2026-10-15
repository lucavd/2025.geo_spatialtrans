"""review/exact/review_sample.py — (1) verifica esatta indipendente di pezzi O2 (2 µm) e coppie O1 campionati a caso contro truth/*.parquet.
Uso: python review_sample.py o2 <key> <seed>   |   python review_sample.py o1 <key> <seed>    (un nullo per processo; parallelo da run_sample.sh)
Uscite: /mnt/micron/geo_spatialtrans/S1.3d/review_exact/sample/<op>_<key>.csv"""
import sys, os, json, random, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np, pandas as pd, pyarrow.parquet as pq, rx
D = "/mnt/micron/geo_spatialtrans/S1.3d"; OUT = f"{D}/review_exact/sample"; os.makedirs(OUT, exist_ok=True)
op, key, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]); rng = random.Random(seed); t0 = time.time()
meta = json.load(open(f"{D}/geom/{key}_meta.json")); gs = meta["gs"]
rows = []
if op == "o2":
    T2 = pq.read_table(f"{D}/truth/{key}_o2_s2.parquet").to_pandas(); T8 = pq.read_table(f"{D}/truth/{key}_o2_s8.parquet").to_pandas()
    pick = sorted(rng.sample(range(len(T2)), 30))                      # 30 pezzi a caso
    full = rng.sample(sorted(set(T2.i)), 3)                            # 3 celle: tutti i pezzi a 2 e 8 µm (pezzi mancanti/spuri)
    cells = rx.read_polys(f"{D}/geom/{key}_cells.parquet", want=set(T2.i.iloc[pick]) | set(full))
    nx2, ny2 = gs["nx8"] * 4, gs["ny8"] * 4
    for q in pick:
        r = T2.iloc[q]; ci, ix, iy = int(r.i), int(r.ix), int(r.iy); s = 2
        a = rx.area_box(cells[ci], F(gs["X0"]) + s * ix, F(gs["X0"]) + s * (ix + 1), F(gs["Y0"]) + s * iy, F(gs["Y0"]) + s * (iy + 1))
        eq, rel = rx.cmpf(a, float(r.area)); rows.append(dict(kind="piece_s2", key=key, i=ci, ix=ix, iy=iy, exact=float(a), truth=float(r.area), eq=eq, rel=rel))
    for ci in full:
        for s, T, nx, ny in ((2, T2, nx2, ny2), (8, T8, gs["nx8"], gs["ny8"])):
            mine = rx.pieces(cells[ci], gs["X0"], gs["Y0"], s, nx, ny); tt = T[T.i == ci]
            td = {(int(a), int(b)): float(v) for a, b, v in zip(tt.ix, tt.iy, tt.area)}
            ks = set(mine) | set(td)
            for k in ks:
                ex = mine.get(k, F(0)); tr = td.get(k, None)
                eq, rel = rx.cmpf(ex, tr) if tr is not None else (ex == 0, float("inf"))
                rows.append(dict(kind=f"cellfull_s{s}", key=key, i=ci, ix=k[0], iy=k[1], exact=float(ex), truth=(np.nan if tr is None else tr), eq=eq, rel=rel,
                                 missing_in_truth=(tr is None), spurious_in_truth=(k not in mine)))
            # somma esatta dei pezzi = area esatta della cella (controllo interno del revisore)
            rows.append(dict(kind=f"cellsum_s{s}", key=key, i=ci, exact=float(sum(mine.values())), truth=float(rx.area(cells[ci])), eq=(sum(mine.values()) == rx.area(cells[ci])), rel=0.0))
elif op == "o1":
    TP = pq.read_table(f"{D}/truth/{key}_o1_pairs.parquet").to_pandas(); TC = pq.read_table(f"{D}/truth/{key}_o1_cells.parquet").to_pandas()
    nz = list(np.nonzero(~TP.zero.values)[0]); z = list(np.nonzero(TP.zero.values)[0])
    pick = rng.sample(nz, min(5, len(nz))); pick += rng.sample(z, 10 - len(pick))
    need = set(TP.i.iloc[pick]) | set(TP.j.iloc[pick]); cells = rx.read_polys(f"{D}/geom/{key}_cells.parquet", want=need)
    for q in pick:
        i, j = int(TP.i.iloc[q]), int(TP.j.iloc[q]); a = rx.area_inter(cells[i], cells[j]); b = rx.area_inter(cells[j], cells[i])
        eq, rel = rx.cmpf(a, float(TP.a.iloc[q]))
        rows.append(dict(kind="pair", key=key, i=i, j=j, exact=float(a), truth=float(TP.a.iloc[q]), eq=eq, rel=rel, sym=(a == b),
                         zero_ok=((a == 0) == bool(TP.zero.iloc[q])), truth_nonzero=not bool(TP.zero.iloc[q])))
    tc = TC.set_index("cell")
    for ci in sorted(need):
        a = rx.area(cells[ci]); sh = abs(rx.shoelace(cells[ci])); eq, rel = rx.cmpf(a, float(tc.a.loc[ci]))
        rows.append(dict(kind="cell_area", key=key, i=ci, exact=float(a), truth=float(tc.a.loc[ci]), eq=eq, rel=rel, sym=(a == sh)))
pd.DataFrame(rows).to_csv(f"{OUT}/{op}_{key}.csv", index=False)
print(op, key, len(rows), "eq", sum(r["eq"] for r in rows), f"{time.time()-t0:.1f}s")
