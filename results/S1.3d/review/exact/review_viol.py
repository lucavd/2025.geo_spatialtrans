"""review/exact/review_viol.py — ricalcolo esatto indipendente (rx) di OGNI violazione registrata in res/*.rds (estratte da dump_viol.R):
coppie e celle di O1 (V0), pezzi O2 a 2 e 8 µm (F2, F3, F5, anche S). Uso: python review_viol.py <key>  -> review_exact/viol/<key>.csv"""
import sys, os, json, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np, pandas as pd, pyarrow.parquet as pq, rx
D = "/mnt/micron/geo_spatialtrans/S1.3d"; OUT = f"{D}/review_exact/viol"; os.makedirs(OUT, exist_ok=True)
key = sys.argv[1]; t0 = time.time()
V = pd.read_parquet(f"{D}/review_exact/viol_all.parquet"); V = V[(V.key == key) & (V.what.isin(["pair_ij", "pair_ji", "cell"]) | V.k.notna())]
meta = json.load(open(f"{D}/geom/{key}_meta.json")); gs = meta["gs"]
o2 = V.k.notna(); ids = set()
for _, r in V[~o2].iterrows(): ids |= {int(r.i)} | ({int(r.j)} if pd.notna(r.j) else set())
for _, r in V[o2].iterrows(): ids.add(int(str(r.k).split()[0]))
cells = rx.read_polys(f"{D}/geom/{key}_cells.parquet", want=ids)
rows = []
for _, r in V.iterrows():
    if pd.isna(r.k):
        i = int(r.i)
        if r.what == "cell": ex = rx.area(cells[i]); ent = dict(i=i, j=None)
        else: j = int(r.j); ex = rx.area_inter(cells[i], cells[j]); ent = dict(i=i, j=j)
        eq, rel = rx.cmpf(ex, float(r.truth))
        rows.append(dict(key=key, vname=r.vname, what=r.what, **ent, s=None, exact=float(ex), truth=float(r.truth), cand=float(r.cand_val), tau=float(r.tau),
                         eq=eq, rel=rel, err_over_tau=abs(float(F(float(r.cand_val)) - ex)) / float(r.tau), real=abs(float(F(float(r.cand_val)) - ex)) > float(r.tau)))
    else:
        ci, ix, iy = map(int, str(r.k).split()); s = 8 if (r.vname.endswith("_8") or r.vname.startswith("o3") or r.file.endswith("s8.rds")) else 2
        ex = rx.area_box(cells[ci], F(gs["X0"]) + s * ix, F(gs["X0"]) + s * (ix + 1), F(gs["Y0"]) + s * iy, F(gs["Y0"]) + s * (iy + 1))
        tr = 0.0 if pd.isna(r.area_t) else float(r.area_t); eq, rel = rx.cmpf(ex, tr)
        tau = 1e-6 * (float(r.perim_t) if pd.notna(r.perim_t) else (float(r.perim) if pd.notna(r.perim) else 4 * s))
        rows.append(dict(key=key, vname=r.vname if r.vname != "cmp" else r.file, what="piece", i=ci, j=None, s=s, ix=ix, iy=iy, exact=float(ex), truth=tr,
                         truth_missing=pd.isna(r.area_t), cand=float(r.area), tau=tau, eq=eq, rel=rel,
                         err_over_tau=abs(float(F(float(r.area)) - ex)) / tau, real=abs(float(F(float(r.area)) - ex)) > tau))
pd.DataFrame(rows).to_csv(f"{OUT}/{key}.csv", index=False)
print(key, len(rows), "eq", sum(x["eq"] for x in rows), "real", sum(x["real"] for x in rows), f"{time.time()-t0:.1f}s")
