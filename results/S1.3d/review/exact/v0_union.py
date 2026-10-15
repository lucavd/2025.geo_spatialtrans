"""review/exact/v0_union.py — punto (c): per un nullo in cui V0 sbaglia l'unione, verità esatta indipendente (rx):
(1) validità esatta di ogni cella (anelli semplici: nessun incrocio/tocco fra lati non adiacenti, test di orientazione in Fraction);
(2) Σ a_i esatta; (3) area esatta di tutte le coppie dell'elenco e di TUTTE le coppie con riquadri sovrapposti non in elenco (completezza);
(4) limiti di Bonferroni Σa − Σa_ij ≤ U ≤ Σa; (5) deficit dell'unione V0 e cella(e) la cui area esatta lo spiega.
Uso: python v0_union.py <key>  -> review_exact/v0/<key>.json"""
import sys, os, json, time, itertools; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np, pandas as pd, pyarrow.parquet as pq, shapely, rx
D = "/mnt/micron/geo_spatialtrans/S1.3d"; O = f"{D}/review_exact/v0"; os.makedirs(O, exist_ok=True)
key = sys.argv[1]; t0 = time.time(); meta = json.load(open(f"{D}/geom/{key}_meta.json")); n = meta["n"]
cells = rx.read_polys(f"{D}/geom/{key}_cells.parquet")
def orient(a, b, c): v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]); return (v > 0) - (v < 0)
def onseg(a, b, p): return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])
def seg_touch(p, q, r, s):
    o1, o2, o3, o4 = orient(p, q, r), orient(p, q, s), orient(r, s, p), orient(r, s, q)
    if o1 != o2 and o3 != o4 and 0 not in (o1, o2, o3, o4): return True
    return (o1 == 0 and onseg(p, q, r)) or (o2 == 0 and onseg(p, q, s)) or (o3 == 0 and onseg(r, s, p)) or (o4 == 0 and onseg(r, s, q))
def ring_simple(rg):
    m = len(rg); E = [(rg[k], rg[(k + 1) % m]) for k in range(m)]
    if len(set(rg)) != m: return False
    for a in range(m):
        for b in range(a + 1, m):
            if b == a + 1 or (a == 0 and b == m - 1):          # adiacenti: non devono essere collineari sovrapposti
                p, q = E[a]; r, s = E[b]; sh = q if b == a + 1 else p; ot = s if b == a + 1 else r; ot0 = p if b == a + 1 else q
                if orient(ot0, sh, ot) == 0 and ((ot[0] - sh[0]) * (ot0[0] - sh[0]) + (ot[1] - sh[1]) * (ot0[1] - sh[1])) > 0: return False
                continue
            if seg_touch(*E[a], *E[b]): return False
    return True
nonsimple = [c for c, R in cells.items() if not all(ring_simple(r) for r in R)]
multi = [c for c, R in cells.items() if len(R) > 1]
a = {c: rx.area(R) for c, R in cells.items()}; Sa = sum(a.values())
pt = pq.read_table(f"{D}/geom/{key}_pairs.parquet").to_pandas(); listed = set(zip(pt.i, pt.j))
pa = {(i, j): rx.area_inter(cells[i], cells[j]) for i, j in listed}; Sov = sum(pa.values())
tp = pq.read_table(f"{D}/truth/{key}_o1_pairs.parquet").to_pandas(); tpd = {(i, j): v for i, j, v in zip(tp.i, tp.j, tp.a)}
pair_eq = sum(float(v) == tpd[k] for k, v in pa.items())
# completezza: tutte le coppie con riquadri (allargati di 1e-9) sovrapposti
ids = sorted(cells); bb = np.array([rx.fbox(cells[c]) for c in ids])
boxes = shapely.box(bb[:, 0] - 1e-9, bb[:, 2] - 1e-9, bb[:, 1] + 1e-9, bb[:, 3] + 1e-9); tree = shapely.STRtree(boxes)
qi, qj = tree.query(boxes, predicate="intersects"); cand = {(min(ids[x], ids[y]), max(ids[x], ids[y])) for x, y in zip(qi, qj) if x != y}
extra = sorted(cand - listed); ea = {k: rx.area_inter(cells[k[0]], cells[k[1]]) for k in extra}
extra_pos = {k: v for k, v in ea.items() if v > 0}
U_lo = Sa - Sov - sum(extra_pos.values()); U_hi = Sa
tr = json.load(open(f"{D}/truth/{key}_o1.json"))["o1"]
V = pd.read_parquet(f"{D}/review_exact/viol_all.parquet"); V = V[(V.key == key) & (V.vname == "o1_V0")]
vu = V[V.what == "union"]; v0_union = float(vu.cand_val.iloc[0]) if len(vu) else None
out = dict(key=key, n=n, n_cells_read=len(cells), n_nonsimple_exact=len(nonsimple), nonsimple=nonsimple[:20], n_multi_ring=len(multi),
           sum_a=float(Sa), sum_ov_listed=float(Sov), n_pairs_listed=len(listed), pairs_eq_truth=pair_eq, n_bbox_pairs_not_listed=len(extra),
           n_not_listed_positive=len(extra_pos), max_not_listed=float(max(ea.values())) if ea else 0.0, U_lo=float(U_lo), U_hi=float(U_hi), U_width=float(U_hi - U_lo),
           truth_union=tr["union"], truth_union_within_bounds=(F(tr["union"]) >= U_lo - F(1, 10**9) and F(tr["union"]) <= U_hi + F(1, 10**9)),
           truth_minus_Ulo=float(F(tr["union"]) - U_lo), truth_is_nearest_double_of_Ulo=(float(U_lo) == tr["union"]), v0_union=v0_union)
if v0_union is not None:
    d = U_lo - F(v0_union); out["deficit"] = float(d); tauU = float(vu.tau.iloc[0]); out["tauU"] = tauU
    best = sorted(a.items(), key=lambda kv: abs(kv[1] - d))[:3]
    out["best_cells"] = [(c, float(v), float(abs(v - d))) for c, v in best]
    out["deficit_is_one_cell"] = bool(abs(best[0][1] - d) <= tauU)
    # coppie di celle (fra le 200 più vicine in area) che spiegano il deficit, se nessuna singola lo fa
    if not out["deficit_is_one_cell"]:
        near = sorted(a.items(), key=lambda kv: abs(kv[1] - d / 2))[:300]
        pr = min(((c1, c2, abs(v1 + v2 - d)) for (c1, v1), (c2, v2) in itertools.combinations(near, 2)), key=lambda z: z[2])
        out["best_pair_cells"] = (pr[0], pr[1], float(pr[2]))
out["pair_viol"] = [dict(what=r.what, i=int(r.i), j=int(r.j), cand=float(r.cand_val), truth=float(r.truth), a_i=float(a[int(r.i)]), a_j=float(a[int(r.j)]))
                    for _, r in V[V.what.isin(["pair_ij", "pair_ji"])].iterrows()]
out["t_s"] = time.time() - t0
json.dump(out, open(f"{O}/{key}.json", "w"), indent=1); print(json.dumps({k: v for k, v in out.items() if k not in ("pair_viol",)}))
