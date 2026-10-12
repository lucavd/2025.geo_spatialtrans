
import sys, csv, json, math, time; sys.path.insert(0, '.')
from fractions import Fraction as F
from exactgeom import *
import shapely
from shapely.strtree import STRtree
exec(open('p1_exact.py').read().split("geo, T, nbad")[0])
E = "/mnt/micron/geo_spatialtrans/S1.3b/review_exact"
def load(fn):
    from collections import defaultdict
    G = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for r in csv.DictReader(open(fn)):
        G[r['role']][int(r['part'])][int(r['ring'])][int(r['k'])] = (float.fromhex(r['x']), float.fromhex(r['y']))
    return {k: [[[p[rg][kk] for kk in sorted(p[rg])] for rg in sorted(p)] for p in [parts[q] for q in sorted(parts)]] for k, parts in G.items()}
g = load(f"{E}/A5_r5_RSA_15_coords_hex.csv"); meta = list(csv.DictReader(open(f"{E}/A5_r5_RSA_15_sf_meta.csv")))[0]
roles = sorted(k for k in g if k != 'region'); n = len(roles)
R = tofrac(g['region']); C = [tofrac(g[k]) for k in roles]
AR = area(R); Ac = [area(c) for c in C]; S = sum(Ac)
print("n celle", n, "vertici regione", sum(len(r) for p in g['region'] for r in p), "parti/anelli", len(g['region']), [len(p) for p in g['region']])
print("area regione esatta %.17g  sf %s  w %s" % (float(AR), meta['area_reg_sf'], meta['area_reg_w']))
print("somma aree esatta %.17g  sf %s" % (float(S), meta['sum_area']))
print("somma - regione (esatto) = %.6g  rel %.3g" % (float(S - AR), float((S - AR) / AR)))
# coppie vicine
Sh = [toshp(g[k]) for k in roles]; tree = STRtree(Sh)
pairs = sorted({(min(a, b), max(a, b)) for a, b in zip(*tree.query(Sh, predicate='intersects')) if a != b})
t0 = time.time(); ov = []
for a, b in pairs:
    x = inter_area(C[a], C[b])
    if x != 0: ov.append((a, b, x))
print("coppie che si toccano", len(pairs), "con area esatta > 0:", len(ov), "somma", float(sum(x for _, _, x in ov)), "max", float(max([x for *_, x in ov], default=0)), "t", round(time.time() - t0, 1))
# contenimento nella regione: area(c ∩ R) = area(c) per le celle vicine al bordo
Rs = toshp(g['region']); bnd = Rs.boundary; near = [k for k in range(n) if Sh[k].distance(bnd) < 1e-3]
t0 = time.time(); out_reg = []
for k in near:
    x = inter_area(C[k], R)
    if x != Ac[k]: out_reg.append((k, float(Ac[k] - x)))
print("celle vicine al bordo", len(near), "con area fuori regione != 0:", len(out_reg), "somma fuori", sum(v for _, v in out_reg), "t", round(time.time() - t0, 1))
gap = AR - (S - sum(x for *_, x in ov) - F(0))  # trascura sovrapposizioni triple (verificate sotto)
print("lacune esatte (regione - somma + sovrapposizioni, al netto di fuori regione) = %.6g" % (float(gap) + sum(v for _, v in out_reg)))
# shapely: unione e celle mancanti
U = shapely.unary_union(Sh); print("shapely unary_union area %.17g deficit %.6g valid %s" % (U.area, sum(s.area for s in Sh) - U.area, U.is_valid))
Ucov = [shapely.intersection(s, U).area for s in Sh]; pos = [shapely.point_on_surface(s) for s in Sh]
miss = [k for k in range(n) if abs(Ucov[k] - Sh[k].area) > 1e-6 or not U.contains(pos[k])]
print("shapely: celle non coperte dall'unione", [(roles[k], Sh[k].area, Ucov[k]) for k in miss])
json.dump(dict(n=n, area_region_exact=str(AR), sum_area_exact=str(S), sum_minus_region=float(S - AR), n_pairs=len(pairs),
    overlaps=[(roles[a], roles[b], str(x), float(x)) for a, b, x in ov], out_region=[(roles[k], v) for k, v in out_reg],
    shapely_union_area=U.area, shapely_missing=[roles[k] for k in miss]), open(f"{E}/A5_exact_results.json", "w"), indent=1)
