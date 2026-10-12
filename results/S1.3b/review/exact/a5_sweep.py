
import sys, csv, math; sys.path.insert(0,'.')
exec(open('a5_exact.py').read().split("g = load(")[0])
from shapely.strtree import STRtree
g = load(f"{E}/A5_r5_RSA_15_coords_hex.csv"); roles = sorted(k for k in g if k != 'region')
Sh = [toshp(g[k]) for k in roles]; tree = STRtree(Sh)
pairs = sorted({(min(a, b), max(a, b)) for a, b in zip(*tree.query(Sh, predicate='intersects')) if a != b})
nd = 0; bad = []
for a, b in pairs:
    Va = {v for p in g[roles[a]] for r in p for v in r}; Vb = {v for p in g[roles[b]] for r in p for v in r}
    diff = [v for v in Va - Vb if min(math.hypot(v[0]-w[0], v[1]-w[1]) for w in Vb) < 1e-6]
    nd += bool(diff)
    rel = shapely.relate(Sh[a], Sh[b]); ab = shapely.intersection(Sh[a], Sh[b]).area; ba = shapely.intersection(Sh[b], Sh[a]).area
    tau = 1e-6 * (Sh[a].length + Sh[b].length)
    if ab > tau or ba > tau: bad.append((roles[a], roles[b], rel, ab, ba, len(diff)))
print("coppie", len(pairs), "con vertice condiviso diverso (<1e-6, !=0):", nd)
print("coppie con intersezione shapely > 1e-6*(Pi+Pj):", len(bad)); [print(" ", x) for x in bad]
