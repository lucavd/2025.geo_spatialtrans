
import sys, csv, math; sys.path.insert(0,'.')
exec(open('a5_exact.py').read().split("g = load(")[0])
g = load(f"{E}/A5_r5_RSA_15_coords_hex.csv"); R = tofrac(g['region'])
a, b = 'c00071', 'c00175'
Pa, Pb = tofrac(g[a]), tofrac(g[b]); Sa, Sb = toshp(g[a]), toshp(g[b])
print("tipi/vertici", len(g[a]), len(g[b]), [len(r) for p in g[a] for r in p], [len(r) for p in g[b] for r in p])
print("area esatta 71 %.12f 175 %.12f" % (float(area(Pa)), float(area(Pb))))
A, pcs = inter_area(Pa, Pb, detail=True); print("intersezione esatta 71∩175 =", A, float(A), "pezzi", len(pcs), "SH", (sh_clip_area(Pa, Pb) if is_convex_ring(Pb) else None), is_convex_ring(Pa), is_convex_ring(Pb))
Va = {v for p in Pa for r in p for v in strip_close(r)}; Vb = {v for p in Pb for r in p for v in strip_close(r)}
print("vertici identici", len(Va & Vb))
for v in Va - Vb:
    w = min(Vb, key=lambda q: (q[0]-v[0])**2 + (q[1]-v[1])**2); d2 = (w[0]-v[0])**2 + (w[1]-v[1])**2
    if d2 < F(1, 10**12): print(" diverso", (float(v[0]), float(v[1])), (float(w[0]), float(w[1])), "ulp", ulp_diff(v[0], w[0]), ulp_diff(v[1], w[1]), "d", math.sqrt(float(d2)),
        "v_in_175", pip(v, edges(Pb)), "w_in_71", pip(w, edges(Pa)), "dreg", math.sqrt(float(bnd_dist2(v, R))), math.sqrt(float(bnd_dist2(w, R))))
# T-junction: vertici di una cella sul bordo dell'altra (non vertici)
for name, V, P in (("71 su bordo 175", Va - Vb, Pb), ("175 su bordo 71", Vb - Va, Pa)):
    print(name, sum(pip(v, edges(P)) == 'on' for v in V))
print("shapely relate", shapely.relate(Sa, Sb), "inter(71,175)", shapely.intersection(Sa, Sb).area, "inter(175,71)", shapely.intersection(Sb, Sa).area, "union deficit", Sa.area + Sb.area - shapely.union(Sa, Sb).area)
print("shapely grid1e-6 inter", shapely.intersection(Sa, Sb, grid_size=1e-6).area)
print("71 coords:", [(x, y) for x, y in g[a][0][0]]); print("175 coords:", [(x, y) for x, y in g[b][0][0]])
# P1: shapely in entrambi gli ordini
exec(open('p1_exact.py').read().split("geo, T, nbad")[0].replace('def load(fn)', 'def load1(fn)'))
geo, T, nb = load1(f"{E}/P1_pairs_coords_hex.csv")
for case in sorted({k[0] for k in geo}):
    Si, Sj = toshp(geo[(case,'i')]), toshp(geo[(case,'j')])
    print(case, "shp inter(i,j) %.6g inter(j,i) %.6g diff(i,j) %.6g union_def %.6g" % (shapely.intersection(Si,Sj).area, shapely.intersection(Sj,Si).area, shapely.difference(Si,Sj).area, Si.area+Sj.area-shapely.union(Si,Sj).area))
