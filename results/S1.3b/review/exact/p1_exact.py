
import sys, csv, json, math; sys.path.insert(0, '.')
from collections import defaultdict
from fractions import Fraction as F
from exactgeom import *
import shapely
from shapely.geometry import Polygon, MultiPolygon, Point
E = "/mnt/micron/geo_spatialtrans/S1.3b/review_exact"
def load(fn):
    G = defaultdict(lambda: defaultdict(lambda: defaultdict(dict))); T = {}; nbad = 0
    for r in csv.DictReader(open(fn)):
        x, y = float.fromhex(r['x']), float.fromhex(r['y'])
        if x != float(r['xd']) or y != float(r['yd']): nbad += 1
        G[(r['case'], r['role'])][int(r['part'])][int(r['ring'])][int(r['k'])] = (x, y); T[(r['case'], r['role'])] = r['type']
    out = {}
    for key, parts in G.items():
        out[key] = [[[parts[p][rg][k] for k in sorted(parts[p][rg])] for rg in sorted(parts[p])] for p in sorted(parts)]
    return out, T, nbad
def tofrac(g): return [[[(F(x), F(y)) for x, y in r] for r in part] for part in g]
def toshp(g):
    ps = [Polygon(part[0], part[1:]) for part in g]
    return ps[0] if len(ps) == 1 else MultiPolygon(ps)
geo, T, nbad = load(f"{E}/P1_pairs_coords_hex.csv"); print("hex != decimale:", nbad)
meta = {r['case']: r for r in csv.DictReader(open(f"{E}/P1_pairs_sf_meta.csv"))}
p1 = list(csv.DictReader(open("../../S1.3b_posthoc_P1_pairs_v0.csv")))
rows = []
for r in p1:
    case = f"{r['archetype']}_{r['roi_id']}_{r['model']}_{int(r['rep']):02d}"; m = meta[case]
    gi, gj, gr = geo[(case, 'i')], geo[(case, 'j')], geo[(case, 'region')]
    Pi, Pj, Rg = tofrac(gi), tofrac(gj), tofrac(gr)
    Ai, Aj = area(Pi), area(Pj)
    Aij, pieces = inter_area(Pi, Pj, detail=True); Aji = inter_area(Pj, Pi)
    cvx_i, cvx_j = is_convex_ring(Pi), is_convex_ring(Pj)
    sh = sh_clip_area(Pi, Pj) if cvx_j else (sh_clip_area(Pj, Pi) if cvx_i else None)
    # vertici condivisi: identici e diversi (< 1e-6)
    Vi = {v for part in Pi for rg in part for v in strip_close(rg)}; Vj = {v for part in Pj for rg in part for v in strip_close(rg)}
    ident = Vi & Vj; diff = []
    for a in Vi - ident:
        b = min(Vj - ident, key=lambda q: (q[0]-a[0])**2 + (q[1]-a[1])**2)
        d2 = (b[0]-a[0])**2 + (b[1]-a[1])**2
        if d2 < F(1, 10**12): diff.append((a, b, d2))
    assert len(diff) >= 1
    Ei, Ej = edges(Pi), edges(Pj)
    dinfo = []
    for a, b, d2 in diff:
        dinfo.append(dict(vi=(float(a[0]).hex(), float(a[1]).hex()), vj=(float(b[0]).hex(), float(b[1]).hex()),
            ulp_x=ulp_diff(a[0], b[0]), ulp_y=ulp_diff(a[1], b[1]), dist=math.sqrt(float(d2)), dist2_exact=str(d2),
            vi_in_j=pip(a, Ej), vj_in_i=pip(b, Ei),
            dreg_vi=math.sqrt(float(bnd_dist2(a, Rg))), dreg_vj=math.sqrt(float(bnd_dist2(b, Rg))),
            vi_on_reg_exact=bnd_dist2(a, Rg) == 0, vj_on_reg_exact=bnd_dist2(b, Rg) == 0))
    # geometria restituita da sf: e' una delle due celle?
    gs = geo.get((case, 'sf_inter')); Ps = tofrac(gs); As = area(Ps)
    Vs = {v for part in Ps for rg in part for v in strip_close(rg)}
    As_in_i, As_in_j = inter_area(Ps, Pi), inter_area(Ps, Pj)
    # shapely (GEOS 3.14.1)
    Si, Sj = toshp(gi), toshp(gj)
    x0 = shapely.intersection(Si, Sj); x6 = shapely.intersection(Si, Sj, grid_size=1e-6)
    Si6, Sj6 = shapely.set_precision(Si, 1e-6), shapely.set_precision(Sj, 1e-6)
    x6b = shapely.intersection(Si6, Sj6)
    u0 = shapely.union(Si, Sj)
    rows.append(dict(case=case, i=int(r['i']), j=int(r['j']), type_i=T[(case,'i')], type_j=T[(case,'j')],
        n_vert_i=len(Vi), n_vert_j=len(Vj), n_shared_ident=len(ident), n_shared_diff=len(diff),
        convex_i=cvx_i, convex_j=cvx_j,
        area_i_exact=float(Ai), area_j_exact=float(Aj), area_i_sf=float(m['area_i_sf']), area_j_sf=float(m['area_j_sf']),
        inter_exact=str(Aij), inter_exact_float=float(Aij), inter_exact_ji_equal=Aij == Aji, inter_SH_equal=(None if sh is None else sh == Aij),
        n_pieces_boundary=len(pieces), interiors_overlap_exact=Aij > 0,
        relate_sf=m['relate'], II_sf_correct=(m['relate'][0] == '2') == (Aij > 0),
        inter_sf=float(m['area_inter_sf']), inter_sf_minus_exact=float(F(float(m['area_inter_sf'])) - Aij),
        sf_inter_area_exact=float(As), sf_inter_vertices_eq_i=Vs == Vi, sf_inter_vertices_eq_j=Vs == Vj,
        sf_inter_vertices_subset_i=Vs <= Vi, sf_inter_minus_i_area=float(As - As_in_i), sf_inter_inside_j_area=float(As_in_j),
        area_v0_csv=float(r['area_v0']), area_i_csv=float(r['area_i']), area_j_csv=float(r['area_j']),
        v0_eq_i_csv=r['area_v0'] == r['area_i'], v0_eq_i_rel=abs(float(r['area_v0'])-float(r['area_i']))/float(r['area_i']),
        shp_area=x0.area, shp_relate=shapely.relate(Si, Sj), shp_overlaps=bool(shapely.overlaps(Si, Sj)), shp_valid=bool(Si.is_valid and Sj.is_valid),
        shp_union_deficit=Si.area + Sj.area - u0.area,
        shp_area_grid1e6=x6.area, shp_area_setprec=x6b.area, shp_relate_setprec=shapely.relate(Si6, Sj6),
        diff=dinfo))
json.dump(rows, open(f"{E}/P1_exact_results.json", "w"), indent=1, default=str)
flat = []
for x in rows:
    d = dict(x); dd = d.pop('diff')
    for k, v in enumerate(dd):
        for kk, vv in v.items(): d[f"diff{k+1}_{kk}"] = vv if not isinstance(vv, tuple) else ";".join(vv)
    flat.append(d)
keys = sorted({k for d in flat for k in d}, key=lambda k: list(flat[0]).index(k) if k in flat[0] else 999)
with open("P1_exact_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); [w.writerow(d) for d in flat]
for x in rows:
    print(x['case'], x['i'], x['j'], "conv", x['convex_i'], x['convex_j'], "id/diff", x['n_shared_ident'], x['n_shared_diff'],
          "| exact", x['inter_exact_float'], x['inter_exact_ji_equal'], x['inter_SH_equal'], "| sf II ok", x['II_sf_correct'],
          "| sf_inter==i", x['sf_inter_vertices_eq_i'], "in_j", x['sf_inter_inside_j_area'],
          "| shp", x['shp_area'], x['shp_relate'], x['shp_area_grid1e6'], x['shp_area_setprec'], x['shp_relate_setprec'])
    for d in x['diff']: print("   diff ulp", d['ulp_x'], d['ulp_y'], "dist", d['dist'], "vi_in_j", d['vi_in_j'], "vj_in_i", d['vj_in_i'], "dreg", d['dreg_vi'], d['dreg_vj'], d['vi_on_reg_exact'], d['vj_on_reg_exact'])
