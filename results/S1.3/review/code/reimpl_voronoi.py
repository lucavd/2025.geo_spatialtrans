#!/usr/bin/env python3
"""results/S1.3/review/code/reimpl_voronoi.py — revisione del codice S1.3 (RA-code).
Reimplementazione INDIPENDENTE di tessellate_voronoi() (cs = 0):
  * tile: scipy.spatial.Voronoi (Qhull), per regione, con 8 generatori fittizi lontani (> diag. della regione)
    -> tutte le celle reali limitate e identiche alla cella vera dentro il riquadro della regione;
  * area(tile ∩ regione) senza GEOS: Sutherland-Hodgman di ogni anello della regione contro i semipiani
    delle bisettrici (coordinate locali al generatore); l'area con segno e' esatta anche per anelli non convessi;
  * pezzi del ritaglio: shapely (GEOS) intersection; regola D-S1.3.2 applicata da noi con lunghezze di confine
    stimate per campionamento (passo H) dei bordi dell'orfano, spostati di EPS verso l'esterno.
Uso: review_venv/bin/python reimpl_voronoi.py <lab> [<lab> ...]   (ingressi in review_code_data, da export_inputs.R)
"""
import sys, json, time
import numpy as np, pandas as pd, shapely
from shapely.geometry import Polygon, MultiPolygon, Point, box
from shapely import STRtree
from scipy.spatial import Voronoi

D = "/mnt/micron/geo_spatialtrans/S1.3/review_code_data"
H = 0.002      # passo di campionamento dei bordi (µm)
EPS = 1e-5     # spostamento verso l'esterno (µm)

def shoelace(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))

def sh_clip(P, hps):
    """Sutherland-Hodgman vettoriale: P (m,2) anello aperto; hps = [(nx, ny, c)] con interno nx*x+ny*y <= c."""
    for nx, ny, cc in hps:
        if len(P) < 3: return None
        f = P[:, 0] * nx + P[:, 1] * ny - cc
        if np.all(f <= 0): continue
        if np.all(f > 0): return None
        Q = np.roll(P, -1, axis=0); fq = np.roll(f, -1)
        keep = f <= 0
        cross = ((f < 0) & (fq > 0)) | ((f > 0) & (fq < 0))
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(cross, f / (f - fq), 0.0)
        I = P + (Q - P) * t[:, None]
        out = np.empty((2 * len(P), 2)); m = np.zeros(2 * len(P), bool)
        out[0::2] = P; out[1::2] = I; m[0::2] = keep; m[1::2] = cross
        P = out[m]
    return P if len(P) >= 3 else None

def region_rings(g):
    """lista di (anello aperto (m,2), segno): +1 esterno, -1 buco; con il bbox di ciascun anello."""
    polys = list(g.geoms) if isinstance(g, MultiPolygon) else [g]
    rings = []
    for p in polys:
        for k, r in enumerate([p.exterior] + list(p.interiors)):
            A = np.asarray(r.coords)[:-1]
            s = np.sign(shoelace(A))
            rings.append((A, s * (1 if k == 0 else -1), A.min(0), A.max(0)))
    return rings

def voronoi_tiles(x, y, bb):
    """tile di Voronoi (Qhull) + semipiani delle bisettrici per ogni generatore reale."""
    n = len(x); xmin, ymin, xmax, ymax = bb
    L = max(xmax - xmin, ymax - ymin); cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2; R = 3 * L + 10
    ang = np.arange(8) * np.pi / 4
    dum = np.c_[cx + R * np.cos(ang), cy + R * np.sin(ang)]
    pts = np.vstack([np.c_[x, y], dum])
    vor = Voronoi(pts)
    nb = [[] for _ in range(n)]
    for a, b in vor.ridge_points:
        if a < n: nb[a].append(b)
        if b < n: nb[b].append(a)
    tiles, hps = [], []
    for i in range(n):
        reg = vor.regions[vor.point_region[i]]
        assert -1 not in reg and len(reg) >= 3, "cella illimitata"
        V = vor.vertices[reg]; cvx = V.mean(0)
        V = V[np.argsort(np.arctan2(V[:, 1] - cvx[1], V[:, 0] - cvx[0]))]
        tiles.append(Polygon(V))
        q = pts[nb[i]] - pts[i]
        hps.append([(qq[0], qq[1], 0.5 * (qq[0] ** 2 + qq[1] ** 2)) for qq in q])
    return tiles, hps

def sh_area(i, x, y, hp, rings, tb):
    """area(tile_i ∩ regione) senza GEOS (coordinate locali al generatore)."""
    tot = 0.0
    for A, s, lo, hi in rings:
        if hi[0] < tb[0] or lo[0] > tb[2] or hi[1] < tb[1] or lo[1] > tb[3]: continue
        P = A - np.array([x, y])
        loc = [(1, 0, tb[2] - x + 1e-9), (-1, 0, -(tb[0] - x) + 1e-9), (0, 1, tb[3] - y + 1e-9), (0, -1, -(tb[1] - y) + 1e-9)]
        P = sh_clip(P, loc)
        if P is None: continue
        P = sh_clip(P, hp)
        if P is None: continue
        tot += s * shoelace(P)
    return tot

def boundary_samples(poly):
    """punti medi di segmenti di lunghezza <= H sul bordo di poly, con la lunghezza di ciascun segmento e la normale."""
    out = []
    for r in [poly.exterior] + list(poly.interiors):
        C = np.asarray(r.coords)
        for a, b in zip(C[:-1], C[1:]):
            d = b - a; l = np.hypot(*d)
            if l == 0: continue
            k = max(1, int(np.ceil(l / H)))
            t = (np.arange(k) + 0.5) / k
            M = a + t[:, None] * d
            nrm = np.array([-d[1], d[0]]) / l
            out.append(np.c_[M, np.full(k, l / k), np.tile(nrm, (k, 1))])
    return np.vstack(out)

def shared_lengths(orph, terrs, tree):
    """lunghezza del bordo dell'orfano adiacente a ciascun territorio (campionamento, passo H)."""
    S = boundary_samples(orph)
    p1 = S[:, :2] + EPS * S[:, 3:5]; p2 = S[:, :2] - EPS * S[:, 3:5]
    in1 = shapely.contains_xy(orph, p1[:, 0], p1[:, 1])
    P = np.where(in1[:, None], p2, p1)                     # il lato fuori dall'orfano
    pts = shapely.points(P)
    qi, ti = tree.query(pts, predicate="within")
    lens = {}
    for a, b in zip(qi, ti):
        lens[b] = lens.get(b, 0.0) + S[a, 2]
    return lens

def pieces(g):
    if g.is_empty: return []
    if isinstance(g, Polygon): return [g]
    if isinstance(g, MultiPolygon): return list(g.geoms)
    return [q for q in getattr(g, "geoms", []) if isinstance(q, Polygon) and q.area > 0] + \
           [z for q in getattr(g, "geoms", []) if isinstance(q, MultiPolygon) for z in q.geoms]

def run(lab):
    t0 = time.time()
    cen = pd.read_csv(f"{D}/{lab}_cen.csv"); reg = pd.read_csv(f"{D}/{lab}_reg.csv"); tv = pd.read_csv(f"{D}/{lab}_tv.csv")
    assert (cen.cell_id.values == np.sort(cen.cell_id.values)).all() and (tv.cell_id.values == cen.cell_id.values).all()
    regg = dict(zip(reg.region_id, shapely.from_wkb(reg.wkb.values)))
    Rterr = shapely.from_wkb(tv.terr_wkb.values); Rtile = shapely.from_wkb(tv.tile_wkb.values)
    n = len(cen)
    my_area = np.full(n, np.nan); sh_tot = np.full(n, np.nan); npieces = np.zeros(n, int)
    tile_symd = np.full(n, np.nan); terr_symd = np.full(n, np.nan); my_terr = [None] * n
    orph_rows = []; region_rows = []
    for rid, g in regg.items():
        k = np.where(cen.region_id.values == rid)[0]
        if len(k) == 0: continue
        x = cen.x.values[k]; y = cen.y.values[k]
        rw = reg.loc[reg.region_id == rid, ["xmin", "xmax", "ymin", "ymax"]].values[0]
        rings = region_rings(g); bb = g.bounds
        if len(k) == 1:
            tiles = [box(bb[0] - 1e3, bb[1] - 1e3, bb[2] + 1e3, bb[3] + 1e3)]; hps = [[]]
        else:
            tiles, hps = voronoi_tiles(x, y, bb)
        rect = box(rw[0], rw[2], rw[1], rw[3])
        main = []; orphans = []
        for j, i in enumerate(k):
            tb = tiles[j].bounds
            sh_tot[i] = sh_area(i, x[j], y[j], hps[j], rings, tb)
            tile_symd[i] = tiles[j].intersection(rect).symmetric_difference(Rtile[i]).area
            pcs = pieces(tiles[j].intersection(g))
            npieces[i] = len(pcs)
            gp = Point(x[j], y[j])
            own = [q for q, pc in enumerate(pcs) if pc.intersects(gp)]
            own = own[0] if own else int(np.argmax([pc.area for pc in pcs]))
            main.append(pcs[own])
            for q, pc in enumerate(pcs):
                if q != own: orphans.append((i, pc))
        # D-S1.3.2 (nostra applicazione): passata 1 contro i soli pezzi principali; orfani senza contatto -> passate successive
        cur = {i: [main[j]] for j, i in enumerate(k)}
        pending = list(range(len(orphans))); assign = {}
        for ps in range(10):
            terrs_i = list(k); geoms = [shapely.union_all(cur[i]) if len(cur[i]) > 1 else cur[i][0] for i in terrs_i]
            tree = STRtree(geoms); newly = {}
            for o in pending:
                don, pc = orphans[o]
                L = shared_lengths(pc, geoms, tree)
                if not L: continue
                mx = max(L.values())
                best = sorted([terrs_i[b] for b, v in L.items() if v >= mx - 1e-12])[0]
                second = sorted(L.values())[-2] if len(L) > 1 else 0.0
                newly[o] = (best, mx, second, terrs_i[max(L, key=L.get)])
            for o, (best, mx, second, _) in newly.items():
                cur[best].append(orphans[o][1]); assign[o] = (best, mx, second, ps + 1)
            pending = [o for o in pending if o not in newly]
            if not pending or not newly: break
        for o in pending:                         # nessun confinante: resta al generatore (deviazione 3)
            don = orphans[o][0]; cur[don].append(orphans[o][1]); assign[o] = (don, 0.0, 0.0, -1)
        for j, i in enumerate(k):
            my_terr[i] = shapely.union_all(cur[i]) if len(cur[i]) > 1 else cur[i][0]
            my_area[i] = sum(q.area for q in cur[i])
            terr_symd[i] = my_terr[i].symmetric_difference(Rterr[i]).area
        # destinatario secondo l'output R: territorio R che contiene un punto interno dell'orfano
        Rtree = STRtree(Rterr[k])
        for o, (don, pc) in enumerate(orphans):
            rp = pc.point_on_surface()
            hit = Rtree.query(rp, predicate="within")
            r_rec = int(cen.cell_id.values[k[hit[0]]]) if len(hit) else -1
            best, mx, second, ps = assign[o]
            orph_rows.append(dict(case=lab, region_id=rid, donor=int(cen.cell_id.values[don]), orphan_area=pc.area,
                                  my_recipient=int(cen.cell_id.values[best]), R_recipient=r_rec, shared_len_best=mx,
                                  shared_len_second=second, pass_=ps))
        A = g.area
        region_rows.append(dict(case=lab, region_id=rid, n_cells=len(k), area_shapely=A,
                                area_sf=float(reg.loc[reg.region_id == rid, "area_sf"].values[0]),
                                sum_sh_tiles=float(np.nansum(sh_tot[k])), sum_my=float(np.nansum(my_area[k])),
                                sum_R=float(tv.territory_area.values[k].sum()), n_orphans=sum(1 for o in orphans)))
    Ra = tv.territory_area.values
    mono = (tv.n_pieces_lost.values == 0) & (tv.n_pieces_gained.values == 0)
    cell = pd.DataFrame(dict(case=lab, cell_id=cen.cell_id.values, region_id=cen.region_id.values, R_area=Ra,
                             my_area=my_area, sh_area=sh_tot, n_pieces_my=npieces, R_lost=tv.n_pieces_lost.values,
                             R_gained=tv.n_pieces_gained.values, R_gtype=tv.gtype.values, tile_symdiff=tile_symd, terr_symdiff=terr_symd,
                             R_valid=shapely.is_valid(Rterr),
                             gen_in_R=shapely.intersects_xy(Rterr, cen.x.values, cen.y.values)))
    cell.to_csv(f"{D}/review_{lab}_cells.csv", index=False)
    pd.DataFrame(orph_rows).to_csv(f"{D}/review_{lab}_orphans.csv", index=False)
    pd.DataFrame(region_rows).to_csv(f"{D}/review_{lab}_regions.csv", index=False)
    rel = lambda a, b: float(np.max(np.abs(a - b) / b)) if len(a) else float("nan")
    inv = ~mono
    s = dict(case=lab, n_cells=n, n_mono_R=int(mono.sum()), n_frag_cells_R=int(inv.sum()),
             n_cells_multipiece_my=int((npieces > 1).sum()), n_orphans=len(orph_rows),
             mono_maxrel_R_vs_SH=rel(Ra[mono], sh_tot[mono]), mono_maxrel_R_vs_my=rel(Ra[mono], my_area[mono]),
             mono_n_pieces_my_gt1=int((npieces[mono] > 1).sum()),
             frag_maxrel_R_vs_my=rel(Ra[inv], my_area[inv]),
             n_recipient_disagree=int(sum(r["my_recipient"] != r["R_recipient"] for r in orph_rows)),
             n_orphans_near_tie=int(sum(r["shared_len_best"] - r["shared_len_second"] < 5 * H for r in orph_rows)),
             max_tile_symdiff=float(np.nanmax(tile_symd)), max_terr_symdiff=float(np.nanmax(terr_symd)),
             sum_terr_symdiff=float(np.nansum(terr_symd)),
             max_rel_region_sumSH=max(abs(r["sum_sh_tiles"] - r["area_shapely"]) / r["area_shapely"] for r in region_rows),
             max_rel_region_sumR=max(abs(r["sum_R"] - r["area_shapely"]) / r["area_shapely"] for r in region_rows),
             n_R_invalid=int((~cell.R_valid).sum()), n_R_not_polygon=int((cell.R_gtype != "POLYGON").sum()),
             n_gen_outside_R=int((~cell.gen_in_R).sum()), elapsed_s=round(time.time() - t0, 1))
    pd.DataFrame([s]).to_csv(f"{D}/review_{lab}_summary.csv", index=False)
    print(json.dumps(s))

if __name__ == "__main__":
    from multiprocessing import Pool
    labs = sys.argv[1:]
    with Pool(min(len(labs), 12)) as p:
        p.map(run, labs)
