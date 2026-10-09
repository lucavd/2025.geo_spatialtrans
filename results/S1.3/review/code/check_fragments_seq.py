#!/usr/bin/env python3
"""results/S1.3/review/code/check_fragments_seq.py — revisione S1.3 (RA-code), compito 2(ii) seconda parte.
Riapplica D-S1.3.2 in modo SEQUENZIALE come il codice R (orfani in ordine di cell_id del donatore; ogni unione
aggiorna subito il territorio ricevente; fino a 10 passate), con la NOSTRA misura delle lunghezze condivise
(campionamento dei bordi, reimpl_voronoi.shared_lengths). Se i disaccordi con R spariscono, la differenza con la
versione "simultanea" di reimpl_voronoi.py e' dovuta all'ordine di elaborazione (regola ambigua sulle catene).
Inoltre: territori R non POLYGON (pezzi, contatto fra i pezzi) e celle di arbiter_systematic.csv discordanti (Qhull).
"""
import sys, numpy as np, pandas as pd, shapely
from shapely.geometry import Point, MultiPolygon, Polygon
from shapely import STRtree
sys.path.insert(0, "/home/user/2025.geo_spatialtrans/results/S1.3/review/code")
import reimpl_voronoi as rv
D = rv.D; OUT = "/home/user/2025.geo_spatialtrans/results/S1.3/review/code"

def run(lab):
    cen = pd.read_csv(f"{D}/{lab}_cen.csv"); reg = pd.read_csv(f"{D}/{lab}_reg.csv"); tv = pd.read_csv(f"{D}/{lab}_tv.csv")
    regg = dict(zip(reg.region_id, shapely.from_wkb(reg.wkb.values))); Rterr = shapely.from_wkb(tv.terr_wkb.values)
    rows = []; mp_rows = []; area_rows = []
    for rid, g in regg.items():
        k = np.where(cen.region_id.values == rid)[0]
        if len(k) < 2: continue
        x = cen.x.values[k]; y = cen.y.values[k]
        tiles, hps = rv.voronoi_tiles(x, y, g.bounds)
        geoms = []; orphans = []
        for j in range(len(k)):
            pcs = rv.pieces(tiles[j].intersection(g)); gp = Point(x[j], y[j])
            own = [q for q, pc in enumerate(pcs) if pc.intersects(gp)]; own = own[0] if own else int(np.argmax([p.area for p in pcs]))
            geoms.append(pcs[own]); orphans += [(j, pc) for q, pc in enumerate(pcs) if q != own]
        pending = list(range(len(orphans))); rec = {}
        for ps in range(10):
            prog = False
            for o in list(pending):
                don, pc = orphans[o]
                tree = STRtree(geoms)
                L = rv.shared_lengths(pc, geoms, tree)
                if not L or max(L.values()) <= 0: continue
                mx = max(L.values()); best = min(b for b, v in L.items() if v >= mx - 1e-12)
                sec = sorted(L.values())[-2] if len(L) > 1 else 0.0
                geoms[best] = shapely.union_all([geoms[best], pc]); rec[o] = (best, mx, sec, ps + 1); pending.remove(o); prog = True
            if not pending or not prog: break
        for o in pending:
            don, pc = orphans[o]; geoms[don] = shapely.union_all([geoms[don], pc]); rec[o] = (don, 0.0, 0.0, -1)
        Rtree = STRtree(Rterr[k])
        for o, (don, pc) in enumerate(orphans):
            hit = Rtree.query(pc.point_on_surface(), predicate="within"); best, mx, sec, ps = rec[o]
            rows.append(dict(case=lab, donor=int(cen.cell_id.values[k[don]]), orphan_area=pc.area, seq_recipient=int(cen.cell_id.values[k[best]]),
                             R_recipient=int(cen.cell_id.values[k[hit[0]]]) if len(hit) else -1, shared_len_best=mx, shared_len_second=sec, pass_=ps))
        for j, i in enumerate(k):
            my = geoms[j]; a = my.area
            rows[-1:]  # noqa
            if not isinstance(Rterr[i], Polygon):
                pcs = list(Rterr[i].geoms) if hasattr(Rterr[i], "geoms") else [Rterr[i]]
                inter = [pcs[u].intersection(pcs[v]).geom_type for u in range(len(pcs)) for v in range(u + 1, len(pcs))]
                dist = min(pcs[u].distance(pcs[v]) for u in range(len(pcs)) for v in range(u + 1, len(pcs)))
                mp_rows.append(dict(case=lab, cell_id=int(cen.cell_id.values[i]), R_type=Rterr[i].geom_type, n_parts=len(pcs),
                                    part_areas=";".join(f"{p.area:.6g}" for p in pcs), parts_contact=";".join(inter), min_dist_parts=dist,
                                    R_lost=int(tv.n_pieces_lost.values[i]), R_gained=int(tv.n_pieces_gained.values[i]),
                                    seq_area=a, R_area=Rterr[i].area, symdiff_seq_R=my.symmetric_difference(Rterr[i]).area))
        for j, i in enumerate(k):
            area_rows.append(dict(case=lab, cell_id=int(cen.cell_id.values[i]), seq_area=geoms[j].area, R_area=float(tv.territory_area.values[i]),
                                  involved=bool(tv.n_pieces_lost.values[i] > 0 or tv.n_pieces_gained.values[i] > 0)))
    o = pd.DataFrame(rows); m = pd.DataFrame(mp_rows); ar_ = pd.DataFrame(area_rows)
    ar_['rel'] = (ar_.seq_area - ar_.R_area).abs() / ar_.R_area
    ar_.to_csv(f'{D}/review_{lab}_seq_areas.csv', index=False)
    iv = ar_[ar_.involved]
    print(lab, 'celle coinvolte (R)', len(iv), 'max rel |seq - R| coinvolte:', iv.rel.max() if len(iv) else float('nan'), ' tutte:', ar_.rel.max(), ' n rel>1e-9:', int((ar_.rel > 1e-9).sum()))
    # aree per cella: sequenziale vs R
    return o, m

if __name__ == "__main__":
    allo, allm = [], []
    for lab in sys.argv[1:]:
        o, m = run(lab); allo.append(o); allm.append(m)
        n_dis = int((o.seq_recipient != o.R_recipient).sum()) if len(o) else 0
        print(lab, "orfani", len(o), "disaccordi sequenziale vs R:", n_dis, " territori R non POLYGON:", len(m))
    O = pd.concat(allo); M = pd.concat(allm)
    O.to_csv(f"{OUT}/fragments_sequential_vs_R.csv", index=False); M.to_csv(f"{OUT}/R_multipart_territories.csv", index=False)
    print(O[O.seq_recipient != O.R_recipient].to_string()); print(M.to_string())
    # celle di arbiter_systematic.csv con hp != GEOS (rel > 1e-9): cella vera con Qhull sull'insieme completo
    import pyarrow.parquet as pq
    from scipy.spatial import Voronoi
    def qcell(x, y, i):
        L = max(np.ptp(x), np.ptp(y)); cx, cy = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2; R = 3 * L + 10
        ang = np.arange(8) * np.pi / 4
        pts = np.vstack([np.c_[x, y], np.c_[cx + R * np.cos(ang), cy + R * np.sin(ang)]]); vor = Voronoi(pts)
        V = vor.vertices[vor.regions[vor.point_region[i]]]; m = V.mean(0); V = V[np.argsort(np.arctan2(V[:, 1] - m[1], V[:, 0] - m[0]))]
        return 0.5 * abs(np.dot(V[:, 0], np.roll(V[:, 1], -1)) - np.dot(np.roll(V[:, 0], -1), V[:, 1])), len(V)
    st = pd.read_csv(f"{D}/arbiter_systematic.csv"); st = st[st.rel > 1e-9]; rows = []
    for _, r in st.iterrows():
        A, roi, meth = r.case.split("_", 2)
        t = pq.read_table(f"/mnt/micron/geo_spatialtrans/R3/real/{A}_{roi}_{meth}_gen.parquet").to_pandas()
        a, nv = qcell(t.x.values, t.y.values, int(r.cell) - 1)
        rows.append(dict(case=r.case, cell=int(r.cell), a_geos=r.a_geos, a_hp=r.a_hp, a_qhull=a, nv_qhull=nv, ns_geos=r.ns_geos, ns_hp=r.ns_hp,
                         K=r.K, rel_geos_qhull=abs(r.a_geos - a) / a, rel_hp_qhull=abs(r.a_hp - a) / a))
    # reticolo esagonale con buco (stessa costruzione di arbiter_check.R), celle 677, 1120, 1121
    a0 = 10; I, J = np.meshgrid(np.arange(41), np.arange(47), indexing="ij"); I = I.ravel(order="F"); J = J.ravel(order="F")
    x = I * a0 + (J % 2) * a0 / 2; y = J * a0 * np.sqrt(3) / 2; keep = (x - 200) ** 2 + (y - 200) ** 2 > 60 ** 2; x = x[keep]; y = y[keep]
    k = pd.read_csv(f"{OUT}/arbiter_known_cases.csv"); kb = k[k.case.str.contains("buco")]
    for _, r in kb[kb.rel_err > 1e-9].iterrows():
        a, nv = qcell(x, y, int(r.cell) - 1)
        rows.append(dict(case="buco r=60 (esagonale)", cell=int(r.cell), a_geos=r.a_ref, a_hp=r.a_hp, a_qhull=a, nv_qhull=nv, ns_geos=r.ns_ref, ns_hp=r.ns_hp,
                         K=r.K, rel_geos_qhull=abs(r.a_ref - a) / a, rel_hp_qhull=abs(r.a_hp - a) / a))
    Q = pd.DataFrame(rows); Q.to_csv(f"{OUT}/arbiter_failures_qhull.csv", index=False); print(Q.to_string())
