#!/usr/bin/env python
"""results/S1.4c/review/recalc/r03_obs.py — revisione avversariale S1.4c (recalc).
Ri-deriva dai grezzi sim_out (misure raster dei simulati) le distanze per seme B1, B2, B7, B8, S1, S2, S3 con codice
proprio: territori osservati = Voronoi dei centroidi raster con shapely.voronoi_polygons (non tessellate_voronoi),
ritagliato alla finestra del ROI, orfani riassegnati con la regola D-S1.3.2 ricodificata (RECALC_FRAG=0 la disattiva); d = distanza dal bordo del proprio territorio; B-8 con cKDTree (d_ij < r_i + r_j).
Prima valida l'operatore sui generatori reali (a_t, d_gen di S1.4) e ricalcola b7_real/b8_real.
Uso: .venv/bin/python r03_obs.py [ncores]
"""
import sys, os, glob, json, numpy as np, pandas as pd, shapely
from scipy.spatial import cKDTree
from scipy.stats import ks_2samp
from multiprocessing import Pool
B = "/mnt/micron/geo_spatialtrans/S1.4c/"; D = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/dump/"
FRAG = os.environ.get("RECALC_FRAG", "1") == "1"   # 1: riassegnazione degli orfani (D-S1.3.2); 0: celle multiparte intere
OUTD = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/obs_rederived" + ("" if FRAG else "_nofrag") + "/"; os.makedirs(OUTD, exist_ok=True)
SEG = dict(A1="cellpose_rgb", A2="cellpose_rgb", A3="cellpose_rgb", A4="stardist_he", A5="cellpose_rgb", A6="stardist_he")
ks = lambda a, b: float(ks_2samp(np.asarray(a, float), np.asarray(b, float), method="asymp").statistic)

def window(A, roi):
    w = pd.read_parquet(D + f"win_{A}_{roi}.parquet"); return shapely.from_wkb(w.wkb.iloc[0])

def _reassign(cell, pts):
    """D-S1.3.2 ricodificata: pezzo col generatore resta; orfani al territorio col confine condiviso piu' lungo
    (> 1e-6 µm; decisioni simultanee per passata; pari merito -> indice minore); orfani isolati tornano al donatore."""
    cell = cell.copy(); multi = np.nonzero(shapely.get_num_geometries(cell) > 1)[0]
    orph = []
    for i in multi:
        pcs = shapely.get_parts(cell[i]); own = np.nonzero(shapely.intersects(pcs, pts[i]))[0]
        own = own[0] if len(own) else int(np.argmax(shapely.area(pcs)))
        cell[i] = pcs[own]; orph += [(pcs[j], i) for j in range(len(pcs)) if j != own]
    if not orph: return cell, 0
    tree = shapely.STRtree(cell); pending = list(range(len(orph)))
    for _ in range(50):
        dec = {}
        for o in pending:
            g, fr = orph[o]
            cand = [j for j in tree.query(shapely.buffer(g, 1e-7)) if j != fr]
            if not cand: continue
            sh = shapely.length(shapely.intersection(shapely.boundary(g), shapely.boundary(cell[cand])))
            ok = sh > 1e-6
            if not ok.any(): continue
            cand = np.array(cand)[ok]; sh = sh[ok]
            best = cand[sh == sh.max()].min(); dec[o] = int(best)
        if not dec: break
        for o, j in dec.items():
            cell[j] = shapely.union(cell[j], orph[o][0])
        tree = shapely.STRtree(cell); pending = [o for o in pending if o not in dec]
        if not pending: break
    for o in pending:
        g, fr = orph[o]; cell[fr] = shapely.union(cell[fr], g)
    return cell, len(orph)

def vor(x, y, region):
    """a_t, d per punto; NaN per i punti fuori dalla regione o mancanti. Restituisce anche n esclusi."""
    n = len(x); a_t = np.full(n, np.nan); d = np.full(n, np.nan)
    ok = np.isfinite(x) & np.isfinite(y)
    pts_all = shapely.points(np.where(ok, x, 0), np.where(ok, y, 0))
    shapely.prepare(region)
    ok[ok] = shapely.intersects(region, pts_all[ok])
    id_ = np.nonzero(ok)[0]; pts = pts_all[id_]
    xmin, ymin, xmax, ymax = region.bounds; pad = 10 * max(xmax - xmin, ymax - ymin)
    env = shapely.box(xmin - pad, ymin - pad, xmax + pad, ymax + pad)
    vc = shapely.get_parts(shapely.voronoi_polygons(shapely.multipoints(pts), extend_to=env, ordered=True))
    assert len(vc) == len(pts)
    assert shapely.contains_xy(vc, x[id_], y[id_]).all(), "ordine delle celle di Voronoi"
    inside = shapely.within(vc, region)
    cell = vc.copy(); cell[~inside] = shapely.intersection(vc[~inside], region)
    if FRAG: cell, _ = _reassign(cell, pts)
    a_t[id_] = shapely.area(cell); d[id_] = shapely.distance(pts, shapely.boundary(cell))
    return a_t, d, int(n - ok.sum())

def b8(x, y, area, ev):
    ok = np.isfinite(x) & np.isfinite(area) & (area > 0); id_ = np.nonzero(ok)[0]
    r = np.sqrt(area[id_] / np.pi); P = np.column_stack([x[id_], y[id_]])
    pr = cKDTree(P).query_pairs(2 * r.max(), output_type="ndarray")
    dd = np.hypot(*(P[pr[:, 0]] - P[pr[:, 1]]).T); h = dd < r[pr[:, 0]] + r[pr[:, 1]]
    hit = np.zeros(len(x), bool); hit[id_[pr[h, 0]]] = True; hit[id_[pr[h, 1]]] = True
    return float(hit[ev & ok].mean())

def real_ctx(A, roi):
    cc = pd.read_parquet(D + f"cells_{A}_{roi}.parquet")
    assert (cc.idx.to_numpy() == np.arange(1, len(cc) + 1)).all()
    ev = cc.ev_prep.to_numpy().astype(bool)
    rs = pd.read_parquet(B + f"real_shape/{A}_{roi}_{SEG[A]}_shape.parquet").set_index("label").loc[cc.shape_label.to_numpy()[ev]]
    return cc, ev, dict(nc=(cc.a_n / cc.a_t).to_numpy()[ev], an=cc.a_n.to_numpy()[ev], s1=float((rs.n_c2 > 0).mean()),
                         s2=(rs.n_c2 / rs.n_bnd).to_numpy(), s3=rs.sol.to_numpy())

def validate_real(args):
    A, roi = args
    cc, ev, _ = real_ctx(A, roi); reg = window(A, roi)
    a_t, d, nex = vor(cc.x.to_numpy(), cc.y.to_numpy(), reg)
    an = cc.a_n.to_numpy(); matched = cc.matched.to_numpy().astype(bool)
    b7 = float((np.sqrt(an[ev] / np.pi) > cc.d_gen.to_numpy()[ev]).mean())
    ar = np.where(matched & np.isfinite(an), an, cc.a_own.to_numpy())
    b8r = b8(cc.x.to_numpy(), cc.y.to_numpy(), ar, ev)
    return dict(archetype=A, roi=roi, n=len(cc), excluded=nex, at_rel_max=float(np.nanmax(np.abs(a_t - cc.a_t) / cc.a_t)),
                d_abs_max=float(np.nanmax(np.abs(d - cc.d_gen))), b7_real_mine=b7, b8_real_mine=b8r,
                win_area_check=float(reg.area))

def obs_one(f):
    b = os.path.basename(f)[:-8]; A, roi, fold, s, V = b.split("_")
    o = OUTD + b + ".json"
    if os.path.exists(o): return json.load(open(o))
    cc, ev, R = real_ctx(A, roi); n = len(cc)
    ps = pd.read_csv(D + "prep_scalars.csv"); pr = ps[(ps.archetype == A) & (ps.roi == roi)].iloc[0]
    m = pd.read_parquet(B + "sim_out/" + os.path.basename(f))
    def v(col):
        z = np.full(n, np.nan); z[m.label.to_numpy() - 1] = m[col].to_numpy(float); return z
    cx, cy, area, sol, nb, c2 = v("cx"), v("cy"), v("area"), v("sol"), v("n_bnd"), v("n_c2")
    a_t, d, nex = vor(cx, cy, window(A, roi))
    e = ev & np.isfinite(cx) & np.isfinite(a_t)
    b7 = float((np.sqrt(area[e] / np.pi) > d[e]).mean())
    out = dict(archetype=A, roi=roi, fold=int(fold[4:]), seed_idx=int(s[1:]), variant={"Kc": "K-c", "K1": "K-1", "K2": "K-2"}[V],
               n_eval=int(e.sum()), n_excluded=nex, b7=b7, b8=b8(cx, cy, area, ev), s1=float((c2[e] > 0).mean()),
               B1=ks(area[e] / a_t[e], R["nc"]), B2=ks(area[e], R["an"]), S1=abs(float((c2[e] > 0).mean()) - R["s1"]),
               S2=ks(c2[e] / nb[e], R["s2"]), S3=ks(sol[e], R["s3"]))
    out["B7"] = abs(b7 - pr.b7_real); out["B8"] = abs(out["b8"] - pr.b8_real)
    json.dump(out, open(o, "w")); return out

if __name__ == "__main__":
    nc = int(sys.argv[1]) if len(sys.argv) > 1 else 28
    rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    with Pool(nc) as p:
        val = pd.DataFrame(p.map(validate_real, list(zip(rois.archetype, rois.roi_id))))
    ps = pd.read_csv(D + "prep_scalars.csv"); val = val.merge(ps, on=["archetype", "roi"])
    val["d_b7"] = val.b7_real_mine - val.b7_real; val["d_b8"] = val.b8_real_mine - val.b8_real
    val.to_csv(OUTD.rstrip("/") + "_validate_real.csv", index=False)
    print("VALIDA operatore sui generatori reali: max rel a_t", val.at_rel_max.max(), "max |d|", val.d_abs_max.max(), "esclusi", val.excluded.sum())
    print("b7_real max|mio-prep|", np.abs(val.d_b7).max(), "b8_real max|mio-prep|", np.abs(val.d_b8).max())
    print("finestra: max|area poly - area_um2|", np.abs(val.win_area_check - val.win_area).max())
    fs = sorted(glob.glob(B + "sim_out/*_s[01][0-9]_K[c12].parquet")); fs = [f for f in fs if "_s00_" not in f]
    with Pool(nc) as p:
        res = pd.DataFrame(p.map(obs_one, fs, chunksize=1))
    res.to_csv(OUTD.rstrip("/") + ".csv", index=False)
    ob = pd.read_csv(D + "obs_all.csv")
    mm = res.merge(ob, on=["archetype", "roi", "fold", "seed_idx", "variant"], suffixes=("", "_rep"))
    print("righe", len(res), "unite", len(mm))
    for col in ["n_eval", "n_excluded", "b7", "b8", "s1", "B1", "B2", "B7", "B8", "S1", "S2", "S3"]:
        dif = np.abs(mm[col] - mm[col + "_rep"])
        print(f"{col}: max|mio - obs| {dif.max():.3e}  n>1e-9 {int((dif > 1e-9).sum())}")
