#!/usr/bin/env python
"""review/code/04_raster_checks.py — prove eseguibili sull'operatore raster di S1.4c.
 (a) rasterizzazione indipendente (shapely.contains_xy sui centri dei pixel) contro skimage.draw.polygon del banco;
 (b) bordo/contatto con codice proprio (spostamenti espliciti) contro sim_out e real_shape;
 (c) fessura di 2e-5 um fra nuclei K-1/K-2: pixel persi e effetto su area, S-1, S-2, S-3 (chiudendo la fessura);
 (d) nuclei multiparte (K-1/K-2) e etichette sconnesse (reali e simulati);
 (e) solidita' di skimage su oggetti piccoli: dischi a sfasamento casuale, valori distinti, dipendenza dalla taglia.
Uso (root del repo): .venv/bin/python results/S1.4c/review/code/04_raster_checks.py [ncores]
"""
import os, sys, glob, json, importlib.util
import numpy as np, pandas as pd, shapely
from scipy import ndimage as ndi
from multiprocessing import Pool
spec = importlib.util.spec_from_file_location("rb", "tools/S1.4c_raster.py"); rb = importlib.util.module_from_spec(spec); spec.loader.exec_module(rb)
S = "/mnt/micron/geo_spatialtrans/S1.4c"; SC = "/mnt/micron/geo_spatialtrans/S1.4c_review/code"; OUTD = "results/S1.4c/review/code"
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
SEG = rb.SHAPE_SEG

def own_raster(geoms, ids, shape, upp):
    lab = np.zeros(shape, np.int64); H, W = shape
    for g, i in zip(geoms, ids):
        if g is None or shapely.is_empty(g): continue
        x0, y0, x1, y1 = shapely.bounds(g)
        c0, c1 = max(int(np.ceil(x0 / upp)), 0), min(int(np.floor(x1 / upp)), W - 1)
        r0, r1 = max(int(np.ceil(y0 / upp)), 0), min(int(np.floor(y1 / upp)), H - 1)
        if c1 < c0 or r1 < r0: continue
        rr, cc = np.mgrid[r0:r1 + 1, c0:c1 + 1]
        ins = shapely.contains_xy(g, cc.ravel() * upp, rr.ravel() * upp)
        lab[rr.ravel()[ins], cc.ravel()[ins]] = i
    return lab

def own_contact(lab, d):
    """n_bnd e n_c<d> per etichetta con spostamenti espliciti (fuori immagine = 0)."""
    H, W = lab.shape; P = np.pad(lab, d + 1)
    core = P[d + 1:d + 1 + H, d + 1:d + 1 + W]
    bnd = np.zeros_like(lab, bool)
    for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        bnd |= P[d + 1 + dr:d + 1 + dr + H, d + 1 + dc:d + 1 + dc + W] != core
    bnd &= core > 0
    oth = np.zeros_like(lab, bool)
    for dr in range(-d, d + 1):
        for dc in range(-d, d + 1):
            nb = P[d + 1 + dr:d + 1 + dr + H, d + 1 + dc:d + 1 + dc + W]; oth |= (nb > 0) & (nb != core)
    M = int(lab.max()) + 1
    return np.bincount(lab[bnd], minlength=M), np.bincount(lab[bnd & oth], minlength=M)

def roi_info(A, roi):
    rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]; return float(rw.um_per_px), int(rw.side_px)

def ev_of(A, roi):
    cc = pd.read_parquet(f"{SC}/exp/cells_{A}_{roi}.parquet"); return cc.ev_prep.to_numpy()

def sim_case(f):
    b = os.path.basename(f)[:-8]; A, roi, fold, s, V = b.split("_"); upp, side = roi_info(A, roi)
    t = pd.read_parquet(f); geoms = shapely.from_wkb(t.wkb.to_numpy()); ids = t.idx.to_numpy()
    labB, ncf = rb.rasterize(geoms, ids, (side, side), upp)
    labO = own_raster(geoms, ids, (side, side), upp)
    diff = int(np.sum(labB != labO))
    so = pd.read_parquet(f"{S}/sim_out/{os.path.basename(f)}").set_index("label")
    nb, n2 = own_contact(labB, 2); _, n1 = own_contact(labB, 1)
    L = so.index.to_numpy()
    dcont = int(np.sum(nb[L] != so.n_bnd.to_numpy()) + np.sum(n2[L] != so.n_c2.to_numpy()) + np.sum(n1[L] != so.n_c1.to_numpy()))
    nparts = shapely.get_num_geometries(geoms); multi = int(np.sum(nparts > 1))
    ncomp = 0
    sl = ndi.find_objects(labB)
    for i, sli in enumerate(sl, start=1):
        if sli is None: continue
        _, k = ndi.label(labB[sli] == i); ncomp += k > 1
    res = dict(file=b, archetype=A, variant=V, n=len(ids), pix_diff_vs_own_raster=diff, contact_mismatch=dcont, n_conflict=ncf,
               n_multipart=multi, n_label_disconnected=int(ncomp))
    # (c) fessura: pixel non etichettati il cui centro dista <= 5e-5 um da un nucleo
    if V in ("K1", "K2"):
        un = labB == 0
        halo = un & (ndi.maximum_filter(labB, size=3) > 0)
        rr, cc_ = np.nonzero(halo); pts = shapely.points(cc_ * upp, rr * upp)
        tree = shapely.STRtree(geoms); qi, gi = tree.query(pts, predicate="dwithin", distance=5e-5)
        gap_r, gap_c, gap_l = rr[qi], cc_[qi], ids[gi]
        # un pixel puo' stare entro 5e-5 da due nuclei (e' nella fessura fra i due): conta i pixel, assegna al primo
        _, first = np.unique(qi, return_index=True)
        labG = labB.copy(); labG[gap_r[first], gap_c[first]] = gap_l[first]
        n_two = int(np.sum(np.bincount(qi) >= 2))
        ev = ev_of(A, roi)
        mB = rb.measure(labB, upp, frac_out=False).set_index("label"); mG = rb.measure(labG, upp, frac_out=False).set_index("label")
        e = np.nonzero(ev)[0] + 1; e = e[np.isin(e, mB.index) & np.isin(e, mG.index)]
        def summ(m):
            z = m.loc[e]; return dict(S1=float(np.mean(z.n_c2 > 0)), S2_med=float(np.median(z.n_c2 / z.n_bnd)), S3_med=float(np.median(z.sol)),
                                      area_mean=float(z.area.mean()), S2_mean=float(np.mean(z.n_c2 / z.n_bnd)), S3_mean=float(np.mean(z.sol)))
        sB, sG = summ(mB), summ(mG)
        res.update(gap_pixels=int(len(first)), gap_pixels_between_two=n_two, labeled_pixels=int(np.sum(labB > 0)),
                   n_label_changed=int(np.sum(np.abs(mB.loc[e].n_px.to_numpy() - mG.loc[e].n_px.to_numpy()) > 0)),
                   **{f"d_{k}": sG[k] - sB[k] for k in sB},
                   ks_S2=float(rb_ks(mB.loc[e].n_c2 / mB.loc[e].n_bnd, mG.loc[e].n_c2 / mG.loc[e].n_bnd)),
                   ks_S3=float(rb_ks(mB.loc[e].sol, mG.loc[e].sol)))
    return res

def rb_ks(a, b):
    a = np.sort(np.asarray(a, float)); b = np.sort(np.asarray(b, float)); v = np.unique(np.concatenate([a, b]))
    return np.max(np.abs(np.searchsorted(a, v, "right") / len(a) - np.searchsorted(b, v, "right") / len(b)))

def real_case(A, roi):
    seg = SEG[A]; gen = pd.read_parquet(f"{rb.R3}/real/{A}_{roi}_{seg}_gen.parquet")
    lab = np.load(f"{rb.MASKS}/{A}_{roi}_{seg}_native.npz")["labels"].astype(np.int64)
    keep = np.zeros(int(lab.max()) + 1, bool); keep[gen.label.to_numpy().astype(np.int64)] = True; lab[~keep[lab]] = 0
    rs = pd.read_parquet(f"{S}/real_shape/{A}_{roi}_{seg}_shape.parquet").set_index("label"); L = rs.index.to_numpy()
    nb, n2 = own_contact(lab, 2); _, n1 = own_contact(lab, 1)
    dcont = int(np.sum(nb[L] != rs.n_bnd.to_numpy()) + np.sum(n2[L] != rs.n_c2.to_numpy()) + np.sum(n1[L] != rs.n_c1.to_numpy()))
    ncomp = 0
    for i, sli in enumerate(ndi.find_objects(lab), start=1):
        if sli is None: continue
        _, k = ndi.label(lab[sli] == i); ncomp += k > 1
    # contatto con l'arbitro sulle coppie (contact_pairs) contro own: etichette con n_c2>0 == etichette in qualche coppia
    pr = pd.read_parquet(f"{S}/real_shape/{A}_{roi}_{seg}_pairs.parquet"); inpair = np.zeros(int(lab.max()) + 1, bool)
    inpair[pr.label_i.to_numpy()] = True; inpair[pr.label_j.to_numpy()] = True
    n_pair_mismatch = int(np.sum(inpair[L] != (n2[L] > 0)))
    return dict(file=f"{A}_{roi}_{seg}_real", archetype=A, n=len(L), contact_mismatch=dcont, n_label_disconnected=int(ncomp),
                pairs_vs_nc2_mismatch=n_pair_mismatch, area_px_median=float(np.median(rs.n_px)), n_sol_distinct=int(rs.sol.nunique()))

def solidity_study():
    from skimage.measure import regionprops
    rng = np.random.default_rng(1); rows = []
    for area_px in (60, 120, 190, 300, 500, 800):
        r = np.sqrt(area_px / np.pi)
        for t in range(200):
            x0, y0 = rng.uniform(20, 21, 2); yy, xx = np.mgrid[0:60, 0:60]
            L = ((xx - x0) ** 2 + (yy - y0) ** 2 <= r ** 2).astype(np.int64)
            p = regionprops(L)[0]; rows.append(dict(area_px=area_px, sol=p.solidity))
    d = pd.DataFrame(rows)
    return d.groupby("area_px").sol.agg(["mean", "std", "min", "max", "nunique"]).reset_index()

if __name__ == "__main__":
    nc = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    fs = sorted(glob.glob(f"{S}/gen_in/*_fold1_s01_K[c12].parquet")) + sorted(glob.glob(f"{S}/gen_in/*_fold3_s07_K[12].parquet"))
    with Pool(nc) as p:
        sim = p.map(sim_case, fs, chunksize=1)
        real = p.starmap(real_case, [(r.archetype, r.roi_id) for r in rois.itertuples()])
    pd.DataFrame(sim).to_csv(f"{OUTD}/04_raster_sim.csv", index=False); pd.DataFrame(real).to_csv(f"{OUTD}/04_raster_real.csv", index=False)
    sol = solidity_study(); sol.to_csv(f"{OUTD}/04_solidity_discs.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(sim).to_string()); print(pd.DataFrame(real).to_string()); print(sol.to_string())
