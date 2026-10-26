#!/usr/bin/env python
"""tools/S1.5e_explore.py - S1.5e, esplorazione dichiarata PRIMA della pre-registrazione (docs/steps/S1.5e.md).
Domande (solo dati veri e unita' gia' simulate di S1.5d, nessuna simulazione nuova):
(E1) la densita' di riferimento (conteggio manuale di R2b) e' riferita all'area con i vuoti? Le finestre di R2b, poste a caso nel ROI,
     contengono tanto spazio vuoto e tanti nuclei segmentati per area quanto il ROI intero?
(E2) l'affollamento dei centri nelle finestre (vicini entro r, normalizzati per la densita' della finestra; distanza dal vicino):
     punti manuali di Luca contro nuclei segmentati e bracci U, N, RV di S1.5d. Se RV e' gia' come il manuale, togliere cellule
     abbasserebbe l'affollamento sotto il vero.
(E3) da dove viene l'eccesso di contatti nei bracci con i vuoti: distanza dal vicino e area del nucleo dei focali, N contro RV.
(E4) forma dei vuoti: componenti connesse dei pixel della maschera lontani dal centro piu' vicino piu' di 3 distanze mediane dal
     vicino (stessa definizione di H_rel_3); allungamento (radice del rapporto degli autovalori dei momenti), spessore
     (2 x massimo della trasformata di distanza), coerenza dell'orientamento; vero contro RV e VF.
Ingressi in sola lettura: maschere R2, nuclei veri (R6, come S1.5b-S1.5d), finestre e annotazioni R2b, celle e raster di S1.5d.
Uscite: results/S1.5e/S1.5e_explore_*.csv. Uso: .venv/bin/python tools/S1.5e_explore.py [n_proc]
"""
import sys, json, importlib.util, numpy as np, pandas as pd
from multiprocessing import Pool
from types import SimpleNamespace
from pathlib import Path
from PIL import Image
from scipy.spatial import cKDTree
from scipy import ndimage as ndi
sys.path.insert(0, "tools")

def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py")
ROIS = pd.read_csv("results/R2/R2_rois_checked.csv")
WINS = pd.read_csv("results/R2b/R2b_windows.csv")
WJ = {w["win_id"]: w for w in json.load(open("results/R2b/annotations/R2b_annotations_Luca.json"))["windows"]}
S5D = Path("/mnt/micron/geo_spatialtrans/S1.5d"); RES = Path("results/S1.5e")
ARMS = ["U", "N", "RV"]; ARMS_VOID = ["RV", "VF"]; SEEDS = range(1, 11); SEEDS_VOID = range(1, 4)
RADII = [4.0, 6.0, 8.0, 10.0, 15.0]; FOCAL_B = 15.0; FRAME_B = 30.0; U_VOID = 3.0; MIN_COMP_UM2 = 50.0

def mask_of(A, r):
    m = np.array(Image.open(f"results/R2/tissue_masks/{A}_{r}_valid_ds4.png"))
    if m.ndim == 3: m = m[..., 0]
    row = ROIS[(ROIS.archetype == A) & (ROIS.roi_id == r)].iloc[0]
    return m > 127, row.side_px * row.um_per_px / m.shape[1], row.side_px * row.um_per_px

def real_centres(A, r):
    _, d = cal.load(A); e = d[d.roi_id == r]
    return np.c_[e.x.values, e.y.values], e.interior.values.astype(bool)

def manual_points(w):
    upp = w["um_per_px"]; mg = w["margin_px"]; p = np.array(w["points"], float).reshape(-1, 2)
    return np.c_[w["x0_um"] + (p[:, 0] - mg) * upp, w["y0_um"] + (p[:, 1] - mg) * upp]

def nn_median(P, I=None):
    dd, _ = cKDTree(P).query(P, k=2); return float(np.median(dd[I, 1] if I is not None else dd[:, 1]))

def void_pixels(P, m, px, side, nnm):
    """Pixel della maschera (centro del pixel) a distanza > U_VOID * nnm dal centro piu' vicino, ad almeno FRAME_B dal bordo del ROI."""
    rr, cc = np.nonzero(m); x = (cc + 0.5) * px; y = (rr + 0.5) * px
    k = (x > FRAME_B) & (y > FRAME_B) & (x < side - FRAME_B) & (y < side - FRAME_B)
    d, _ = cKDTree(P).query(np.c_[x[k], y[k]])
    v = np.zeros_like(m, bool); v[rr[k][d > U_VOID * nnm], cc[k][d > U_VOID * nnm]] = True
    valid = np.zeros_like(m, bool); valid[rr[k], cc[k]] = True
    return v, valid

# ---- E1: le finestre di R2b sono rappresentative del ROI per spazio vuoto e densita'? ----
def e1_window(w):
    A, r = w.archetype, w.roi_id; m, px, side = mask_of(A, r); P, I = real_centres(A, r); nnm = nn_median(P, I)
    v, valid = void_pixels(P, m, px, side, nnm)
    s = w.side_px * w.um_per_px; c0, c1 = int(w.x0_um / px), int((w.x0_um + s) / px); r0, r1 = int(w.y0_um / px), int((w.y0_um + s) / px)
    vw = v[r0:r1, c0:c1]; mw = m[r0:r1, c0:c1]; valw = valid[r0:r1, c0:c1]
    inw = (P[:, 0] >= w.x0_um) & (P[:, 0] < w.x0_um + s) & (P[:, 1] >= w.y0_um) & (P[:, 1] < w.y0_um + s)
    inmask_roi = m[np.minimum((P[:, 1] / px).astype(int), m.shape[0] - 1), np.minimum((P[:, 0] / px).astype(int), m.shape[1] - 1)]
    Pm = manual_points(WJ[w.win_id]); inm = (Pm[:, 0] >= w.x0_um) & (Pm[:, 0] < w.x0_um + s) & (Pm[:, 1] >= w.y0_um) & (Pm[:, 1] < w.y0_um + s)
    return dict(win_id=w.win_id, archetype=A, roi_id=r, nn_median_roi=nnm,
                void_frac_window=vw[valw].mean() if valw.any() else np.nan, void_frac_roi=v[valid].mean(),
                mask_frac_window=mw.mean(), dens_seg_window=inw.sum() / (mw.sum() * px * px) * 1e6,
                dens_seg_roi=inmask_roi.sum() / (m.sum() * px * px) * 1e6, dens_manual_window=inm.sum() / (s * s) * 1e6)

# ---- E2, E3: affollamento dei centri e contatti nelle finestre ----
def crowd(P, x0, y0, s, extra=None):
    F = (P[:, 0] >= x0 + FOCAL_B) & (P[:, 0] < x0 + s - FOCAL_B) & (P[:, 1] >= y0 + FOCAL_B) & (P[:, 1] < y0 + s - FOCAL_B)
    inw = (P[:, 0] >= x0) & (P[:, 0] < x0 + s) & (P[:, 1] >= y0) & (P[:, 1] < y0 + s)
    lam = inw.sum() / (s * s); tr = cKDTree(P); o = dict(n_window=int(inw.sum()), n_focal=int(F.sum()), dens_window=lam * 1e6)
    if F.sum() < 3: return o
    for rad in RADII:
        k = np.array([len(x) - 1 for x in tr.query_ball_point(P[F], rad)], float)
        o[f"k{int(rad)}"] = k.mean(); o[f"knorm{int(rad)}"] = k.mean() / (lam * np.pi * rad * rad)
    dd, _ = tr.query(P[F], k=2); o["nn_median"] = float(np.median(dd[:, 1])); o["nn_q10"] = float(np.quantile(dd[:, 1], 0.1))
    o["frac_nn_lt6"] = float(np.mean(dd[:, 1] < 6.0))
    if extra is not None:
        a, c = extra; o["nuc_area_median"] = float(np.median(a[F])); o["s1_all_focal"] = float(np.mean(c[F] > 0))
        # contatto per classe di distanza dal vicino (in unita' della mediana del tessuto non serve: riportiamo in um)
        for lo, hi in [(0, 6), (6, 8), (8, 10), (10, 14), (14, np.inf)]:
            k = (dd[:, 1] >= lo) & (dd[:, 1] < hi); cf = c[F] > 0
            o[f"frac_nn_{lo}_{hi}"] = float(k.mean()); o[f"s1_nn_{lo}_{hi}"] = float(cf[k].mean()) if k.any() else np.nan
    return o

def e2_window(w):
    A, r = w.archetype, w.roi_id; s = w.side_px * w.um_per_px; x0, y0 = w.x0_um, w.y0_um; rows = []
    rows.append(dict(win_id=w.win_id, archetype=A, roi_id=r, source="manual", seed=0, **crowd(manual_points(WJ[w.win_id]), x0, y0, s)))
    P, _ = real_centres(A, r); rows.append(dict(win_id=w.win_id, archetype=A, roi_id=r, source="segmented", seed=0, **crowd(P, x0, y0, s)))
    for arm in ARMS:
        for sd in SEEDS:
            t = pd.read_parquet(S5D / "raster" / arm / f"{A}_{r}_s{sd:02d}.parquet", columns=["cx", "cy", "area", "n_c2"])
            rows.append(dict(win_id=w.win_id, archetype=A, roi_id=r, source=arm, seed=sd,
                             **crowd(np.c_[t.cx.values, t.cy.values], x0, y0, s, extra=(t.area.values, t.n_c2.values))))
    return rows

# ---- E4: forma dei vuoti ----
def components(v, px):
    lab, n = ndi.label(v, structure=np.ones((3, 3), bool)); out = []
    if n == 0: return out
    dt = ndi.distance_transform_edt(v) * px
    objs = ndi.find_objects(lab)
    for i, sl in enumerate(objs, start=1):
        cm = lab[sl] == i; a = cm.sum() * px * px
        if a < MIN_COMP_UM2: continue
        yy, xx = np.nonzero(cm); X = np.c_[xx, yy] * px; C = np.cov(X.T)
        ev, evec = np.linalg.eigh(C); ev = np.maximum(ev, 1e-12)
        out.append(dict(area_um2=a, elong=float(np.sqrt(ev[1] / ev[0])), thick_um=float(2 * dt[sl][cm].max()),
                        theta=float(np.arctan2(evec[1, 1], evec[0, 1]))))
    return out

def summarise_comps(cs):
    if not cs: return dict(n_comp=0, void_area_um2=0.0, elong_wmed=np.nan, frac_area_elong_gt3=np.nan, thick_wmed=np.nan, orient_coherence=np.nan)
    d = pd.DataFrame(cs); wgt = d.area_um2 / d.area_um2.sum()
    def wmed(x):
        o = np.argsort(x.values); cw = np.cumsum(wgt.values[o]); return float(x.values[o][np.searchsorted(cw, 0.5)])
    big = d.elong > 1.5
    coh = float(np.abs(np.sum(wgt[big] * np.exp(2j * d.theta[big]))) / max(wgt[big].sum(), 1e-12)) if big.any() else np.nan
    return dict(n_comp=len(d), void_area_um2=float(d.area_um2.sum()), elong_wmed=wmed(d.elong), frac_area_elong_gt3=float(wgt[d.elong > 3].sum()),
                thick_wmed=wmed(d.thick_um), orient_coherence=coh)

def e4_roi(args):
    A, r = args; m, px, side = mask_of(A, r); P, I = real_centres(A, r); nnm = nn_median(P, I)
    v, _ = void_pixels(P, m, px, side, nnm); rows = [dict(archetype=A, roi_id=r, source="real", seed=0, **summarise_comps(components(v, px)))]
    for arm in ARMS_VOID:
        for sd in SEEDS_VOID:
            t = pd.read_parquet(S5D / "cells" / arm / f"{A}_{r}_s{sd:02d}.parquet", columns=["x", "y", "interior"])
            Ps = np.c_[t.x.values, t.y.values]; vs, _ = void_pixels(Ps, m, px, side, nn_median(Ps, t.interior.values.astype(bool)))
            rows.append(dict(archetype=A, roi_id=r, source=arm, seed=sd, **summarise_comps(components(vs, px))))
    return rows

if __name__ == "__main__":
    RES.mkdir(parents=True, exist_ok=True); nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    wl = [SimpleNamespace(**d) for d in WINS.to_dict("records")]   # namedtuple di itertuples non si serializza fra processi
    with Pool(nproc) as p:
        e1 = p.map(e1_window, wl); e2 = p.map(e2_window, wl)
        e4 = p.map(e4_roi, list(zip(ROIS.archetype, ROIS.roi_id)))
    pd.DataFrame(e1).to_csv(RES / "S1.5e_explore_density_ref.csv", index=False)
    pd.DataFrame([x for o in e2 for x in o]).to_csv(RES / "S1.5e_explore_crowding.csv", index=False)
    pd.DataFrame([x for o in e4 for x in o]).to_csv(RES / "S1.5e_explore_void_shape.csv", index=False)
    print("ok")
