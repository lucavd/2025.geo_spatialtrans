#!/usr/bin/env python
"""tools/S1.4f_explore.py - S1.4f, esplorazione dichiarata prima della pre-registrazione (docs/steps/S1.4f.md). Solo dati veri.
F1. Angolo fra l'asse maggiore del nucleo vero e l'asse maggiore del suo territorio di Voronoi (centri = centroidi del segmentatore di
    forma, come E5 di S1.4e), per tessuto e per classe di distanza dal vicino; forma della distribuzione: von Mises assiale (2*delta ~ vM(0, k))
    contro miscela (quota w allineata esattamente, il resto uniforme) e contro Cauchy avvolta (2*delta ~ WC(0, rho)), stimate per massima verosimiglianza, con l'errore di misura dell'asse
    del territorio non modellato (limite dichiarato).
Cellule usate: nucleo allungato (asse maggiore / minore >= 1.2), territorio allungato (>= 1.1), centro a >= 20 um dal bordo del ROI,
territorio <= 4 volte la mediana del ROI (esclude le cellule al bordo di spazi vuoti). Convenzione dell'asse del nucleo come E5:
direzione (x, y) = (sin t, cos t) con t = `orientation` di skimage.
Uscite: results/S1.4f/explore/S1.4f_F1_cells.csv.gz, S1.4f_F1_summary.csv"""
import numpy as np, pandas as pd, shapely
from pathlib import Path
from scipy.spatial import cKDTree
from scipy.special import i0e, i1e
from scipy.optimize import brentq, minimize_scalar
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/S1.4f/explore"; OUT.mkdir(parents=True, exist_ok=True)
SHAPE_SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}
BINS = [0, 5, 6, 8, 10, np.inf]; LAB = ["<5", "5-6", "6-8", "8-10", ">=10"]

def poly_axis(g):
    """Angolo dell'asse maggiore (rad, da +x) e allungamento sqrt(l1/l2) di un poligono semplice (momenti con la formula di shoelace)."""
    xy = np.asarray(shapely.geometry.polygon.orient(g, 1.0).exterior.coords); xy = xy - xy[0]
    x, y = xy[:-1, 0], xy[:-1, 1]; x2, y2 = xy[1:, 0], xy[1:, 1]; cr = x * y2 - x2 * y; a = cr.sum() / 2
    cx = ((x + x2) * cr).sum() / 6 / a; cy = ((y + y2) * cr).sum() / 6 / a
    mxx = ((x * x + x * x2 + x2 * x2) * cr).sum() / 12 / a - cx * cx; myy = ((y * y + y * y2 + y2 * y2) * cr).sum() / 12 / a - cy * cy
    mxy = ((x * y2 + 2 * x * y + 2 * x2 * y2 + x2 * y) * cr).sum() / 24 / a - cx * cy
    t = np.sqrt(((mxx - myy) / 2) ** 2 + mxy ** 2); l1 = (mxx + myy) / 2 + t; l2 = max((mxx + myy) / 2 - t, 1e-15)
    return 0.5 * np.arctan2(2 * mxy, mxx - myy), np.sqrt(l1 / l2), a

def axial(d):
    """Differenza assiale in (-pi/2, pi/2]."""
    return (d + np.pi / 2) % np.pi - np.pi / 2

def fit_vm(phi):
    """phi = 2*delta in (-pi, pi]; MLE di k per vM(0, k): A(k) = I1/I0 = media di cos(phi)."""
    c = np.mean(np.cos(phi))
    if c <= 0: return 0.0, -len(phi) * np.log(2 * np.pi)
    f = lambda k: i1e(k) / i0e(k) - c
    k = 1e4 if f(1e4) <= 0 else brentq(f, 1e-9, 1e4)
    ll = np.sum(k * np.cos(phi)) - len(phi) * (np.log(2 * np.pi) + np.log(i0e(k)) + k)
    return float(k), float(ll)

def fit_mix(delta, h):
    """Miscela: con prob w delta ~ N(0, h) (allineato, h = larghezza fissa che rappresenta l'errore di misura), altrimenti uniforme su
    (-pi/2, pi/2]. MLE di w a h fissato; restituisce w e log-verosimiglianza (densita' di delta)."""
    g = np.exp(-0.5 * (delta / h) ** 2) / (h * np.sqrt(2 * np.pi)); u = 1 / np.pi
    nll = lambda w: -np.sum(np.log(w * g + (1 - w) * u))
    r = minimize_scalar(nll, bounds=(1e-6, 1 - 1e-6), method="bounded"); return float(r.x), float(-r.fun)

def fit_wc(phi):
    """phi = 2*delta; MLE di rho per la Cauchy avvolta: f(phi) = (1 - rho^2) / (2 pi (1 + rho^2 - 2 rho cos phi)). Log-verosimiglianza
    come densita' di delta (+ n log 2)."""
    nll = lambda r: -np.sum(np.log(1 - r * r) - np.log(2 * np.pi) - np.log(1 + r * r - 2 * r * np.cos(phi)))
    o = minimize_scalar(nll, bounds=(0.0, 0.999), method="bounded"); return float(o.x), float(-o.fun + len(phi) * np.log(2))

def main():
    rois = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv"); rows = []
    for r in rois.itertuples():
        A = r.archetype; seg = SHAPE_SEG[A]; side = r.side_px * r.um_per_px
        d = pd.read_parquet(ROOT / f"results/R2/nuclei/{A}_{r.roi_id}_{seg}_native.parquet")
        d = d[(d.scale == 1) & d.keep].reset_index(drop=True)
        P = d[["x_um", "y_um"]].to_numpy()
        vor = shapely.voronoi_polygons(shapely.multipoints(P), extend_to=shapely.box(-50, -50, side + 50, side + 50), ordered=True)
        V = list(shapely.get_parts(vor)); assert len(V) == len(P)
        V = [v.intersection(shapely.box(0, 0, side, side)) for v in V]
        ax = np.array([poly_axis(v) if v.geom_type == "Polygon" and v.area > 0 else (np.nan, np.nan, np.nan) for v in V])
        dist, nn = cKDTree(P).query(P, k=2); dist = dist[:, 1]
        t = d.orientation.to_numpy(); thN = np.arctan2(np.cos(t), np.sin(t))           # direzione (sin t, cos t) in (x, y)
        ecc = d.eccentricity.to_numpy(); asp = 1 / np.sqrt(np.clip(1 - ecc ** 2, 1e-9, None))
        de = np.minimum.reduce([P[:, 0], side - P[:, 0], P[:, 1], side - P[:, 1]])
        rows.append(pd.DataFrame(dict(archetype=A, roi_id=r.roi_id, nn_um=dist, aspect=asp, terr_aspect=ax[:, 1], terr_area=ax[:, 2],
                                      delta=axial(thN - ax[:, 0]), edge_um=de, terr_rel=ax[:, 2] / np.nanmedian(ax[:, 2]))))
    E = pd.concat(rows, ignore_index=True); E["bin"] = pd.cut(E.nn_um, BINS, right=False, labels=LAB)
    E["use"] = (E.aspect >= 1.2) & (E.terr_aspect >= 1.1) & (E.edge_um >= 20) & (E.terr_rel <= 4)
    E.to_csv(OUT / "S1.4f_F1_cells.csv.gz", index=False)
    S = []
    U = E[E.use]; groups = [(A, "all", g) for A, g in U.groupby("archetype")] + [(A, str(b), g) for (A, b), g in U.groupby(["archetype", "bin"], observed=True)]
    for A, b, g in groups:
        dl = g.delta.to_numpy(); phi = 2 * dl
        k, llv = fit_vm(phi); llv_d = llv + len(dl) * np.log(2)                          # densita' di delta = 2 x densita' di phi
        w, llm = fit_mix(dl, np.radians(10)); rho, llw = fit_wc(phi)
        per = g.groupby("roi_id").delta.apply(lambda s: float(np.mean(np.abs(s) < np.pi / 4)))
        S.append(dict(archetype=A, bin=b, n=len(g), n_roi=g.roi_id.nunique(),
                      R2=float(np.hypot(np.mean(np.cos(phi)), np.mean(np.sin(phi)))), mean_cos2=float(np.mean(np.cos(phi))),
                      frac_lt45=float(np.mean(np.abs(dl) < np.pi / 4)), frac_lt45_roi_min=float(per.min()), frac_lt45_roi_max=float(per.max()),
                      median_abs_deg=float(np.degrees(np.median(np.abs(dl)))), kappa_vm=k, ll_vm=llv_d, w_mix=w, ll_mix=llm, dll_vm_minus_mix=llv_d - llm,
                      rho_wc=rho, ll_wc=llw, dll_vm_minus_wc=llv_d - llw))
    S = pd.DataFrame(S); S.to_csv(OUT / "S1.4f_F1_summary.csv", index=False); print(S.round(3).to_string())

if __name__ == "__main__":
    main()
