#!/usr/bin/env python
"""tools/S1.4f_posthoc.py - S1.4f, analisi post hoc PH-1 (dichiarata dopo il risultato di P-5, docs/steps/S1.4f.md).
Domanda: a 5-8 um dal vicino l'asse del territorio e' di traverso al vicino anche nei tessuti veri, o solo nei simulati?
Quota dei territori allungati (>= 1.1) con l'asse maggiore di traverso al vicino piu' vicino (> 45 gradi), per classe di distanza:
  veri: territori di Voronoi dei centroidi del segmentatore di forma (stessa selezione di F1: centro a >= 20 um dal bordo, territorio <= 4 x mediana);
  simulati: territori delle unita' del banco B, semi 1-5, cellule interne (asse e allungamento dai momenti del poligono, come G-4).
Accanto: la stessa quota per i nuclei veri allungati (>= 1.2), dalla stessa selezione di cellule. Unita' = ROI (5 per tessuto).
Uscita: results/S1.4f/S1.4f_PH1_territory_nn.csv. Uso: .venv/bin/python tools/S1.4f_posthoc.py [ncores]"""
import importlib.util, glob, os, sys
import numpy as np, pandas as pd, shapely
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
_s = importlib.util.spec_from_file_location("m14e", ROOT / "tools/S1.4e_measure.py"); ME = importlib.util.module_from_spec(_s); sys.modules["m14e"] = ME; _s.loader.exec_module(ME)
M = Path("/mnt/micron/geo_spatialtrans/S1.4f"); SHAPE_SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}
CLS = {"nn<5": (0, 5), "nn5-6": (5, 6), "nn6-8": (6, 8)}

def frac_by_cls(d, psi, keep):
    return {c: (int((keep & (d >= lo) & (d < hi)).sum()), float(np.mean(psi[keep & (d >= lo) & (d < hi)] > 45)) if (keep & (d >= lo) & (d < hi)).any() else np.nan)
            for c, (lo, hi) in CLS.items()}

def real_one(r):
    A = r.archetype; side = r.side_px * r.um_per_px
    d = pd.read_parquet(ROOT / f"results/R2/nuclei/{A}_{r.roi_id}_{SHAPE_SEG[A]}_native.parquet"); d = d[(d.scale == 1) & d.keep].reset_index(drop=True)
    P = d[["x_um", "y_um"]].to_numpy()
    V = list(shapely.get_parts(shapely.voronoi_polygons(shapely.multipoints(P), extend_to=shapely.box(-50, -50, side + 50, side + 50), ordered=True)))
    V = [v.intersection(shapely.box(0, 0, side, side)) for v in V]
    S = np.array([ME.poly_shape(v) if v.area > 0 else (np.nan,) * 5 for v in V]); area = np.array([v.area for v in V])
    dist, nn = cKDTree(P).query(P, k=2); dist, nn = dist[:, 1], nn[:, 1]; u = (P[nn] - P) / dist[:, None]
    de = np.minimum.reduce([P[:, 0], side - P[:, 0], P[:, 1], side - P[:, 1]]); base = (de >= 20) & (area <= 4 * np.nanmedian(area))
    psi_t = ME.psi_deg(S[:, 2], u[:, 0], u[:, 1])
    t = d.orientation.to_numpy(); psi_n = np.degrees(np.arccos(np.clip(np.abs(np.sin(t) * u[:, 0] + np.cos(t) * u[:, 1]), 0, 1)))
    asp_n = 1 / np.sqrt(np.clip(1 - d.eccentricity.to_numpy() ** 2, 1e-9, None))
    out = []
    for what, psi, keep in [("territory", psi_t, base & (S[:, 3] >= 1.1)), ("nucleus", psi_n, base & (asp_n >= 1.2))]:
        for c, (n, f) in frac_by_cls(dist, psi, keep).items(): out.append(dict(source="real", what=what, archetype=A, roi_id=r.roi_id, seed="", cls=c, n=n, frac=f))
    return out

def sim_one(f):
    t = pd.read_parquet(f); T = shapely.from_wkb(t.terr_wkb.to_numpy()); S = np.array([ME.poly_shape(g) for g in T])
    key = os.path.basename(f)[:-8]; cells = pd.read_parquet(M / "B" / "cells" / f"{key}.parquet", columns=["interior"])
    P = t[["x", "y"]].to_numpy(); dist, nn = cKDTree(P).query(P, k=2); dist, nn = dist[:, 1], nn[:, 1]; u = (P[nn] - P) / dist[:, None]
    psi = ME.psi_deg(S[:, 2], u[:, 0], u[:, 1]); A, roi, sd = key.split("_")
    return [dict(source="sim", what="territory", archetype=A, roi_id=roi, seed=sd, cls=c, n=n, frac=fr)
            for c, (n, fr) in frac_by_cls(dist, psi, cells.interior.to_numpy() & (S[:, 3] >= 1.1)).items()]

def main(nc):
    from multiprocessing import Pool
    rois = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")
    with Pool(nc) as p:
        a = p.map(real_one, [pd.Series(x) for x in rois.to_dict('records')]); b = p.map(sim_one, sorted(glob.glob(str(M / "B" / "g4" / "V2" / "*.parquet"))))
    D = pd.DataFrame([x for y in a + b for x in y]); D.to_csv(ROOT / "results/S1.4f/S1.4f_PH1_territory_nn.csv", index=False)
    R = D.groupby(["source", "what", "archetype", "roi_id", "cls"]).frac.mean().groupby(["source", "what", "archetype", "cls"]).agg(["mean", "std", "count"])
    print(R.round(3).unstack("cls").to_string())

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 15)
