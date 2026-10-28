#!/usr/bin/env python
"""tools/S1.4e_explore_e5.py - S1.4e, esplorazione dichiarata (docs/steps/S1.4e.md), E5: orientamento dei nuclei veri rispetto al
vicino piu' vicino, su tutti i nuclei del segmentatore di forma (SHAPE_SEG) nei 5 ROI per tessuto (scala nativa, keep).
psi = angolo fra asse maggiore e direzione del vicino (0 = il vicino sta sull'asse lungo, 90 = di fianco). Asse maggiore dalla colonna
`orientation` di skimage (angolo dall'asse delle righe): direzione (x, y) = (sin t, cos t). Controllo della convenzione: psi ricalcolato
dai pixel della maschera sulle coppie di E2 (S1.4e_E12_pairs.csv) contro psi da `orientation` degli stessi nuclei.
Nullo: orientazioni permutate fra i nuclei dello stesso ROI (stessa distribuzione di allungamenti, nessun legame col vicino).
Uscite: results/S1.4e/explore/S1.4e_E5_*.csv."""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/S1.4e/explore"
SHAPE_SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}
BINS = [0, 5, 6, 7, 8, 10, np.inf]; LAB = ["<5", "5-6", "6-7", "7-8", "8-10", ">=10"]
rng = np.random.default_rng(20261028)

def psi_of(theta, ux, uy):
    vx, vy = np.sin(theta), np.cos(theta)
    return np.degrees(np.arccos(np.clip(np.abs(vx * ux + vy * uy), 0, 1)))

def main():
    rois = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")
    rows = []
    for r in rois.itertuples():
        A = r.archetype; seg = SHAPE_SEG[A]
        d = pd.read_parquet(ROOT / f"results/R2/nuclei/{A}_{r.roi_id}_{seg}_native.parquet")
        d = d[(d.scale == 1) & d.keep].reset_index(drop=True)
        P = d[["x_um", "y_um"]].to_numpy(); dist, nn = cKDTree(P).query(P, k=2); dist, nn = dist[:, 1], nn[:, 1]
        u = (P[nn] - P) / dist[:, None]
        psi = psi_of(d.orientation.to_numpy(), u[:, 0], u[:, 1])
        perm = rng.permutation(len(d)); psi0 = psi_of(d.orientation.to_numpy()[perm], u[:, 0], u[:, 1])
        ecc = d.eccentricity.to_numpy(); asp = 1 / np.sqrt(np.clip(1 - ecc ** 2, 1e-9, None))
        rows.append(pd.DataFrame(dict(archetype=A, roi_id=r.roi_id, label=d.label, nn_um=dist, psi=psi, psi_null=psi0, aspect=asp,
                                      aspect_nn=asp[nn], area=d.area_um2)))
    E = pd.concat(rows); E["bin"] = pd.cut(E.nn_um, BINS, right=False, labels=LAB); E.to_csv(OUT / "S1.4e_E5_nuclei.csv.gz", index=False)
    # convenzione: psi dai pixel (E2) contro psi da orientation
    pr = pd.read_csv(OUT / "S1.4e_E12_pairs.csv"); pr = pr[pr.status == "ok"]
    pts = pd.read_csv(ROOT / "results/R2b/R2b_points.csv"); pts = pts[pts.rater == "Luca"].set_index(["win_id", "idx"])
    wins = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv").set_index("win_id"); chk = []
    for x in pr.itertuples():
        seg = SHAPE_SEG[x.archetype]; roi = wins.loc[x.win_id].roi_id
        d = pd.read_parquet(ROOT / f"results/R2/nuclei/{x.archetype}_{roi}_{seg}_native.parquet", columns=["label", "x_um", "y_um", "orientation", "scale"])
        d = d[d.scale == 1].set_index("label")
        li, lj = int(pts.loc[(x.win_id, x.i)][f"label_{seg}"]), int(pts.loc[(x.win_id, x.j)][f"label_{seg}"])
        gi, gj = d.loc[li], d.loc[lj]; v = np.array([gj.x_um - gi.x_um, gj.y_um - gi.y_um]); v /= np.hypot(*v)
        chk.append(dict(psi_pix=x.psi_i, psi_orient=float(psi_of(gi.orientation, v[0], v[1]))))
    C = pd.DataFrame(chk); C["diff"] = (C.psi_pix - C.psi_orient).abs()
    print("convenzione: coppie", len(C), "| mediana |diff| gradi", round(C["diff"].median(), 2), "| quota < 5 gradi", round((C["diff"] < 5).mean(), 3))
    # riassunto per tessuto e classe, incertezza per ROI (bootstrap dei ROI, B = 2000)
    S = []
    for (A, b), g in E.groupby(["archetype", "bin"], observed=True):
        el = g[g.aspect >= 1.2]
        per = el.groupby("roi_id").agg(f=("psi", lambda s: np.mean(s > 45)), f0=("psi_null", lambda s: np.mean(s > 45)), n=("psi", "size"))
        bs = [];
        for _ in range(2000):
            k = rng.integers(0, len(per), len(per)); w = per.n.to_numpy()[k]
            bs.append(np.sum((per.f.to_numpy()[k] - per.f0.to_numpy()[k]) * w) / w.sum())
        S.append(dict(archetype=A, bin=b, n=len(g), n_elong=len(el), frac_psi_gt45=float(np.mean(el.psi > 45)),
                      frac_psi_gt45_null=float(np.mean(el.psi_null > 45)), diff=float(np.mean(el.psi > 45) - np.mean(el.psi_null > 45)),
                      diff_lo=float(np.quantile(bs, .025)), diff_hi=float(np.quantile(bs, .975)), aspect_med=float(g.aspect.median())))
    S = pd.DataFrame(S); S.to_csv(OUT / "S1.4e_E5_summary.csv", index=False); print(S.round(3).to_string())

if __name__ == "__main__":
    main()
