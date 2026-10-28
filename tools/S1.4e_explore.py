#!/usr/bin/env python
"""tools/S1.4e_explore.py - S1.4e, esplorazione dichiarata PRIMA della pre-registrazione (docs/steps/S1.4e.md).
Solo dati veri e unita' gia' simulate (S1.5d, braccio RV = semina dei preset di oggi), nessuna simulazione nuova.
E1  coppie giudicate da Luca in R2c (kind == 'pair'): risposta (touch / separate / unsure) contro distanza fra i suoi punti (R2b).
E2  stesse coppie, nuclei veri dalla maschera del segmentatore di forma (SHAPE_SEG di tools/S1.4c_raster.py) agganciati ai punti
    (colonne label_<seg> di R2b_points.csv): aree, estensione di ciascun nucleo lungo l'asse della coppia, angolo fra asse maggiore
    (momenti dei pixel) e asse della coppia, distanza di Chebyshev fra le due maschere (contatto dell'operatore = <= 2 px).
E3  taglia contro distanza dal vicino: punti di Luca agganciati (area della maschera, distanza dal vicino fra i punti) contro nuclei
    simulati del raster RV di S1.5d nella stessa finestra (10 semi; area del raster, distanza dal vicino fra i centri del raster).
Uscite: results/S1.4e/explore/S1.4e_E*.csv.  Uso (root del repo): .venv/bin/python tools/S1.4e_explore.py"""
import json
import numpy as np, pandas as pd
from pathlib import Path
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
R2B = ROOT / "results/R2b"; R2C = ROOT / "results/R2c"; OUT = ROOT / "results/S1.4e/explore"
MASKS = Path("/mnt/micron/geo_spatialtrans/R2/masks"); RAS = Path("/mnt/micron/geo_spatialtrans/S1.5d/raster/RV")
SHAPE_SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}
BINS = [0, 5, 6, 7, 8, 10, np.inf]; LAB = ["<5", "5-6", "6-7", "7-8", "8-10", ">=10"]
SEEDS = [f"s{k:02d}" for k in range(1, 11)]
_mask = {}

def mask(A, roi):
    k = (A, roi)
    if k not in _mask:
        _mask.clear(); _mask[k] = np.load(MASKS / f"{A}_{roi}_{SHAPE_SEG[A]}_native.npz")["labels"]
    return _mask[k]

def dbin(d):
    return pd.cut(d, BINS, right=False, labels=LAB)

def major_axis(pts):
    c = np.cov(pts.T); w, v = np.linalg.eigh(c); return v[:, 1], float(np.sqrt(w[1] / max(w[0], 1e-12)))

def pair_geom(M, upp, li, lj, pi, pj, pad_um=14.0):
    """Geometria della coppia di etichette li, lj intorno ai punti pi, pj (um). None se un'etichetta manca nel ritaglio."""
    lo = np.floor((np.minimum(pi, pj) - pad_um) / upp).astype(int); hi = np.ceil((np.maximum(pi, pj) + pad_um) / upp).astype(int)
    c0, r0 = max(lo[0], 0), max(lo[1], 0); c1, r1 = min(hi[0], M.shape[1]), min(hi[1], M.shape[0])
    sub = M[r0:r1, c0:c1]; mi = sub == li; mj = sub == lj
    if not mi.any() or not mj.any():
        return None
    ri, ci_ = np.nonzero(mi); rj, cj_ = np.nonzero(mj)
    Pi = np.c_[(ci_ + c0) * upp, (ri + r0) * upp]; Pj = np.c_[(cj_ + c0) * upp, (rj + r0) * upp]
    gi, gj = Pi.mean(0), Pj.mean(0); d = float(np.hypot(*(gj - gi))); u = (gj - gi) / d
    ei = float(((Pi - gi) @ u).max() + upp / 2); ej = float(((Pj - gj) @ (-u)).max() + upp / 2)
    vi, ari = major_axis(Pi); vj, arj = major_axis(Pj)
    psi = lambda v: float(np.degrees(np.arccos(min(1.0, abs(v @ u)))))     # 0 = asse maggiore lungo la coppia, 90 = di traverso
    cheb = ndi.distance_transform_cdt(~mi, metric="chessboard")[mj].min()
    ai, aj = mi.sum() * upp ** 2, mj.sum() * upp ** 2
    return dict(d_seg=d, ext_i=ei, ext_j=ej, gap_axis=d - ei - ej, psi_i=psi(vi), psi_j=psi(vj), aspect_i=ari, aspect_j=arj,
                area_i=ai, area_j=aj, rsum_circle=float(np.sqrt(ai / np.pi) + np.sqrt(aj / np.pi)), cheb_px=int(cheb),
                contact_op=bool(cheb <= 2))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wins = pd.read_csv(R2B / "R2b_windows.csv").set_index("win_id")
    pts = pd.read_csv(R2B / "R2b_points.csv"); pts = pts[pts.rater == "Luca"].set_index(["win_id", "idx"])
    it = pd.read_csv(R2C / "R2c_items_Luca.csv"); a = pd.DataFrame(json.load(open(R2C / "annotations/R2c_answers_Luca.json"))["answers"])
    J = it.merge(a, on="item", validate="1:1"); J = J[J.kind == "pair"].copy(); J["archetype"] = J.win_id.str[:2]
    rows = []
    for r in J.sort_values(["win_id"]).itertuples():
        w = wins.loc[r.win_id]; A = r.archetype; seg = SHAPE_SEG[A]; upp = float(w.um_per_px)
        p_i, p_j = pts.loc[(r.win_id, int(r.i))], pts.loc[(r.win_id, int(r.j))]
        P_i = np.array([p_i.x_um, p_i.y_um]); P_j = np.array([p_j.x_um, p_j.y_um])
        base = dict(archetype=A, win_id=r.win_id, i=int(r.i), j=int(r.j), d_manual=float(np.hypot(*(P_j - P_i))), d_um_items=float(r.d_um),
                    answer=r.answer, seg=seg, matched_i=bool(p_i[f"match_{seg}"]), matched_j=bool(p_j[f"match_{seg}"]))
        li, lj = p_i[f"label_{seg}"], p_j[f"label_{seg}"]
        base["status"] = "unmatched"
        if base["matched_i"] and base["matched_j"]:
            li, lj = int(li), int(lj)
            if li == lj:
                base["status"] = "same_label"
            else:
                g = pair_geom(mask(A, w.roi_id), upp, li, lj, P_i, P_j)
                if g is not None:
                    base.update(g); base["status"] = "ok"
        rows.append(base)
    E = pd.DataFrame(rows); E["bin"] = dbin(E.d_manual); E.to_csv(OUT / "S1.4e_E12_pairs.csv", index=False)
    print("d_manual vs d_um di R2c_items, max |diff|:", float((E.d_manual - E.d_um_items).abs().max()))
    S = []
    for (A, b), g in E.groupby(["archetype", "bin"], observed=True):
        jd = g[g.answer != "unsure"]; ok = g[g.status == "ok"]; okj = ok[ok.answer != "unsure"]
        S.append(dict(archetype=A, bin=b, n_pairs=len(g), n_judged=len(jd), p_touch=float((jd.answer == "touch").mean()) if len(jd) else np.nan,
                      p_unsure=float((g.answer == "unsure").mean()), n_ok=len(ok), p_same_label=float((g.status == "same_label").mean()),
                      p_contact_op=float(ok.contact_op.mean()) if len(ok) else np.nan,
                      p_touch_ok=float((okj.answer == "touch").mean()) if len(okj) else np.nan,
                      sens_op=float(okj[okj.answer == "touch"].contact_op.mean()) if (okj.answer == "touch").any() else np.nan,
                      spec_op=float(1 - okj[okj.answer == "separate"].contact_op.mean()) if (okj.answer == "separate").any() else np.nan,
                      psi_med=float(np.median(np.r_[ok.psi_i, ok.psi_j])) if len(ok) else np.nan,
                      frac_psi_gt45=float(np.mean(np.r_[ok.psi_i, ok.psi_j] > 45)) if len(ok) else np.nan,
                      ext_over_d=float(np.median((ok.ext_i + ok.ext_j) / ok.d_seg)) if len(ok) else np.nan,
                      rcirc_over_d=float(np.median(ok.rsum_circle / ok.d_seg)) if len(ok) else np.nan,
                      area_med=float(np.median(np.r_[ok.area_i, ok.area_j])) if len(ok) else np.nan))
    pd.DataFrame(S).to_csv(OUT / "S1.4e_E12_summary.csv", index=False)
    # E3: taglia contro distanza dal vicino, veri (punti di Luca agganciati) e simulati (RV), stesse finestre
    allp = pd.read_csv(R2B / "R2b_points.csv"); allp = allp[allp.rater == "Luca"]
    T = []
    for wid, g in allp.groupby("win_id"):
        w = wins.loc[wid]; A = wid[:2]; seg = SHAPE_SEG[A]; upp = float(w.um_per_px); s = float(w.side_um)
        x0, y0 = float(w.x0_um), float(w.y0_um); P = g[["x_um", "y_um"]].to_numpy()
        nn = cKDTree(P).query(P, k=2)[0][:, 1]; de = np.minimum.reduce([P[:, 0] - x0, x0 + s - P[:, 0], P[:, 1] - y0, y0 + s - P[:, 1]])
        M = mask(A, w.roi_id); lab = g[f"label_{seg}"].to_numpy(); mt = g[f"match_{seg}"].to_numpy().astype(bool)
        cnt = np.bincount(M.ravel()) if True else None
        for k in np.nonzero(mt & (de >= 10))[0]:
            T.append(dict(source="manual_seg", archetype=A, win_id=wid, seed="", nn_um=nn[k], area=cnt[int(lab[k])] * upp ** 2))
        for sd in SEEDS:
            f = RAS / f"{A}_{w.roi_id}_{sd}.parquet"
            t = pd.read_parquet(f, columns=["cx", "cy", "area"]); Q = t[["cx", "cy"]].to_numpy()
            inw = (Q[:, 0] >= x0 + 10) & (Q[:, 0] < x0 + s - 10) & (Q[:, 1] >= y0 + 10) & (Q[:, 1] < y0 + s - 10)
            nq = cKDTree(Q).query(Q[inw], k=2)[0][:, 1]
            T += [dict(source="RV", archetype=A, win_id=wid, seed=sd, nn_um=a_, area=b_) for a_, b_ in zip(nq, t.area.to_numpy()[inw])]
    T = pd.DataFrame(T); T["bin"] = dbin(T.nn_um); T.to_csv(OUT / "S1.4e_E3_size_nn.csv", index=False)
    S3 = T.groupby(["archetype", "source", "bin"], observed=True).agg(n=("area", "size"), area_med=("area", "median"),
         area_q25=("area", lambda x: x.quantile(.25)), area_q75=("area", lambda x: x.quantile(.75))).reset_index()
    S3.to_csv(OUT / "S1.4e_E3_summary.csv", index=False)
    print("fatto:", len(E), "coppie,", len(T), "nuclei (E3)")

if __name__ == "__main__":
    main()
