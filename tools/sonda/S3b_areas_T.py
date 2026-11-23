"""tools/sonda/S3b_areas_T.py - S3b (docs/sonda/REGOLA_S3b.md): aree di T omologhe alle 10 aree di V, cancelli G-R, G-O, G-C, controprova CP-2.
Catena di coordinate: µm di array di H2 -> pixel di piena risoluzione (affine per minimi quadrati sui quadratini da 8 µm in tissue) -> hires di H2
-> hires di H1 (registrazione globale di results/sonda/S3b/overlap.json, invariata) -> µm di array di H1. Area di T = quadrato di 750 µm allineato
agli assi di H1, centrato nel centro riportato, angolo arrotondato a 2 µm. G-O: frazione del quadrato vero riportato (4 angoli) coperta dal quadrato
di T (shapely). G-C: copertura in_tissue (quadratini da 8 µm) del quadrato di T. CP-2: punteggio r di REGOLA_S3 (matrice filtrata da 8 µm).
Uscite: results/sonda/S3b/areas_H1_hom.csv, gates_areas.json."""
import os, sys, json, numpy as np, pandas as pd
from scipy.stats import spearmanr; from shapely.geometry import Polygon, box
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
OUT = f"{C.REPO}/results/sonda/S3b"; I = json.load(open(f"{OUT}/overlap.json")); Wh = np.array(I["H2hires_to_H1hires_colrow"]); LE = json.load(open(f"{OUT}/local_ecc.json"))
co = {s: np.array(I[s]["array_to_fullres"]).T for s in ("H1", "H2")}; sf = {s: I[s]["hires_scalef"] for s in ("H1", "H2")}
E = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"]; S = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A3"]
def um_to_hires(sec, xy):  # xy: (x, y) µm di array -> (col, row) hires
    r, c = xy[:, 1] / 8 - 0.5, xy[:, 0] / 8 - 0.5; P = np.c_[r, c, np.ones(len(r))] @ co[sec]; return np.c_[P[:, 1], P[:, 0]] * sf[sec]
def hires_to_um(sec, xy):
    A = co[sec][:2].T; b = co[sec][2]; rc = np.linalg.solve(A, (xy[:, ::-1] / sf[sec] - b).T).T; return np.c_[(rc[:, 1] + 0.5) * 8, (rc[:, 0] + 0.5) * 8]
def h2_to_h1(xy): return hires_to_um("H1", (Wh @ np.c_[um_to_hires("H2", xy), np.ones(len(xy))].T).T[:, :2])
def score(sec, sq):
    tp = C.positions(sec, "008"); h5 = f"{C.bdir(sec, '008')}/filtered_feature_bc_matrix.h5"; ids, names = C.features(h5); bc = C.barcodes(h5)
    out = []
    for x0, y0, side in sq:
        x = (tp.array_col.to_numpy() + 0.5) * 8; y = (tp.array_row.to_numpy() + 0.5) * 8
        inside = (x >= x0) & (x < x0 + side) & (y >= y0) & (y < y0 + side); n = int(round(side / 8)) ** 2
        cov = float(tp.in_tissue.to_numpy()[inside].sum() / n); b = tp[inside & (tp.in_tissue.to_numpy() == 1)]
        if sec == "H2": C.guard((b.array_col.to_numpy() + 0.5) * 8, C.split_m("H2"))
        col = pd.Series(np.arange(len(bc)), index=bc).reindex(b.barcode.to_numpy()).to_numpy(); X, _ = C.read_cols(h5, col.astype(np.int64))
        nm = pd.Series(np.arange(len(names)), index=names); gE = [g for g in E if g in nm.index]; gS = [g for g in S if g in nm.index]
        t = X.sum(); e = X[:, nm[gE].to_numpy()].sum() / t; s = X[:, nm[gS].to_numpy()].sum() / t
        out.append(dict(coverage=cov, e=float(e), s=float(s), r=float(np.log2((e + 1e-4) / (s + 1e-4)))))
    return pd.DataFrame(out)
V = C.areas("H2").reset_index(drop=True); rows = []
for a in V.itertuples():
    sq = np.array([[a.x0, a.y0], [a.x0 + a.side, a.y0], [a.x0 + a.side, a.y0 + a.side], [a.x0, a.y0 + a.side]]); ctr = np.array([[a.x0 + a.side / 2, a.y0 + a.side / 2]])
    q = h2_to_h1(sq); cm = h2_to_h1(ctr)[0]; x0 = 2 * np.round((cm[0] - a.side / 2) / 2); y0 = 2 * np.round((cm[1] - a.side / 2) / 2)
    true = Polygon(q); T = box(x0, y0, x0 + a.side, y0 + a.side)
    rows.append(dict(section="H1", tissue=a.tissue, k=int(a.k), x0=float(x0), y0=float(y0), side=float(a.side), cx_H1=float(cm[0]), cy_H1=float(cm[1]),
                     frac_true_covered=float(true.intersection(T).area / true.area), true_area_ratio=float(true.area / a.side ** 2), V_x0=a.x0, V_y0=a.y0, V_r=a.r))
A = pd.DataFrame(rows)
sc = score("H1", A[["x0", "y0", "side"]].to_numpy()); A = pd.concat([A, sc.add_prefix("")], axis=1)
Vre = score("H2", V[["x0", "y0", "side"]].to_numpy()); A["V_r_recomputed"] = Vre.r.to_numpy()
To = C.areas("H1").sort_values(["tissue", "k"]); A["m_um"] = C.split_m("H1"); A["xmax_area_um"] = A.x0 + A.side; A["in_old_H1_dev_half"] = A.xmax_area_um < A.m_um - 50
A["flag_low_coverage"] = A.coverage < 0.9; A.to_csv(f"{OUT}/areas_H1_hom.csv", index=False)
g = {}
g["G-R"] = dict(p95_um=LE["sigma2"]["p95_um"], pos_ctrl_um=LE["positive_control_40um"]["median_dx_minus_base_um"], pass_=bool(LE["sigma2"]["p95_um"] <= 75 and abs(LE["positive_control_40um"]["median_dx_minus_base_um"] - 40) <= 5))
g["G-O"] = dict(min_frac=float(A.frac_true_covered.min()), pass_=bool((A.frac_true_covered >= 0.95).all()))
g["G-C"] = dict(min_cov=float(A.coverage.min()), n_below_0_9=int(A.flag_low_coverage.sum()), pass_=bool((A.coverage >= 0.8).all()))
g["r_identity_V"] = dict(max_abs_dev=float(np.abs(A.V_r_recomputed - A.V_r).max()))
rho_hom = spearmanr(A.V_r, A.r).correlation; sign_hom = int(((A.r > 0) == (A.V_r > 0)).sum())
Vs = V.sort_values(["tissue", "k"]).reset_index(drop=True); Ts = To.reset_index(drop=True); rho_old = spearmanr(Vs.r, Ts.r).correlation
g["CP-2"] = dict(spearman_hom=float(rho_hom), sign_agree_hom=sign_hom, spearman_S3_T_vs_V_same_order=float(rho_old), pass_=bool(rho_hom >= 0.7 and sign_hom >= 8))
json.dump(g, open(f"{OUT}/gates_areas.json", "w"), indent=1)
print(A[["tissue", "k", "x0", "y0", "frac_true_covered", "coverage", "V_r", "r", "in_old_H1_dev_half"]].round(3).to_string(index=False)); print(json.dumps(g, indent=0))
