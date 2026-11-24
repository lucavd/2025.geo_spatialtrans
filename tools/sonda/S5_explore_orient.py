"""tools/sonda/S5_explore_orient.py - S5, esplorazione prima della regola (nessuna espressione letta): orientamento delle griglie di X1 e X2
rispetto a H2 e posizione delle 10 aree di V nelle sezioni X. Come S3b_explore_overlap.py, con quattro rotazioni iniziali dell'ECC (0, 90, 180,
270 gradi) e quella con la correlazione piu' alta. Per ogni sezione X: affine array -> pixel, ECC euclidea poi affine (H2 fissa), residui locali
su tessere da 1 mm, rotazione della griglia di array di X rispetto a quella di H2, dove cadono le aree di V (meta' di sviluppo o tenuta da parte di X).
Uscita: /mnt/micron/geo_spatialtrans/S5_sonda/explore/orient_<X>.{png,json}, results/sonda/S5/orient.json"""
import os, sys, json, numpy as np, pandas as pd, cv2, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image; from skimage.registration import phase_cross_correlation; from shapely.geometry import Polygon, box
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
for s in ("X1", "X2"): C.SEC[s] = f"{C.M}/data_lung/{s}/visium"
OUT = f"{C.M}/S5_sonda/explore"; RES = f"{C.REPO}/results/sonda/S5"; os.makedirs(OUT, exist_ok=True); os.makedirs(RES, exist_ok=True)
def load(sec):
    tp = C.positions(sec, "008"); t = tp[tp.in_tissue == 1]
    A = np.c_[t.array_row, t.array_col, np.ones(len(t))]; Y = np.c_[t.pxl_row_in_fullres, t.pxl_col_in_fullres]
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    sf = json.load(open(f"{C.bdir(sec, '008')}/spatial/scalefactors_json.json"))["tissue_hires_scalef"]
    g = cv2.cvtColor(np.asarray(Image.open(f"{C.bdir(sec, '008')}/spatial/tissue_hires_image.png"))[..., :3], cv2.COLOR_RGB2GRAY)
    R, Cc = int(tp.array_row.max()) + 1, int(tp.array_col.max()) + 1; M = np.zeros((R, Cc), np.uint8); M[tp.array_row, tp.array_col] = tp.in_tissue
    return dict(coef=coef, sf=sf, g=g, M=M, m=C.split_m(sec), R=R, C=Cc, resid=float(np.abs(A @ coef - Y).max()))
def prep(g, shape=None):
    s = cv2.resize(g, (g.shape[1] // 2, g.shape[0] // 2), interpolation=cv2.INTER_AREA).astype(np.float32); s = cv2.GaussianBlur(255 - s, (0, 0), 2)
    if shape is None: return s
    o = np.zeros(shape, np.float32); h, w = min(shape[0], s.shape[0]), min(shape[1], s.shape[1]); o[:h, :w] = s[:h, :w]; return o
def um_to_hires(S, xy):
    r, c = xy[:, 1] / 8 - 0.5, xy[:, 0] / 8 - 0.5; P = np.c_[r, c, np.ones(len(r))] @ S["coef"]; return np.c_[P[:, 1], P[:, 0]] * S["sf"]
def hires_to_um(S, xy):
    A = S["coef"][:2].T; b = S["coef"][2]; rc = np.linalg.solve(A, (xy[:, ::-1] / S["sf"] - b).T).T; return np.c_[(rc[:, 1] + 0.5) * 8, (rc[:, 0] + 0.5) * 8]
H2 = load("H2"); f2 = prep(H2["g"]); Hh, Ww = f2.shape; crit = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 500, 1e-7)
V = C.areas("H2").reset_index(drop=True); summary = {}
for sec in sys.argv[1:]:
    S = load(sec); f1 = prep(S["g"], f2.shape); info = dict(array_to_fullres=S["coef"].T.tolist(), max_abs_resid_px=S["resid"], hires_scalef=S["sf"], m_um=S["m"])
    ctr = np.array([Ww / 2, Hh / 2]); best = None; tries = {}
    for ang in (0, 90, 180, 270):
        R0 = cv2.getRotationMatrix2D((float(ctr[0]), float(ctr[1])), ang, 1.0).astype(np.float32)
        try:
            cc, w = cv2.findTransformECC(f2, f1, R0.copy(), cv2.MOTION_EUCLIDEAN, crit, None, 5)
        except cv2.error as e:
            tries[ang] = dict(error=str(e)[:80]); continue
        tries[ang] = dict(cc=float(cc), angle_deg=float(np.degrees(np.arctan2(w[1, 0], w[0, 0]))))
        if best is None or cc > best[0]: best = (cc, w, ang)
    info["ecc_init_tries"] = tries
    cc_a, wa = cv2.findTransformECC(f2, f1, best[1].copy(), cv2.MOTION_AFFINE, crit, None, 5)
    info["ecc"] = dict(init_deg=best[2], euclid_cc=float(best[0]), affine_cc=float(cc_a), angle_deg=float(np.degrees(np.arctan2(wa[1, 0], wa[0, 0]))), warp_halfhires=wa.astype(float).tolist())
    Wh = np.r_[np.c_[wa[:, :2], wa[:, 2] * 2], [[0, 0, 1]]]; info["H2hires_to_Xhires_colrow"] = Wh.tolist()  # punto di H2 (hires) -> punto di X (hires)
    def h2_to_x(xy): return hires_to_um(S, (Wh @ np.c_[um_to_hires(H2, xy), np.ones(len(xy))].T).T[:, :2])
    # rotazione della griglia di array di X rispetto a quella di H2 (parte lineare della mappa µm H2 -> µm X)
    P0 = h2_to_x(np.array([[1000.0, 1000.0], [2000.0, 1000.0], [1000.0, 2000.0]])); L = np.c_[P0[1] - P0[0], P0[2] - P0[0]] / 1000.0
    info["array_map_H2_to_X"] = dict(linear=L.tolist(), det=float(np.linalg.det(L)), angle_deg=float(np.degrees(np.arctan2(L[1, 0], L[0, 0]))))
    # residui locali su tessere da 1 mm
    mpp2 = C.upp("H2") / H2["sf"]; tile = int(round(1000 / mpp2)); g2 = 255 - H2["g"].astype(np.float32)
    g1w = cv2.warpAffine(255 - S["g"].astype(np.float32), Wh[:2], (g2.shape[1], g2.shape[0]), flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP); loc = []
    for r0 in range(0, g2.shape[0] - tile, tile):
        for c0 in range(0, g2.shape[1] - tile, tile):
            a, b = g2[r0:r0 + tile, c0:c0 + tile], g1w[r0:r0 + tile, c0:c0 + tile]
            if (a > 40).mean() < 0.8 or (b > 40).mean() < 0.8: continue
            sh, err, _ = phase_cross_correlation(a, b, upsample_factor=10); loc.append(dict(r0=r0, c0=c0, dy_um=float(sh[0] * mpp2), dx_um=float(sh[1] * mpp2)))
    Lt = pd.DataFrame(loc); Lt["norm_um"] = np.hypot(Lt.dx_um, Lt.dy_um) if len(Lt) else []
    info["local_tiles"] = dict(n=len(Lt), median_um=float(Lt.norm_um.median()) if len(Lt) else None, p95_um=float(Lt.norm_um.quantile(0.95)) if len(Lt) else None, max_um=float(Lt.norm_um.max()) if len(Lt) else None)
    Lt.to_csv(f"{OUT}/local_tiles_{sec}.csv", index=False)
    # maschere e meta' nello spazio di array di H2
    R2, C2 = np.meshgrid(np.arange(H2["R"]), np.arange(H2["C"]), indexing="ij"); rr, cc = R2.ravel(), C2.ravel()
    xy = np.c_[(cc + 0.5) * 8, (rr + 0.5) * 8]; q = h2_to_x(xy); ci, ri = np.rint(q[:, 0] / 8 - 0.5).astype(int), np.rint(q[:, 1] / 8 - 0.5).astype(int)
    ok = (ri >= 0) & (ri < S["R"]) & (ci >= 0) & (ci < S["C"]); tX = np.zeros(len(rr), bool); tX[ok] = S["M"][ri[ok], ci[ok]] == 1
    t2 = H2["M"][rr, cc] == 1; d2 = xy[:, 0] < H2["m"] - 50; dX = q[:, 0] < S["m"] - 50; hX = q[:, 0] > S["m"] + 50; a8 = 64e-6
    info["overlap_mm2"] = dict(tissue_H2=float(t2.sum() * a8), tissue_both=float((t2 & tX).sum() * a8), mask_dice=float(2 * (t2 & tX).sum() / (t2.sum() + S["M"].sum())),
                               devH2_and_devX=float((t2 & tX & d2 & dX).sum() * a8), devH2_and_holdoutX=float((t2 & tX & d2 & hX).sum() * a8))
    rows = []
    for a in V.itertuples():
        sq = np.array([[a.x0, a.y0], [a.x0 + a.side, a.y0], [a.x0 + a.side, a.y0 + a.side], [a.x0, a.y0 + a.side]]); qq = h2_to_x(sq); cm = h2_to_x(np.array([[a.x0 + a.side / 2, a.y0 + a.side / 2]]))[0]
        x0 = 2 * np.round((cm[0] - a.side / 2) / 2); y0 = 2 * np.round((cm[1] - a.side / 2) / 2); tr = Polygon(qq); Tb = box(x0, y0, x0 + a.side, y0 + a.side)
        ins = (xy[:, 0] >= a.x0) & (xy[:, 0] < a.x0 + a.side) & (xy[:, 1] >= a.y0) & (xy[:, 1] < a.y0 + a.side) & t2
        rows.append(dict(section=sec, tissue=a.tissue, k=int(a.k), x0=float(x0), y0=float(y0), side=float(a.side), frac_true_covered=float(tr.intersection(Tb).area / tr.area),
                         tissue_in_X=float(tX[ins].mean()), xmax_um=float(x0 + a.side), xmin_um=float(x0), in_X_dev=bool(x0 + a.side < S["m"] - 50), in_X_holdout=bool(x0 > S["m"] + 50)))
    A = pd.DataFrame(rows); A.to_csv(f"{OUT}/areas_hom_{sec}.csv", index=False); info["areas"] = dict(n_in_dev=int(A.in_X_dev.sum()), n_in_holdout=int(A.in_X_holdout.sum()), n_straddle=int((~A.in_X_dev & ~A.in_X_holdout).sum()), min_tissue_in_X=float(A.tissue_in_X.min()))
    Z = np.zeros(len(rr)); Z[t2] = 1; Z[t2 & tX] = 2; Z[t2 & tX & d2 & dX] = 3; Z[t2 & tX & d2 & hX] = 4
    fig, ax = plt.subplots(1, 2, figsize=(14, 7)); ax[0].imshow(Z.reshape(R2.shape), cmap="tab10", vmin=0, vmax=9, extent=(0, H2["C"] * 8, H2["R"] * 8, 0))
    ax[0].set_title(f"H2 array space: 1 H2 only, 2 both, 3 dev H2 & dev {sec}, 4 dev H2 & holdout {sec}")
    for r_ in V.itertuples(): ax[0].add_patch(plt.Rectangle((r_.x0, r_.y0), r_.side, r_.side, fill=False, color="k"))
    ax[1].imshow(np.dstack([g2 / 255, g1w / 255, g2 / 255])); ax[1].set_title(f"H2 (magenta) vs {sec} warped (green)"); fig.savefig(f"{OUT}/orient_{sec}.png", dpi=80)
    json.dump(info, open(f"{OUT}/orient_{sec}.json", "w"), indent=1)
    summary[sec] = {k: info[k] for k in ("ecc_init_tries", "ecc", "array_map_H2_to_X", "local_tiles", "overlap_mm2", "areas")}; summary[sec]["ecc"].pop("warp_halfhires")
json.dump(summary, open(f"{RES}/orient.json", "w"), indent=1); print(json.dumps(summary, indent=0))
