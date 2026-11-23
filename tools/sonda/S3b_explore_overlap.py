"""tools/sonda/S3b_explore_overlap.py - S3b, esplorazione prima della regola (nessuna espressione letta):
1. affine array (riga, colonna, quadratini da 8 um) -> pixel di piena risoluzione per sezione (minimi quadrati sui quadratini in tissue);
2. registrazione dell'H&E hires di H1 su quella di H2 (grigio, ECC euclidea di OpenCV, poi affine) e residui locali su tessere da 1 mm
   (correlazione di fase di skimage dopo la trasformazione globale);
3. metà di sviluppo di H1 (x < m - 50 um) e di H2 portate nello spazio array di H2: area di tessuto in ciascuna e nell'intersezione.
Uscita: /mnt/micron/geo_spatialtrans/S3b_sonda/explore/overlap.{png,json}"""
import os, sys, json, numpy as np, pandas as pd, cv2, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image; from skimage.registration import phase_cross_correlation
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
OUT = "/mnt/micron/geo_spatialtrans/S3b_sonda/explore"; os.makedirs(OUT, exist_ok=True); info = {}
S = {}
for sec in ("H1", "H2"):
    tp = C.positions(sec, "008"); t = tp[tp.in_tissue == 1]
    A = np.c_[t.array_row, t.array_col, np.ones(len(t))]; Y = np.c_[t.pxl_row_in_fullres, t.pxl_col_in_fullres]
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None); res = np.abs(A @ coef - Y).max()
    sf = json.load(open(f"{C.bdir(sec, '008')}/spatial/scalefactors_json.json"))["tissue_hires_scalef"]
    g = cv2.cvtColor(np.asarray(Image.open(f"{C.bdir(sec, '008')}/spatial/tissue_hires_image.png"))[..., :3], cv2.COLOR_RGB2GRAY)
    R, Cc = int(tp.array_row.max()) + 1, int(tp.array_col.max()) + 1
    M = np.zeros((R, Cc), np.uint8); M[tp.array_row, tp.array_col] = tp.in_tissue
    S[sec] = dict(coef=coef, sf=sf, g=g, M=M, m=C.split_m(sec), R=R, C=Cc)
    info[sec] = dict(array_to_fullres=coef.T.tolist(), max_abs_resid_px=float(res), hires_scalef=sf, m_um=S[sec]["m"])
# ECC: H1 (moving) -> H2 (fixed), immagini a meta' risoluzione hires, grigio invertito (tessuto chiaro su fondo scuro)
def prep(g):
    s = cv2.resize(g, (g.shape[1] // 2, g.shape[0] // 2), interpolation=cv2.INTER_AREA).astype(np.float32); return cv2.GaussianBlur(255 - s, (0, 0), 2)
f2, f1 = prep(S["H2"]["g"]), prep(S["H1"]["g"]); H, W = f2.shape
f1p = np.zeros_like(f2); f1p[:min(H, f1.shape[0]), :min(W, f1.shape[1])] = f1[:H, :W]
crit = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 500, 1e-7)
out = {}
for name, mode in (("euclidean", cv2.MOTION_EUCLIDEAN), ("affine", cv2.MOTION_AFFINE)):
    w = np.eye(2, 3, dtype=np.float32); cc, w = cv2.findTransformECC(f2, f1p, w, mode, crit, None, 5); out[name] = (float(cc), w.astype(float))
info["ecc"] = {k: dict(cc=v[0], warp_halfhires=v[1].tolist(), angle_deg=float(np.degrees(np.arctan2(v[1][1, 0], v[1][0, 0])))) for k, v in out.items()}
w = out["affine"][1]
# warp in coordinate hires piene: x2 = w x1 (stessa scala per entrambe a meta' risoluzione) -> hires piena: traslazione * 2
Wh = np.r_[np.c_[w[:, :2], w[:, 2] * 2], [[0, 0, 1]]]  # hires H2 (x=col, y=row) = Wh @ hires H1? findTransformECC: template f2, input f1 warped by w -> f1(w x) ~ f2(x)
# quindi un punto x di H2 (hires, col,row) corrisponde a w x in H1. Wh: H2 -> H1.
info["H2hires_to_H1hires_colrow"] = Wh.tolist()
# residui locali: H1 riportata su H2 con la trasformazione globale, tessere da 1 mm (in hires H2), solo tessere con >= 80% tessuto
mpp2 = C.upp("H2") / S["H2"]["sf"]; tile = int(round(1000 / mpp2)); g2 = 255 - S["H2"]["g"].astype(np.float32)
g1w = cv2.warpAffine(255 - S["H1"]["g"].astype(np.float32), Wh[:2], (g2.shape[1], g2.shape[0]), flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP)
loc = []
for r0 in range(0, g2.shape[0] - tile, tile):
    for c0 in range(0, g2.shape[1] - tile, tile):
        a, b = g2[r0:r0 + tile, c0:c0 + tile], g1w[r0:r0 + tile, c0:c0 + tile]
        if (a > 40).mean() < 0.8 or (b > 40).mean() < 0.8: continue
        sh, err, _ = phase_cross_correlation(a, b, upsample_factor=10); loc.append(dict(r0=r0, c0=c0, dy_um=float(sh[0] * mpp2), dx_um=float(sh[1] * mpp2)))
L = pd.DataFrame(loc); L["norm_um"] = np.hypot(L.dx_um, L.dy_um)
info["local_tiles"] = dict(n=len(L), tile_um=1000, um_per_hires_px=mpp2, median_um=float(L.norm_um.median()), p95_um=float(L.norm_um.quantile(0.95)), max_um=float(L.norm_um.max()))
L.to_csv(f"{OUT}/local_tiles.csv", index=False)
# spazio array H2 (8 um): per ogni quadratino di H2 la posizione in H1 (array) -> in tissue in H1? in metà di sviluppo di H1?
def arr_to_hires(sec, rr, cc):
    co = S[sec]["coef"]; P = np.c_[rr, cc, np.ones(len(rr))] @ co; return P[:, 1] * S[sec]["sf"], P[:, 0] * S[sec]["sf"]  # col, row hires
def hires_to_arr(sec, xc, yr):
    co = S[sec]["coef"]; A2 = co[:2].T; b = co[2]; F = np.c_[yr, xc] / S[sec]["sf"] - b; rc = np.linalg.solve(A2, F.T).T; return rc[:, 0], rc[:, 1]
R2, C2 = np.meshgrid(np.arange(S["H2"]["R"]), np.arange(S["H2"]["C"]), indexing="ij"); rr, cc = R2.ravel(), C2.ravel()
xh, yh = arr_to_hires("H2", rr, cc); P1 = Wh @ np.c_[xh, yh, np.ones(len(xh))].T; r1, c1 = hires_to_arr("H1", P1[0], P1[1])
ri, ci = np.rint(r1).astype(int), np.rint(c1).astype(int); ok = (ri >= 0) & (ri < S["H1"]["R"]) & (ci >= 0) & (ci < S["H1"]["C"])
t1 = np.zeros(len(rr), bool); t1[ok] = S["H1"]["M"][ri[ok], ci[ok]] == 1
t2 = S["H2"]["M"][rr, cc] == 1
d2 = (cc + 0.5) * 8 < S["H2"]["m"] - 50; d1 = (c1 + 0.5) * 8 < S["H1"]["m"] - 50
a = 64e-6  # mm^2 per quadratino
info["overlap_in_H2_array_space_mm2"] = dict(tissue_H2=float(t2.sum() * a), tissue_both=float((t2 & t1).sum() * a), dev_H2=float((t2 & d2).sum() * a),
    dev_H1_mapped=float((t2 & t1 & d1).sum() * a), dev_both=float((t2 & t1 & d1 & d2).sum() * a), holdout_both=float((t2 & t1 & ~d1 & ~d2).sum() * a),
    mask_dice=float(2 * (t1 & t2).sum() / (t1.sum() + t2.sum())))
Z = np.zeros(len(rr)); Z[t2] = 1; Z[t2 & d2] = 2; Z[t2 & t1 & d1] = 3; Z[t2 & t1 & d1 & d2] = 4
fig, ax = plt.subplots(1, 2, figsize=(14, 7))
ax[0].imshow(Z.reshape(R2.shape), cmap="tab10", vmin=0, vmax=9, extent=(0, S["H2"]["C"] * 8, S["H2"]["R"] * 8, 0))
ax[0].set_title("H2 array space: 1 tissue, 2 dev H2 only, 3 dev H1 only, 4 dev both")
for r_ in C.areas("H2").itertuples(): ax[0].add_patch(plt.Rectangle((r_.x0, r_.y0), r_.side, r_.side, fill=False, color="k"))
ax[1].imshow(np.dstack([g2 / 255, g1w / 255, g2 / 255])); ax[1].set_title("H2 (magenta) vs H1 warped (green), hires")
fig.savefig(f"{OUT}/overlap.png", dpi=80); json.dump(info, open(f"{OUT}/overlap.json", "w"), indent=1)
print(json.dumps({k: info[k] for k in ("ecc", "local_tiles", "overlap_in_H2_array_space_mm2")}, indent=0)); print({s: info[s]["max_abs_resid_px"] for s in ("H1", "H2")})
