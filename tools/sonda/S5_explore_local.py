"""tools/sonda/S5_explore_local.py - S5, esplorazione: residuo locale della registrazione globale X -> H2 (orient_<X>.json), come
S3b_explore_local.py: tessere da 1 mm a passo 500 um dell'H&E hires di H2 con >= 80 % di tessuto in entrambe, traslazione ECC (MOTION_TRANSLATION)
di X gia' riportata su H2, sfocature 2 e 6 px; controllo positivo: X spostata di 20 px hires in x. Uso: python S5_explore_local.py X1 [X2]"""
import os, sys, json, numpy as np, pandas as pd, cv2
from PIL import Image
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
for s in ("X1", "X2"): C.SEC[s] = f"{C.M}/data_lung/{s}/visium"
E = f"{C.M}/S5_sonda/explore"; RES = f"{C.REPO}/results/sonda/S5"; out = {}
crit = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 200, 1e-6)
def gray(s): return 255 - cv2.cvtColor(np.asarray(Image.open(f"{C.bdir(s, '008')}/spatial/tissue_hires_image.png"))[..., :3], cv2.COLOR_RGB2GRAY).astype(np.float32)
g2 = gray("H2"); sf2 = json.load(open(f"{C.bdir('H2', '008')}/spatial/scalefactors_json.json"))["tissue_hires_scalef"]; mpp = C.upp("H2") / sf2; T = int(round(1000 / mpp))
def ecc(a, b, sig):
    a = cv2.GaussianBlur(a, (0, 0), sig); b = cv2.GaussianBlur(b, (0, 0), sig); w = np.eye(2, 3, dtype=np.float32)
    try: cc, w = cv2.findTransformECC(a, b, w, cv2.MOTION_TRANSLATION, crit, None, 5); return float(cc), float(w[0, 2] * mpp), float(w[1, 2] * mpp)
    except cv2.error: return np.nan, np.nan, np.nan
for sec in sys.argv[1:]:
    I = json.load(open(f"{E}/orient_{sec}.json")); Wh = np.array(I["H2hires_to_Xhires_colrow"])
    g1 = cv2.warpAffine(gray(sec), Wh[:2], (g2.shape[1], g2.shape[0]), flags=cv2.INTER_LINEAR | cv2.WARP_INVERSE_MAP); rows = []
    for r0 in range(0, g2.shape[0] - T + 1, T // 2):
        for c0 in range(0, g2.shape[1] - T + 1, T // 2):
            a, b = g2[r0:r0 + T, c0:c0 + T], g1[r0:r0 + T, c0:c0 + T]
            if (a > 40).mean() < 0.8 or (b > 40).mean() < 0.8: continue
            row = dict(r0=r0, c0=c0)
            for sig in (2, 6):
                cc, dx, dy = ecc(a, b, sig); row.update({f"cc_{sig}": cc, f"dx_{sig}": dx, f"dy_{sig}": dy, f"norm_{sig}": float(np.hypot(dx, dy))})
            bs = np.roll(b, 20, axis=1); cc, dx, dy = ecc(a, bs, 2); row.update(pos_dx=dx, pos_dy=dy); rows.append(row)
    L = pd.DataFrame(rows); L.to_csv(f"{E}/local_ecc_{sec}.csv", index=False)
    s = {f"sigma{sig}": dict(n=int(L[f"norm_{sig}"].notna().sum()), median_um=float(L[f"norm_{sig}"].median()), p95_um=float(L[f"norm_{sig}"].quantile(0.95)), max_um=float(L[f"norm_{sig}"].max())) for sig in (2, 6)}
    s["positive_control_40um"] = dict(median_dx_minus_base_um=float((L.pos_dx - L.dx_2).median()), expected_abs_um=float(20 * mpp)); s["um_per_hires_px"] = mpp; out[sec] = s
json.dump(out, open(f"{RES}/local_ecc.json", "w"), indent=1); print(json.dumps(out, indent=0))
