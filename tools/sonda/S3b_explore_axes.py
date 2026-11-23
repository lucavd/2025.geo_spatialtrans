"""tools/sonda/S3b_explore_axes.py - S3b, controprova della complementarità delle metà di sviluppo (nessuna espressione letta).
(a) direzione dell'asse x di array (colonne) di H1 e H2 nell'immagine hires, dalla sola affine array -> pixel (senza registrazione);
(b) i centri delle 20 aree di S3 portati nell'altra sezione con la registrazione globale (overlap.json): in quale metà cadono, e a che distanza dalla linea m - 50 um."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
I = json.load(open("/mnt/micron/geo_spatialtrans/S3b_sonda/explore/overlap.json")); Wh = np.array(I["H2hires_to_H1hires_colrow"])
co = {s: np.array(I[s]["array_to_fullres"]).T for s in ("H1", "H2")}; sf = {s: I[s]["hires_scalef"] for s in ("H1", "H2")}; m = {s: I[s]["m_um"] for s in ("H1", "H2")}
for s in ("H1", "H2"):
    d = co[s][1] * sf[s]  # d(pixel row, pixel col)/d(array_col) in hires
    print(s, "array_col axis in hires (d row, d col) per 8 um:", np.round(d, 3), " array_row axis:", np.round(co[s][0] * sf[s], 3))
def a2h(s, r, c): p = np.array([r, c, 1.0]) @ co[s]; return np.array([p[1], p[0]]) * sf[s]          # (col, row) hires
def h2a(s, xy): A = co[s][:2].T; b = co[s][2]; rc = np.linalg.solve(A, np.array([xy[1], xy[0]]) / sf[s] - b); return rc
Winv = np.linalg.inv(Wh)
for s, o in (("H2", "H1"), ("H1", "H2")):
    for a in C.areas(s).itertuples():
        cx, cy = a.x0 + a.side / 2, a.y0 + a.side / 2; r, c = cy / 8 - 0.5, cx / 8 - 0.5
        xy = a2h(s, r, c); xy2 = (Wh if s == "H2" else Winv) @ np.r_[xy, 1]; ro, co_ = h2a(o, xy2[:2]); xo = (co_ + 0.5) * 8
        print(f"{s} {a.tissue} {a.k}: centre x {cx:.0f} -> {o} x {xo:.0f}  ({'dev' if xo < m[o] - 50 else 'HOLDOUT'}, margin {m[o] - 50 - xo:+.0f} um)")
