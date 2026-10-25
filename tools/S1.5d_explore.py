"""tools/S1.5d_explore.py - S1.5d, esplorazione dichiarata PRIMA della pre-registrazione (docs/steps/S1.5d.md).
Domande: (E1) quanto spazio vuoto c'e' dentro la maschera di tessuto, reale contro simulato (bracci U e N-a di S1.5c)?
(E2) gli spazi vuoti reali sono tessuto acellulare o nuclei persi dal segmentatore (colore H&E dei pixel lontani dai nuclei)?
(E3) nelle finestre contate a mano da Luca (R2b), i vuoti dei punti manuali sono come quelli dei nuclei segmentati?
(E4) la densita' locale varia piu' del Poisson, e piu' del simulato, a 25, 50, 100 um (CV^2 in eccesso, invariante al diradamento indipendente)?
Ingressi in sola lettura: maschere R2 (results/R2/tissue_masks), nuclei veri con territorio (R6, come S1.5b/S1.5c), ROI H&E (R2), annotazioni R2b,
celle simulate di S1.5c (/mnt/micron/geo_spatialtrans/S1.5c/cells). Uscite: results/S1.5d/S1.5d_explore_*.csv.
Uso: .venv/bin/python tools/S1.5d_explore.py [n_proc]"""
import sys, json, glob, importlib.util, numpy as np, pandas as pd
from multiprocessing import Pool
from PIL import Image
from scipy.spatial import cKDTree
sys.path.insert(0, "tools")
import R2_common as C
from skimage import color
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py")
OUT5C = "/mnt/micron/geo_spatialtrans/S1.5c/cells"; RES = "results/S1.5d"
ARMS = ["U", "N-a"]; SEEDS = range(1, 11)
U_REL = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]; D_ABS = [10.0, 15.0, 20.0, 30.0]; QSIDE = [25.0, 50.0, 100.0]
FRAME_B = 30.0          # pixel e quadrati ad almeno 30 um dal bordo del ROI (i centri fuori dal ROI non esistono)
ROIS = pd.read_csv("results/R2/R2_rois_checked.csv")

def mask_of(A, r):
    m = np.array(Image.open(f"results/R2/tissue_masks/{A}_{r}_valid_ds4.png"))
    if m.ndim == 3: m = m[..., 0]
    m = m > 127; row = ROIS[(ROIS.archetype == A) & (ROIS.roi_id == r)].iloc[0]
    return m, row.side_px * row.um_per_px / m.shape[1], row.side_px * row.um_per_px

def pix_xy(m, px, side):
    rr, cc = np.nonzero(m); x = (cc + 0.5) * px; y = (rr + 0.5) * px
    k = (x > FRAME_B) & (y > FRAME_B) & (x < side - FRAME_B) & (y < side - FRAME_B)
    return np.c_[x[k], y[k]], rr[k], cc[k]

def empty_space(P, Q, interior):
    tr = cKDTree(P); dd, ii = tr.query(P, k=2); m_all = np.median(dd[interior, 1]) if interior is not None else np.median(dd[:, 1])
    d, idx = tr.query(Q)
    o = {"nn_median": m_all, "n_centers": len(P)}
    for u in U_REL: o[f"H_rel_{u}"] = np.mean(d / m_all > u)
    for a in D_ABS: o[f"H_abs_{int(a)}"] = np.mean(d > a)
    o["q99_rel"] = np.quantile(d / m_all, 0.99)
    return o, d / m_all

def quadrats(P, m, px, side):
    o = {}
    for q in QSIDE:
        n = int((side - 2 * FRAME_B) // q); cnt = []
        for i in range(n):
            for j in range(n):
                x0 = FRAME_B + i * q; y0 = FRAME_B + j * q
                c0, c1, r0, r1 = int(x0 / px), int((x0 + q) / px), int(y0 / px), int((y0 + q) / px)
                if m[r0:r1, c0:c1].mean() < 0.99: continue
                cnt.append(np.sum((P[:, 0] >= x0) & (P[:, 0] < x0 + q) & (P[:, 1] >= y0) & (P[:, 1] < y0 + q)))
        cnt = np.asarray(cnt, float); o[f"nq_{int(q)}"] = len(cnt)
        o[f"xcv2_{int(q)}"] = (cnt.var(ddof=1) - cnt.mean()) / cnt.mean() ** 2 if len(cnt) >= 5 and cnt.mean() > 0 else np.nan
    return o

def he_classes(A, r, m, px, Q, rr, cc, drel):
    rgb = C.load_roi(A, r); f = 4
    h, w = (rgb.shape[0] // f) * f, (rgb.shape[1] // f) * f
    hed = color.rgb2hed(rgb[:h, :w]); H = hed[..., 0].reshape(h // f, f, w // f, f).mean(axis=(1, 3))
    E = hed[..., 1].reshape(h // f, f, w // f, f).mean(axis=(1, 3))
    V = color.rgb2hsv(rgb[:h, :w])[..., 2].reshape(h // f, f, w // f, f).mean(axis=(1, 3))
    ok = (rr < H.shape[0]) & (cc < H.shape[1]); hv = H[rr[ok], cc[ok]]; ev = E[rr[ok], cc[ok]]; vv = V[rr[ok], cc[ok]]; dv = drel[ok]
    rows = []
    for lo, hi in [(0, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, np.inf)]:
        k = (dv >= lo) & (dv < hi)
        rows.append(dict(d_lo=lo, d_hi=hi, n_px=int(k.sum()), frac_px=k.mean(), H_mean=hv[k].mean() if k.any() else np.nan,
                         H_q90=np.quantile(hv[k], 0.9) if k.any() else np.nan, E_mean=ev[k].mean() if k.any() else np.nan,
                         V_mean=vv[k].mean() if k.any() else np.nan))
    return rows

def one_roi(args):
    A, r = args
    m, px, side = mask_of(A, r); Q, rr, cc = pix_xy(m, px, side)
    _, d = cal.load(A); e = d[d.roi_id == r]
    P = np.c_[e.x.values, e.y.values]; I = e.interior.values.astype(bool)
    es, drel = empty_space(P, Q, I)
    rows_es = [dict(archetype=A, roi_id=r, source="real", seed=0, **es, **quadrats(P, m, px, side))]
    for arm in ARMS:
        for s in SEEDS:
            cdf = pd.read_parquet(f"{OUT5C}/{arm}/{A}_{r}_s{s:02d}.parquet")
            Ps = np.c_[cdf.x.values, cdf.y.values]; o, _ = empty_space(Ps, Q, cdf.interior.values.astype(bool))
            rows_es.append(dict(archetype=A, roi_id=r, source=arm, seed=s, **o, **quadrats(Ps, m, px, side)))
    he = [dict(archetype=A, roi_id=r, **x) for x in he_classes(A, r, m, px, Q, rr, cc, drel)]
    return rows_es, he

def manual_windows():
    j = json.load(open("results/R2b/annotations/R2b_annotations_Luca.json"))["windows"]
    rows = []
    for w in j:
        A, r, upp = w["archetype"], w["roi_id"], w["um_per_px"]; side = w["side_px"] * upp; mg = w["margin_px"]
        pts = np.array(w["points"], float).reshape(-1, 2)
        Pm = np.c_[w["x0_um"] + (pts[:, 0] - mg) * upp, w["y0_um"] + (pts[:, 1] - mg) * upp]
        _, d = cal.load(A); e = d[d.roi_id == r]; Ps = np.c_[e.x.values, e.y.values]
        x0, y0 = w["x0_um"], w["y0_um"]; mgu = mg * upp
        near = (Ps[:, 0] > x0 - mgu) & (Ps[:, 0] < x0 + side + mgu) & (Ps[:, 1] > y0 - mgu) & (Ps[:, 1] < y0 + side + mgu)
        Ps = Ps[near]
        g = np.arange(0.25, side, 0.5); X, Y = np.meshgrid(x0 + g, y0 + g); Q = np.c_[X.ravel(), Y.ravel()]
        m, px, _ = mask_of(A, r); ci = np.minimum((Q[:, 0] / px).astype(int), m.shape[1] - 1); ri = np.minimum((Q[:, 1] / px).astype(int), m.shape[0] - 1)
        Q = Q[m[ri, ci]]
        for lab, P in [("manual_Luca", Pm), ("segmented", Ps)]:
            inwin = (P[:, 0] >= x0) & (P[:, 0] < x0 + side) & (P[:, 1] >= y0) & (P[:, 1] < y0 + side)
            tr = cKDTree(P); dd, _ = tr.query(P[inwin], k=2); m_all = np.median(dd[:, 1])
            dq, _ = tr.query(Q); edge = np.minimum.reduce([Q[:, 0] - (x0 - mgu), x0 + side + mgu - Q[:, 0], Q[:, 1] - (y0 - mgu), y0 + side + mgu - Q[:, 1]])
            k = dq < edge                                  # correzione di bordo: il vicino piu' vicino non puo' stare fuori dall'area annotata
            o = dict(win_id=w["win_id"], archetype=A, roi_id=r, source=lab, n_in_window=int(inwin.sum()), n_total=len(P), nn_median=m_all,
                     frac_px_kept=k.mean())
            for u in U_REL: o[f"H_rel_{u}"] = np.mean(dq[k] / m_all > u)
            for a in D_ABS: o[f"H_abs_{int(a)}"] = np.mean(dq[k] > a)
            rows.append(o)
    return pd.DataFrame(rows)

if __name__ == "__main__":
    import os; os.makedirs(RES, exist_ok=True)
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    jobs = [(a, r) for a, r in zip(ROIS.archetype, ROIS.roi_id)]
    with Pool(nproc) as p: out = p.map(one_roi, jobs)
    pd.DataFrame([x for o in out for x in o[0]]).to_csv(f"{RES}/S1.5d_explore_space.csv", index=False)
    pd.DataFrame([x for o in out for x in o[1]]).to_csv(f"{RES}/S1.5d_explore_he.csv", index=False)
    manual_windows().to_csv(f"{RES}/S1.5d_explore_manual.csv", index=False)
    print("ok")
