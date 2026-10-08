# rv_task1_real_pcf.py — revisione avversariale S1.2 / checkB, compito 1.
# Stimatore g(r) indipendente (numpy/scipy, nessuna funzione del progetto), su TUTTI i 30 ROI.
# Finestra = maschera results/R2/tissue_masks/<A>_<roi>_valid_ds4.png (riga 1 in alto, y verso il basso),
# pixel = side_um/ncol, side_um = side_px*um_per_px. Nuclei: R2_nuclei_all.parquet, scale==1 & keep.
# Stime:  (a) "anelli": conteggio coppie per anello [r-0.25, r+0.25), correzione di traslazione e = |W|/|W ∩ W+v|
#         (covarianza d'insieme via FFT della maschera, interpolazione bilineare), nessun lisciamento;
#         (b) "kernel": Epanechnikov con sd = bw del progetto (semi-ampiezza h = sqrt(5) bw), divisore d (1/(2 pi d_ij)),
#             normalizzazione n(n-1)/|W|; tre varianti al bordo r=0: nessuna, rinormalizzazione della massa del kernel
#             su [0, inf) ("renorm"), riflessione.
#         (c) controprove di orientamento/unita': y NON invertita (specchio verticale), e coordinate x0.5 / x2.
# Uso: cd ~/2025.geo_spatialtrans && .venv/bin/python results/S1.2/review/checkB/rv_task1_real_pcf.py
import numpy as np, pandas as pd, pyarrow.parquet as pq
from scipy.spatial import cKDTree
from scipy.signal import fftconvolve
from scipy.ndimage import map_coordinates
from skimage import io
RV = "results/S1.2/review/checkB"
R = np.arange(0, 30.0001, 0.5)
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
prim = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A5": "cellpose_rgb", "A4": "spaceranger", "A6": "spaceranger"}
ref = pd.read_csv(f"{RV}/data/real_pcf_export.csv")
nuc = pq.read_table("results/R2/R2_nuclei_all.parquet", columns=["archetype", "roi_id", "method", "scale", "keep", "x_um", "y_um"]).to_pandas()
nuc = nuc[(nuc.scale == 1) & nuc.keep]

def window(A, roi):
    rr = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]
    side = rr.side_px * rr.um_per_px
    m = io.imread(f"results/R2/tissue_masks/{A}_{roi}_valid_ds4.png")
    if m.ndim == 3: m = m[..., 0]
    M = (m > 127).astype(float)            # riga 0 = alto
    px = side / M.shape[1]
    # covarianza d'insieme: C[dy, dx] = area(W ∩ (W + v)), v in pixel; centro all'indice (ny-1, nx-1)
    C = fftconvolve(M, M[::-1, ::-1], mode="full") * px * px
    return dict(M=M.astype(bool), px=px, side=side, area=M.sum() * px * px, C=C, ny=M.shape[0], nx=M.shape[1])

def inside(w, x, y_down):
    col = np.floor(x / w["px"]).astype(int); row = np.floor(y_down / w["px"]).astype(int)
    ok = (col >= 0) & (col < w["nx"]) & (row >= 0) & (row < w["ny"])
    ok2 = ok.copy(); ok2[ok] = w["M"][row[ok], col[ok]]
    return ok2

def pairs(w, x, y, rmax):
    P = np.c_[x, y]; t = cKDTree(P)
    ij = t.query_pairs(rmax, output_type="ndarray")
    dx = P[ij[:, 1], 0] - P[ij[:, 0], 0]; dy = P[ij[:, 1], 1] - P[ij[:, 0], 1]
    d = np.hypot(dx, dy)
    cy = w["ny"] - 1 + dy / w["px"]; cx = w["nx"] - 1 + dx / w["px"]
    Cv = map_coordinates(w["C"], [cy, cx], order=1)
    Cv2 = map_coordinates(w["C"], [w["ny"] - 1 - dy / w["px"], w["nx"] - 1 - dx / w["px"]], order=1)  # simmetria (controllo)
    e = w["area"] / Cv
    return d, e, np.max(np.abs(Cv - Cv2) / Cv)

def g_ring(d, e, n, area):
    l2a = n * (n - 1) / area
    out = np.full(R.size, np.nan)
    for k, r in enumerate(R):
        if r == 0: continue
        a, b = r - 0.25, r + 0.25
        s = (d >= a) & (d < b)
        out[k] = 2 * e[s].sum() / (np.pi * (b * b - a * a)) / l2a     # fattore 2: coppie ordinate
    return out

def g_kernel(d, e, n, area, bw, zc="none"):
    h = np.sqrt(5) * bw; l2a = n * (n - 1) / area
    out = np.zeros(R.size)
    for k, r in enumerate(R):
        u = (r - d) / h; s = np.abs(u) < 1
        kv = 0.75 / h * (1 - u[s] ** 2)
        if zc == "reflect":
            u2 = (r + d) / h; s2 = np.abs(u2) < 1
            out[k] += 2 * np.sum(0.75 / h * (1 - u2[s2] ** 2) * e[s2] / (2 * np.pi * d[s2])) / l2a
        out[k] += 2 * np.sum(kv * e[s] / (2 * np.pi * d[s])) / l2a
        if zc == "renorm":  # massa del kernel centrato in r che cade in [0, inf)
            lo = max(-1.0, -r / h); mass = 0.75 * ((1 - lo) - (1 - lo ** 3) / 3) if lo > -1 else 1.0
            out[k] /= mass
    return out

rows, curves = [], []
for _, rr in rois.iterrows():
    A, roi = rr.archetype, rr.roi_id
    w = window(A, roi)
    for role, meth in (("primary", prim[A]), ("secondary", "stardist_he")):
        rf = ref[(ref.archetype == A) & (ref.roi_id == roi) & (ref.role == role)].sort_values("r")
        assert rf.method.iloc[0] == meth, (A, roi, role, rf.method.iloc[0])
        bw = rf.bw.iloc[0]
        dd = nuc[(nuc.archetype == A) & (nuc.roi_id == roi) & (nuc.method == meth)]
        x, y = dd.x_um.to_numpy(), dd.y_um.to_numpy()
        ok = inside(w, x, y)
        ok_flip = inside(w, x, w["side"] - y)       # orientazione sbagliata (y non invertita rispetto alla maschera)
        xs, ys = x[ok], y[ok]; n = ok.sum()
        d, e, asym = pairs(w, xs, ys, 30 + np.sqrt(5) * bw + 1)
        gr = g_ring(d, e, n, w["area"])
        gk = {z: g_kernel(d, e, n, w["area"], bw, z) for z in ("none", "renorm", "reflect")}
        # controprova: pattern con y specchiata ma finestra non specchiata (errore tipico), stessi punti tenuti
        yf = w["side"] - y[ok_flip]; xf = x[ok_flip]
        df_, ef_, _ = pairs(w, xf, yf, 30 + np.sqrt(5) * bw + 1)
        gk_flip = g_kernel(df_, ef_, ok_flip.sum(), w["area"], bw, "renorm")
        gref = rf.g.to_numpy()
        I = R > 0
        trap = lambda f: np.trapezoid(f[I], R[I])
        rec = dict(archetype=A, roi_id=roi, role=role, method=meth, n_ref=int(rf.n.iloc[0]), n_lost_ref=int(rf.n_lost.iloc[0]),
                   n_rv=int(n), n_lost_rv=int((~ok).sum()), n_lost_if_yflip=int((~ok_flip).sum()),
                   area_mm2_ref=rf.area_mm2.iloc[0], area_mm2_rv=w["area"] / 1e6, px_um=w["px"], bw=bw, setcov_asym=asym,
                   lambda_rv=n / w["area"] * 1e6, lambda_ref=rf["lambda"].iloc[0])
        for z, g in gk.items():
            rec[f"maxabs_kernel_{z}"] = np.nanmax(np.abs(g[I] - gref[I])); rec[f"ISE_kernel_{z}"] = trap((g - gref) ** 2)
        rec["maxabs_kernel_renorm_r_ge_5"] = np.max(np.abs(gk["renorm"] - gref)[R >= 5])
        rec["ISE_ring_vs_ref"] = trap((gr - gref) ** 2)
        rec["ISE_yflip_vs_ref"] = trap((gk_flip - gref) ** 2)
        for rv in (2, 4, 6, 30):
            k = np.where(R == rv)[0][0]; rec[f"g{rv}_ref"] = gref[k]; rec[f"g{rv}_ring"] = gr[k]; rec[f"g{rv}_kern"] = gk["renorm"][k]
        rows.append(rec)
        curves.append(pd.DataFrame(dict(archetype=A, roi_id=roi, role=role, r=R, g_ref=gref, g_ring=gr, g_kern_none=gk["none"],
                                        g_kern_renorm=gk["renorm"], g_kern_reflect=gk["reflect"], g_kern_yflip=gk_flip)))
    print(A, roi, flush=True)
out = pd.DataFrame(rows)
out.to_csv(f"{RV}/rv_task1_real_pcf_summary.csv", index=False)
pd.concat(curves).to_csv(f"{RV}/rv_task1_real_pcf_curves.csv", index=False)
