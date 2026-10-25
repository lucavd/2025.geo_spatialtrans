"""tools/S1.5d_metrics.py - S1.5d: metriche della disposizione dei centri, stesse funzioni per i nuclei veri e per le simulazioni
(definizioni: docs/steps/S1.5d.md). Ingresso per ROI: centri (um), indicatore di cella interna, area del territorio, maschera di tessuto valido
di R2 (results/R2/tissue_masks/<A>_<roi>_valid_ds4.png, pixel di circa 1.1 um).
- H_rel_u: frazione dei pixel della maschera (ad almeno 30 um dal bordo del ROI) la cui distanza dal centro piu' vicino supera u volte la
  mediana della distanza dal vicino piu' vicino delle celle interne (u = 1.5, 2, 3, 4). Adimensionale; invariante al diradamento
  indipendente solo per un processo di Poisson (non per centri con distanza minima: nota D-3 di docs/steps/S1.5d.md).
- xcv2_q: CV^2 in eccesso dei conteggi in quadrati di lato q um interamente nella maschera ((varianza - media) / media^2; q = 25, 50, 100).
  0 per un processo di Poisson; esattamente invariante al diradamento indipendente ((p^2 V + p(1-p) m - p m) / (p m)^2 = (V - m) / m^2).
- nn_cv: CV della distanza dal vicino piu' vicino delle celle interne.
- frac_gt3: frazione dei territori delle celle interne sopra 3 volte la mediana (come PH-1 di S1.5c)."""
import numpy as np, pandas as pd
from PIL import Image
from scipy.spatial import cKDTree
U_REL = [1.5, 2.0, 3.0, 4.0]; QSIDE = [25.0, 50.0, 100.0]; FRAME_B = 30.0
_ROIS = None
def rois_table():
    global _ROIS
    if _ROIS is None: _ROIS = pd.read_csv("results/R2/R2_rois_checked.csv")
    return _ROIS
def mask_of(A, r):
    m = np.array(Image.open(f"results/R2/tissue_masks/{A}_{r}_valid_ds4.png"))
    if m.ndim == 3: m = m[..., 0]
    m = m > 127; row = rois_table(); row = row[(row.archetype == A) & (row.roi_id == r)].iloc[0]
    side = row.side_px * row.um_per_px
    return m, side / m.shape[1], side
def mask_pixels(m, px, side):
    rr, cc = np.nonzero(m); x = (cc + 0.5) * px; y = (rr + 0.5) * px
    k = (x > FRAME_B) & (y > FRAME_B) & (x < side - FRAME_B) & (y < side - FRAME_B)
    return np.c_[x[k], y[k]]
def quadrat_boxes(m, px, side):
    out = {}
    for q in QSIDE:
        n = int((side - 2 * FRAME_B) // q); bx = []
        for i in range(n):
            for j in range(n):
                x0 = FRAME_B + i * q; y0 = FRAME_B + j * q
                c0, c1, r0, r1 = int(x0 / px), int((x0 + q) / px), int(y0 / px), int((y0 + q) / px)
                if m[r0:r1, c0:c1].mean() >= 0.99: bx.append((x0, y0))
        out[q] = np.array(bx, float).reshape(-1, 2)
    return out
class RoiGeom:
    """Pixel della maschera e quadrati di un ROI, calcolati una volta e usati per il reale e per ogni simulazione."""
    def __init__(self, A, r):
        self.m, self.px, self.side = mask_of(A, r); self.Q = mask_pixels(self.m, self.px, self.side); self.boxes = quadrat_boxes(self.m, self.px, self.side)
def position_metrics(x, y, interior, a_t, g):
    P = np.c_[x, y]; I = np.asarray(interior, bool); tr = cKDTree(P)
    dd, _ = tr.query(P, k=2); nn = dd[I, 1]; m_all = np.median(nn)
    d, _ = tr.query(g.Q)
    o = {"nn_median": m_all, "nn_cv": nn.std(ddof=1) / nn.mean()}
    for u in U_REL: o[f"H_rel_{u}"] = float(np.mean(d / m_all > u))
    for q in QSIDE:
        bx = g.boxes[q]
        if len(bx) < 5: o[f"xcv2_{int(q)}"] = np.nan; continue
        ix = np.floor((P[:, 0] - bx[:, 0].min()) / q).astype(int); iy = np.floor((P[:, 1] - bx[:, 1].min()) / q).astype(int)
        key = pd.Series(1, index=pd.MultiIndex.from_arrays([ix, iy])).groupby(level=[0, 1]).size()
        bi = np.round((bx[:, 0] - bx[:, 0].min()) / q).astype(int); bj = np.round((bx[:, 1] - bx[:, 1].min()) / q).astype(int)
        cnt = key.reindex(pd.MultiIndex.from_arrays([bi, bj]), fill_value=0).values.astype(float)
        o[f"xcv2_{int(q)}"] = (cnt.var(ddof=1) - cnt.mean()) / cnt.mean() ** 2 if cnt.mean() > 0 else np.nan
    a = np.asarray(a_t, float)[I]; o["frac_gt3"] = float(np.mean(a / np.median(a) > 3))
    return o

def interior_by_distance(x, y, g, edge_um=10.0):
    """Regola della calibrazione dell'eterogeneita' (docs/steps/S1.5d.md): centro ad almeno edge_um dal bordo della maschera e FRAME_B dal bordo
    del ROI. Distanza dal bordo = trasformata di distanza euclidea della maschera (pixel di circa 1.1 um), stessa regola per vero e simulato."""
    from scipy.ndimage import distance_transform_edt
    if not hasattr(g, "edt"): g.edt = distance_transform_edt(g.m) * g.px
    ci = np.clip((np.asarray(x) / g.px).astype(int), 0, g.m.shape[1] - 1); ri = np.clip((np.asarray(y) / g.px).astype(int), 0, g.m.shape[0] - 1)
    x = np.asarray(x); y = np.asarray(y)
    return (g.edt[ri, ci] >= edge_um) & (x > FRAME_B) & (y > FRAME_B) & (x < g.side - FRAME_B) & (y < g.side - FRAME_B)
CAL_TARGETS = ["H_rel_1.5", "H_rel_2.0", "H_rel_3.0", "H_rel_4.0", "xcv2_25", "xcv2_50", "xcv2_100", "nn_cv"]
def calibration_metrics(x, y, g):
    I = interior_by_distance(x, y, g); o = position_metrics(x, y, I, np.ones(len(x)), g); o.pop("frac_gt3"); o["n_interior_cal"] = int(I.sum()); return o
