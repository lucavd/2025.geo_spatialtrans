"""tools/S1.5e_metrics.py - S1.5e: forma dei vuoti e densita' del ROI (definizioni: docs/steps/S1.5e.md, esplorazione E4), stesse funzioni per
nuclei veri e simulazioni. Riusa la geometria del ROI di tools/S1.5d_metrics.py (maschera di tessuto valido di R2, pixel di circa 1.1 um,
cornice di 30 um dal bordo del ROI).
- void_shape(x, y, interior, g): pixel vuoti = pixel della maschera (centro del pixel, nella cornice) a distanza > 3 volte la mediana della distanza
  dal vicino delle celle interne dal centro piu' vicino (come H_rel_3); componenti connesse a 8 vicini di almeno 50 um^2; per componente
  allungamento = sqrt(autovalore massimo / minimo) della covarianza delle coordinate dei pixel (um); uscite: elong_wmed (mediana pesata per area;
  1.0 se nessuna componente), frac_elong_gt3 (quota dell'area vuota in componenti con allungamento > 3; NaN se nessuna), n_comp.
  Identica a tools/S1.5e_explore.py (components, summarise_comps) per costruzione: G-4 la ricalcola con skimage.measure.regionprops.
- calibration_metrics_e: gli 8 bersagli di S1.5d + elong_wmed, con le celle interne della regola della calibrazione di S1.5d."""
import importlib.util, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy import ndimage as ndi
def _mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
m5d = _mod("m5d_e", "tools/S1.5d_metrics.py")
U_VOID = 3.0; MIN_COMP_UM2 = 50.0
class RoiGeomE(m5d.RoiGeom):
    def __init__(self, A, r):
        super().__init__(A, r)
        rr, cc = np.nonzero(self.m); x = (cc + 0.5) * self.px; y = (rr + 0.5) * self.px; F = m5d.FRAME_B
        k = (x > F) & (y > F) & (x < self.side - F) & (y < self.side - F); self.rr, self.cc = rr[k], cc[k]
        self.mask_area_um2 = float(self.m.sum() * self.px * self.px)
def void_image(x, y, interior, g):
    P = np.c_[x, y]; I = np.asarray(interior, bool); tr = cKDTree(P); dd, _ = tr.query(P, k=2); nnm = float(np.median(dd[I, 1]))
    d, _ = tr.query(np.c_[(g.cc + 0.5) * g.px, (g.rr + 0.5) * g.px])
    v = np.zeros_like(g.m, bool); k = d > U_VOID * nnm; v[g.rr[k], g.cc[k]] = True
    return v
def components(v, px):
    lab, n = ndi.label(v, structure=np.ones((3, 3), bool)); out = []
    if n == 0: return out
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        cm = lab[sl] == i; a = cm.sum() * px * px
        if a < MIN_COMP_UM2: continue
        yy, xx = np.nonzero(cm); X = np.c_[xx, yy] * px; ev = np.maximum(np.linalg.eigh(np.cov(X.T))[0], 1e-12)
        out.append((a, float(np.sqrt(ev[1] / ev[0]))))
    return out
def summarise(cs):
    if not cs: return dict(elong_wmed=1.0, frac_elong_gt3=np.nan, n_comp=0)
    a = np.array([c[0] for c in cs]); e = np.array([c[1] for c in cs]); w = a / a.sum(); o = np.argsort(e); cw = np.cumsum(w[o])
    return dict(elong_wmed=float(e[o][np.searchsorted(cw, 0.5)]), frac_elong_gt3=float(w[e > 3].sum()), n_comp=int(len(cs)))
def void_shape(x, y, interior, g):
    return summarise(components(void_image(x, y, interior, g), g.px))
def roi_density(n_cells, g):
    return n_cells / g.mask_area_um2 * 1e6
CAL_TARGETS_E = m5d.CAL_TARGETS + ["elong_wmed"]
def calibration_metrics_e(x, y, g):
    I = m5d.interior_by_distance(x, y, g); o = m5d.position_metrics(x, y, I, np.ones(len(x)), g); o.pop("frac_gt3"); o["n_interior_cal"] = int(I.sum())
    o.update(void_shape(x, y, I, g)); return o
