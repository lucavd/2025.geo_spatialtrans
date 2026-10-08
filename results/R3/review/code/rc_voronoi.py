#!/usr/bin/env python
"""results/R3/review/code/rc_voronoi.py — revisione avversariale R3 (revisore 'code').
Ricalcolo INDIPENDENTE dai soli ingressi grezzi (R2_nuclei_all.parquet, maschere ds4, label native),
senza usare tools/R3_*. Voronoi: shapely.voronoi_polygons(ordered=True); momenti del territorio per
decomposizione in triangoli + autovettori (numpy.linalg.eigh); contenimento per punto-in-poligono sul
territorio (non KD-tree); orientazione nucleare per autovettori della covarianza dei pixel.
Uso: .venv/bin/python results/R3/review/code/rc_voronoi.py <A> <roi> <method> [jitter_um] [seed]
"""
import sys, json, os, time
import numpy as np, pandas as pd, shapely
from shapely.geometry import box
from PIL import Image
from scipy.signal import fftconvolve

OUTD = "results/R3/review/code"
A, ROI, METHOD = sys.argv[1:4]
JIT = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
JSEED = int(sys.argv[5]) if len(sys.argv) > 5 else 1
tag = f"{A}_{ROI}_{METHOD}" + (f"_jit{JIT:g}_s{JSEED}" if JIT > 0 else "")
t0 = time.time()

rois = pd.read_csv("results/R2/R2_rois_checked.csv")
rr = rois[(rois.archetype == A) & (rois.roi_id == ROI)].iloc[0]
S = float(rr.side_px) * float(rr.um_per_px)

# ---- maschera: pixel (i,j) copre x in [j*px,(j+1)*px], y in [i*px,(i+1)*px] (y verso il basso)
m = np.array(Image.open(f"results/R2/tissue_masks/{A}_{ROI}_valid_ds4.png"))
if m.ndim == 3: m = m[..., 0]
m = m > 127
H, W = m.shape; assert H == W
px = S / W
boxes = []
for i in range(H):
    row = np.concatenate([[False], m[i], [False]]).astype(np.int8)
    d = np.diff(row); st = np.where(d == 1)[0]; en = np.where(d == -1)[0]
    for a, b in zip(st, en):
        boxes.append(box(a * px, i * px, b * px, (i + 1) * px))
mask = shapely.union_all(boxes)
shapely.prepare(mask)
area_mask = mask.area; area_px = m.sum() * px * px

# ---- nuclei
nuc = pd.read_parquet("results/R2/R2_nuclei_all.parquet",
                      filters=[("archetype", "==", A), ("roi_id", "==", ROI), ("method", "==", METHOD)])
nuc = nuc[(nuc["scale"] == 1) & (nuc["keep"])].reset_index(drop=True)
n_keep = len(nuc)
ci = np.floor(nuc.x_um.values / px).astype(int); ri = np.floor(nuc.y_um.values / px).astype(int)
ok = (ci >= 0) & (ri >= 0) & (ci < W) & (ri < H)
inm = np.zeros(len(nuc), bool); inm[ok] = m[ri[ok], ci[ok]]
upp = float(nuc.um_per_px.iloc[0])
# convenzione: centroidi in centro-pixel (pixel 0 -> 0); con lati-pixel sarebbero +upp/2
ci2 = np.floor((nuc.x_um.values + upp / 2) / px).astype(int); ri2 = np.floor((nuc.y_um.values + upp / 2) / px).astype(int)
ok2 = (ci2 >= 0) & (ri2 >= 0) & (ci2 < W) & (ri2 < H)
inm2 = np.zeros(len(nuc), bool); inm2[ok2] = m[ri2[ok2], ci2[ok2]]
g = nuc[inm].reset_index(drop=True)
n = len(g)
x = g.x_um.values.copy(); y = g.y_um.values.copy()
if JIT > 0:
    rng = np.random.default_rng(JSEED)
    x = x + rng.uniform(-JIT, JIT, n); y = y + rng.uniform(-JIT, JIT, n)
n_dup = n - len(np.unique(np.round(np.column_stack([x, y]), 9), axis=0))

# ---- Voronoi (ordinato come i punti) ritagliato al riquadro allargato di 50 um
rw = box(-50, -50, S + 50, S + 50)
vor = shapely.voronoi_polygons(shapely.MultiPoint(np.column_stack([x, y])), extend_to=rw, ordered=True)
tiles = np.array(shapely.get_parts(vor))
assert len(tiles) == n, (len(tiles), n)
tiles = shapely.intersection(tiles, rw)
pts = shapely.points(x, y)
gen_in_own = shapely.contains(tiles, pts)
area_tile = shapely.area(tiles)
within = shapely.contains(mask, tiles)
frame = box(0, 0, S, S); shapely.prepare(frame)
in_frame = shapely.contains(frame, tiles)
area = area_tile.copy()
for k in np.where(~within)[0]:
    loc = shapely.clip_by_rect(mask, *tiles[k].bounds)
    area[k] = 0.0 if loc.is_empty else shapely.intersection(tiles[k], loc).area
interior = within | (in_frame & (area / area_tile > 0.999))

# ---- momenti del territorio: triangoli a ventaglio dal primo vertice + autovettori
def poly_shape(poly):
    xy = np.asarray(poly.exterior.coords)[:-1]
    p0 = xy[0]; q = xy - p0
    a, b = q[1:-1], q[2:]
    Ar = 0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])          # aree con segno
    At = Ar.sum()
    # integrali su triangolo (0,a,b): ∫x = A(ax+bx)/3; ∫x^2 = A/6 (ax^2+bx^2+ax bx); ∫xy = A/12 (2ax ay + 2bx by + ax by + bx ay)
    Sx = (Ar * (a[:, 0] + b[:, 0]) / 3).sum(); Sy = (Ar * (a[:, 1] + b[:, 1]) / 3).sum()
    Sxx = (Ar / 6 * (a[:, 0]**2 + b[:, 0]**2 + a[:, 0] * b[:, 0])).sum()
    Syy = (Ar / 6 * (a[:, 1]**2 + b[:, 1]**2 + a[:, 1] * b[:, 1])).sum()
    Sxy = (Ar / 12 * (2 * a[:, 0] * a[:, 1] + 2 * b[:, 0] * b[:, 1] + a[:, 0] * b[:, 1] + b[:, 0] * a[:, 1])).sum()
    cx, cy = Sx / At, Sy / At
    C = np.array([[Sxx / At - cx**2, Sxy / At - cx * cy], [Sxy / At - cx * cy, Syy / At - cy**2]])
    w, v = np.linalg.eigh(C)
    ecc = np.sqrt(max(1 - w[0] / w[1], 0.0)) if w[1] > 0 else 0.0
    th = np.arctan2(v[1, 1], v[0, 1])
    return abs(At), ecc, th, len(xy)
sh = np.array([poly_shape(t) for t in tiles])
ecc_T = sh[:, 1]; theta_T = sh[:, 2]; nsides_vert = sh[:, 3].astype(int)

# ---- nuclei: pixel nativi (x = c*upp, y = r*upp), momenti e contenimento
lab = np.load(f"/mnt/micron/geo_spatialtrans/R2/masks/{A}_{ROI}_{METHOD}_native.npz")["labels"]
lab2gen = np.full(int(lab.max()) + 1, -1, np.int64); lab2gen[g.label.values.astype(np.int64)] = np.arange(n)
rr_, cc_ = np.nonzero(lab); L = lab[rr_, cc_]; own = lab2gen[L]; kk = own >= 0
rr_, cc_, own = rr_[kk], cc_[kk], own[kk]
X = cc_ * upp; Y = rr_ * upp
npx = np.bincount(own, minlength=n).astype(float)
mx = np.bincount(own, X, n) / npx; my = np.bincount(own, Y, n) / npx
dx = X - mx[own]; dy = Y - my[own]
cxx = np.bincount(own, dx * dx, n) / npx; cyy = np.bincount(own, dy * dy, n) / npx; cxy = np.bincount(own, dx * dy, n) / npx
Cn = np.stack([np.stack([cxx, cxy], -1), np.stack([cxy, cyy], -1)], -2)
wN, vN = np.linalg.eigh(Cn)
with np.errstate(invalid="ignore", divide="ignore"):
    ecc_N = np.sqrt(np.clip(1 - wN[:, 0] / wN[:, 1], 0, None))
theta_N = np.arctan2(vN[:, 1, 1], vN[:, 0, 1])
cent_err = np.hypot(mx - g.x_um.values, my - g.y_um.values)  # label <-> generatore (senza jitter)
shapely.prepare(tiles)
inside = shapely.contains_xy(tiles[own], X, Y)
frac_out = np.bincount(own[~inside], minlength=n) / npx

# ---- CV_loc indipendente: kernel gaussiano con correzione di bordo uniforme e(u) = (k * W)(u)
lam = n / area_px; sig = 5 / np.sqrt(lam)
spx = sig / px; R = int(np.ceil(5 * spx))
gg = np.arange(-R, R + 1) * px; K = np.exp(-(gg[:, None]**2 + gg[None, :]**2) / (2 * sig**2)); K *= px * px / (2 * np.pi * sig**2)
e = fftconvolve(m.astype(float), K, mode="same")
def kern_sum(loo):
    out = np.zeros(n)
    from scipy.spatial import cKDTree
    T = cKDTree(np.column_stack([x, y]))
    nb = T.query_ball_point(np.column_stack([x, y]), r=6 * sig, workers=4)
    for i, js in enumerate(nb):
        js = np.asarray(js); d2 = (x[js] - x[i])**2 + (y[js] - y[i])**2
        v = np.exp(-d2 / (2 * sig**2)); out[i] = v.sum() - (1.0 if loo else 0.0)
    return out / (2 * np.pi * sig**2)
eu = e[np.clip(np.floor(y / px).astype(int), 0, H - 1), np.clip(np.floor(x / px).astype(int), 0, W - 1)]
lam_loc = kern_sum(False) / eu; lam_loc_loo = kern_sum(True) / eu

cells = pd.DataFrame(dict(idx=np.arange(1, n + 1), label=g.label.values, x=x, y=y, area_tile=area_tile, area=area,
                          interior=interior, within=within, gen_in_own=gen_in_own, nsides=nsides_vert, ecc_T=ecc_T,
                          theta_T=theta_T, area_nuc=g.area_um2.values, ecc_N=ecc_N, theta_N=theta_N, n_px=npx,
                          frac_out=frac_out, cent_err=cent_err, ecc_pq=g.eccentricity.values,
                          theta_pq=np.pi / 2 - g.orientation.values, lam_loc=lam_loc, lam_loc_loo=lam_loc_loo))
cells.to_parquet(f"{OUTD}/cells/{tag}.parquet", index=False)

def axd(a, b):
    d = np.abs(a - b) % np.pi; return np.minimum(d, np.pi - d)
I = interior
a_i = area[I]
el = I & (ecc_N >= 0.8) & (ecc_T >= 0.5) & np.isfinite(theta_N)
dth = np.degrees(axd(theta_T[el], theta_N[el]))
rng = np.random.default_rng(20261008)
perm = np.array([np.median(np.degrees(axd(theta_T[el], rng.permutation(theta_N[el])))) for _ in range(2000)])
al = a_i * lam_loc[I]; al2 = a_i * lam_loc_loo[I]
res = dict(archetype=A, roi_id=ROI, method=METHOD, jitter_um=JIT, jseed=JSEED, side_um=S, px_mask=px, upp=upp,
           n_keep=n_keep, n=n, n_inmask_edgeconv=int(inm2.sum()), n_dup=int(n_dup),
           area_mask_poly=area_mask, area_mask_px=area_px, sum_area=float(area.sum()),
           rel_err_sum=float(area.sum() / area_mask - 1), gen_in_own=float(gen_in_own.mean()),
           n_interior=int(I.sum()), frac_interior=float(I.mean()),
           median_area=float(np.median(a_i)), mean_area=float(a_i.mean()), median_eq_r=float(np.median(np.sqrt(a_i / np.pi))),
           cv=float(a_i.std(ddof=1) / a_i.mean()), median_ecc_T=float(np.median(ecc_T[I])),
           mean_nsides=float(nsides_vert[I].mean()),
           median_nc=float(np.median(g.area_um2.values[I] / a_i)), median_ratio=float(np.median(np.sqrt(g.area_um2.values[I] / a_i))),
           frac_cut=float(np.mean(frac_out[I] > 0.05)), median_frac_out=float(np.median(frac_out[I])),
           n_eligible=int(el.sum()), median_dtheta=float(np.median(dth)) if el.any() else None,
           perm_p=float(np.mean(perm <= np.median(dth))) if el.any() else None, perm_median=float(np.median(perm)) if el.any() else None,
           median_dtheta_mirror=float(np.median(np.degrees(axd(theta_T[el], -theta_N[el])))) if el.any() else None,
           sigma_loc=float(sig), cv_loc=float(al.std(ddof=1) / al.mean()), cv_loc_loo=float(al2.std(ddof=1) / al2.mean()),
           cent_err_max=float(np.nanmax(cent_err)), cent_err_median=float(np.nanmedian(cent_err)),
           n_gen_without_px=int((npx == 0).sum()), secs=time.time() - t0)
json.dump(res, open(f"{OUTD}/roi/{tag}.json", "w"), indent=1)
print(json.dumps(res))
