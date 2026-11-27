"""Diagnostic of chain leg 1 (Xenium -> H&E, 10x alignment matrix), independent of expression.
Xenium nucleus centroids (nucleus_boundaries vertex means) carried through the chain vs the
haematoxylin channel of the H&E, both on a 1 um grid of the 2 um-array frame. Shift h maximises
corr(N(x + h), Hem(x)) (same sign convention as s). Also checks that the Xenium H&E and the Visium HD
H&E are the same image (pixel equality on the region crops)."""
import json, time, numpy as np, pandas as pd, pyarrow.parquet as pq, tifffile, zarr
from scipy.ndimage import gaussian_filter
from skimage.color import rgb2hed
from g1_common import *
chain = load_chain(); H = np.array(chain["H_px2arr"]); Hinv = np.linalg.inv(H)
BTF = f"{D}/visium/Visium_HD_Human_Lung_Cancer_post_Xenium_v1_Experiment1_tissue_image.btf"
nb = pq.read_table(f"{XEN}_nucleus_boundaries.parquet", columns=["cell_id", "vertex_x", "vertex_y"]).to_pandas()
cen = nb.groupby("cell_id", sort=False)[["vertex_x", "vertex_y"]].mean().to_numpy(float)
na = xenium_to_array(cen, chain)
def open_level0(path):
    tf = tifffile.TiffFile(path); z = zarr.open(tf.series[0].levels[0].aszarr(), mode="r")
    if not hasattr(z, "shape"):          # multiscale group: take the full-resolution array
        keys = sorted(z.array_keys()) if hasattr(z, "array_keys") else sorted(z.keys())
        z = z[keys[0]]
    return tf, z
tfx, zx = open_level0(f"{XEN}_he_image.ome.tif"); tfv, zv = open_level0(BTF)
print("shapes", zx.shape, zv.shape, flush=True)
def crop(z, r0, r1, c0, c1):
    a = np.asarray(z[r0:r1, c0:c1]) if z.shape[-1] == 3 else np.moveaxis(np.asarray(z[:, r0:r1, c0:c1]), 0, -1)
    return a
def xcorr(A, B, win):
    """offset (dx, dy) in grid units maximising corr(A shifted, B) over the inner part; A(x+d) vs B(x)."""
    n0, n1 = B.shape; m = win
    b = B[m:n0 - m, m:n1 - m]; b = (b - b.mean()) / b.std()
    best = (-9, 0, 0); C = np.full((2 * win + 1, 2 * win + 1), np.nan)
    for dy in range(-win, win + 1):
        for dx in range(-win, win + 1):
            a = A[m + dy:n0 - m + dy, m + dx:n1 - m + dx]
            C[dy + win, dx + win] = ((a - a.mean()) * b).mean() / (a.std() + 1e-12)
    iy, ix = np.unravel_index(np.nanargmax(C), C.shape); off = np.zeros(2)
    if 0 < iy < 2 * win and 0 < ix < 2 * win:   # separable parabolic refinement (diagnostic only)
        for k, (cm, c0, cp) in enumerate([(C[iy, ix - 1], C[iy, ix], C[iy, ix + 1]), (C[iy - 1, ix], C[iy, ix], C[iy + 1, ix])]):
            den = cm - 2 * c0 + cp; off[k] = 0.5 * (cm - cp) / den if den < 0 else 0
    return ix - win + off[0], iy - win + off[1], float(np.nanmax(C)), bool(iy in (0, 2 * win) or ix in (0, 2 * win))
rows = []; same = []
for R in chain["regions"]:
    t0 = time.time(); side = chain["region_side_um"]; mg = 40.0   # um margin
    a0 = np.array([R["x0"] / 2, R["y0"] / 2]) - mg / 2; a1 = a0 + (side + 2 * mg) / 2
    corners = np.array([[a0[0], a0[1]], [a1[0], a0[1]], [a0[0], a1[1]], [a1[0], a1[1]]])
    pc = proj_apply(Hinv, corners); c0, r0 = np.floor(pc.min(0)).astype(int); c1, r1 = np.ceil(pc.max(0)).astype(int)
    img = crop(zx, r0, r1, c0, c1); imgv = crop(zv, r0, r1, c0, c1)
    same.append(dict(region=R["name"], n_px=int(img.size), n_diff=int((img != imgv).sum()), max_absdiff=int(np.abs(img.astype(int) - imgv.astype(int)).max())))
    hem = rgb2hed(img)[..., 0]
    rr, cc = np.mgrid[r0:r1, c0:c1]
    aa = proj_apply(H, np.c_[cc.ravel(), rr.ravel()].astype(float))
    u = 2 * (aa - a0)                                  # local um, 1 um cells
    n = int(side + 2 * mg); k = np.floor(u).astype(int); ok = (k[:, 0] >= 0) & (k[:, 0] < n) & (k[:, 1] >= 0) & (k[:, 1] < n)
    idx = k[ok, 1] * n + k[ok, 0]
    hsum = np.bincount(idx, weights=hem.ravel()[ok], minlength=n * n); hcnt = np.bincount(idx, minlength=n * n)
    Hm = (hsum / np.maximum(hcnt, 1)).reshape(n, n)
    un = 2 * (na - a0); kn = np.floor(un).astype(int); okn = (kn[:, 0] >= 0) & (kn[:, 0] < n) & (kn[:, 1] >= 0) & (kn[:, 1] < n)
    Nm = np.bincount(kn[okn, 1] * n + kn[okn, 0], minlength=n * n).reshape(n, n).astype(float)
    Hs = gaussian_filter(Hm, 1.5); Ns = gaussian_filter(Nm, 2.5)
    hx, hy, pk, edge = xcorr(Ns, Hs, 30)
    rows.append(dict(region=R["name"], scope="region", tx=-1, ty=-1, hx_um=hx, hy_um=hy, peak=pk, edge=edge, n_nuclei=int(okn.sum())))
    # per 200 um tile of the core (core starts at mg um), tile + 30 um margin each side
    for ty in range(5):
        for tx in range(5):
            x0 = int(mg + 200 * tx) - 30; y0 = int(mg + 200 * ty) - 30
            hx_, hy_, pk_, e_ = xcorr(Ns[y0:y0 + 260, x0:x0 + 260], Hs[y0:y0 + 260, x0:x0 + 260], 30)
            rows.append(dict(region=R["name"], scope="tile", tx=tx, ty=ty, hx_um=hx_, hy_um=hy_, peak=pk_, edge=e_, n_nuclei=-1))
    print(R["name"], rows[-26], same[-1], round(time.time() - t0), "s", flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/he_leg1.csv", index=False)
pd.DataFrame(same).to_csv(f"{OUT}/he_identity.csv", index=False)
print(pd.DataFrame(same).to_string())
