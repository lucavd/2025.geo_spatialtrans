"""S1 (sonda): registration gate Xenium -> H&E -> Visium HD 2 um bins.
Implements docs/sonda/REGOLA_S1.md (commit 29602d4). Stages:
  setup   : hold-out split, bin<->pixel projective, H&E identity, variant choice, regions
  se      : split-half precision per tile size (chooses L, never prints residuals)
  main    : residuals with chosen L, positive and negative controls, local affine
Usage: python S1_registration.py <X1|X2> <stage> [L]
"""
import sys, json, glob, os
import numpy as np, pandas as pd, pyarrow.parquet as pq, h5py

SEED = 20261121
D = "/mnt/micron/geo_spatialtrans/data_lung"
OUT = "/mnt/micron/geo_spatialtrans/S1_sonda"
XPIX = 0.2125
REG_UM = 1000.0

def log(*a):
    print(*a, flush=True)

def paths(tag):
    v = f"{D}/{tag}/visium"; x = f"{D}/{tag}/xenium"
    b2 = glob.glob(f"{v}/binned/**/square_002um", recursive=True)[0]
    b8 = glob.glob(f"{v}/binned/**/square_008um", recursive=True)[0]
    return dict(v=v, x=x, b2=b2, b8=b8,
                btf=glob.glob(f"{v}/*tissue_image.btf")[0],
                he=glob.glob(f"{x}/*he_image.ome.tif")[0],
                align=glob.glob(f"{x}/*he_imagealignment.csv")[0],
                panel=glob.glob(f"{x}/*gene_panel.json")[0],
                cells=glob.glob(f"{x}/*_cells.parquet")[0],
                tx=glob.glob(f"{x}/*transcripts.parquet")[0],
                h5=f"{b2}/filtered_feature_bc_matrix.h5")

def positions(bdir):
    tp = pq.read_table(f"{bdir}/spatial/tissue_positions.parquet").to_pandas()
    return tp

def fit_homography(src, dst):
    # dst ~ H src (DLT, normalised)
    def norm(p):
        m = p.mean(0); s = np.sqrt(2) / np.mean(np.linalg.norm(p - m, axis=1))
        T = np.array([[s, 0, -s * m[0]], [0, s, -s * m[1]], [0, 0, 1]])
        return T, (np.c_[p, np.ones(len(p))] @ T.T)[:, :2]
    Ts, s = norm(src); Td, d = norm(dst)
    A = []
    for (x, y), (u, v) in zip(s, d):
        A.append([-x, -y, -1, 0, 0, 0, u * x, u * y, u])
        A.append([0, 0, 0, -x, -y, -1, v * x, v * y, v])
    _, _, Vt = np.linalg.svd(np.asarray(A), full_matrices=False)
    Hn = Vt[-1].reshape(3, 3)
    H = np.linalg.inv(Td) @ Hn @ Ts
    return H / H[2, 2]

def apply_h(H, p):
    q = np.c_[p, np.ones(len(p))] @ H.T
    return q[:, :2] / q[:, 2:3]

# ---------------------------------------------------------------- setup
def stage_setup(tag):
    P = paths(tag); os.makedirs(f"{OUT}/{tag}", exist_ok=True); res = {"tag": tag}
    tp8 = positions(P["b8"]); tp2 = positions(P["b2"])
    # hold-out split (8 um in-tissue mask only)
    x8 = (tp8.array_col.values + 0.5) * 8.0
    m = float(np.median(x8[tp8.in_tissue.values == 1]))
    res["split_m_um"] = m
    # bin -> pixel projective on 2 um bins
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(tp2), size=min(200000, len(tp2)), replace=False)
    src = tp2[["array_col", "array_row"]].values[idx].astype(float)
    dst = tp2[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values[idx].astype(float)
    H = fit_homography(src, dst)
    allsrc = tp2[["array_col", "array_row"]].values.astype(float)
    alldst = tp2[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    err = np.linalg.norm(apply_h(H, allsrc) - alldst, axis=1)
    res["bin2px_max_resid_px"] = float(err.max()); res["bin2px_p99_resid_px"] = float(np.percentile(err, 99))
    # 8 um centres vs mean of 2 um centres (aggregation floor(i/4))
    g = tp2.assign(c8=tp2.array_col // 4, r8=tp2.array_row // 4).groupby(["c8", "r8"])[["pxl_col_in_fullres", "pxl_row_in_fullres"]].mean()
    g.index.names = ["array_col", "array_row"]
    j = tp8.set_index(["array_col", "array_row"])[["pxl_col_in_fullres", "pxl_row_in_fullres"]].join(g, how="inner", rsuffix="_2")
    res["agg8_max_px"] = float(np.abs(j.iloc[:, :2].values - j.iloc[:, 2:].values).max())
    sf = json.load(open(f"{P['b2']}/spatial/scalefactors_json.json"))
    res["microns_per_pixel"] = sf.get("microns_per_pixel")
    # H&E identity
    res["he"] = he_identity(P)
    # Xenium -> H&E variants
    M = np.loadtxt(P["align"], delimiter=",")
    if M.shape != (3, 3):
        M = M.reshape(3, 3)
    cells = pq.read_table(P["cells"], columns=["x_centroid", "y_centroid"]).to_pandas().values
    intiss = set(map(tuple, tp8.loc[tp8.in_tissue == 1, ["array_col", "array_row"]].values))
    Hinv = np.linalg.inv(H)
    A_he = np.array(res["he"]["affine_x_to_v"]) if not res["he"]["identical"] else np.eye(3)
    var = {}
    for dname, MM in [("direct", M), ("inverse", np.linalg.inv(M))]:
        for uname, sc in [("px0.2125", 1 / XPIX), ("um", 1.0)]:
            px = apply_h(A_he @ MM, cells * sc)
            a = apply_h(Hinv, px)
            k8 = np.floor(np.round(a) / 4).astype(int)
            f = np.mean([tuple(t) in intiss for t in k8])
            var[f"{dname}|{uname}"] = float(f)
    res["variants"] = var
    order = sorted(var, key=var.get, reverse=True)
    best, second = order[0], order[1]
    res["variant"] = best
    res["variant_ok"] = bool(var[best] >= 0.8 and var[best] - var[second] >= 0.2)
    T = composite(M, H, A_he, best)
    np.save(f"{OUT}/{tag}/composite_parts.npy", np.stack([M, H, A_he]))
    # regions
    a = apply_t(T, cells)
    res.update(choose_regions(tp8, a, m))
    json.dump(res, open(f"{OUT}/{tag}/setup.json", "w"), indent=1)
    log(json.dumps({k: v for k, v in res.items() if k != "candidates"}, indent=1))

def composite(M, H, A_he, variant):
    d, u = variant.split("|")
    MM = M if d == "direct" else np.linalg.inv(M)
    sc = 1 / XPIX if u.startswith("px") else 1.0
    S = np.diag([sc, sc, 1.0])
    return np.linalg.inv(H) @ A_he @ MM @ S   # Xenium um -> continuous 2um array (col,row)

def apply_t(T, p):
    return apply_h(T, p)

def load_T(tag):
    M, H, A = np.load(f"{OUT}/{tag}/composite_parts.npy")
    s = json.load(open(f"{OUT}/{tag}/setup.json"))
    return composite(M, H, A, s["variant"]), s

def he_identity(P):
    import tifffile
    from skimage.registration import phase_cross_correlation
    out = {}
    with tifffile.TiffFile(P["btf"]) as tv, tifffile.TiffFile(P["he"]) as tx:
        sv = tv.series[0].levels[0].shape; sx = tx.series[0].levels[0].shape
        out["shape_v"] = list(sv); out["shape_x"] = list(sx)
        def opener(tf):
            try:
                import zarr
                a = zarr.open(tf.series[0].levels[0].aszarr(), mode="r")
                _ = np.asarray(a[(slice(0, 2),) * (len(a.shape) - 1)])
                return a, "zarr"
            except Exception as e:
                log("zarr read failed, full read:", repr(e)[:200])
                return tf.series[0].levels[0].asarray(), "full"
        av, mv = opener(tv); ax, mx = opener(tx); out["read_mode"] = [mv, mx]
        def gray(a, r0, c0, n):
            sh = a.shape
            if len(sh) == 3 and sh[0] == 3:
                t = np.asarray(a[:, r0:r0 + n, c0:c0 + n]).astype(float).mean(0)
            elif len(sh) == 3:
                t = np.asarray(a[r0:r0 + n, c0:c0 + n, :]).astype(float).mean(2)
            else:
                t = np.asarray(a[r0:r0 + n, c0:c0 + n]).astype(float)
            return t
        hv, wv = (sv[0], sv[1]) if sv[-1] == 3 else (sv[-2], sv[-1])
        hx, wx = (sx[0], sx[1]) if sx[-1] == 3 else (sx[-2], sx[-1])
        n = 2048; pts = []; shifts = []
        for fr in (0.3, 0.5, 0.7):
            for fc in (0.3, 0.5, 0.7):
                r0 = int(min(hv, hx) * fr) - n // 2; c0 = int(min(wv, wx) * fc) - n // 2
                gv = gray(av, r0, c0, n); gx = gray(ax, r0, c0, n)
                if gv.std() < 1 or gx.std() < 1:
                    continue
                sh, _, _ = phase_cross_correlation(gv, gx, upsample_factor=10)
                shifts.append([float(sh[1]), float(sh[0])]); pts.append([c0 + n / 2, r0 + n / 2])
        out["crop_shifts_colrow"] = shifts
        same = (hv, wv) == (hx, wx) and len(shifts) >= 5 and max(np.hypot(*s) for s in shifts) <= 1.0
        out["identical"] = bool(same)
        if not same and len(shifts) >= 3:
            p = np.asarray(pts); q = p + np.asarray(shifts)   # x-image point p maps to v-image q
            Aaff = np.linalg.lstsq(np.c_[p, np.ones(len(p))], q, rcond=None)[0].T
            out["affine_x_to_v"] = np.vstack([Aaff, [0, 0, 1]]).tolist()
        else:
            out["affine_x_to_v"] = np.eye(3).tolist()
    return out

def choose_regions(tp8, a, m):
    # covered: in tissue + >=1 Xenium cell centroid in 3x3 neighbourhood of 8um bins
    nc = int(tp8.array_col.max()) + 1; nr = int(tp8.array_row.max()) + 1
    occ = np.zeros((nr, nc), bool)
    k = np.floor(np.round(a) / 4).astype(int)
    ok = (k[:, 0] >= 0) & (k[:, 0] < nc) & (k[:, 1] >= 0) & (k[:, 1] < nr)
    occ[k[ok, 1], k[ok, 0]] = True
    from scipy.ndimage import binary_dilation
    occ3 = binary_dilation(occ, structure=np.ones((3, 3), bool))
    tis = np.zeros((nr, nc), bool)
    t = tp8[tp8.in_tissue == 1]; tis[t.array_row.values, t.array_col.values] = True
    cov = (tis & occ3).astype(float)
    cs = np.pad(cov.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    xmax_dev = m - 50.0
    out = {}
    for side, thr in [(REG_UM, 0.9), (REG_UM, 0.8), (750.0, 0.8)]:
        nb = int(round(side / 8)); cands = []
        for x0 in np.arange(0, nc * 8 - side + 1e-9, 100.0):
            if x0 + side > xmax_dev: break
            c0 = int(round(x0 / 8))
            for y0 in np.arange(0, nr * 8 - side + 1e-9, 100.0):
                r0 = int(round(y0 / 8))
                s = cs[r0 + nb, c0 + nb] - cs[r0, c0 + nb] - cs[r0 + nb, c0] + cs[r0, c0]
                f = s / nb ** 2
                if f >= thr:
                    cands.append((x0, y0, f))
        if len(cands) < 3:
            continue
        C = np.array(cands); ctr = C[:, :2] + side / 2
        rr, cc = np.nonzero(cov > 0)
        dev = (cc + 0.5) * 8 < xmax_dev
        bary = np.array([((cc[dev] + 0.5) * 8).mean(), ((rr[dev] + 0.5) * 8).mean()])
        def disjoint(i, js):
            return all(abs(C[i, 0] - C[j, 0]) >= side or abs(C[i, 1] - C[j, 1]) >= side for j in js)
        def tiebreak(keys):
            return np.lexsort((C[:, 1], C[:, 0], keys))
        d0 = np.linalg.norm(ctr - bary, axis=1)
        r1 = int(tiebreak(d0)[0])
        d1 = -np.linalg.norm(ctr - ctr[r1], axis=1)
        r2 = next((int(i) for i in tiebreak(d1) if disjoint(i, [r1])), None)
        if r2 is None: continue
        d2 = -np.minimum(np.linalg.norm(ctr - ctr[r1], axis=1), np.linalg.norm(ctr - ctr[r2], axis=1))
        r3 = next((int(i) for i in tiebreak(d2) if disjoint(i, [r1, r2])), None)
        if r3 is None: continue
        out["region_side_um"] = side; out["coverage_threshold"] = thr; out["n_candidates"] = len(cands)
        out["regions"] = [dict(name=f"R{i+1}", x0=float(C[j, 0]), y0=float(C[j, 1]), cover=float(C[j, 2])) for i, j in enumerate([r1, r2, r3])]
        return out
    out["regions"] = []; out["region_failure"] = "no 3 disjoint candidates"
    return out

# ---------------------------------------------------------------- data for regions
def load_region_data(tag, T, s, margin_um=240.0):
    P = paths(tag)
    tp2 = positions(P["b2"])
    panel = json.load(open(P["panel"]))["payload"]["targets"]
    sym2ens = {}
    for t in panel:
        dd = t.get("type", {}).get("data", {})
        if str(dd.get("id", "")).startswith("ENSG"):
            sym2ens[dd["name"]] = dd["id"]
    with h5py.File(P["h5"], "r") as f:
        g = f["matrix"]
        fid = [x.decode() for x in g["features/id"][:]]
        bcs = [x.decode() for x in g["barcodes"][:]]
        indptr = g["indptr"][:]
        all_ind = g["indices"][:]; all_dat = g["data"][:]
        shared = sorted(set(sym2ens.values()) & set(fid))
        gi = {e: i for i, e in enumerate(shared)}
        frow = np.array([gi.get(e, -1) for e in fid])
        bcpos = tp2.set_index("barcode").loc[bcs, ["array_col", "array_row"]].values.astype(np.int64)
        regs = {}
        Tinv = np.linalg.inv(T)
        for R in s["regions"]:
            side = s["region_side_um"]
            c0 = int(round((R["x0"] - margin_um) / 2)); c1 = int(round((R["x0"] + side + margin_um) / 2))
            r0 = int(round((R["y0"] - margin_um) / 2)); r1 = int(round((R["y0"] + side + margin_um) / 2))
            sel = np.nonzero((bcpos[:, 0] >= c0) & (bcpos[:, 0] < c1) & (bcpos[:, 1] >= r0) & (bcpos[:, 1] < r1))[0]
            lens = (indptr[sel + 1] - indptr[sel]).astype(np.int64)
            if lens.sum():
                starts = np.repeat(indptr[sel].astype(np.int64), lens)
                within = np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens)
                pos = starts + within
                fr = frow[all_ind[pos]]; cols = np.repeat(sel, lens); vals = all_dat[pos]
                k = fr >= 0; rows = fr[k]; cols = cols[k]; vals = vals[k]
            else:
                rows = cols = vals = np.array([], int)
            V = pd.DataFrame(dict(g=rows, c=bcpos[cols, 0] - c0, r=bcpos[cols, 1] - r0, n=vals))
            corners = np.array([[c0, r0], [c1, r0], [c0, r1], [c1, r1]], float)
            xb = apply_h(Tinv, corners)
            regs[R["name"]] = dict(V=V, c0=c0, r0=r0, nc=c1 - c0, nr=r1 - r0, xbox=(xb.min(0) - 30, xb.max(0) + 30))
    # Xenium transcripts, streamed in batches, filtered on qv and on the union of region boxes
    pf = pq.ParquetFile(P["tx"]); acc = {k: [] for k in regs}
    for bt in pf.iter_batches(columns=["feature_name", "x_location", "y_location", "qv"], batch_size=5_000_000):
        df = bt.to_pandas()
        if len(df) and isinstance(df.feature_name.iloc[0], bytes):
            df["feature_name"] = df.feature_name.str.decode("utf-8")
        df = df[df.qv.values >= 20]
        for name, Rg in regs.items():
            lo, hi = Rg["xbox"]
            k = (df.x_location.values >= lo[0]) & (df.x_location.values <= hi[0]) & (df.y_location.values >= lo[1]) & (df.y_location.values <= hi[1])
            if k.any():
                acc[name].append(df.loc[k, ["feature_name", "x_location", "y_location"]])
    for name, Rg in regs.items():
        sub = pd.concat(acc[name]) if acc[name] else pd.DataFrame(columns=["feature_name", "x_location", "y_location"])
        ens = sub.feature_name.astype(str).map(sym2ens)
        k = ens.isin(gi.keys()).values
        sub = sub[k]; ens = ens[k]
        a = apply_h(T, sub[["x_location", "y_location"]].values.astype(float))
        Rg["X"] = pd.DataFrame(dict(g=ens.map(gi).values.astype(int), c=a[:, 0] - Rg["c0"], r=a[:, 1] - Rg["r0"]))
        log(name, "xenium tx", len(Rg["X"]), "visium umi", int(Rg["V"].n.sum()) if len(Rg["V"]) else 0)
    return regs, len(shared)

# ---------------------------------------------------------------- maps and CC
class Pts:
    """points of one layer in a region, sorted by row; coordinates in 2um bins (centres at integers)"""
    def __init__(self, g, c, r, w):
        o = np.argsort(r, kind="stable")
        self.g = np.asarray(g)[o]; self.c = np.asarray(c, float)[o]; self.r = np.asarray(r, float)[o]; self.w = np.asarray(w, np.float32)[o]
    def shifted(self, dc, dr):
        q = Pts.__new__(Pts); q.g, q.w = self.g, self.w; q.c = self.c + dc; q.r = self.r + dr; return q
    def window(self, c_lo, r_lo, n):
        # bins c_lo..c_lo+n-1 ; a point belongs to bin floor(x+0.5)
        i0 = np.searchsorted(self.r, r_lo - 0.5, "left"); i1 = np.searchsorted(self.r, r_lo + n - 0.5, "left")
        c = self.c[i0:i1]; k = (c >= c_lo - 0.5) & (c < c_lo + n - 0.5)
        return self.g[i0:i1][k], (np.floor(c[k] + 0.5) - c_lo).astype(int), (np.floor(self.r[i0:i1][k] + 0.5) - r_lo).astype(int), self.w[i0:i1][k]

def total_map(P, nc, nr):
    ci = np.floor(P.c + 0.5).astype(int); ri = np.floor(P.r + 0.5).astype(int)
    k = (ci >= 0) & (ci < nc) & (ri >= 0) & (ri < nr)
    A = np.zeros((nr, nc), np.float64); np.add.at(A, (ri[k], ci[k]), P.w[k]); return A

def coarse_shift(Xp, Vp, nc, nr, side_b, margin_b):
    X2 = total_map(Xp, nc, nr); V2 = total_map(Vp, nc, nr)
    def agg(t):
        h, w = t.shape; return t[: h // 4 * 4, : w // 4 * 4].reshape(h // 4, 4, w // 4, 4).sum((1, 3))
    X8 = agg(X2); V8 = agg(V2)
    m8 = margin_b // 4; s8 = side_b // 4; win = 25   # +-200 um
    v = V8[m8:m8 + s8, m8:m8 + s8]; v = (v - v.mean()) / (v.std() + 1e-9)
    best = (-np.inf, 0, 0)
    for dy in range(-win, win + 1):
        for dx in range(-win, win + 1):
            x = X8[m8 + dy:m8 + dy + s8, m8 + dx:m8 + dx + s8]
            if x.shape != v.shape: continue
            x = (x - x.mean()) / (x.std() + 1e-9)
            cc = float((x * v).mean())
            if cc > best[0]: best = (cc, dx, dy)
    return best[1] * 4, best[2] * 4, best[0]

W = 6   # +-12 um in 2um bins

def tile_cc(Xp, Vp, c0, r0, nb, cx, cy, min_gene=10, min_tot=500):
    vg, vc, vr, vw = Vp.window(c0, r0, nb)
    xg, xc, xr, xw = Xp.window(c0 + cx - W, r0 + cy - W, nb + 2 * W)
    core = (xc >= W) & (xc < W + nb) & (xr >= W) & (xr < W + nb)
    ng = int(max(vg.max(initial=-1), xg.max(initial=-1))) + 1
    if ng <= 0: return None
    nv = np.bincount(vg, weights=vw, minlength=ng); nx = np.bincount(xg[core], weights=xw[core], minlength=ng)
    keep = (nv >= min_gene) & (nx >= min_gene)
    if keep.sum() == 0 or nv[keep].sum() < min_tot or nx[keep].sum() < min_tot: return None
    gmap = -np.ones(ng, int); gmap[np.nonzero(keep)[0]] = np.arange(keep.sum()); G = int(keep.sum())
    kv = keep[vg]; kx = keep[xg]
    Va = np.zeros((G, nb, nb), np.float32); np.add.at(Va, (gmap[vg[kv]], vr[kv], vc[kv]), vw[kv])
    Xa = np.zeros((G, nb + 2 * W, nb + 2 * W), np.float32); np.add.at(Xa, (gmap[xg[kx]], xr[kx], xc[kx]), xw[kx])
    Vz = (Va - Va.mean((1, 2), keepdims=True)) / (Va.std((1, 2), keepdims=True) + 1e-9)
    Xz = (Xa - Xa.mean((1, 2), keepdims=True)) / (Xa.std((1, 2), keepdims=True) + 1e-9)
    Vf = Vz.reshape(G, -1)
    cc = np.empty((2 * W + 1, 2 * W + 1))
    for iy in range(2 * W + 1):
        for ix in range(2 * W + 1):
            cc[iy, ix] = np.einsum("ij,ij->", np.ascontiguousarray(Xz[:, iy:iy + nb, ix:ix + nb]).reshape(G, -1), Vf)
    iy, ix = np.unravel_index(np.argmax(cc), cc.shape)
    if iy in (0, 2 * W) or ix in (0, 2 * W):
        return dict(edge=True, dx=float(ix - W), dy=float(iy - W), genes=G)
    nbh = cc[iy - 1:iy + 2, ix - 1:ix + 2]
    yy, xx = np.mgrid[-1:2, -1:2]
    Aq = np.c_[np.ones(9), xx.ravel(), yy.ravel(), xx.ravel() ** 2, yy.ravel() ** 2, (xx * yy).ravel()]
    b = np.linalg.lstsq(Aq, nbh.ravel(), rcond=None)[0]
    Hq = np.array([[2 * b[3], b[5]], [b[5], 2 * b[4]]])
    try:
        off = -np.linalg.solve(Hq, b[1:3])
    except np.linalg.LinAlgError:
        off = np.zeros(2)
    if np.any(np.abs(off) > 1) or not np.all(np.linalg.eigvalsh(Hq) < 0):
        off = np.zeros(2)
    return dict(edge=False, dx=float(ix - W + off[0]), dy=float(iy - W + off[1]), genes=G)

def tiles_run(Xp, Vp, side_um, L_um, margin_b, coarse):
    nb = int(round(L_um / 2)); nt = int(round(side_um / L_um)); out = []
    cx, cy = coarse
    for ty in range(nt):
        for tx in range(nt):
            r0 = margin_b + ty * nb; c0 = margin_b + tx * nb
            rr = tile_cc(Xp, Vp, c0, r0, nb, cx, cy)
            d = dict(tx=tx, ty=ty, cx_um=(tx + 0.5) * L_um, cy_um=(ty + 0.5) * L_um)
            if rr is None:
                d.update(valid=False, edge=np.nan)
            else:
                d.update(valid=True, edge=rr["edge"], genes=rr["genes"], sx_um=2 * (cx + rr["dx"]), sy_um=2 * (cy + rr["dy"]),
                         lx_um=2 * rr["dx"], ly_um=2 * rr["dy"])
            out.append(d)
    return pd.DataFrame(out)

def layer_pts(Rg, xmask=None, vw=None):
    X = Rg["X"] if xmask is None else Rg["X"][xmask]
    Xp = Pts(X.g.values, X.c.values, X.r.values, np.ones(len(X), np.float32))
    V = Rg["V"]; w = V.n.values if vw is None else vw
    k = np.asarray(w) > 0
    Vp = Pts(V.g.values[k], V.c.values[k].astype(float), V.r.values[k].astype(float), np.asarray(w)[k])
    return Xp, Vp

def prep(tag):
    T, s = load_T(tag)
    regs, ng = load_region_data(tag, T, s)
    margin_b = int(round(240.0 / 2)); side_b = int(round(s["region_side_um"] / 2))
    return T, s, regs, ng, margin_b, side_b

def okmask(a, b):
    k = a.valid.values.astype(bool) & b.valid.values.astype(bool)
    ea = a.edge.astype(float).fillna(1).values.astype(bool); eb = b.edge.astype(float).fillna(1).values.astype(bool)
    return k & ~ea & ~eb

def stage_se(tag):
    T, s, regs, ng, mb, sb = prep(tag)
    rng = np.random.default_rng(SEED); rows = []; coarse = {}
    for name, Rg in regs.items():
        Xp, Vp = layer_pts(Rg)
        cx, cy, ccv = coarse_shift(Xp, Vp, Rg["nc"], Rg["nr"], sb, mb); coarse[name] = (int(cx), int(cy), float(ccv))
        log(name, "coarse shift (bins)", cx, cy, "cc", round(ccv, 3))
        xm = rng.random(len(Rg["X"])) < 0.5
        vA = rng.binomial(Rg["V"].n.values.astype(int), 0.5); vB = Rg["V"].n.values - vA
        XA, VA = layer_pts(Rg, xmask=xm, vw=vA); XB, VB = layer_pts(Rg, xmask=~xm, vw=vB)
        for L in (100.0, 200.0, 250.0):
            a = tiles_run(XA, VA, s["region_side_um"], L, mb, (cx, cy)); b = tiles_run(XB, VB, s["region_side_um"], L, mb, (cx, cy))
            k = okmask(a, b)
            d = np.hypot(a.sx_um[k].values - b.sx_um[k].values, a.sy_um[k].values - b.sy_um[k].values) / 2 if k.sum() else np.array([])
            rows.append(dict(region=name, L=L, n_tiles=len(a), n_used=int(k.sum()), rms_se_um=float(np.sqrt(np.mean(d ** 2))) if k.sum() else np.nan))
            log(rows[-1])
    se = pd.DataFrame(rows); se.to_csv(f"{OUT}/{tag}/se_by_L.csv", index=False)
    Lc = None
    for L in (100.0, 200.0, 250.0):
        q = se[se.L == L]
        if len(q) == 3 and (q.rms_se_um <= 0.5).all():
            Lc = L; break
    json.dump(dict(L=Lc, coarse=coarse), open(f"{OUT}/{tag}/L_choice.json", "w"))
    log(se.to_string(index=False)); log("chosen L:", Lc)

def p95(v):
    return float(np.percentile(v, 95)) if len(v) else np.nan

def stage_main(tag):
    T, s, regs, ng, mb, sb = prep(tag)
    ch = json.load(open(f"{OUT}/{tag}/L_choice.json")); L = ch["L"]
    assert L is not None, "no tile size meets the precision rule"
    names = list(regs); summary = []; alltiles = []; pos = []; neg = []
    pts = {name: layer_pts(regs[name]) for name in names}
    side = s["region_side_um"]
    for i, name in enumerate(names):
        cx, cy, _ = ch["coarse"][name]
        Xp, Vp = pts[name]
        t = tiles_run(Xp, Vp, side, L, mb, (cx, cy)); t["region"] = name
        alltiles.append(t)
        v = t[t.valid.astype(bool)].copy(); e = v.edge.astype(bool).values
        r = np.where(e, 1e6, np.hypot(v.sx_um.values, v.sy_um.values))
        res_loc = np.full(len(v), 1e6); par = ((v.tx + v.ty) % 2).values; good = ~e
        A = np.c_[v.cx_um.values, v.cy_um.values, np.ones(len(v))]
        for f in (0, 1):
            fit = (par == f) & good; ev = par != f
            if fit.sum() >= 6:
                cxy = np.linalg.lstsq(A[fit], v[["sx_um", "sy_um"]].values[fit], rcond=None)[0]
                pr = A @ cxy
                rl = np.hypot(v.sx_um.values - pr[:, 0], v.sy_um.values - pr[:, 1])
                res_loc[ev] = np.where(good[ev], rl[ev], 1e6)
        for inj in [(1.0, 0.0), (0.0, 2.5), (1.5, -1.5)]:
            ti = tiles_run(Xp.shifted(inj[0] / 2, inj[1] / 2), Vp, side, L, mb, (cx, cy))
            k = okmask(t, ti)
            dx = float(np.median(ti.sx_um[k] - t.sx_um[k])); dy = float(np.median(ti.sy_um[k] - t.sy_um[k]))
            pos.append(dict(region=name, inj_x=inj[0], inj_y=inj[1], rec_x=dx, rec_y=dy, n=int(k.sum()),
                            ok=bool(abs(dx - inj[0]) <= 0.3 and abs(dy - inj[1]) <= 0.3)))
        other = names[(i + 1) % len(names)]
        Vo = pts[other][1]
        # Visium of the other region, re-expressed in this region's frame (same tile indices)
        Ro, Rs = regs[other], regs[name]
        Vo2 = Vo.shifted(0, 0)
        tn = tiles_run(Xp, Vo2, side, L, mb, (cx, cy))
        vn = tn[tn.valid.astype(bool)]
        near = np.where(vn.edge.astype(bool).values, False, np.hypot(vn.lx_um.values, vn.ly_um.values) <= 2.0) if len(vn) else np.array([])
        neg.append(dict(region=name, vs=other, n=int(len(vn)), frac_within2=float(near.mean()) if len(vn) else np.nan,
                        ok=bool(len(vn) > 0 and near.mean() <= 0.10)))
        summary.append(dict(region=name, L=L, coarse_x_um=2 * cx, coarse_y_um=2 * cy, n_tiles=len(t), n_valid=len(v), n_edge=int(e.sum()),
                            p95_global_um=p95(r), p95_local_um=p95(res_loc), median_global_um=float(np.median(r)) if len(r) else np.nan,
                            within1=float(np.mean(r <= 1)), within2=float(np.mean(r <= 2)), within5=float(np.mean(r <= 5)),
                            loc_within2=float(np.mean(res_loc <= 2)), pass_region=bool(p95(r) <= 2 or p95(res_loc) <= 2)))
        log(summary[-1])
    S = pd.DataFrame(summary); Pc = pd.DataFrame(pos); Nc = pd.DataFrame(neg)
    pd.concat(alltiles).to_csv(f"{OUT}/{tag}/tiles.csv", index=False)
    S.to_csv(f"{OUT}/{tag}/summary.csv", index=False); Pc.to_csv(f"{OUT}/{tag}/positive.csv", index=False); Nc.to_csv(f"{OUT}/{tag}/negative.csv", index=False)
    meter_ok = bool(Pc.ok.all() and Nc.ok.all())
    verdict = dict(tag=tag, L=L, meter_ok=meter_ok, section_pass=bool(meter_ok and S.pass_region.all()))
    json.dump(verdict, open(f"{OUT}/{tag}/verdict.json", "w"))
    pd.set_option("display.width", 220)
    log(S.round(3).to_string(index=False)); log(Pc.round(3).to_string(index=False)); log(Nc.round(3).to_string(index=False)); log(verdict)

if __name__ == "__main__":
    tag, stage = sys.argv[1], sys.argv[2]
    {"setup": stage_setup, "se": stage_se, "main": stage_main}[stage](tag)
