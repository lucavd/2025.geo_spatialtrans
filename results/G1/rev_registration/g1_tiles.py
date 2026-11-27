"""Per-tile shift estimator (own design, differs from S1):
 - Xenium transcripts shifted continuously (0.25 um raster + summed-area table), re-binned on the 2 um grid;
 - both maps smoothed with a Gaussian (sigma = 2 um); similarity = count-weighted mean over genes of the
   Pearson correlation (weights sqrt(nX_g * nV_g)); S1 instead sums unweighted z-score products on raw bins;
 - two-stage grid argmax (1 um over +-12 um around the region coarse shift, then 0.25 um over +-1 um);
   no quadric refinement. Edge: coarse argmax at +-12 um on either axis.
Convention as in the rule: s maximises sim(X(x + s), V(x)); s in um, x = col, y = row of the 2 um array.
Usage: python g1_tiles.py <R1|R2|R3> <pass> [L_um]
pass: main | pos:<dx>:<dy> | warp | neg_swap | neg_shuffle | halfA | halfB"""
import sys, json, time, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
from g1_common import OUT, MARGIN_BINS, load_chain
SEED_HALF = 4242; SEED_SHUF = 777
W_UM = 12; MB = 3; SIG = 1.0          # window, smoothing margin (bins), sigma in bins (2 um)
MIN_GENE, MIN_TOT = 10, 500

def load_region(name):
    z = np.load(f"{OUT}/region_{name}.npz")
    return {k: z[k] for k in z.files}

def total_map(c, r, w, n):
    ci = np.floor(c + 0.5).astype(int); ri = np.floor(r + 0.5).astype(int)
    k = (ci >= 0) & (ci < n) & (ri >= 0) & (ri < n)
    return np.bincount(ri[k] * n + ci[k], weights=w[k], minlength=n * n).reshape(n, n)

def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum() + 1e-12))

def coarse_shift(R, n_core):
    """region-level coarse shift (integer bins) from smoothed total maps, +-100 um search."""
    n = int(R["n"]); m = MARGIN_BINS
    X = gaussian_filter(total_map(R["xc"], R["xr"], np.ones(len(R["xc"])), n), 2.0)
    V = gaussian_filter(total_map(R["vc"].astype(float), R["vr"].astype(float), R["vn"].astype(float), n), 2.0)
    v = V[m:m + n_core, m:m + n_core]
    def cc(dx, dy):
        return pearson(X[m + dy:m + dy + n_core, m + dx:m + dx + n_core], v)
    best = max(((cc(dx, dy), dx, dy) for dy in range(-48, 49, 4) for dx in range(-48, 49, 4)))
    _, bx, by = best
    best = max(((cc(dx, dy), dx, dy) for dy in range(by - 4, by + 5) for dx in range(bx - 4, bx + 5)))
    return best[1], best[2], best[0]

class Tile:
    def __init__(self, xg, xu, yu, vg, vc, vr, vn, bc0, br0, nb, sc, G_all, shuffle_seed=None):
        self.nb = nb; nm = nb + 2 * MB; self.sc = sc
        # Visium dense map over tile + smoothing margin
        jc = vc - (bc0 - MB); jr = vr - (br0 - MB)
        k = (jc >= 0) & (jc < nm) & (jr >= 0) & (jr < nm)
        Vfull = np.bincount((vg[k] * nm + jr[k]) * nm + jc[k], weights=vn[k], minlength=G_all * nm * nm).reshape(G_all, nm, nm)
        if shuffle_seed is not None:
            p = np.random.default_rng(shuffle_seed).permutation(nm * nm)
            Vfull = Vfull.reshape(G_all, -1)[:, p].reshape(G_all, nm, nm)
        nV = Vfull[:, MB:MB + nb, MB:MB + nb].sum((1, 2))
        # Xenium raster (0.25 um) covering every block requested by offsets in [-13, 13] um
        self.U0x = 2 * (bc0 - MB) - 1 + sc[0] - 13; self.U0y = 2 * (br0 - MB) - 1 + sc[1] - 13
        N = 8 * nm + 4 * 26; self.N = N
        kx = np.floor((xu - self.U0x) / 0.25).astype(np.int64); ky = np.floor((yu - self.U0y) / 0.25).astype(np.int64)
        kk = (kx >= 0) & (kx < N) & (ky >= 0) & (ky < N)
        # Xenium core counts at the coarse shift (gene filter)
        ex0, ey0 = 2 * bc0 - 1 + sc[0], 2 * br0 - 1 + sc[1]
        core = kk & (xu >= ex0) & (xu < ex0 + 2 * nb) & (yu >= ey0) & (yu < ey0 + 2 * nb)
        nX = np.bincount(xg[core], minlength=G_all)
        keep = (nV >= MIN_GENE) & (nX >= MIN_GENE)
        self.share_v = float(nV[keep].sum() / max(nV.sum(), 1)); self.share_x = float(nX[keep].sum() / max(nX.sum(), 1))
        self.nV_tot = float(nV.sum()); self.nX_tot = float(nX.sum())
        self.G = int(keep.sum()); self.nVk = float(nV[keep].sum()); self.nXk = float(nX[keep].sum())
        self.valid = self.G > 0 and self.nVk >= MIN_TOT and self.nXk >= MIN_TOT
        if not self.valid: return
        gi = -np.ones(G_all, np.int64); gi[np.nonzero(keep)[0]] = np.arange(self.G)
        sel = kk & keep[xg]
        flat = (gi[xg[sel]] * N + ky[sel]) * N + kx[sel]
        A = np.bincount(flat, minlength=self.G * N * N).reshape(self.G, N, N).astype(np.int32)
        S = np.zeros((self.G, N + 1, N + 1), np.int32); S[:, 1:, 1:] = A.cumsum(1).cumsum(2); self.S = S
        Vs = gaussian_filter(Vfull[keep], (0, SIG, SIG), mode="constant", truncate=3.0)[:, MB:MB + nb, MB:MB + nb]
        Vs = Vs - Vs.mean((1, 2), keepdims=True)
        self.Vz = Vs / np.sqrt((Vs ** 2).sum((1, 2), keepdims=True) + 1e-12)
        w = np.sqrt(nV[keep] * nX[keep]); self.w = w / w.sum()
    def xmap(self, ox, oy):
        """Xenium counts on the tile+margin grid for offset (ox, oy) um from the coarse shift."""
        qx = int(round(4 * (ox + 13))); qy = int(round(4 * (oy + 13))); nm = self.nb + 2 * MB; S = self.S
        sy0 = slice(qy, qy + 8 * nm, 8); sy1 = slice(qy + 8, qy + 8 + 8 * nm, 8)
        sx0 = slice(qx, qx + 8 * nm, 8); sx1 = slice(qx + 8, qx + 8 + 8 * nm, 8)
        return (S[:, sy1, sx1] - S[:, sy0, sx1] - S[:, sy1, sx0] + S[:, sy0, sx0]).astype(np.float32)
    def sim(self, ox, oy):
        X = gaussian_filter(self.xmap(ox, oy), (0, SIG, SIG), mode="constant", truncate=3.0)[:, MB:MB + self.nb, MB:MB + self.nb]
        X = X - X.mean((1, 2), keepdims=True)
        r = (X * self.Vz).sum((1, 2)) / np.sqrt((X ** 2).sum((1, 2)) + 1e-12)
        return float((self.w * r).sum())
    def search(self):
        g = np.arange(-W_UM, W_UM + 1)
        C = np.array([[self.sim(ox, oy) for ox in g] for oy in g])
        iy, ix = np.unravel_index(np.argmax(C), C.shape); ox, oy = g[ix], g[iy]
        edge = abs(ox) == W_UM or abs(oy) == W_UM
        f = np.arange(-1, 1.0001, 0.25)
        F = np.array([[self.sim(ox + a, oy + b) for a in f] for b in f])
        jy, jx = np.unravel_index(np.argmax(F), F.shape)
        return dict(edge=bool(edge), ox=float(ox + f[jx]), oy=float(oy + f[jy]), peak=float(F.max()),
                    coarse_peak=float(C.max()), median_sim=float(np.median(C)))

def run(name, pas, L):
    chain = load_chain(); regs = [r["name"] for r in chain["regions"]]
    R = load_region(name); side_b = int(round(chain["region_side_um"] / 2))
    cxb, cyb, ccv = coarse_shift(R, side_b); sc = (2 * cxb, 2 * cyb)
    xg = R["xg"].astype(np.int64); xu = 2 * R["xc"]; yu = 2 * R["xr"]
    vg, vc, vr, vn = R["vg"].astype(np.int64), R["vc"].astype(np.int64), R["vr"].astype(np.int64), R["vn"].astype(float)
    G_all = len(json.load(open(f"{OUT}/genes.json"))["genes"]); shuffle = False
    if pas.startswith("pos:"):
        _, dx, dy = pas.split(":"); xu = xu + float(dx); yu = yu + float(dy)
    elif pas == "warp":                       # non-affine: y += 4 sin(2 pi x / 400 um), phase 0 at core edge
        yu = yu + 4.0 * np.sin(2 * np.pi * (xu - (2 * MARGIN_BINS - 1)) / 400.0)
    elif pas == "neg_swap":
        Ro = load_region(regs[(regs.index(name) + 1) % len(regs)])
        vg, vc, vr, vn = Ro["vg"].astype(np.int64), Ro["vc"].astype(np.int64), Ro["vr"].astype(np.int64), Ro["vn"].astype(float)
    elif pas == "neg_shuffle":
        shuffle = True
    elif pas in ("halfA", "halfB"):
        rng = np.random.default_rng(SEED_HALF + regs.index(name))
        hx = rng.random(len(xg)) < 0.5; vA = rng.binomial(vn.astype(np.int64), 0.5).astype(float)
        if pas == "halfA": m = hx; vn2 = vA
        else: m = ~hx; vn2 = vn - vA
        xg, xu, yu = xg[m], xu[m], yu[m]; k = vn2 > 0; vg, vc, vr, vn = vg[k], vc[k], vr[k], vn2[k]
    nb = int(round(L / 2)); nt = int(round(chain["region_side_um"] / L)); rows = []; t0 = time.time()
    for ty in range(nt):
        for tx in range(nt):
            bc0 = MARGIN_BINS + tx * nb; br0 = MARGIN_BINS + ty * nb
            T = Tile(xg, xu, yu, vg, vc, vr, vn, bc0, br0, nb, sc, G_all,
                     shuffle_seed=(SEED_SHUF + 100 * regs.index(name) + ty * nt + tx) if shuffle else None)
            d = dict(region=name, pas=pas, L=L, tx=tx, ty=ty, cx_um=(tx + 0.5) * L, cy_um=(ty + 0.5) * L,
                     coarse_x_um=sc[0], coarse_y_um=sc[1], valid=T.valid, genes=T.G, nV=T.nV_tot, nX=T.nX_tot,
                     nV_kept=getattr(T, "nVk", np.nan), nX_kept=getattr(T, "nXk", np.nan), share_v=T.share_v, share_x=T.share_x)
            if T.valid:
                s = T.search(); d.update(s); d["sx_um"] = sc[0] + s["ox"]; d["sy_um"] = sc[1] + s["oy"]
            rows.append(d)
    df = pd.DataFrame(rows); tag = pas.replace(":", "_")
    df.to_csv(f"{OUT}/tiles/{name}_{tag}_L{int(L)}.csv", index=False)
    print(name, pas, L, "coarse um", sc, "cc", round(ccv, 3), "valid", int(df.valid.sum()), "of", len(df), "time", round(time.time() - t0), "s", flush=True)

if __name__ == "__main__":
    name, pas = sys.argv[1], sys.argv[2]; L = float(sys.argv[3]) if len(sys.argv) > 3 else 200.0
    run(name, pas, L)
