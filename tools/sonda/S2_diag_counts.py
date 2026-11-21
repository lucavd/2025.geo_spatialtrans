"""S2 diagnostic (counts only, no shifts): why are X2 tiles invalid? Per tile and split half, genes passing the
per-gene filter (>=10 in both maps), UMI share on those genes, and validity with and without the per-gene filter."""
import sys, os, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import S1_registration as m, S2_registration as s2
OUT = s2.OUT
rows = []
for tag in ("X1", "X2"):
    regs, ng = s2.load(tag)
    rng = np.random.default_rng(m.SEED)
    for name, Rg in regs.items():
        xm = rng.random(len(Rg["X"])) < 0.5
        vA = rng.binomial(Rg["V"].n.values.astype(int), 0.5)
        Xp, Vp = m.layer_pts(Rg, xmask=xm, vw=vA)
        cx, cy = s2.G_["coarse"][name]; mb = s2.G_["mb"]; side = s2.G_["side"]
        for L in (100.0, 200.0, 250.0):
            nb = int(round(L / 2)); nt = int(round(side / L))
            for ty in range(nt):
                for tx in range(nt):
                    r0 = mb + ty * nb; c0 = mb + tx * nb
                    vg, vc, vr, vw = Vp.window(c0, r0, nb)
                    xg, xc, xr, xw = Xp.window(c0 + cx - s2.W, r0 + cy - s2.W, nb + 2 * s2.W)
                    core = (xc >= s2.W) & (xc < s2.W + nb) & (xr >= s2.W) & (xr < s2.W + nb)
                    n_ = int(max(vg.max(initial=-1), xg.max(initial=-1))) + 1
                    nv = np.bincount(vg, weights=vw, minlength=n_); nx = np.bincount(xg[core], weights=xw[core], minlength=n_)
                    keep = (nv >= 10) & (nx >= 10)
                    rows.append(dict(tag=tag, region=name, L=L, tx=tx, ty=ty, v_umi=nv.sum(), x_tx=nx.sum(), genes_any=int(((nv > 0) & (nx > 0)).sum()),
                                     genes_kept=int(keep.sum()), v_umi_kept=nv[keep].sum(), x_tx_kept=nx[keep].sum(),
                                     valid_s1=bool(keep.sum() > 0 and nv[keep].sum() >= 500 and nx[keep].sum() >= 500),
                                     valid_nofilter=bool(nv.sum() >= 500 and nx.sum() >= 500)))
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/diag_counts_half.csv", index=False)
g = D.groupby(["tag", "region", "L"]).agg(tiles=("tx", "size"), valid_s1=("valid_s1", "sum"), valid_nofilter=("valid_nofilter", "sum"),
     v_umi_med=("v_umi", "median"), genes_kept_med=("genes_kept", "median"), share_kept_med=("v_umi_kept", "median"))
g["share_kept_med"] = (D.groupby(["tag", "region", "L"]).apply(lambda q: float(np.median(q.v_umi_kept / q.v_umi.clip(lower=1))))).values
g.to_csv(f"{OUT}/diag_counts_half_summary.csv"); pd.set_option("display.width", 250); print(g.round(2).to_string())
