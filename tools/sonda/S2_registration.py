"""S2 (sonda): registration gate on X2 with count-weighted meters. Implements docs/sonda/REGOLA_S2.md.
Meters: S1 (reference, per-gene z-scored CC), A<s> (count-weighted CC), B<s> (Poisson likelihood); s = Gaussian sigma in um.
Stages:
  synth                 : S1 synthetic test for every meter + X2-like sparse diagnostic
  se <TAG>              : split-half precision per meter and L (only SEs and tile counts are written)
  main <TAG> <M> <L>    : residuals, controls, local affine, verdict for one meter (X1: comparison with S1 tiles)
  select                : admission and choice rule -> choice.json
Workers: env S2_WORKERS (default 6).
"""
import sys, os, json
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.ndimage import gaussian_filter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import S1_registration as m

OUT = "/mnt/micron/geo_spatialtrans/S2_sonda"
S1OUT = m.OUT
W = m.W
CANDS = ["A2", "A4", "B2", "B4"]
METERS = ["S1"] + CANDS
EPS = 0.1
NW = int(os.environ.get("S2_WORKERS", "6"))
log = m.log

# ------------------------------------------------------------ similarity surfaces
def tile_data(Xp, Vp, c0, r0, nb, cx, cy, pad, min_gene=10, min_tot=500):
    vg, vc, vr, vw = Vp.window(c0, r0, nb)
    n = nb + 2 * (W + pad)
    xg, xc, xr, xw = Xp.window(c0 + cx - W - pad, r0 + cy - W - pad, n)
    o = W + pad
    core = (xc >= o) & (xc < o + nb) & (xr >= o) & (xr < o + nb)
    ng = int(max(vg.max(initial=-1), xg.max(initial=-1))) + 1
    if ng <= 0: return None
    nv = np.bincount(vg, weights=vw, minlength=ng); nx = np.bincount(xg[core], weights=xw[core], minlength=ng)
    keep = (nv >= min_gene) & (nx >= min_gene)
    if keep.sum() == 0 or nv[keep].sum() < min_tot or nx[keep].sum() < min_tot: return None
    gmap = -np.ones(ng, int); gmap[np.nonzero(keep)[0]] = np.arange(keep.sum()); G = int(keep.sum())
    kv = keep[vg]; kx = keep[xg]
    Xa = np.zeros((G, n, n), np.float64); np.add.at(Xa, (gmap[xg[kx]], xr[kx], xc[kx]), xw[kx])
    return dict(G=G, vg=gmap[vg[kv]], vr=vr[kv], vc=vc[kv], vw=vw[kv].astype(np.float64), Xa=Xa, nb=nb, pad=pad)

def surface(d, kind, sig_b):
    Xs = gaussian_filter(d["Xa"], sigma=(0, sig_b, sig_b), mode="constant", truncate=3.0)
    p = d["pad"]
    if p: Xs = Xs[:, p:-p, p:-p]
    nb, G = d["nb"], d["G"]
    if kind == "B":
        Xs = Xs + EPS * Xs.mean((1, 2), keepdims=True)
        F = np.log(Xs)
    else:
        F = Xs
    I = np.pad(Xs.cumsum(1).cumsum(2), ((0, 0), (1, 0), (1, 0)))
    Ng = np.bincount(d["vg"], weights=d["vw"], minlength=G)
    vg, vr, vc, vw = d["vg"], d["vr"], d["vc"], d["vw"]
    cc = np.empty((2 * W + 1, 2 * W + 1))
    for iy in range(2 * W + 1):
        for ix in range(2 * W + 1):
            t1 = float(np.dot(vw, F[vg, vr + iy, vc + ix]))
            S = I[:, iy + nb, ix + nb] - I[:, iy, ix + nb] - I[:, iy + nb, ix] + I[:, iy, ix]
            t2 = float(np.dot(Ng, np.log(S))) if kind == "B" else float(np.dot(Ng / nb ** 2, S))
            cc[iy, ix] = t1 - t2
    return cc

def refine(cc):
    # identical to the peak step of S1_registration.tile_cc
    iy, ix = np.unravel_index(np.argmax(cc), cc.shape)
    if iy in (0, 2 * W) or ix in (0, 2 * W):
        return dict(edge=True, dx=float(ix - W), dy=float(iy - W))
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
    return dict(edge=False, dx=float(ix - W + off[0]), dy=float(iy - W + off[1]))

def tile_shift(meter, Xp, Vp, c0, r0, nb, cx, cy):
    if meter == "S1":
        return m.tile_cc(Xp, Vp, c0, r0, nb, cx, cy)
    kind, sb = meter[0], float(meter[1:]) / 2.0
    pad = int(np.ceil(3 * sb))
    d = tile_data(Xp, Vp, c0, r0, nb, cx, cy, pad)
    if d is None: return None
    r = refine(surface(d, kind, sb)); r["genes"] = d["G"]; return r

def tiles_run(meter, Xp, Vp, side_um, L_um, margin_b, coarse):
    nb = int(round(L_um / 2)); nt = int(round(side_um / L_um)); out = []
    cx, cy = coarse
    for ty in range(nt):
        for tx in range(nt):
            r0 = margin_b + ty * nb; c0 = margin_b + tx * nb
            rr = tile_shift(meter, Xp, Vp, c0, r0, nb, cx, cy)
            d = dict(tx=tx, ty=ty, cx_um=(tx + 0.5) * L_um, cy_um=(ty + 0.5) * L_um)
            if rr is None:
                d.update(valid=False, edge=np.nan)
            else:
                d.update(valid=True, edge=rr["edge"], genes=rr["genes"], sx_um=2 * (cx + rr["dx"]), sy_um=2 * (cy + rr["dy"]),
                         lx_um=2 * rr["dx"], ly_um=2 * rr["dy"])
            out.append(d)
    return pd.DataFrame(out)

# ------------------------------------------------------------ synthetic
def synth_data(sparse=False):
    rng = np.random.default_rng(1)
    side = 400
    ctr = rng.uniform(0, side, (1800, 2)); typ = rng.integers(0, 5, len(ctr))
    ngen = 2000 if sparse else 30
    prof = rng.gamma(0.3 if sparse else 0.5, 1, (5, ngen)); prof /= prof.sum(1, keepdims=True)
    def sample(nmol, shift, diff, ctr, typ):
        ci = rng.integers(0, len(ctr), nmol)
        p = ctr[ci] + rng.normal(0, 2.0, (nmol, 2)) + rng.normal(0, diff, (nmol, 2))
        cum = prof[typ[ci]].cumsum(1); g = np.minimum((rng.random(nmol)[:, None] > cum).sum(1), ngen - 1)
        return g, p[:, 0] + shift[0], p[:, 1] + shift[1]
    nx, nv = (1500000, 700000) if sparse else (400000, 150000)
    gx, cx_, rx_ = sample(nx, (0.65, -0.35), 0.0, ctr, typ)
    gv, cv_, rv_ = sample(nv, (0, 0), 1.0, ctr, typ)
    key = np.c_[gv, np.floor(cv_ + 0.5), np.floor(rv_ + 0.5)].astype(int)
    u, n = np.unique(key, axis=0, return_counts=True)
    Xp = m.Pts(gx, cx_, rx_, np.ones(len(gx))); Vp = m.Pts(u[:, 0], u[:, 1].astype(float), u[:, 2].astype(float), n)
    ctr2 = rng.uniform(0, side, (1800, 2)); typ2 = rng.integers(0, 5, len(ctr2))
    gz, cz, rz = sample(nx, (0, 0), 0.0, ctr2, typ2)
    return Xp, Vp, m.Pts(gz, cz, rz, np.ones(len(gz)))

def synth_one(meter):
    Xp, Vp, Xz = synth_data(False); ok = True; rec = dict(meter=meter)
    for L in (100.0, 200.0):
        t = tiles_run(meter, Xp, Vp, 600.0, L, 60, (0, 0)); v = t[t.valid.astype(bool)]
        mx, my = v.sx_um.median(), v.sy_um.median()
        rec[f"L{int(L)}_med_x"] = mx; rec[f"L{int(L)}_med_y"] = my; rec[f"L{int(L)}_valid"] = len(v)
        ok &= bool(abs(mx - 1.3) <= 0.3 and abs(my + 0.7) <= 0.3)
    t2 = tiles_run(meter, Xp.shifted(0.75, 0.0), Vp, 600.0, 100.0, 60, (0, 0)); v2 = t2[t2.valid.astype(bool)]
    rec["inj_med_x"] = v2.sx_um.median(); ok &= bool(abs(v2.sx_um.median() - 2.8) <= 0.3)
    tn = tiles_run(meter, Xz, Vp, 600.0, 100.0, 60, (0, 0)); vn = tn[tn.valid.astype(bool)]
    fr = float(np.mean(~vn.edge.astype(bool).values & (np.hypot(vn.lx_um, vn.ly_um) <= 2)))
    rec["neg_within2"] = fr; ok &= fr <= 0.10
    rec["pass"] = bool(ok)
    # X2-like sparse diagnostic (no decision depends on it)
    Xs, Vs, _ = synth_data(True)
    for L in (100.0, 200.0):
        t = tiles_run(meter, Xs, Vs, 600.0, L, 60, (0, 0)); v = t[t.valid.astype(bool) & ~t.edge.astype(float).fillna(1).astype(bool)]
        rec[f"sparse_L{int(L)}_n"] = len(v)
        rec[f"sparse_L{int(L)}_rmse"] = float(np.sqrt(np.mean((v.sx_um - 1.3) ** 2 + (v.sy_um + 0.7) ** 2))) if len(v) else np.nan
    return rec

def stage_synth(meters=None):
    meters = meters or METERS
    os.makedirs(OUT, exist_ok=True)
    with Pool(min(NW, len(meters))) as p:
        recs = p.map(synth_one, meters)
    df = pd.DataFrame(recs)
    if meters == METERS: df.to_csv(f"{OUT}/synth.csv", index=False)
    pd.set_option("display.width", 250); log(df.round(3).to_string(index=False))
    return df

# ------------------------------------------------------------ real data
G_ = {}

def _task(args):
    kind, meter, name, L, which = args
    D = G_
    side, mb = D["side"], D["mb"]; cx, cy = D["coarse"][name]
    if kind == "half":
        Xp, Vp = D["halves"][name][which]
        return (meter, name, L, which), tiles_run(meter, Xp, Vp, side, L, mb, (cx, cy))
    Xp, Vp = D["pts"][name]
    if which == "main":
        return (meter, name, L, which), tiles_run(meter, Xp, Vp, side, L, mb, (cx, cy))
    if which.startswith("pos"):
        inj = [(1.0, 0.0), (0.0, 2.5), (1.5, -1.5)][int(which[3])]
        return (meter, name, L, which), tiles_run(meter, Xp.shifted(inj[0] / 2, inj[1] / 2), Vp, side, L, mb, (cx, cy))
    if which == "neg":
        names = D["names"]; other = names[(names.index(name) + 1) % len(names)]
        return (meter, name, L, which), tiles_run(meter, Xp, D["pts"][other][1].shifted(0, 0), side, L, mb, (cx, cy))

def load(tag):
    T, s, regs, ng, mb, sb = m.prep(tag)
    ch = json.load(open(f"{S1OUT}/{tag}/L_choice.json"))
    G_.update(side=s["region_side_um"], mb=mb, names=list(regs), coarse={k: tuple(v[:2]) for k, v in ch["coarse"].items()},
              pts={k: m.layer_pts(R) for k, R in regs.items()})
    return regs, ng

def stage_se(tag, meters=None):
    meters = meters or METERS
    regs, ng = load(tag)
    rng = np.random.default_rng(m.SEED); G_["halves"] = {}
    for name, Rg in regs.items():   # same draw order as S1 stage_se -> identical halves
        xm = rng.random(len(Rg["X"])) < 0.5
        vA = rng.binomial(Rg["V"].n.values.astype(int), 0.5); vB = Rg["V"].n.values - vA
        G_["halves"][name] = (m.layer_pts(Rg, xmask=xm, vw=vA), m.layer_pts(Rg, xmask=~xm, vw=vB))
    tasks = [("half", mt, nm, L, w) for mt in meters for nm in regs for L in (100.0, 200.0, 250.0) for w in (0, 1)]
    with Pool(NW) as p:
        res = dict(p.map(_task, tasks, chunksize=1))
    rows = []
    for mt in meters:
        for nm in regs:
            for L in (100.0, 200.0, 250.0):
                a, b = res[(mt, nm, L, 0)], res[(mt, nm, L, 1)]
                k = m.okmask(a, b)
                d = np.hypot(a.sx_um[k].values - b.sx_um[k].values, a.sy_um[k].values - b.sy_um[k].values) / 2 if k.sum() else np.array([])
                rows.append(dict(meter=mt, region=nm, L=L, n_tiles=len(a), n_used=int(k.sum()),
                                 rms_se_um=float(np.sqrt(np.mean(d ** 2))) if k.sum() else np.nan))
    se = pd.DataFrame(rows); os.makedirs(f"{OUT}/{tag}", exist_ok=True)
    se.to_csv(f"{OUT}/{tag}/se_by_L.csv", index=False)
    log(se.pivot_table(index=["meter", "L"], columns="region", values=["rms_se_um", "n_used"]).round(2).to_string())

def stage_main(tag, meter, L):
    L = float(L); regs, ng = load(tag); names = G_["names"]
    tasks = [("full", meter, nm, L, w) for nm in names for w in ("main", "pos0", "pos1", "pos2", "neg")]
    with Pool(min(NW, len(tasks))) as p:
        res = dict(p.map(_task, tasks, chunksize=1))
    summary, alltiles, pos, neg = [], [], [], []
    for nm in names:
        t = res[(meter, nm, L, "main")].copy(); t["region"] = nm; alltiles.append(t)
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
        for j, inj in enumerate([(1.0, 0.0), (0.0, 2.5), (1.5, -1.5)]):
            ti = res[(meter, nm, L, f"pos{j}")]; k = m.okmask(t, ti)
            dx = float(np.median(ti.sx_um[k] - t.sx_um[k])); dy = float(np.median(ti.sy_um[k] - t.sy_um[k]))
            pos.append(dict(region=nm, inj_x=inj[0], inj_y=inj[1], rec_x=dx, rec_y=dy, n=int(k.sum()),
                            ok=bool(abs(dx - inj[0]) <= 0.3 and abs(dy - inj[1]) <= 0.3)))
        tn = res[(meter, nm, L, "neg")]; vn = tn[tn.valid.astype(bool)]
        near = np.where(vn.edge.astype(bool).values, False, np.hypot(vn.lx_um.values, vn.ly_um.values) <= 2.0) if len(vn) else np.array([])
        other = names[(names.index(nm) + 1) % len(names)]
        neg.append(dict(region=nm, vs=other, n=int(len(vn)), frac_within2=float(near.mean()) if len(vn) else np.nan,
                        ok=bool(len(vn) > 0 and near.mean() <= 0.10)))
        cx, cy = G_["coarse"][nm]
        summary.append(dict(region=nm, L=L, coarse_x_um=2 * cx, coarse_y_um=2 * cy, n_tiles=len(t), n_valid=len(v), n_edge=int(e.sum()),
                            p95_global_um=m.p95(r), p95_local_um=m.p95(res_loc), median_global_um=float(np.median(r)) if len(r) else np.nan,
                            within1=float(np.mean(r <= 1)), within2=float(np.mean(r <= 2)), within5=float(np.mean(r <= 5)),
                            loc_within2=float(np.mean(res_loc <= 2)), pass_region=bool(m.p95(r) <= 2 or m.p95(res_loc) <= 2)))
    S = pd.DataFrame(summary); Pc = pd.DataFrame(pos); Nc = pd.DataFrame(neg); Tt = pd.concat(alltiles)
    od = f"{OUT}/{tag}/{meter}_L{int(L)}"; os.makedirs(od, exist_ok=True)
    Tt.to_csv(f"{od}/tiles.csv", index=False); S.to_csv(f"{od}/summary.csv", index=False)
    Pc.to_csv(f"{od}/positive.csv", index=False); Nc.to_csv(f"{od}/negative.csv", index=False)
    meter_ok = bool(Pc.ok.all() and Nc.ok.all())
    verdict = dict(tag=tag, meter=meter, L=L, meter_ok=meter_ok, section_pass=bool(meter_ok and S.pass_region.all()))
    s1t = f"{S1OUT}/{tag}/tiles.csv"
    if os.path.exists(s1t):
        R = pd.read_csv(s1t); R = R[R.L == L] if "L" in R.columns else R
        J = Tt.merge(R, on=["region", "tx", "ty"], suffixes=("", "_s1"))
        ok = J.valid.astype(bool) & J.valid_s1.astype(bool) & ~J.edge.astype(float).fillna(1).astype(bool) & ~J.edge_s1.astype(float).fillna(1).astype(bool)
        cmp_rows = []
        for nm in names:
            q = J[ok & (J.region == nm)]
            ddx, ddy = q.sx_um - q.sx_um_s1, q.sy_um - q.sy_um_s1
            cmp_rows.append(dict(region=nm, n=len(q), med_dx=float(ddx.median()), med_dy=float(ddy.median()),
                                 rms_d=float(np.sqrt(np.mean(ddx ** 2 + ddy ** 2))), ok=bool(abs(ddx.median()) <= 0.3 and abs(ddy.median()) <= 0.3)))
        Cc = pd.DataFrame(cmp_rows); Cc.to_csv(f"{od}/vs_S1.csv", index=False)
        verdict["reproduces_S1"] = bool(Cc.ok.all()); log(Cc.round(3).to_string(index=False))
    json.dump(verdict, open(f"{od}/verdict.json", "w"))
    pd.set_option("display.width", 250)
    log(S.round(3).to_string(index=False)); log(Pc.round(3).to_string(index=False)); log(Nc.round(3).to_string(index=False)); log(verdict)

def stage_select():
    sy = pd.read_csv(f"{OUT}/synth.csv").set_index("meter")
    s1v = json.load(open(f"{S1OUT}/X1/verdict.json"))
    se = pd.read_csv(f"{OUT}/X2/se_by_L.csv")
    adm, why, rows = [], {}, []
    for c in CANDS:
        vp = f"{OUT}/X1/{c}_L200/verdict.json"
        if not os.path.exists(vp): why[c] = "X1 main missing"; continue
        v = json.load(open(vp))
        conds = dict(synth=bool(sy.loc[c, "pass"]), x1_controls=v["meter_ok"], x1_reproduces=v.get("reproduces_S1", False),
                     x1_same_verdict=v["section_pass"] == s1v["section_pass"])
        why[c] = conds
        if all(conds.values()): adm.append(c)
    for c in adm:
        q = se[se.meter == c]; Ls = None
        for L in (100.0, 200.0, 250.0):
            z = q[q.L == L]
            if len(z) == 3 and (z.rms_se_um <= 0.5).all() and (z.n_used >= 0.5 * z.n_tiles).all():
                Ls = L; break
        mx = float(q[q.L == Ls].rms_se_um.max()) if Ls else np.nan
        rows.append(dict(meter=c, L_star=Ls, max_se=mx))
    R = pd.DataFrame(rows, columns=["meter", "L_star", "max_se"]); R2 = R.dropna(subset=["L_star"])
    choice = None
    if len(R2):
        R2 = R2.assign(kb=[0 if c[0] == "B" else 1 for c in R2.meter], sg=[float(c[1:]) for c in R2.meter])
        choice = R2.sort_values(["L_star", "max_se", "kb", "sg"]).iloc[0][["meter", "L_star"]].to_dict()
    out = dict(admission=why, admitted=adm, ranking=R.to_dict("records"), choice=choice)
    json.dump(out, open(f"{OUT}/choice.json", "w"), indent=1); log(json.dumps(out, indent=1))

if __name__ == "__main__":
    st = sys.argv[1]
    if st == "synth": stage_synth(sys.argv[2].split(",") if len(sys.argv) > 2 else None)
    elif st == "se": stage_se(sys.argv[2], sys.argv[3].split(",") if len(sys.argv) > 3 else None)
    elif st == "main": stage_main(sys.argv[2], sys.argv[3], sys.argv[4])
    elif st == "select": stage_select()
