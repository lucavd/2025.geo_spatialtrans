"""Synthetic check of the S1 tile cross-correlation: known shift must be recovered."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import S1_registration as m
rng = np.random.default_rng(1)
side = 400
ctr = rng.uniform(0, side, (1800, 2)); typ = rng.integers(0, 5, len(ctr))
prof = rng.gamma(0.5, 1, (5, 30)); prof /= prof.sum(1, keepdims=True)
def sample(nmol, shift, diff):
    ci = rng.integers(0, len(ctr), nmol)
    p = ctr[ci] + rng.normal(0, 2.0, (nmol, 2)) + rng.normal(0, diff, (nmol, 2))
    cum = prof[typ[ci]].cumsum(1); g = (rng.random(nmol)[:, None] > cum).sum(1)
    return g, p[:, 0] + shift[0], p[:, 1] + shift[1]
gx, cx_, rx_ = sample(400000, (0.65, -0.35), 0.0)
gv, cv_, rv_ = sample(150000, (0, 0), 1.0)
key = np.c_[gv, np.floor(cv_ + 0.5), np.floor(rv_ + 0.5)].astype(int)
u, n = np.unique(key, axis=0, return_counts=True)
Xp = m.Pts(gx, cx_, rx_, np.ones(len(gx))); Vp = m.Pts(u[:, 0], u[:, 1].astype(float), u[:, 2].astype(float), n)
ok = True
for L in (100.0, 200.0):
    t = m.tiles_run(Xp, Vp, 600.0, L, 60, (0, 0)); v = t[t.valid.astype(bool)]
    mx, my = v.sx_um.median(), v.sy_um.median()
    print(f"L={L} tiles={len(t)} valid={len(v)} edge={int(v.edge.astype(bool).sum())} median=({mx:.3f},{my:.3f}) um true=(1.300,-0.700) sd=({v.sx_um.std():.3f},{v.sy_um.std():.3f})")
    ok &= abs(mx - 1.3) <= 0.3 and abs(my + 0.7) <= 0.3
t2 = m.tiles_run(Xp.shifted(0.75, 0.0), Vp, 600.0, 100.0, 60, (0, 0)); v2 = t2[t2.valid.astype(bool)]
print(f"injected +1.5 um in x: median sx={v2.sx_um.median():.3f} (expected 2.800)")
ok &= abs(v2.sx_um.median() - 2.8) <= 0.3
# negative: Xenium of an unrelated pattern
ctr = rng.uniform(0, side, (1800, 2)); typ = rng.integers(0, 5, len(ctr))
gz, cz, rz = sample(400000, (0, 0), 0.0)
tn = m.tiles_run(m.Pts(gz, cz, rz, np.ones(len(gz))), Vp, 600.0, 100.0, 60, (0, 0)); vn = tn[tn.valid.astype(bool)]
fr = float(np.mean(~vn.edge.astype(bool).values & (np.hypot(vn.lx_um, vn.ly_um) <= 2)))
print(f"negative: frac within 2 um = {fr:.3f} (n={len(vn)})")
ok &= fr <= 0.10
print("RESULT", "PASS" if ok else "FAIL")
