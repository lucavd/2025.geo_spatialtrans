"""R6 gate G-1: tests of footprint_labels / aggregate_counts on a synthetic label image, plus 3 mutants.
Run: .venv/bin/python tools/R6_test.py  -> PASS/FAIL per assertion, exit 1 on any FAIL."""
import sys, numpy as np, scipy.sparse as sp
sys.path.insert(0, "tools"); import R6_common as R
res = []
def check(name, cond): res.append((name, bool(cond))); print(("PASS" if cond else "FAIL"), name)

# synthetic ROI: 200x200 px at 0.5 um/px; two discs (labels 1, 2) and one far disc (label 3)
um = 0.5; H = W = 200
yy, xx = np.mgrid[0:H, 0:W]; mask = np.zeros((H, W), dtype=np.uint32)
mask[(yy-50)**2 + (xx-50)**2 <= 10**2] = 1      # radius 10 px = 5 um
mask[(yy-50)**2 + (xx-75)**2 <= 8**2] = 2       # neighbour, edge gap 7 px = 3.5 um
mask[(yy-150)**2 + (xx-150)**2 <= 6**2] = 3
# bins every 4 px (2 um) at pixel centres
rows, cols = [a.ravel() for a in np.mgrid[2:H:4, 2:W:4]]
d_exp = 3.0
lab_nuc = R.footprint_labels(mask, rows, cols, um, d_exp, "nuc")
lab_exp = R.footprint_labels(mask, rows, cols, um, d_exp, "exp")
check("T1 nuc: label equals pixel label", np.array_equal(lab_nuc, mask[rows, cols]))
check("T2 nuc subset of exp", np.all((lab_nuc == 0) | (lab_nuc == lab_exp)))
from scipy import ndimage
dist, (ir, ic) = ndimage.distance_transform_edt(mask == 0, return_indices=True)
d_um = dist[rows, cols] * um
check("T3 exp: ambient bins beyond d_exp stay 0", np.all(lab_exp[(lab_nuc == 0) & (d_um > d_exp)] == 0))
check("T4 exp: ambient bins within d_exp get nearest label", np.all(lab_exp[(lab_nuc == 0) & (d_um <= d_exp)] == mask[ir[rows, cols], ic[rows, cols]][(lab_nuc == 0) & (d_um <= d_exp)]))
check("T5 exp: more bins assigned than nuc", (lab_exp > 0).sum() > (lab_nuc > 0).sum())
# expected counts by brute force: bin inside disc i, or within d_exp of disc i and nearer to i than j
def brute(cx, cy, r):
    d = np.hypot(cols - cx, rows - cy); return d
d1, d2, d3 = brute(50, 50, 10), brute(75, 50, 8), brute(150, 150, 6)
e1 = (d1 - 10) * um; e2 = (d2 - 8) * um; e3 = (d3 - 6) * um      # edge distance in um (approx. for circles)
exp_brute = np.zeros_like(lab_exp)
near = np.argmin(np.c_[e1, e2, e3], axis=1); emin = np.min(np.c_[e1, e2, e3], axis=1)
exp_brute[emin <= d_exp] = near[emin <= d_exp] + 1
agree = np.mean(exp_brute == lab_exp)
check(f"T6 exp vs brute-force circles agree on >= 98% of bins ({agree:.3f})", agree >= 0.98)
# out-of-image bins -> 0
lab_oob = R.footprint_labels(mask, np.array([-1, 0, H]), np.array([0, -5, 0]), um, d_exp, "exp")
check("T7 out-of-image bins are ambient", np.all(lab_oob == 0))
# shift mode: labels differ from exp and are a rolled version
lab_sh = R.footprint_labels(mask, rows, cols, um, d_exp, "shift", shift_um=20.0)
s = int(round(20.0 / um)); m2 = np.roll(np.roll(mask, s, 0), s, 1)
check("T8 shift equals exp on rolled mask", np.array_equal(lab_sh, R.footprint_labels(m2, rows, cols, um, d_exp, "exp")))
# aggregation: random counts, labels 1..3 kept, check sums
rng = np.random.default_rng(1); G = 7
Xb = sp.random(G, len(rows), density=0.3, random_state=1, data_rvs=lambda n: rng.integers(1, 5, n)).tocsc()
C, nb, umi = R.aggregate_counts(Xb, lab_exp, np.array([1, 2, 3]))
check("T9 aggregate shape nuclei x genes", C.shape == (3, G))
check("T10 n_bins per nucleus equals label counts", all(nb[i] == (lab_exp == i + 1).sum() for i in range(3)))
check("T11 UMI conserved: sum over nuclei == sum over assigned bins", abs(umi.sum() - Xb[:, lab_exp > 0].sum()) < 1e-9)
check("T12 per-nucleus counts equal brute-force sums", all(np.allclose(C[i].toarray().ravel(), np.asarray(Xb[:, lab_exp == i + 1].sum(axis=1)).ravel()) for i in range(3)))
# keep subset: label 2 excluded -> its bins are not reassigned to 1
C2, nb2, umi2 = R.aggregate_counts(Xb, lab_exp, np.array([1, 3]))
check("T13 excluded nucleus bins are dropped, not merged", nb2[0] == nb[0] and nb2[1] == nb[2] and abs(umi2.sum() - (umi[0] + umi[2])) < 1e-9)

# mutants: each must be caught by the suite above (re-run the relevant assertions)
def run_core(fl, ag):
    ok = True
    le = fl(mask, rows, cols, um, d_exp, "exp"); ln = fl(mask, rows, cols, um, d_exp, "nuc")
    ok &= np.all(le[(ln == 0) & (d_um > d_exp)] == 0)
    ok &= np.array_equal(ln, mask[rows, cols])
    Cm, nbm, umim = ag(Xb, le, np.array([1, 2, 3]))
    ok &= all(nbm[i] == (le == i + 1).sum() for i in range(3))
    ok &= abs(umim.sum() - Xb[:, le > 0].sum()) < 1e-9
    ok &= Cm.shape == (3, G)
    return bool(ok)
def M1(mask, rows, cols, um, d, mode="exp", shift_um=20.0):   # ignores d_exp (infinite expansion)
    return R.footprint_labels(mask, rows, cols, um, 1e9, mode, shift_um)
def M2(mask, rows, cols, um, d, mode="exp", shift_um=20.0):   # swaps rows/cols
    return R.footprint_labels(mask, cols, rows, um, d, mode, shift_um)
def M3(Xb, labels, keep):                                      # drops the last nucleus
    C, nb, umi = R.aggregate_counts(Xb, labels, keep[:-1]); return sp.vstack([C, sp.csr_matrix((1, C.shape[1]))]).tocsr(), np.r_[nb, 0], np.r_[umi, 0]
check("M1 (no expansion limit) detected", not run_core(M1, R.aggregate_counts))
check("M2 (rows/cols swapped) detected", not run_core(M2, R.aggregate_counts))
check("M3 (last nucleus dropped) detected", not run_core(R.footprint_labels, M3))
check("baseline passes core", run_core(R.footprint_labels, R.aggregate_counts))
n_fail = sum(not ok for _, ok in res); print(f"{len(res) - n_fail}/{len(res)} PASS"); sys.exit(1 if n_fail else 0)
