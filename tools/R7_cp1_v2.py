"""tools/R7_cp1_v2.py — correzioni C-1 e C-2 dopo la revisione avversariale: kappa F_box da R7_kappa_v3_all.csv (10^6 repliche) e banda bootstrap con 20 000 estrazioni comuni (prima 4 000). Originale: tools/R7_cp1.py.
 — R7: correlazione dentro le aree (M-4) e previsione P-2 (controprova CP-1).
Campi su tessere di 48 um di una griglia globale nell'unione delle 5 aree: UMI medi per bin di 2 um in tessuto ("umi"),
nuclei Cellpose-RGB per mm^2 di tessuto ("nuclei", ogni nucleo preso dal pezzo disgiunto della sua area).
Semivariogramma fra 48 e 720 um (15 classi) sulle coppie di tessere della stessa area; modello c0 + c1 (1 - exp(-h/L)).
Uso: .venv/bin/python tools/R7_cp1.py [workers]  -> results/R7/cp1/*.csv"""
import os, sys, json, itertools
import numpy as np, pandas as pd, pyarrow.parquet as pq
from scipy.optimize import least_squares
from scipy.spatial.distance import cdist
sys.path.insert(0, "tools")
import R7_common as C
from R4_common import verdict
OUT = "results/R7/cp1_v2"; os.makedirs(OUT, exist_ok=True)
R4RUN = "/mnt/micron/geo_spatialtrans/R4/run1"
ARCH5 = ["A1", "A2", "A3", "A4", "A5"]; AIDX = {a: i for i, a in enumerate(["A1", "A2", "A3", "A4", "A5", "A6"])}
TILE_UM = 48.0; MIN_BINS = 288; NCLS = 15; NREP_B = 20000
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
K3 = pd.read_csv("results/R7/coverage/R7_kappa_v3_all.csv"); KAPPA = K3[K3.matrix == "F_box"].set_index("archetype").kappa_boot.to_dict()

def tiles_for(A):
    g = rois[rois.archetype == A].reset_index(drop=True); ids = list(g.roi_id)
    boxes = g[["c0", "r0", "c1", "r1"]].to_numpy(float); upp = float(g.um_per_px.mean())
    gx, gy, u = [], [], []
    for k, rid in enumerate(ids):
        b = np.load(f"{R4RUN}/roi/{A}_{rid}_bins.npz")
        gx.append(boxes[k, 0] + np.rint(b["x_um"] / g.um_per_px[k])); gy.append(boxes[k, 1] + np.rint(b["y_um"] / g.um_per_px[k])); u.append(b["umi"])
    gx, gy, u = map(np.concatenate, (gx, gy, u))
    key = gx.astype(np.int64) * 10**7 + gy.astype(np.int64); _, first = np.unique(key, return_index=True)
    gx, gy, u = gx[first], gy[first], u[first].astype(float)                     # bin unici
    tpx = TILE_UM / upp
    tx = np.floor(gx / tpx).astype(np.int64); ty = np.floor(gy / tpx).astype(np.int64)
    T = pd.DataFrame(dict(tx=tx, ty=ty, u=u)).groupby(["tx", "ty"]).agg(nb=("u", "size"), us=("u", "sum")).reset_index()
    # nuclei dal pezzo disgiunto
    nuc = pq.read_table("results/R2/R2_nuclei_all.parquet").to_pandas()
    nuc = nuc[(nuc.method == "cellpose_rgb") & (nuc.scale == 1) & nuc.keep & (nuc.archetype == A)]
    nx, ny = [], []
    for k, rid in enumerate(ids):
        d = nuc[nuc.roi_id == rid]; x = boxes[k, 0] + d.x_um.to_numpy() / g.um_per_px[k]; y = boxes[k, 1] + d.y_um.to_numpy() / g.um_per_px[k]
        keep = C.assign_pieces(x, y, boxes) == k; nx.append(x[keep]); ny.append(y[keep])
    nx, ny = np.concatenate(nx), np.concatenate(ny)
    N = pd.DataFrame(dict(tx=np.floor(nx / tpx).astype(np.int64), ty=np.floor(ny / tpx).astype(np.int64))).groupby(["tx", "ty"]).size().rename("nn").reset_index()
    T = T.merge(N, on=["tx", "ty"], how="left").fillna({"nn": 0})
    T = T[T.nb >= MIN_BINS].copy()
    T["umi"] = T.us / T.nb; T["nuclei"] = T.nn / (T.nb * 4.0 * 1e-6)            # nuclei per mm^2 di tessuto in bin
    T["cx"] = (T.tx + 0.5) * TILE_UM; T["cy"] = (T.ty + 0.5) * TILE_UM           # centro in um
    cpx = (T.tx.to_numpy() + 0.5) * tpx; cpy = (T.ty.to_numpy() + 0.5) * tpx
    member = (cpx[:, None] >= boxes[None, :, 0]) & (cpx[:, None] < boxes[None, :, 2]) & (cpy[:, None] >= boxes[None, :, 1]) & (cpy[:, None] < boxes[None, :, 3])
    return T.reset_index(drop=True), member, ids, boxes

EDGES = TILE_UM * (np.arange(NCLS + 1) + 0.5)          # 24, 72, ..., 744 um -> centri 48 ... 720
def area_variograms(T, member, field):
    """Per ogni area: somme dei quadrati delle differenze e numero di coppie per classe (solo coppie della stessa area)."""
    v = T[field].to_numpy(); P = T[["cx", "cy"]].to_numpy(); out = []
    for k in range(member.shape[1]):
        idx = np.where(member[:, k])[0]
        D = cdist(P[idx], P[idx]); dv = (v[idx][:, None] - v[idx][None, :]) ** 2
        iu = np.triu_indices(len(idx), 1); d = D[iu]; s = dv[iu]
        cls = np.digitize(d, EDGES) - 1; ok = (cls >= 0) & (cls < NCLS)
        out.append((np.bincount(cls[ok], weights=s[ok], minlength=NCLS), np.bincount(cls[ok], minlength=NCLS).astype(float)))
    return out

H = TILE_UM * np.arange(1, NCLS + 1)
def fit_vgm(S, Nn):
    gam = S / (2 * Nn)
    def res(p):
        c0, c1, L = p; m = c0 + c1 * (1 - np.exp(-H / L)); return np.sqrt(Nn) * (gam / m - 1)
    p0 = [max(gam[0] * 0.5, 1e-9), max(gam[-1] - gam[0] * 0.5, 1e-9), 150.0]
    r = least_squares(res, p0, bounds=([0, 0, 10], [np.inf, np.inf, 2000]))
    return r.x, gam

def implied_R(T, member, p):
    c0, c1, L = p; P = T[["cx", "cy"]].to_numpy(); n = member.shape[1]; Cv = np.zeros((n, n))
    idx = [np.where(member[:, k])[0] for k in range(n)]
    for i in range(n):
        for j in range(i, n):
            D = cdist(P[idx[i]], P[idx[j]]); M = c1 * np.exp(-D / L); M[D == 0] = c0 + c1
            Cv[i, j] = Cv[j, i] = M.mean()
    d = np.sqrt(np.diag(Cv)); return Cv / np.outer(d, d)

def coverage_under(R, kappa, Z):
    L = np.linalg.cholesky(R + 1e-12 * np.eye(len(R))); X = Z @ L.T; lo, hi = C.exact_boot_ci(X)
    return C.coverage_widened(X, 0.0, kappa, lo, hi)

def run(task):
    A, field = task
    T, member, ids, boxes = tiles_for(A)
    vg = area_variograms(T, member, field)
    S = sum(x[0] for x in vg); Nn = sum(x[1] for x in vg)
    p, gam = fit_vgm(S, Nn); Rp = implied_R(T, member, p)
    kappa = KAPPA[A]
    rng = np.random.default_rng([C.SEED, AIDX[A], 100 + (field == "nuclei")])
    cov_pt = coverage_under(Rp, kappa, rng.standard_normal((C.NREP, 5)))
    Z = rng.standard_normal((NREP_B, 5)); covs, pars = [], []
    for bi in range(C.B):
        pick = rng.integers(0, 5, 5)
        Sb = sum(vg[k][0] for k in pick); Nb = sum(vg[k][1] for k in pick)
        pb, _ = fit_vgm(Sb, Nb); covs.append(coverage_under(implied_R(T, member, pb), kappa, Z)); pars.append(pb)
    covs = np.array(covs); lo, hi = np.percentile(covs, [2.5, 97.5])
    F = C.box_overlap_matrix(boxes)
    pairs = [dict(archetype=A, field=field, i=ids[i], j=ids[j], F_box=F[i, j], R_implied=Rp[i, j]) for i, j in itertools.combinations(range(5), 2)]
    pars = np.array(pars)
    row = dict(archetype=A, field=field, n_tiles=len(T), c0=p[0], c1=p[1], L_um=p[2], L_lo=float(np.percentile(pars[:, 2], 2.5)), L_hi=float(np.percentile(pars[:, 2], 97.5)),
               kappa_Fbox=kappa, coverage=cov_pt, cov_lo=lo, cov_hi=hi, verdict=verdict(lo, hi, *C.BRADLEY),
               rho_bar_implied=C.rho_bar(Rp), rho_bar_F=C.rho_bar(F), kappa_needed_implied=C.calibrate_kappa(Rp, np.random.default_rng([C.SEED, AIDX[A], 200 + (field == "nuclei")])))
    vrows = [dict(archetype=A, field=field, h_um=H[q], gamma=gam[q], n_pairs=Nn[q], model=p[0] + p[1] * (1 - np.exp(-H[q] / p[2]))) for q in range(NCLS)]
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    return row, pairs, vrows

if __name__ == "__main__":
    from multiprocessing import Pool
    tasks = [(a, f) for a in ARCH5 for f in ("umi", "nuclei")]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        res = pool.map(run, tasks)
    pd.DataFrame([r[0] for r in res]).to_csv(f"{OUT}/R7_P2.csv", index=False)
    pd.DataFrame([x for r in res for x in r[1]]).to_csv(f"{OUT}/R7_cp1_pairs.csv", index=False)
    pd.DataFrame([x for r in res for x in r[2]]).to_csv(f"{OUT}/R7_cp1_variograms.csv", index=False)
