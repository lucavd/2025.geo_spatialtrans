
# results/R3/review/claims/rv_interp.py — controprove per le interpretazioni (a) P3, (b) P2, (d) interne, (e) d*.
import numpy as np, pandas as pd, glob
from scipy.spatial import cKDTree
from scipy.stats import gamma
R = "results/R3"; OUT = "/mnt/micron/geo_spatialtrans/R3"; RV = "results/R3/review/claims"
PRIM = dict(A1="cellpose_rgb", A2="cellpose_rgb", A3="cellpose_rgb", A4="spaceranger", A5="cellpose_rgb", A6="spaceranger")
DSTAR = dict(A1=2.5, A2=3.5, A3=3.5, A4=4, A5=7, A6=2.5)
cells = pd.read_parquet(f"{OUT}/cells_all.parquet")
nuc = pd.read_parquet("results/R2/R2_nuclei_all.parquet", columns=["archetype", "roi_id", "method", "scale", "keep", "label", "major_um", "minor_um"])
nuc = nuc[(nuc.scale == 1) & nuc.keep]
rr = pd.read_csv(f"{R}/R3_roi_summary.csv"); ns = pd.read_csv(f"{R}/R3_null_summary.csv")
def adiff(a, b): d = np.abs(a - b) % np.pi; return np.minimum(d, np.pi - d)
rng = np.random.default_rng(777)
G = []; P2 = []; NN = []
for (A, roi, me), g in cells.groupby(["archetype", "roi_id", "method"]):
    gen = pd.read_parquet(f"{OUT}/real/{A}_{roi}_{me}_gen.parquet")
    g = g.merge(gen[["idx", "label"]], on="idx").merge(nuc[(nuc.archetype == A) & (nuc.roi_id == roi) & (nuc.method == me)][["label", "major_um"]], on="label", how="left")
    xy = g[["x", "y"]].values; tr = cKDTree(xy); dd, ii = tr.query(xy, k=2); dnn = dd[:, 1]; jnn = ii[:, 1]
    r = rr[(rr.archetype == A) & (rr.roi_id == roi) & (rr.method == me)].iloc[0]
    lam = r.n / r.area_poly_um2
    NN.append(dict(archetype=A, roi_id=roi, method=me, frac_nn_lt_dstar=(dnn < DSTAR[A]).mean(), nn_q01=np.quantile(dnn, .01), nn_q05=np.quantile(dnn, .05),
                   nn_median=np.median(dnn), dstar=DSTAR[A], packing_dstar=lam * np.pi * (DSTAR[A] / 2) ** 2,
                   # (a) P3 / (d)
                   eqr_of_lambda=np.sqrt(1 / (np.pi * lam)), eqr_med_csr_theory=np.sqrt(gamma(3.5, scale=1 / 3.5).median() / (np.pi * lam)),
                   eqr_med_interior=np.median(np.sqrt(g.area[g.interior] / np.pi)), mean_area_int_x_lambda=g.area[g.interior].mean() * lam,
                   frac_interior=g.interior.mean(), lambda_loc_int_over_all=g.lambda_loc[g.interior].mean() / g.lambda_loc.mean()))
    if me != PRIM[A] and me != "stardist_he": continue
    g["dnn"] = dnn; g["jnn"] = jnn
    e = g[g.interior & (g.ecc_N >= .8) & (g.ecc_T >= .5) & np.isfinite(g.theta_N) & np.isfinite(g.major_um)].copy()
    e["dth"] = np.degrees(adiff(e.theta_T.values, e.theta_N.values))
    e["q"] = e.major_um / e.dnn                                  # P2 (direzionale)
    e["q_loc"] = e.major_um * np.sqrt(e.lambda_loc)              # non direzionale: asse / spaziatura media locale
    vx = xy[e.jnn.values, 0] - e.x.values; vy = xy[e.jnn.values, 1] - e.y.values
    e["phi_nn"] = np.degrees(adiff(np.arctan2(vy, vx), e.theta_N.values))  # angolo fra asse nucleo e direzione del vicino
    thn_all = g.theta_N.values
    e["dth_nb"] = np.degrees(adiff(e.theta_N.values, thn_all[e.jnn.values]))  # co-allineamento con il nucleo vicino
    ok = np.isfinite(e.dth_nb)
    perm = [np.median(np.degrees(adiff(e.theta_N.values[ok], rng.permutation(thn_all[np.isfinite(thn_all)])[:ok.sum()]))) for _ in range(500)]
    G.append(dict(archetype=A, roi_id=roi, method=me, n=len(e), med_dth=e.dth.median(), med_phi_nn=e.phi_nn.median(),
                  med_dth_neighbour=np.median(e.dth_nb[ok]), p_neighbour=(np.array(perm) <= np.median(e.dth_nb[ok])).mean(),
                  q_med=e.q.median(), qloc_med=e.q_loc.median()))
    e["archetype"] = A; e["method"] = me; e["roi_id"] = roi
    P2.append(e[["archetype", "method", "roi_id", "dth", "q", "q_loc", "phi_nn", "dth_nb", "ecc_N"]])
pd.DataFrame(NN).to_csv(f"{RV}/rv_nn_dstar_interior.csv", index=False)
pd.DataFrame(G).to_csv(f"{RV}/rv_cp2_alternatives_roi.csv", index=False)
P2 = pd.concat(P2)
P2["qloc_bin"] = pd.cut(P2.q_loc, [0, .5, .75, 1, 1.25, 1.5, 2, np.inf], right=False)
P2["q_bin"] = pd.cut(P2.q, [0, .25, .5, .75, 1, 1.5, np.inf], right=False)
P2["ecc_bin"] = pd.cut(P2.ecc_N, [.8, .9, .95, 1.0001], right=False)
a = P2.groupby(["archetype", "method", "qloc_bin"], observed=True).agg(n=("dth", "size"), med_dth=("dth", "median"), med_phi_nn=("phi_nn", "median")).reset_index()
a.to_csv(f"{RV}/rv_P2_by_qloc.csv", index=False)
b = P2.groupby(["archetype", "method", "q_bin"], observed=True).agg(n=("dth", "size"), med_dth=("dth", "median"), med_phi_nn=("phi_nn", "median")).reset_index()
b.to_csv(f"{RV}/rv_P2_by_q_phi.csv", index=False)
c = P2.groupby(["archetype", "method", "ecc_bin"], observed=True).agg(n=("dth", "size"), med_dth=("dth", "median"), med_q=("q", "median")).reset_index()
c.to_csv(f"{RV}/rv_P2_by_ecc.csv", index=False)
# per-ROI replication of the q effect (P2 was pooled): median dth for q>=1.5 vs q<1
rep = P2.groupby(["archetype", "method", "roi_id"]).apply(lambda z: pd.Series(dict(n_hi=(z.q >= 1.5).sum(), med_hi=z.dth[z.q >= 1.5].median(),
        n_lo=(z.q < 1).sum(), med_lo=z.dth[z.q < 1].median()))).reset_index()
rep.to_csv(f"{RV}/rv_P2_roi_replication.csv", index=False)
pd.set_option("display.width", 250)
print(pd.DataFrame(NN).round(4).to_string()); print(pd.DataFrame(G).round(4).to_string())
print(a[a.n >= 30].round(2).to_string()); print(rep.round(2).to_string())

