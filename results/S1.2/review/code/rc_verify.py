# rc_verify.py — verifiche indipendenti (numpy/scipy, nessuna funzione del progetto) per la revisione S1.2 track "code"
# Uso (root del repo): .venv/bin/python results/S1.2/review/code/rc_verify.py
import os, json, glob, math
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy import stats
OUT = "results/S1.2/review/code/out"; VER = "results/S1.2/review/code/verify"; os.makedirs(VER, exist_ok=True)
RES = {}

def rd(case, f):
    p = os.path.join(OUT, case, f)
    return pd.read_csv(p, dtype=str, keep_default_na=False) if os.path.exists(p) else None

def num(s):
    return pd.to_numeric(s.replace({"NA": np.nan, "": np.nan}), errors="coerce")

def rings_by_region(poly):
    out = {}
    for (rid,), g in poly.groupby(["region_id"], sort=False):
        rr = []
        for _, h in g.groupby(["part", "ring"], sort=False):
            rr.append(np.c_[h.x.astype(float).values, h.y.astype(float).values])
        out[rid] = rr
    return out

def pip(x, y, rings):  # even-odd ray casting su tutti gli anelli (bordo: misura nulla)
    ins = np.zeros(len(x), bool)
    for R in rings:
        x1, y1 = R[:-1, 0], R[:-1, 1]; x2, y2 = R[1:, 0], R[1:, 1]
        for a, b, c_, d in zip(x1, y1, x2, y2):
            if b == d: continue
            cr = ((b > y) != (d > y)) & (x < (c_ - a) * (y - b) / (d - b) + a)
            ins ^= cr
    return ins

def exp_dmin(ct, factor=2/3):
    dens = num(ct["density"]).values
    md = num(ct["min_dist_um"]).values if "min_dist_um" in ct else np.full(len(ct), np.nan)
    with np.errstate(divide="ignore"):
        eq = np.sqrt(1e6 / (np.pi * dens))
    return dict(zip(ct["cell_type"], np.where(np.isnan(md), factor * eq, md))), dict(zip(ct["cell_type"], dens))

def verify_case(case):
    meta = rd(case, "meta.csv").iloc[0].to_dict()
    r = {"case": case, "error": meta["error"], "warnings": meta["warnings"][:200], "elapsed_s": float(meta["elapsed_s"])}
    cen = rd(case, "centroids.csv")
    if cen is None: return r
    x = cen.x.astype(float).values; y = cen.y.astype(float).values
    ct = rd(case, "in_cell_types.csv"); comp = rd(case, "in_comp.csv"); ird = rd(case, "in_region_df.csv")
    rdf = rd(case, "region_df.csv"); tdf = rd(case, "type_df.csv"); poly = rd(case, "in_polygons.csv")
    dmap, dens = exp_dmin(ct)
    # 1) distanze minime globali (cKDTree)
    dt = np.array([dmap[t] for t in cen.cell_type]); r["n"] = len(cen)
    if len(cen) > 1 and np.nanmax(dt) > 0:
        tr = cKDTree(np.c_[x, y]); P = tr.query_pairs(r=float(np.nanmax(dt)) * 1.5, output_type="ndarray")
        dd = np.hypot(x[P[:, 0]] - x[P[:, 1]], y[P[:, 0]] - y[P[:, 1]]); thr = (dt[P[:, 0]] + dt[P[:, 1]]) / 2
        v = dd < thr - 1e-9; cross = cen.region_id.values[P[:, 0]] != cen.region_id.values[P[:, 1]]
        r.update(pairs=len(P), viol=int(v.sum()), viol_cross=int((v & cross).sum()), cross_pairs=int(cross.sum()),
                 min_ratio=float(np.min(dd / thr)) if len(P) else np.nan)
        tp = cen.cell_type.values
        mr = {}
        for i, j, q in zip(tp[P[:, 0]], tp[P[:, 1]], dd / thr):
            k = "|".join(sorted([i, j])); mr[k] = min(mr.get(k, np.inf), q)
        r["min_ratio_by_pair"] = json.dumps({k: round(v_, 4) for k, v_ in mr.items()})
    # 2) appartenenza al poligono (ray casting proprio)
    if poly is not None:
        RB = rings_by_region(poly); ok = np.zeros(len(cen), bool)
        for rid, idx in cen.groupby("region_id").indices.items():
            if rid in RB: ok[idx] = pip(x[idx], y[idx], RB[rid])
        r["outside"] = int((~ok).sum())
    # 3) conteggi per regione e coerenza con region_df
    cnt = cen.region_id.value_counts()
    n_cells = num(rdf.n_cells); n_target = num(rdf.n_target); n_failed = num(rdf.n_failed)
    obs = np.array([cnt.get(k, 0) for k in rdf.region_id])
    r["regions_count_mismatch"] = int((obs != n_cells.values).sum())
    r["sum_obs"] = int(obs.sum()); r["sum_n_cells_rdf"] = int(n_cells.sum())
    r["regions_ncells_ne_target"] = int((n_cells.values != n_target.values).sum()); r["n_failed"] = int(n_failed.sum())
    r["centroid_region_not_in_input"] = int((~cen.region_id.isin(ird.region_id)).sum())
    # 4) densita' armonica ricalcolata dagli input dichiarati
    comp = comp.assign(fraction=num(comp.fraction))
    comp = comp[comp.fraction > 0]
    rho = {}
    for cl, g in comp.groupby("cluster_id"):
        f = g.fraction.values / g.fraction.sum(); rho[cl] = 1.0 / np.sum(f / np.array([dens[t] for t in g.cell_type]))
    exp_rho = np.array([rho.get(c, np.nan) for c in ird.cluster_id]); got = num(rdf.target_density_weighted).values
    r["max_rel_err_rho"] = float(np.nanmax(np.abs(got - exp_rho) / exp_rho)) if np.isfinite(exp_rho).any() else np.nan
    lam = exp_rho * num(ird.area_um2).values / 1e6
    r["n_target_outside_floor_ceil"] = int(np.sum(~((n_target.values == np.floor(lam)) | (n_target.values == np.ceil(lam)))))
    # 5) |n_i - n f_i| < 1 per regione e tipo (dove n_failed = 0)
    worst = 0.0
    tab = cen.groupby(["region_id", "cell_type"]).size()
    for k, row in rdf.iterrows():
        if float(row.n_failed) > 0 or row.cluster_id not in set(comp.cluster_id): continue
        g = comp[comp.cluster_id == row.cluster_id]; f = g.fraction.values / g.fraction.sum()
        for t, fi in zip(g.cell_type, f):
            worst = max(worst, abs(tab.get((row.region_id, t), 0) - float(row.n_cells) * fi))
    r["max_dev_type_alloc"] = worst
    # 6) ordine dei tipi atteso (densita' decrescente, parita' -> nome in ordine di code point)
    used = sorted(set(comp.cell_type), key=lambda t: (-dens[t], t))
    r["type_order_got"] = "|".join(tdf.cell_type); r["type_order_codepoint"] = "|".join(used)
    r["dmin_typedf_vs_expected_maxabs"] = float(max(abs(float(a) - dmap[t]) for t, a in zip(tdf.cell_type, tdf.min_dist_um)))
    return r

cases = sorted(d for d in os.listdir(OUT) if os.path.isdir(os.path.join(OUT, d)) and os.path.exists(os.path.join(OUT, d, "meta.csv")))
rows = [verify_case(c) for c in cases]
V = pd.DataFrame(rows); V.to_csv(os.path.join(VER, "V_cases.csv"), index=False)

# ---- A02 invarianza alla griglia ----
a, b = rd("A02_near", "centroids.csv"), rd("A02_far", "centroids.csv")
if a is not None and b is not None:
    A = a[a.region_id.isin(["1", "2"])]; B = b[b.region_id.isin(["1", "2"])]
    same12 = len(A) == len(B) and np.array_equal(A.x.astype(float).values, B.x.astype(float).values) and np.array_equal(A.y.astype(float).values, B.y.astype(float).values)
    A3 = a[a.region_id == "3"]; B3 = b[b.region_id == "3"]
    d3 = float(np.max(np.abs((B3.x.astype(float).values - 199000) - A3.x.astype(float).values))) if len(A3) == len(B3) and len(A3) else np.nan
    g = float(rd("A02_far", "meta.csv").grid_cell_um[0]); xs = b.x.astype(float).values; ys = b.y.astype(float).values
    occ = pd.Series(list(zip(np.floor(xs / g).astype(int), np.floor(ys / g).astype(int)))).value_counts()
    RES["A02"] = dict(identical_regions12=bool(same12), n12=len(A), n12_far=len(B), max_dx_region3_shifted=d3, grid_cell_far=g,
                      grid_cell_near=float(rd("A02_near", "meta.csv").grid_cell_um[0]), max_pts_per_grid_cell_far=int(occ.max()), K_initial=8)

# ---- A11 effetti di bordo ----
def a11():
    out = {}
    sqc = pd.read_csv(os.path.join(OUT, "A11_square", "centroids_20seeds.csv"))
    dist = np.minimum.reduce([sqc.x.values, 1000 - sqc.x.values, sqc.y.values, 1000 - sqc.y.values])
    edges = np.array([0, 1, 2, 3, 4, 5, 7.5, 10, 15, 20, 30, 50])
    nseed = sqc.seed.nunique()
    def band_area(a_, b_): return (1000 - 2 * a_) ** 2 - (1000 - 2 * b_) ** 2
    inner = dist >= 50; rho_in = inner.sum() / ((1000 - 100) ** 2 * nseed)
    prof = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = ((dist >= lo) & (dist < hi)).sum(); ar = band_area(lo, hi) * nseed
        prof.append(dict(band=f"[{lo},{hi})", n=int(k), rel_density=k / ar / rho_in, se=math.sqrt(k) / ar / rho_in))
    out["square_profile"] = prof
    for nm in ["A11_adjLR", "A11_adjRL"]:
        d = pd.read_csv(os.path.join(OUT, nm, "centroids_20seeds.csv")); d = d[(d.y > 50) & (d.y < 950)]
        res = []
        for lo, hi in [(0, 1), (1, 2), (2, 3.5), (3.5, 5), (5, 10), (10, 20)]:
            L = ((d.x < 500 - lo) & (d.x >= 500 - hi)).sum(); R = ((d.x >= 500 + lo) & (d.x < 500 + hi)).sum()
            ar = (hi - lo) * 900 * d.seed.nunique()
            Lr = d[(d.x < 500 - lo) & (d.x >= 500 - hi)].region_row.unique().tolist()
            res.append(dict(band=f"[{lo},{hi})", left_rel=L / ar / rho_in, right_rel=R / ar / rho_in, left_row=Lr, ratio_left_right=L / max(R, 1)))
        out[nm] = res
    return out
try: RES["A11"] = a11()
except Exception as e: RES["A11"] = {"error": str(e)}

# ---- A16 Madow ----
def a16():
    out = {}
    D = pd.read_csv(os.path.join(OUT, "A16a_madow_direct.csv"))
    specs = {"n7": (7, [.37, .29, .21, .13]), "n1q": (1, [.25] * 4), "n10": (10, [.1, .2, .3, .4])}
    for k, (n, f) in specs.items():
        M = D[D.case == k].iloc[:, 1:5].values; e = n * np.array(f); fl = np.floor(e + 1e-9); fr = e - fl
        p_hat = (M > fl).mean(0); se = np.sqrt(np.maximum(fr * (1 - fr), 1e-300) / len(M))
        out[k] = dict(reps=len(M), sum_always_n=bool((M.sum(1) == n).all()), p_incl=p_hat.round(5).tolist(), frac=fr.round(5).tolist(),
                      max_z=float(np.max(np.abs(p_hat - fr) / np.where(fr * (1 - fr) > 0, se, np.inf))) if (fr * (1 - fr) > 0).any() else 0.0,
                      never_outside_floor_ceil=bool(((M == fl) | (M == fl + 1)).all()),
                      mean=M.mean(0).round(4).tolist(), expected_mean=e.round(4).tolist())
    # implementazione indipendente del campionamento sistematico di Madow e confronto delle distribuzioni congiunte
    rng = np.random.default_rng(1); n, f = 7, np.array([.37, .29, .21, .13]); e = n * f; fl = np.floor(e + 1e-9); fr = e - fl; r_ = int(n - fl.sum())
    cs = np.cumsum(fr); cs = cs * r_ / cs[-1]
    U = rng.random(200000)[:, None] + np.arange(r_)[None, :]
    hit = np.stack([((U > np.r_[0, cs][i]) & (U <= cs[i])).sum(1) for i in range(4)], 1)
    pat_py = pd.Series(["".join(map(str, h)) for h in hit]).value_counts(normalize=True)
    M7 = D[D.case == "n7"].iloc[:, 1:5].values - fl.astype(int)
    pat_r = pd.Series(["".join(map(str, h)) for h in M7]).value_counts(normalize=True)
    out["joint_patterns_R"] = pat_r.round(4).to_dict(); out["joint_patterns_python"] = pat_py.round(4).to_dict()
    B = pd.read_csv(os.path.join(OUT, "A16b_madow_full.csv"))
    Mb = B[["w", "x", "y", "z"]].values
    out["full_fn"] = dict(regions_x_seeds=len(B), n_target_always_7=bool((B.n_target == 7).all()), p_incl=(Mb > fl).mean(0).round(4).tolist(),
                          frac=fr.round(4).tolist(), mean=Mb.mean(0).round(4).tolist(), expected_mean=e.round(4).tolist(),
                          max_z=float(np.max(np.abs((Mb > fl).mean(0) - fr) / np.sqrt(fr * (1 - fr) / len(Mb)))))
    return out
try: RES["A16"] = a16()
except Exception as e: RES["A16"] = {"error": str(e)}

# ---- A17 C1b ----
def a17():
    D = pd.read_csv(os.path.join(OUT, "A17_c1b_1000seeds.csv")); p = 1185 * 200 / 1e6
    s200 = D[D.seed <= 200].n_cells
    sim = np.random.default_rng(7).binomial(1, p, size=(1000, 2000)).sum(1)   # schema floor + Bernoulli simulato da zero (floor = 0)
    sd_th = math.sqrt(2000 * p * (1 - p))
    chi = (len(D) - 1) * D.n_cells.var() / sd_th ** 2
    return dict(lambda_region=p, expected=2000 * p, sd_theory=sd_th, mean_1000=D.n_cells.mean(), sd_1000=D.n_cells.std(),
                z_mean_1000=(D.n_cells.mean() - 2000 * p) / (sd_th / math.sqrt(len(D))),
                p_var_chi2_two_sided=float(2 * min(stats.chi2.cdf(chi, len(D) - 1), stats.chi2.sf(chi, len(D) - 1))),
                mean_seeds1_200=s200.mean(), sd_seeds1_200=s200.std(), any_failed=int(D.n_failed.sum()), target_eq_cells=bool((D.n_cells == D.n_target).all()),
                python_sim_mean=float(sim.mean()), python_sim_sd=float(sim.std(ddof=1)))
try: RES["A17"] = a17()
except Exception as e: RES["A17"] = {"error": str(e)}

# ---- mutanti M1/M2/M4 del progetto: rilevati per la ragione giusta? ----
def mut():
    out = {}
    c1 = rd("M1_proj_I1", "centroids.csv"); poly = rd("M1_proj_I1", "in_polygons.csv")
    if c1 is not None:
        RB = rings_by_region(poly); x = c1.x.astype(float).values; y = c1.y.astype(float).values; ok = np.zeros(len(c1), bool); inbb = np.zeros(len(c1), bool)
        for rid, idx in c1.groupby("region_id").indices.items():
            ok[idx] = pip(x[idx], y[idx], RB[rid]); allv = np.vstack(RB[rid])
            inbb[idx] = (x[idx] >= allv[:, 0].min()) & (x[idx] <= allv[:, 0].max()) & (y[idx] >= allv[:, 1].min()) & (y[idx] <= allv[:, 1].max())
        out["M1"] = dict(outside=int((~ok).sum()), outside_but_in_own_bbox=int(((~ok) & inbb).sum()))
    c2 = rd("M2_proj_I1", "centroids.csv")
    if c2 is not None:
        ct = rd("M2_proj_I1", "in_cell_types.csv"); dmap, _ = exp_dmin(ct); x = c2.x.astype(float).values; y = c2.y.astype(float).values
        dt = np.array([dmap[t] for t in c2.cell_type]); P = cKDTree(np.c_[x, y]).query_pairs(r=float(dt.max()), output_type="ndarray")
        dd = np.hypot(x[P[:, 0]] - x[P[:, 1]], y[P[:, 0]] - y[P[:, 1]]); v = P[dd < (dt[P[:, 0]] + dt[P[:, 1]]) / 2 - 1e-9]
        g = float(rd("M2_proj_I1", "meta.csv").grid_cell_um[0]); pl = rd("M2_proj_I1", "in_polygons.csv")
        x0, y0 = pl.x.astype(float).min(), pl.y.astype(float).min()
        ix = np.floor((x - x0) / g).astype(int); iy = np.floor((y - y0) / g).astype(int)
        same = (ix[v[:, 0]] == ix[v[:, 1]]) & (iy[v[:, 0]] == iy[v[:, 1]])
        cheb = np.maximum(abs(ix[v[:, 0]] - ix[v[:, 1]]), abs(iy[v[:, 0]] - iy[v[:, 1]]))
        out["M2"] = dict(violations=len(v), in_same_cell=int(same.sum()), in_adjacent_cells=int((cheb == 1).sum()), grid_cell=g)
    for nm in ["M4_proj_unif", "M4_ref_unif"]:
        c4 = rd(nm, "centroids.csv")
        if c4 is None: continue
        x = c4.x.astype(float).values; y = c4.y.astype(float).values
        q = np.bincount(np.minimum((x // 100).astype(int), 9) + 10 * np.minimum((y // 100).astype(int), 9), minlength=100)
        out[nm] = dict(n=len(c4), empty_quadrats=int((q == 0).sum()), chi2_p=float(stats.chisquare(q).pvalue), var_over_mean=float(q.var(ddof=1) / q.mean()))
    return out
try: RES["mutants"] = mut()
except Exception as e: RES["mutants"] = {"error": str(e)}


# ---- extra: A08 (locale), A05e (scambio di regioni), A04b, A12 tentativi, A13b, valori derivati della pre-registrazione ----
ex = {}
c_, e_ = rd("A08_tie_C", "centroids.csv"), rd("A08_tie_enGBUTF8", "centroids.csv")
if c_ is not None and e_ is not None:
    ex["A08"] = dict(type_order_C="|".join(rd("A08_tie_C", "type_df.csv").cell_type), type_order_enGB="|".join(rd("A08_tie_enGBUTF8", "type_df.csv").cell_type),
                     identical_xy=bool(np.array_equal(c_.x.values, e_.x.values) and np.array_equal(c_.y.values, e_.y.values)),
                     n_C=len(c_), n_enGB=len(e_), n_a_C=int((c_.cell_type == "a").sum()), n_a_enGB=int((e_.cell_type == "a").sum()))
for k in ["A05e_rid_factor_revlev", "A05a_rid_factor", "A04b_zero_area_poly", "A13a_area_x10"]:
    cen = rd(k, "centroids.csv"); rdf_ = rd(k, "region_df.csv"); ird = rd(k, "in_region_df.csv"); poly = rd(k, "in_polygons.csv")
    if cen is None: continue
    d = dict(input_region_ids=ird.region_id.tolist(), input_area=ird.area_um2.tolist(), rdf_n_target=rdf_.n_target.tolist(), rdf_n_cells=rdf_.n_cells.tolist(),
             centroid_region_ids=cen.region_id.value_counts().to_dict())
    if poly is not None:
        RB = rings_by_region(poly); x = cen.x.astype(float).values; y = cen.y.astype(float).values
        d["centroids_inside_region_by_label"] = {rid: int(pip(x, y, RB[rid]).sum()) for rid in RB}
        d["poly_shoelace_area"] = {rid: float(sum(0.5 * abs(np.dot(R[:-1, 0], R[1:, 1]) - np.dot(R[1:, 0], R[:-1, 1])) for R in RB[rid][:1])) for rid in RB}
        if k == "A04b_zero_area_poly":
            d["y_of_points_with_y_ge_150"] = sorted(set(np.round(y[y >= 150], 9).tolist()))
    ex[k] = d
for k in ["A12a_strip_1x1000", "A12b_strip_rot45"]:
    m = rd(k, "meta.csv")
    if m is not None: ex[k] = m[["elapsed_s", "n_cells", "n_target", "n_failed", "n_attempts"]].iloc[0].to_dict()
p13 = os.path.join(OUT, "A13b_area_consistency.csv")
if os.path.exists(p13): ex["A13b"] = pd.read_csv(p13).to_dict(orient="records")
tot = np.array([12419, 8924, 3100, 28096, 1185, 988.]); evi = np.array([10291, 7376, 2311, 26844, 953, 906.])
eq = np.sqrt(1e6 / (np.pi * tot))
ex["prereg_derived"] = dict(width=(1 - evi / tot).round(3).tolist(), eq_r=eq.round(2).tolist(), d_rule=(2 / 3 * eq).round(2).tolist(),
                            d_jam=np.sqrt(4 * 0.547 / (np.pi * tot / 1e6)).round(2).tolist(), eta_rule=float(np.pi * (2 / 3) ** 2 / (4 * np.pi)),
                            c1b_expected=2000 * 200 * 1185 / 1e6, c1b_3se=3 * math.sqrt(2000 * .237 * .763) / math.sqrt(200),
                            c4b_rho_mix=1 / (.9 / 8000 + .08 / 3000 + .02 / 20000), c4b_expected=[2000 * 400 / 1e6 / (.9 / 8000 + .08 / 3000 + .02 / 20000) * f for f in (.9, .08, .02)])
RES["extra"] = ex

json.dump(RES, open(os.path.join(VER, "V_special.json"), "w"), indent=1, default=str)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_colwidth", 60)
print(json.dumps({k: RES[k] for k in ["A11", "extra"]}, indent=1, default=str))
