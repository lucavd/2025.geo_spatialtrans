import json, glob, numpy as np, pandas as pd
D = "results/R3/review/code"
J = pd.DataFrame([json.load(open(f)) for f in glob.glob(f"{D}/roi/*.json")])
J.to_csv(f"{D}/rc_roi_all.csv", index=False)
base = J[J.jitter_um == 0].copy()
rs = pd.read_csv("results/R3/R3_roi_summary.csv")
cp2 = pd.read_csv("results/R3/R3_cp2_roi.csv")
ca = pd.read_parquet("/mnt/micron/geo_spatialtrans/R3/cells_all.parquet")
ci = ca[ca.interior]
g = ci.groupby(["archetype", "roi_id", "method"])
theirs_cell = pd.DataFrame({"nc_med": g.nc.median(), "ratio_med": g.ratio.median(),
                            "frac_cut": g.frac_out.apply(lambda v: np.mean(v > 0.05)),
                            "eq_r_med_cells": g.eq_r.median()}).reset_index()
M = base.merge(rs, on=["archetype", "roi_id", "method"], suffixes=("", "_R3"), how="left") \
        .merge(theirs_cell, on=["archetype", "roi_id", "method"], how="left") \
        .merge(cp2[["archetype", "roi_id", "n_eligible", "median_dtheta", "perm_p"]], on=["archetype", "roi_id"], how="left", suffixes=("", "_cp2"))
pairs = [("n", "n_R3"), ("n_keep", "n_keep_R3"), ("sum_area", "sum_area_R3"), ("area_mask_poly", "area_poly_um2"),
         ("frac_interior", "frac_interior_R3"), ("median_area", "median_area_R3"), ("median_eq_r", "median_eq_r_R3"),
         ("cv", "cv_R3"), ("median_ecc_T", "median_ecc_T_R3"), ("mean_nsides", "mean_nsides_R3"), ("cv_loc", "cv_loc_R3"),
         ("median_nc", "nc_med"), ("median_ratio", "ratio_med"), ("frac_cut", "frac_cut"),
         ("n_eligible", "n_eligible_cp2"), ("median_dtheta", "median_dtheta_cp2"), ("perm_p", "perm_p_cp2"), ("sigma_loc", "sigma_loc_um")]
# ri-nomina per evitare collisioni (colonne R3 senza suffisso dove non c'e' clash)
for a, b in pairs:
    if b not in M.columns and b.replace("_R3", "") in rs.columns: M[b] = M[b.replace("_R3", "")]
rows = []
for _, r in M.iterrows():
    for a, b in pairs:
        va, vb = r.get(a), r.get(b)
        try: rel = abs(float(va) - float(vb)) / max(abs(float(vb)), 1e-12)
        except Exception: rel = np.nan
        rows.append(dict(archetype=r.archetype, roi_id=r.roi_id, method=r.method, metric=a, recomputed=va, declared=vb, rel_diff=rel))
C = pd.DataFrame(rows); C.to_csv(f"{D}/rc_compare_roi.csv", index=False)
print(C.groupby("metric").rel_diff.agg(["count", "max", "median"]).to_string())
print(base[["n_keep", "n", "n_inmask_edgeconv", "n_dup", "gen_in_own", "cent_err_max", "n_gen_without_px", "rel_err_sum"]].describe().loc[["min", "max"]].to_string())
print("edge-conv differences:", (base.n_inmask_edgeconv - base.n).abs().sum(), "max per ROI", (base.n_inmask_edgeconv - base.n).abs().max())
print("cv_loc_loo vs cv_loc max rel:", (base.cv_loc_loo / base.cv_loc - 1).abs().max())
