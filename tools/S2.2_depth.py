"""S2.2, regola R-3: fattore di profondità per tessuto, appaiato per area.
f_k = (UMI veri medi per bin di 2 um in tessuto dell'area k, R4_roi.csv mean_umi2, - fondo b) /
      (somma sulle cellule simulate dell'area k del totale medio del loro tipo / (area dei territori / 4 um^2)).
Totale medio del tipo: media dei totali dei nuclei veri classificabili (impronta intera, almeno 20 UMI) di R6 sulle 5 aree.
Fattore del tessuto = media delle 5 aree; intervallo bootstrap sulle aree (B = 2000, seme 20261107).
Uso: .venv/bin/python tools/S2.2_depth.py -> results/S2.2/S2.2_depth_factor.csv, S2.2_depth_per_roi.csv, S2.2_type_totals.csv"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "tools"); import S22_common as S
M = "/mnt/micron/geo_spatialtrans"
B_GLASS = {"A1": 0.161, "A2": 0.197, "A3": 0.197, "A4": 0.042, "A5": 0.309, "A6": 0.117}   # B-342, B-348, B-354, B-360, B-366 (design sezione 9)
KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
tt = []
for A in S.PRESET_TYPES:
    lab = pd.read_parquet(f"{M}/R6/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet")
    lab["ptype"] = lab.cell_type.replace(S.MERGE)
    lab = lab[lab.ptype.isin(S.PRESET_TYPES[A]) & lab.classifiable]
    g = lab.groupby("ptype").umi
    t = pd.DataFrame({"archetype": A, "cell_type": g.mean().index, "n_nuclei": g.size().values, "mean_total": g.mean().values,
                      "median_total": g.median().values, "q01": g.quantile(0.01).values, "q99": g.quantile(0.99).values})
    tt.append(t)
tt = pd.concat(tt, ignore_index=True); tt.to_csv("results/S2.2/S2.2_type_totals.csv", index=False)
lay = pd.read_csv("results/S2.2/S2.2_depth_layers.csv")
lay = lay.merge(tt[["archetype", "cell_type", "mean_total"]], on=["archetype", "cell_type"], how="left")
assert lay.mean_total.notna().all(), lay[lay.mean_total.isna()]
r4 = pd.read_csv("results/R4/R4_roi.csv")[["archetype", "roi_id", "mean_umi2", "n_bins2"]]
per = lay.assign(exp_tot=lay.n_cells * lay.mean_total).groupby(["archetype", "roi_id"]).agg(exp_tot=("exp_tot", "sum"), n_cells=("n_cells", "sum"),
                                                                                            area=("territory_area_um2", "first")).reset_index()
per = per.merge(r4, on=["archetype", "roi_id"])
per["b"] = per.archetype.map(B_GLASS)
per["sim_umi2_f1"] = per.exp_tot / (per.area / 4.0)
per["f"] = (per.mean_umi2 - per.b) / per.sim_umi2_f1
per["area_ratio_sim_real"] = per.area / (per.n_bins2 * 4.0)
per.to_csv("results/S2.2/S2.2_depth_per_roi.csv", index=False)
rng = np.random.default_rng(20261107); out = []
for A, g in per.groupby("archetype"):
    v = g.f.to_numpy(); bs = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(2000)])
    out.append(dict(archetype=A, f=v.mean(), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5), f_min=v.min(), f_max=v.max(), n_roi=len(v),
                    share_explained=1 / v.mean()))
out = pd.DataFrame(out); out.to_csv("results/S2.2/S2.2_depth_factor.csv", index=False)
print(per[["archetype", "roi_id", "n_cells", "mean_umi2", "sim_umi2_f1", "f", "area_ratio_sim_real"]].round(3).to_string(index=False))
print(out.round(3).to_string(index=False))
