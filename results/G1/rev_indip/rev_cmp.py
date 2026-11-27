import pandas as pd, numpy as np
REPO = "/home/user/2025.geo_spatialtrans/"; OUT = "/mnt/micron/geo_spatialtrans/G1_review/rev_indip/"
M = pd.read_csv(OUT + "rev_M_cells.csv")
print("labels disagree total", int((M.n_labels - M.label_agree).sum()), "of", int(M.n_labels.sum()))
t = pd.concat([pd.read_csv(REPO + "results/sonda/S3b/area_values.csv"), pd.read_csv(REPO + "results/sonda/S4/area_values.csv"),
               pd.read_csv(REPO + "results/sonda/S6/area_values.csv")], ignore_index=True)
t["area"] = t.unit.str.rsplit("_", n=1).str[-1].astype(int); t["seed"] = t["k"].astype(str)
t = t[t.arm.isin(M.arm.unique())]
key = ["arm", "A", "area", "seed"]
m = M.merge(t, on=key, suffixes=("_rev", "_dec"), how="outer", indicator=True)
print("merge", m._merge.value_counts().to_dict())
rows = []
for c in ["M1", "M2", "M3", "M4", "M5", "M6", "n_cells", "n_typed", "n_nuc_blocks", "n_pairs"]:
    d = (m[c + "_rev"] - m[c + "_dec"]).abs(); rel = d / m[c + "_dec"].abs()
    rows.append({"col": c, "max_abs": d.max(), "max_rel": rel.max(), "n_gt_1e9rel": int((rel > 1e-9).sum())})
R = pd.DataFrame(rows); print(R.to_string(index=False))
m["max_rel_M"] = np.max([((m[c + "_rev"] - m[c + "_dec"]).abs() / m[c + "_dec"].abs()).to_numpy() for c in ["M1","M2","M3","M4","M5","M6"]], axis=0)
m[key + ["n_labels", "label_agree", "max_rel_M"] + [f"{c}_{s}" for c in ["M1","M2","M3","M4","M5","M6"] for s in ("rev","dec")]].to_csv(OUT + "rev_M_compare.csv", index=False)
print(m.groupby("arm").max_rel_M.max().to_dict())
bad = m[m.max_rel_M > 1e-9]; print(bad[key + ["max_rel_M", "n_labels", "label_agree"]].to_string(index=False) if len(bad) else "all 80 within 1e-9")
