# Revisione avversariale S1.1 — 08_compare.py
# Confronta i ricalcoli indipendenti (py_out) con i CSV riportati e con r0_df salvato dal test.
import json, glob, numpy as np, pandas as pd
REPO = "/home/user/2025.geo_spatialtrans"; REV = "/mnt/micron/geo_spatialtrans/S1.1/review"
summ = {}
for f in glob.glob(f"{REV}/py_out/components_summary_*.json"): summ.update(json.load(open(f)))
per = pd.read_csv(f"{REPO}/results/S1.1/S1.1_per_input.csv").set_index("id")
cp3 = pd.read_csv(f"{REPO}/results/S1.1/S1.1_cp3_null.csv").set_index("dataset")
out = {}
# compito 1: conteggi + multinsieme (cluster, n_px) vs r0_df salvato
t1 = []
for id_ in ["real_full_A1","real_full_A6","null_full_A2A3","I6_syn6800x6500_c2","roi_A6_r2","I2_syn600_c2","adv_checker","adv_pinch_hole","I5_syn600_c2_labels"]:
    mine = pd.read_parquet(f"{REV}/py_out/{id_}_components.parquet")
    theirs = pd.read_parquet(f"{REV}/export/{id_}_r0df.parquet")
    a = sorted(zip(mine.cluster.astype(str), mine.n_px.astype(int)))
    def norm(c):
        try: return str(int(float(c)))
        except: return str(c)
    b = sorted(zip([norm(c) for c in theirs.cluster_id], theirs.n_px.astype(int)))
    a = sorted((norm(c), n) for c, n in a)
    t1.append(dict(id=id_, n_comp_recomputed=len(mine), n_comp_reported=int(per.loc[id_, "n_components"]),
                   multiset_identical=a == b, size_max=int(mine.n_px.max()), n_single_px=int((mine.n_px == 1).sum()),
                   sum_npx=int(mine.n_px.sum()), n_bins=int(per.loc[id_, "n_px"])))
out["task1"] = t1
# compito 2
t2 = []
for g in ["real_full", "null_full"]:
    for a in ["A1", "A2A3", "A4", "A5", "A6"]:
        id_ = f"{g}_{a}"; s = summ[id_]
        t2.append(dict(id=id_, n_regions_100_re=s["n_regions_100"], n_regions_100_rep=int(per.loc[id_, "n_regions_100"]),
                       n_excl_re=s["n_excluded_100"], n_excl_rep=int(per.loc[id_, "n_excluded_100"]),
                       frac_excl_re=s["frac_area_excluded_100"], frac_excl_rep=float(per.loc[id_, "frac_area_excluded_100"]),
                       n_comp_re=s["n_components"], n_comp_rep=int(per.loc[id_, "n_components"]),
                       n_single_px=s["n_single_px"]))
t2 = pd.DataFrame(t2); out["task2"] = t2.to_dict("records")
for a in ["A1", "A2A3", "A4", "A5", "A6"]:
    r, n = summ[f"real_full_{a}"], summ[f"null_full_{a}"]
    out.setdefault("cp3", []).append(dict(ds=a, ratio_re=n["n_regions_100"] / r["n_regions_100"], ratio_rep=float(cp3.loc[a, "ratio_regions"]),
        excl_real=r["frac_area_excluded_100"], excl_null=n["frac_area_excluded_100"],
        comp_ratio=n["n_components"] / r["n_components"],
        comp_real_per_bin=r["n_components"] / r["n_bins"], comp_null_per_bin=n["n_components"] / n["n_bins"]))
# compito 4: I5
tr = pd.read_csv(f"{REV}/export/I5_syn600_c2_labels_truth.csv")
mine = pd.read_parquet(f"{REV}/py_out/I5_syn600_c2_labels_components.parquet")
g = mine.groupby(mine.cluster.astype(int)).agg(n_px=("n_px", "sum"), n_comp=("n_px", "size"), max_comp=("n_px", "max"))
m = tr.merge(g, left_on="patch", right_index=True, how="outer")
out["task4"] = dict(n_patch=int(len(tr)), sum_truth=int(tr.Freq.sum()), all_equal=bool((m.Freq == m.n_px).all()),
                    all_single_component=bool((m.n_comp == 1).all()), min_patch=int(tr.Freq.min()), max_patch=int(tr.Freq.max()),
                    truth_cols=list(tr.columns))
rep4 = pd.read_csv(f"{REPO}/results/S1.1/S1.1_cp4_truth.csv")
mm = rep4.merge(m, on="patch", suffixes=("_rep", "_re"))
out["task4"]["reported_vs_recomputed_identical"] = bool((mm.n_px_rep == mm.n_px_re).all() and (mm.n_components == mm.n_comp).all())
# compito 9
tt = pd.read_csv(f"{REPO}/results/S1.1/S1.1_test_results.csv")
out["task9"] = dict(n=int(len(tt)), counts=tt.status.value_counts().to_dict(),
                    by_check=tt.groupby(["check", "status"]).size().unstack(fill_value=0).to_dict("index"),
                    fails=tt[tt.status != "PASS"].to_dict("records"))
json.dump(out, open(f"{REV}/py_out/compare.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))

