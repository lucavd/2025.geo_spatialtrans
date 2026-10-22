"""R6 gate G-4: consistency of labels and counts. Every kept nucleus has exactly one label; fractions per ROI sum to 1;
bins assigned <= filtered bins of the ROI; UMI assigned <= UMI of the ROI; counts matrix rows == nuclei. Exit 1 on any FAIL."""
import sys, json, glob, numpy as np, pandas as pd, scipy.sparse as sp
M = "/mnt/micron/geo_spatialtrans/R6"; res = []
def check(name, cond): res.append((name, bool(cond))); print(("PASS" if cond else "FAIL"), name)
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
for _, R in rois.iterrows():
    tag = f"{R.archetype}_{R.roi_id}"; info = json.load(open(f"{M}/counts/{tag}_info.json")); tab = pd.read_parquet(f"{M}/counts/{tag}_nuclei.parquet")
    ok_bins = all(info[f"bins_assigned_{fp}"] <= info["n_bins_filtered"] for fp in ["exp", "nuc", "shift"])
    ok_umi = all(info[f"umi_assigned_{fp}"] <= info["umi_roi"] + 1e-6 for fp in ["exp", "nuc", "shift"])
    ok_nuc_sub = (tab.n_bins_nuc <= tab.n_bins_exp).all() and (tab.umi_nuc <= tab.umi_exp + 1e-9).all()
    ok_rows = all(sp.load_npz(f"{M}/counts/{tag}_{fp}.npz").shape[0] == len(tab) for fp in ["exp", "nuc", "shift"])
    ok_sum = all(abs(sp.load_npz(f"{M}/counts/{tag}_{fp}.npz").sum() - tab[f"umi_{fp}"].sum()) < 1e-6 for fp in ["exp", "nuc", "shift"])
    check(f"{tag}: bins/UMI assigned within ROI totals, nucleus-only subset of expanded, matrix rows = nuclei, matrix sum = UMI table", ok_bins and ok_umi and ok_nuc_sub and ok_rows and ok_sum)
for f in sorted(glob.glob(f"{M}/labels/*_nuclei.parquet")) + sorted(glob.glob(f"{M}/labels/*_nuclei_ph1.parquet")):
    lab = pd.read_parquet(f); name = f.split("/")[-1]
    one = lab.cell_type.notna().all() and (lab.cell_type != "").all()
    fr = lab.groupby("roi_id").cell_type.value_counts(normalize=True).groupby(level=0).sum()
    sums = np.allclose(fr.values, 1.0, atol=1e-9)
    typed_have_cluster = lab.loc[~lab.cell_type.isin(["non classificabile"]), "cluster"].notna().all()
    nocl_unclass = (lab.loc[lab.cluster.isna(), "cell_type"] == "non classificabile").all()
    uniq = not lab.duplicated(["roi_id", "label"]).any()
    check(f"{name}: one label per nucleus, fractions sum to 1 per ROI, typed nuclei have a cluster, cluster-less are 'non classificabile', unique (roi,label)", one and sums and typed_have_cluster and nocl_unclass and uniq)
n_fail = sum(not ok for _, ok in res); print(f"{len(res) - n_fail}/{len(res)} PASS")
pd.DataFrame(res, columns=["check", "pass"]).to_csv("results/R6/R6_G4_check.csv", index=False); sys.exit(1 if n_fail else 0)
