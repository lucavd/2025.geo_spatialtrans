# Revisione avversariale S1.1 — 03_sources.py
# Ricostruisce real_full_A1 dai file sorgente (clusters.csv + tissue_positions.parquet), confronta con l'RDS
# esportato; md5 delle sorgenti di tutti i 5 dataset; verifiche di I6 (stride) e del nullo (composizione).
import hashlib, json, numpy as np, pandas as pd
REPO = "/home/user/2025.geo_spatialtrans"
EXP = "/mnt/micron/geo_spatialtrans/S1.1/review/export"
OUT = "/mnt/micron/geo_spatialtrans/S1.1/review/py_out"
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2A3": "Visium_HD_Human_Colon_Cancer",
      "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
def md5(f):
    h = hashlib.md5()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
res = {"md5": {}}
rep = pd.read_csv(f"{REPO}/results/S1.1/S1.1_real_sources_md5.csv")
for tag, ds in DS.items():
    sq = f"data/real/datasets/{ds}/binned_outputs_x/binned_outputs/square_008um"
    for f in [f"{sq}/analysis/clustering/gene_expression_graphclust/clusters.csv", f"{sq}/spatial/tissue_positions.parquet"]:
        m = md5(f"{REPO}/{f}")
        r = rep.loc[rep.file == f, "md5"]
        res["md5"][f] = {"recomputed": m, "reported": (r.iloc[0] if len(r) else None), "match": bool(len(r) and r.iloc[0] == m)}
# ricostruzione di A1 e (economica) di tutti gli altri
res["recon"] = {}
for tag, ds in DS.items():
    sq = f"{REPO}/data/real/datasets/{ds}/binned_outputs_x/binned_outputs/square_008um"
    cl = pd.read_csv(f"{sq}/analysis/clustering/gene_expression_graphclust/clusters.csv")
    pos = pd.read_parquet(f"{sq}/spatial/tissue_positions.parquet")
    d = pos.merge(cl, left_on="barcode", right_on="Barcode", how="inner")
    rec = pd.DataFrame({"x": d.array_col + 1, "y": d.array_row + 1, "cluster": d.Cluster.astype(str)}).sort_values(["y", "x"]).reset_index(drop=True)
    rds = pd.read_parquet(f"{EXP}/real_full_{tag}.parquet")
    rds = rds.assign(x=rds.x.astype(int), y=rds.y.astype(int)).sort_values(["y", "x"]).reset_index(drop=True)
    same = len(rec) == len(rds) and bool((rec.x.values == rds.x.values).all() and (rec.y.values == rds.y.values).all() and (rec.cluster.values == rds.cluster.values).all())
    clustered_not_in_tissue = int((d.in_tissue != 1).sum())
    in_tissue_not_clustered = int(((pos.in_tissue == 1) & ~pos.barcode.isin(cl.Barcode)).sum())
    res["recon"][tag] = {"n_clusters_csv": int(len(cl)), "n_pos": int(len(pos)), "n_in_tissue": int((pos.in_tissue == 1).sum()),
                         "n_merged": int(len(d)), "n_rds": int(len(rds)), "identical_xy_cluster": same,
                         "clustered_not_in_tissue": clustered_not_in_tissue, "in_tissue_not_clustered": in_tissue_not_clustered,
                         "barcodes_csv_not_in_pos": int((~cl.Barcode.isin(pos.barcode)).sum()),
                         "array_row_range": [int(pos.array_row.min()), int(pos.array_row.max())],
                         "array_col_range": [int(pos.array_col.min()), int(pos.array_col.max())]}
    # nullo: stessa griglia, stessa composizione, etichette diverse
    nul = pd.read_parquet(f"{EXP}/null_full_{tag}.parquet")
    nul = nul.assign(x=nul.x.astype(int), y=nul.y.astype(int)).sort_values(["y", "x"]).reset_index(drop=True)
    res["recon"][tag]["null_same_xy"] = bool(len(nul) == len(rds) and (nul.x.values == rds.x.values).all() and (nul.y.values == rds.y.values).all())
    res["recon"][tag]["null_same_composition"] = bool((nul.cluster.value_counts().sort_index() == rds.cluster.value_counts().sort_index()).all())
    res["recon"][tag]["null_frac_label_unchanged"] = float((nul.cluster.values == rds.cluster.values).mean())
    pr = rds.cluster.value_counts(normalize=True)
    res["recon"][tag]["null_expected_frac_unchanged"] = float((pr ** 2).sum())
    print(tag, res["recon"][tag], flush=True)
# I6
i6 = pd.read_parquet(f"{EXP}/I6_syn6800x6500_c2.parquet")
x = i6.x.astype(int).values; y = i6.y.astype(int).values
res["I6"] = {"n": int(len(i6)), "x_min": int(x.min()), "x_max": int(x.max()), "y_min": int(y.min()), "y_max": int(y.max()),
             "x_mod6": sorted(set((x - 1) % 6)), "y_mod6": sorted(set((y - 1) % 6)),
             "n_unique_x": int(len(np.unique(x))), "n_unique_y": int(len(np.unique(y))),
             "gcd_x": int(np.gcd.reduce(np.diff(np.unique(x)))), "gcd_y": int(np.gcd.reduce(np.diff(np.unique(y)))),
             "full_grid_cells": int(len(np.unique(x)) * len(np.unique(y)))}
print("I6", res["I6"])
json.dump(res, open(f"{OUT}/sources_check.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))

