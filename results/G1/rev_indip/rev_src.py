import pandas as pd, numpy as np, json
G = "/mnt/micron/geo_spatialtrans/"; rows = []
for arm in ("Yalt", "YH1", "Yshuf"):
    for t in ("tumour", "stroma"):
        for k in range(1, 6):
            d = f"{G}S6_sonda/runs/{arm}/H2_{t}_{k}_k0"
            exp_src = {"Yalt": f"{G}S3_sonda/runs/V/H2_{t}_{k % 5 + 1}", "YH1": f"{G}S3b_sonda/runs/T/H1_{t}_{k}", "Yshuf": f"{G}S3_sonda/runs/V/H2_{t}_{k}"}[arm]
            s = pd.read_parquet(d + "/cells_src.parquet"); c = pd.read_parquet(d + "/cells.parquet").set_index("label")
            dirs = s.src_dir.unique()
            tl = pd.read_parquet(exp_src + "/labels_transferred.parquet").set_index("label").cell_type
            src_type = tl.reindex(s.src_real_label.to_numpy()).to_numpy(); sim_type = c.cell_type.reindex(s.label.to_numpy()).to_numpy()
            eq = pd.Series(src_type).fillna("NA").to_numpy() == pd.Series(sim_type).fillna("NA").to_numpy()
            # Y (S4) on same area: same nuclei count?
            y = pd.read_parquet(f"{G}S4_sonda/runs/Y/H2_{t}_{k}_k0/cells.parquet")
            comp_sim = pd.Series(sim_type).value_counts(normalize=True); comp_src = tl.dropna().value_counts(normalize=True)
            rows.append({"arm": arm, "A": t, "k": k, "src_ok": bool(len(dirs) == 1 and dirs[0] == exp_src), "n_cells": len(s), "type_eq_src": float(eq.mean()),
                         "n_cells_Y": len(y), "same_labels_as_Y": bool(np.array_equal(np.sort(y.label.to_numpy()), np.sort(c.index.to_numpy())))})
R = pd.DataFrame(rows); print(R.to_string(index=False)); print("all src_ok", R.src_ok.all(), "min type_eq", R.type_eq_src.min())
