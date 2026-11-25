"""tools/sonda/S8_b2c.py - S8 (docs/sonda/REGOLA_S8.md): metodo b2c, bin2cell 0.3.5 come tools/F3.1_b2c.py (destripe, expand_labels
max_bin_distance = 2) sui nuclei `label` dell'unita'; scrive <out>/exp_b2c.npy (etichetta allargata per quadratino, ordine di bins.parquet).
Uso: .venv-b2c/bin/python -B tools/sonda/S8_b2c.py <area_dir> <out_dir>"""
import sys, os, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, bin2cell as b2c
d, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
if os.path.exists(f"{out}/exp_b2c.npy"): sys.exit(0)
bins = pd.read_parquet(f"{d}/bins.parquet"); genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
T = pd.read_parquet(f"{d}/X.parquet"); X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(genes)))
a = ad.AnnData(X.astype(np.float32)); a.obs_names = bins.bin_id.astype(str).to_numpy(); a.var_names = pd.Index(genes).astype(str); a.var_names_make_unique()
a.obs["array_row"] = bins.row.to_numpy(); a.obs["array_col"] = bins.col.to_numpy(); a.obs["labels"] = bins.label.to_numpy().astype(np.int64)
a.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
b2c.destripe(a); b2c.expand_labels(a, labels_key="labels", expanded_labels_key="labels_expanded", max_bin_distance=2)
np.save(f"{out}/exp_b2c.npy", a.obs["labels_expanded"].to_numpy().astype(np.int64))
