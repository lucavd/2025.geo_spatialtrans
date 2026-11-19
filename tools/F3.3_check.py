"""tools/F3.3_check.py - F3.3, cancello G-4: invarianti per ogni esecuzione della catena (F3.3/runs). (1) UMI in tessuto di export.json =
somma ricalcolata dai conteggi sui quadratini in_reg = umi_reg di b2c_label.json; (2) ogni blocco di 8 um = somma dei suoi 16 quadratini
(blocks_counts.mtx contro i conteggi per quadratino); (3) blocchi passati a RCTD = blocchi con almeno 8 quadratini in tessuto; (4) numero
di geni: 18 024 comuni nel tumore del colon (_cg), il pannello della sezione altrove. Uso: .venv/bin/python tools/F3.3_check.py"""
import glob, json, os, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio, h5py
R = "/mnt/micron/geo_spatialtrans/F3.3/runs"; AR = pd.read_csv("results/F3.3/F3.3_areas.csv").set_index("area"); rows = []
NCG = sum(1 for _ in open("results/F3.3/F3.3_crc_common_genes.txt"))
for d in sorted(glob.glob(f"{R}/*")):
    ex = json.load(open(f"{d}/export.json")); bj = json.load(open(f"{d}/b2c_label.json"))
    bins = pd.read_parquet(f"{d}/bins.parquet"); T = pd.read_parquet(f"{d}/X.parquet"); g = pd.read_csv(f"{d}/genes.tsv", header=None)[0]
    X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(g)))
    umi = np.asarray(X.sum(axis=1)).ravel(); reg = bins.in_reg.to_numpy()
    c1 = int(umi[reg].sum()) == ex["umi_reg"] == int(bj["umi_reg"])
    r8, c8 = bins.row.to_numpy() // 4, bins.col.to_numpy() // 4; key = r8.astype(np.int64) * 100000 + c8
    bk, inv = np.unique(key, return_inverse=True); nreg = np.bincount(inv, weights=reg.astype(float), minlength=len(bk)); keepb = nreg >= 8
    Bm = sp.csr_matrix((np.ones(len(bins)), (inv, np.arange(len(bins)))), shape=(len(bk), len(bins))); XB = (Bm @ X)[keepb]
    MB = sio.mmread(f"{d}/blocks_counts.mtx").tocsc().T.tocsr()
    c2 = MB.shape == XB.shape and (abs(MB - XB)).sum() == 0
    blk = pd.read_csv(f"{d}/blocks.csv"); rj = json.load(open(f"{d}/rctd.json"))
    c3 = len(blk) == int(keepb.sum()) == int(rj["n_blocks"])
    a = AR.loc[ex["area"]]
    if ex["genes"] == "common": c4 = len(g) == NCG
    else:
        with h5py.File(f"{a.sq}/square_002um/filtered_feature_bc_matrix.h5", "r") as h: c4 = len(g) == len(h["matrix"]["features"]["id"])
    rows.append({"run": os.path.basename(d), "umi_reg_identity": c1, "blocks_sum_children": c2, "blocks_to_rctd": c3, "n_genes_expected": c4, "n_genes": len(g)})
D = pd.DataFrame(rows); D.to_csv("results/F3.3/F3.3_G4.csv", index=False)
bad = D[~D[["umi_reg_identity", "blocks_sum_children", "blocks_to_rctd", "n_genes_expected"]].all(axis=1)]
print(f"{len(D)} esecuzioni, {len(bad)} con un controllo fallito"); print(bad.to_string())
