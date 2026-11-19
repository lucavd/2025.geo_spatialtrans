"""tools/F3.3_sim_cg.py - F3.3, ramo condizionale del punto 4 di DEC-097: le 20 simulazioni di F3.1 del tumore e dello stroma (5 aree x 2 semi)
con i soli 18 024 geni comuni del tumore del colon. Copia bins, cells, origin, profiles; conteggi e geni ristretti; poi bin2cell e RCTD
di F3.1 (lato simulato) e misure con le etichette trasferite. Uso: .venv/bin/python tools/F3.3_sim_cg.py prepare; poi tools/F3.3_sim_cg.sh"""
import sys, os, glob, shutil, numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq
M = "/mnt/micron/geo_spatialtrans"; O = f"{M}/F3.3/sim_cg"; keep = set(pd.read_csv("results/F3.3/F3.3_crc_common_genes.txt", header=None)[0].astype(str))
for d in sorted(glob.glob(f"{M}/F3.1/main/sim/A[23]_*_k[01]")):
    o = f"{O}/{os.path.basename(d)}"; os.makedirs(o, exist_ok=True)
    if os.path.exists(f"{o}/X.parquet"): continue
    for f in ("bins.parquet", "cells.parquet", "origin.parquet", "profiles.parquet", "background.parquet"):
        if os.path.exists(f"{d}/{f}"): shutil.copy2(f"{d}/{f}", o)
    g = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy(); ok = np.isin(g, list(keep)); new = -np.ones(len(g), np.int64); new[ok] = np.arange(ok.sum())
    T = pd.read_parquet(f"{d}/X.parquet"); T = T[ok[T.gene.to_numpy()]].copy(); T["gene"] = new[T.gene.to_numpy()]
    pq.write_table(pa.Table.from_pandas(T.reset_index(drop=True)), f"{o}/X.parquet"); pd.Series(g[ok]).to_csv(f"{o}/genes.tsv", index=False, header=False)
    print(os.path.basename(d), int(ok.sum()), flush=True)
