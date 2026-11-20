"""tools/riesame_panel_restrict.py - riesame, Parte B: restringe un'esportazione di un'area (formato di F3.1/F3.3: bins.parquet,
X.parquet, genes.tsv, cells.parquet) a un insieme di geni, conservando l'ordine dei geni e i quadratini. Insiemi: panel = i 408 geni
del pannello Xenium dentro i 18 024 comuni (riesame/xenium/panel_ids_in_common.txt); common = i 18 024 geni comuni del tumore del
colon (results/F3.3/F3.3_crc_common_genes.txt); all = nessuna restrizione. Per il simulato con genes = all copia anche origin,
profiles e background (verita' usata da tools/F3.1_b2c.py con side = sim); con una restrizione non li copia, perche' origin non ha i
geni e non si puo' restringere (si usa side = simr: le misure M-1...M-7 non cambiano, cambiano solo le misure di verita').
Uso: .venv/bin/python tools/riesame_panel_restrict.py <src_dir> <out_dir> <panel|common|all> <real|sim>
Scrive <out_dir>/{bins,X,cells}.parquet, genes.tsv, chain_restrict.json (md5 delle sorgenti, geni, UMI)."""
import sys, os, json, time, shutil, hashlib, numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq
GENESETS = {"panel": "/mnt/micron/geo_spatialtrans/riesame/xenium/panel_ids_in_common.txt",
            "common": os.path.expanduser("~/2025.geo_spatialtrans/results/F3.3/F3.3_crc_common_genes.txt")}
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def restrict(src, out, mode, side):
    t0 = time.time(); os.makedirs(out, exist_ok=True)
    srcs = [f for f in ("bins.parquet", "X.parquet", "genes.tsv", "cells.parquet") if os.path.exists(f"{src}/{f}")]
    m_before = {f: md5(f"{src}/{f}") for f in srcs}
    g = pd.read_csv(f"{src}/genes.tsv", header=None)[0].astype(str).to_numpy()
    if mode == "all": ok = np.ones(len(g), bool)
    else:
        keep = set(pd.read_csv(GENESETS[mode], header=None)[0].astype(str)); ok = np.isin(g, list(keep))
        miss = len(keep - set(g))
        if miss: raise SystemExit(f"{miss} genes of set '{mode}' missing from {src}/genes.tsv")
    new = -np.ones(len(g), np.int64); new[ok] = np.arange(ok.sum())
    T = pd.read_parquet(f"{src}/X.parquet"); umi_all = int(T.n.sum())
    T = T[ok[T.gene.to_numpy()]].copy(); T["gene"] = new[T.gene.to_numpy()]
    pq.write_table(pa.Table.from_pandas(T.reset_index(drop=True), preserve_index=False), f"{out}/X.parquet")
    pd.Series(g[ok]).to_csv(f"{out}/genes.tsv", index=False, header=False)
    for f in ("bins.parquet", "cells.parquet"): shutil.copy2(f"{src}/{f}", out)
    if side == "sim" and mode == "all":
        for f in ("origin.parquet", "profiles.parquet", "background.parquet"):
            if os.path.exists(f"{src}/{f}"): shutil.copy2(f"{src}/{f}", out); m_before[f] = md5(f"{src}/{f}")
    bins = pd.read_parquet(f"{out}/bins.parquet"); umi = np.bincount(T.bin.to_numpy(), weights=T.n.to_numpy(), minlength=len(bins))
    info = {"src": src, "genes": mode, "side": side, "n_genes_src": int(len(g)), "n_genes": int(ok.sum()), "umi_src": umi_all,
            "umi": int(T.n.sum()), "umi_reg": int(umi[bins.in_reg.to_numpy()].sum()), "n_bins": int(len(bins)),
            "md5_src": m_before, "md5_src_after": {f: md5(f"{src}/{f}") for f in m_before}, "sec": round(time.time() - t0, 2)}
    info["src_unchanged"] = info["md5_src"] == info["md5_src_after"]
    json.dump(info, open(f"{out}/chain_restrict.json", "w"), indent=1)
    return info
if __name__ == "__main__":
    src, out, mode, side = sys.argv[1:5]
    i = restrict(src, out, mode, side); print(json.dumps({k: v for k, v in i.items() if not k.startswith("md5")}), flush=True)
