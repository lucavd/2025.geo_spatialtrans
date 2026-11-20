"""tools/riesame_xen_targets.py - Riesame Parte B, bracci X1/XC: bersaglio di V (UMI per quadratino in_reg e in_tissue sui 408 geni del pannello)
dalle esportazioni di F3.1 delle 10 aree di taratura di P2, e quota dei 408 geni negli UMI delle simulazioni S di F3.1 (per il fondo di XC, XD-19).
Solo lettura; md5 degli ingressi prima e dopo in results/riesame/xen_targets_md5.tsv. Uso: .venv/bin/python -B tools/riesame_xen_targets.py"""
import hashlib, json, numpy as np, pandas as pd
M = "/mnt/micron/geo_spatialtrans"; F = f"{M}/F3.1/main"
AREAS = [f"A2_r{i}" for i in range(1, 6)] + [f"A3_r{i}" for i in range(1, 6)]
PANEL = pd.read_csv(f"{M}/riesame/xenium/panel_ids_in_common.txt", header=None)[0].astype(str).tolist()
def md5(p):
    h = hashlib.md5(); f = open(p, "rb")
    for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()
files = []
for a in AREAS:
    files += [f"{F}/real/{a}/{x}" for x in ("X.parquet", "bins.parquet", "genes.tsv")]
    for k in (0, 1): files += [f"{F}/sim/{a}_k{k}/{x}" for x in ("X.parquet", "genes.tsv")]
before = {p: md5(p) for p in files}
rows = []
for a in AREAS:
    d = f"{F}/real/{a}"; g = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy(); b = pd.read_parquet(f"{d}/bins.parquet")
    gi = np.flatnonzero(np.isin(g, PANEL)); assert len(gi) == 408, (a, len(gi))
    T = pd.read_parquet(f"{d}/X.parquet"); Tp = T[T.gene.isin(gi)]
    u = np.bincount(Tp.bin.to_numpy(), weights=Tp.n.to_numpy(), minlength=len(b)); ua = np.bincount(T.bin.to_numpy(), weights=T.n.to_numpy(), minlength=len(b))
    reg = b.in_reg.to_numpy(); tis = b.in_tissue.to_numpy()
    r = {"area": a, "tissue": a[:2], "n_bins": len(b), "n_reg": int(reg.sum()), "n_in_tissue": int(tis.sum()),
         "v_umi408_reg": float(u[reg].sum()), "v_mean408_reg": float(u[reg].mean()), "v_mean408_in_tissue": float(u[tis].mean()),
         "v_meanall_reg": float(ua[reg].mean()), "v_panel_share_reg": float(u[reg].sum() / ua[reg].sum())}
    for k in (0, 1):
        s = f"{F}/sim/{a}_k{k}"; gs = pd.read_csv(f"{s}/genes.tsv", header=None)[0].astype(str).to_numpy()
        gsi = np.flatnonzero(np.isin(gs, PANEL)); S = pd.read_parquet(f"{s}/X.parquet")
        r[f"s_panel_share_k{k}"] = float(S.n[S.gene.isin(gsi)].sum() / S.n.sum()); r[f"s_n_panel_genes_k{k}"] = int(len(gsi))
    rows.append(r)
df = pd.DataFrame(rows); df.to_csv("results/riesame/xen_targets.csv", index=False)
after = {p: md5(p) for p in files}
pd.DataFrame({"file": files, "md5_before": [before[p] for p in files], "md5_after": [after[p] for p in files]}).to_csv("results/riesame/xen_targets_md5.tsv", sep="\t", index=False)
print("md5 unchanged:", all(before[p] == after[p] for p in files), len(files))
tt = df.groupby("tissue").agg(v_target=("v_mean408_reg", "mean"), v_target_in_tissue=("v_mean408_in_tissue", "mean"),
                              s_share=("s_panel_share_k0", "mean"), s_share1=("s_panel_share_k1", "mean"), v_share=("v_panel_share_reg", "mean"))
print(df[["area", "n_reg", "v_mean408_reg", "v_mean408_in_tissue", "v_panel_share_reg", "s_panel_share_k0", "s_panel_share_k1", "s_n_panel_genes_k0"]].to_string())
print(tt.to_string())
