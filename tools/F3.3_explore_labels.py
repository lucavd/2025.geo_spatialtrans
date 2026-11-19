"""tools/F3.3_explore_labels.py - F3.3, esplorazione E-2 (solo dati di taratura e simulazioni di F3.1): etichette dei nuclei per
trasferimento dal riferimento di F3.1 invece del raggruppamento di R6 (vero) o della verita' (simulato), e loro effetto su M-3 e M-5.
Etichettatore: profilo medio multinomiale per tipo dalle cellule del riferimento (somma dei conteggi + 0.5 per gene), tipo di massima
verosimiglianza sui conteggi della cellula di bin2cell (nucleo allargato di 4 um), cellule con almeno 20 UMI; sul vero di taratura i
profili escludono le cellule del riferimento della stessa area (fuori area). bin2cell rifatto come tools/F3.1_b2c.py (controllo di
identita' con b2c_label.json). Uso: .venv-b2c/bin/python tools/F3.3_explore_labels.py <area_dir> <A> <sim|real> <out_json>"""
import sys, os, json, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
import anndata as ad, bin2cell as b2c
d, A, side, outp = sys.argv[1:5]
REF = "/mnt/micron/geo_spatialtrans/F3.1/main/ref"
bins = pd.read_parquet(f"{d}/bins.parquet"); genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
T = pd.read_parquet(f"{d}/X.parquet"); X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(genes)))
cells = pd.read_parquet(f"{d}/cells.parquet"); mk = pd.read_csv(f"{REF}/markers_{A}.csv"); types = sorted(mk.cell_type.unique())
reg = bins.in_reg.to_numpy(); nuc = bins["label"].to_numpy().astype(np.int64); umi = np.asarray(X.sum(axis=1)).ravel()
a = ad.AnnData(X.astype(np.float32)); a.obs_names = bins.bin_id.astype(str).to_numpy(); a.var_names = pd.Index(genes).astype(str); a.var_names_make_unique()
a.obs["array_row"] = bins.row.to_numpy(); a.obs["array_col"] = bins.col.to_numpy(); a.obs["labels"] = nuc; a.obs["n_counts"] = umi
b2c.destripe(a); b2c.expand_labels(a, labels_key="labels", expanded_labels_key="labels_expanded", max_bin_distance=2)
exp = a.obs["labels_expanded"].to_numpy().astype(np.int64)
lab_ids, inv = np.unique(nuc[nuc > 0], return_inverse=True)
nb_all = np.bincount(inv, minlength=len(lab_ids)); nb_reg = np.bincount(inv, weights=reg[nuc > 0].astype(float), minlength=len(lab_ids))
cell_set = lab_ids[nb_reg >= 0.5 * nb_all]
pos = pd.Series(np.arange(len(cell_set)), index=cell_set); sel = np.isin(exp, cell_set)
I = sp.csr_matrix((np.ones(sel.sum()), (pos.reindex(exp[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(cell_set), len(bins)))
C = (I @ X).tocsr(); u = np.asarray(C.sum(axis=1)).ravel()
gidx = pd.Series(np.arange(len(genes)), index=genes)
Mmat = np.zeros((len(genes), len(types)))
for k, ty in enumerate(types): Mmat[gidx.reindex(mk.gene[mk.cell_type == ty]).dropna().astype(int).to_numpy(), k] = 1
MU = np.asarray(C @ Mmat); tot = MU.sum(axis=1)
def m3(ct):
    own = np.full(len(cell_set), np.nan)
    for k, ty in enumerate(types): own[ct == ty] = MU[ct == ty, k]
    typed = np.isin(ct, types) & (tot >= 5); f = 1 - own / np.where(tot > 0, tot, np.nan)
    return float(np.nanmean(f[typed])), int(typed.sum())
orig = cells.set_index("label").cell_type.reindex(cell_set).to_numpy()
# riferimento
Xr = sio.mmread(f"{REF}/ref_{A}_counts.mtx").tocsc(); rg = pd.read_csv(f"{REF}/ref_{A}_genes.tsv", header=None)[0].astype(str).to_numpy()
rc = pd.read_csv(f"{REF}/ref_{A}_cells.csv"); roi = os.path.basename(d).split("_")[1]
use = np.ones(len(rc), bool)
if side == "real": use = (rc.roi_id.astype(str) != roi).to_numpy()
gi = gidx.reindex(rg).to_numpy(); okg = ~np.isnan(gi); gi = gi[okg].astype(int)
P = np.zeros((len(types), len(gi)))
for k, ty in enumerate(types):
    cols = np.flatnonzero(use & (rc.cell_type == ty).to_numpy())
    P[k] = np.asarray(Xr[:, cols][okg].sum(axis=1)).ravel() + 0.5
L = np.log(P / P.sum(axis=1, keepdims=True))
Cg = C[:, gi]
ll = np.asarray(Cg @ L.T); tl = np.array(types, dtype=object)[ll.argmax(axis=1)]; tl[u < 20] = None
res = {"area": os.path.basename(d), "A": A, "side": side, "n_cells": int(len(cell_set)), "n_ref_used": int(use.sum())}
m3o, no = m3(orig); m3t, nt = m3(tl)
res.update(M3_orig=m3o, n_typed_orig=no, M3_transfer=m3t, n_typed_transfer=nt)
bj = json.load(open(f"{d}/b2c_label.json")); res["M3_f31"] = bj["foreign_exp"]; res["b2c_identical"] = abs(bj["foreign_exp"] - m3o) < 1e-12 and bj["n_cells"] == len(cell_set)
cmp = np.isin(orig, types) & (u >= 20)
res["agree_with_orig"] = float((tl[cmp] == orig[cmp]).mean()); res["n_cmp"] = int(cmp.sum())
for ty in types:
    k = cmp & (orig == ty); res[f"recall_{ty}"] = float((tl[k] == ty).mean()) if k.any() else None
    res[f"share_orig_{ty}"] = float((orig[cmp] == ty).mean()); res[f"share_transfer_{ty}"] = float((tl[cmp] == ty).mean())
# M-5 con il tipo trasferito del nucleo dominante del blocco
rb = pd.read_csv(f"{d}/rctd_blocks.csv"); tmap = pd.Series(tl, index=cell_set)
run = rb.spot_class.notna(); ok = run & (rb.spot_class != "reject")
def m5(nt):
    hn = ok & nt.notna() & nt.isin(types) & (rb.nuc_bins >= 4)
    return float((rb.first[hn] == nt[hn]).mean()), int(hn.sum())
res["M5_orig"], res["n5_orig"] = m5(rb.nuc_type.where(rb.nuc_type.isin(types)))
res["M5_transfer"], res["n5_transfer"] = m5(pd.Series(tmap.reindex(rb.nuc_label).to_numpy(), index=rb.index))
res["M5_f31"] = json.load(open(f"{d}/rctd.json"))["agree_nucleus"]
os.makedirs(os.path.dirname(outp), exist_ok=True); json.dump(res, open(outp, "w"))
print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items() if not k.startswith(("recall", "share"))}), flush=True)
