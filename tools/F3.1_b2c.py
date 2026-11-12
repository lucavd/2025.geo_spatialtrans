"""tools/F3.1_b2c.py - F3.1: assegnazione dei quadratini alle cellule con bin2cell e misure uguali sui due lati; blocchi di 8 um per RCTD.
Metodo: bin2cell 0.3.5 come nel flusso standard: destripe(), poi expand_labels(max_bin_distance = 2, cioe' 4 um, pareggi decisi
con la PCA dell'espressione), con i nuclei dati (vero: maschera adottata in R6; simulato: nucleo vero, copertura >= 1/2).
Le misure usano i conteggi grezzi (UMI) sommati per cellula, come bin_to_cell() sui conteggi grezzi.
Uso: .venv-b2c/bin/python tools/F3.1_b2c.py <area_dir> <ref_dir> <archetype> <side: sim|real> [label | label_dimg | label_shift]
Scrive in <area_dir>/b2c_<label_col>.json, cells_<label_col>.parquet; e (solo con label) blocks.csv, blocks_counts.mtx."""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
import anndata as ad, bin2cell as b2c
d, ref_dir, A, side = sys.argv[1:5]; lcol = sys.argv[5] if len(sys.argv) > 5 else "label"
t0 = time.time()
bins = pd.read_parquet(f"{d}/bins.parquet"); genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
T = pd.read_parquet(f"{d}/X.parquet")
X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(genes)))
cells = pd.read_parquet(f"{d}/cells.parquet")
mk = pd.read_csv(f"{ref_dir}/markers_{A}.csv"); types = sorted(mk.cell_type.unique())
reg = bins.in_reg.to_numpy()
if lcol == "label_shift":   # controprova: nuclei spostati di 10 quadratini (20 um) in diagonale sulla griglia, come il controllo 'shift' di R6
    SH = 10; key = pd.Series(bins.label.to_numpy(), index=bins.row.to_numpy().astype(np.int64) * 100000 + bins.col.to_numpy())
    nuc = key.reindex((bins.row.to_numpy().astype(np.int64) - SH) * 100000 + (bins.col.to_numpy() - SH)).fillna(0).to_numpy().astype(np.int64)
else:
    nuc = bins[lcol].to_numpy().astype(np.int64)
umi = np.asarray(X.sum(axis=1)).ravel()
# --- bin2cell ---
a = ad.AnnData(X.astype(np.float32)); a.obs_names = bins.bin_id.astype(str).to_numpy(); a.var_names = pd.Index(genes).astype(str)
a.var_names_make_unique()
a.obs["array_row"] = bins.row.to_numpy(); a.obs["array_col"] = bins.col.to_numpy(); a.obs["labels"] = nuc; a.obs["n_counts"] = umi
b2c.destripe(a)
b2c.expand_labels(a, labels_key="labels", expanded_labels_key="labels_expanded", max_bin_distance=2)
exp = a.obs["labels_expanded"].to_numpy().astype(np.int64)
# --- cellule dell'area: nucleo con almeno meta' dei suoi quadratini in reg ---
lab_ids, inv = np.unique(nuc[nuc > 0], return_inverse=True)
nb_all = np.bincount(inv, minlength=len(lab_ids)); nb_reg = np.bincount(inv, weights=reg[nuc > 0].astype(float), minlength=len(lab_ids))
cell_set = lab_ids[nb_reg >= 0.5 * nb_all]
ctype = cells.set_index("label").cell_type.reindex(cell_set).to_numpy()
# marcatori: matrice geni x tipi (indicatrici)
gidx = pd.Series(np.arange(len(genes)), index=genes)
assert gidx.index.is_unique, "gene ids not unique"
mk = mk[mk.gene.isin(gidx.index)]
Mmat = np.zeros((len(genes), len(types)))
for k, ty in enumerate(types): Mmat[gidx.reindex(mk.gene[mk.cell_type == ty]).dropna().astype(int).to_numpy(), k] = 1
def per_cell(lab):
    sel = np.isin(lab, cell_set)
    pos = pd.Series(np.arange(len(cell_set)), index=cell_set)
    I = sp.csr_matrix((np.ones(sel.sum()), (pos.reindex(lab[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(cell_set), len(bins)))
    C = I @ X; u = np.asarray(C.sum(axis=1)).ravel(); nb = np.asarray(I.sum(axis=1)).ravel()
    MU = np.asarray(C @ Mmat)  # UMI su marcatori di ciascun tipo
    return C, u, nb, MU
res = {"area": os.path.basename(d), "archetype": A, "side": side, "label_col": lcol, "n_bins": int(len(bins)), "n_reg": int(reg.sum()),
       "umi_reg": float(umi[reg].sum()), "n_cells": int(len(cell_set)), "types": types}
tab = pd.DataFrame({"label": cell_set, "cell_type": ctype})
for nm, lab in [("nuc", nuc), ("exp", exp)]:
    res[f"share_assigned_{nm}"] = float(umi[reg & (lab > 0)].sum() / umi[reg].sum())
    res[f"bins_assigned_{nm}"] = float(((lab > 0) & reg).mean() / reg.mean())
    C, u, nb, MU = per_cell(lab)
    res[f"umi_per_cell_{nm}"] = float(u.mean()); res[f"bins_per_cell_{nm}"] = float(nb.mean())
    tot = MU.sum(axis=1); own = np.full(len(cell_set), np.nan)
    for k, ty in enumerate(types): own[ctype == ty] = MU[ctype == ty, k]
    typed = np.isin(ctype, types) & (tot >= 5)
    foreign = 1 - own / np.where(tot > 0, tot, np.nan)
    res[f"foreign_{nm}"] = float(np.nanmean(foreign[typed])); res[f"n_typed_{nm}"] = int(typed.sum())
    tab[f"umi_{nm}"] = u; tab[f"bins_{nm}"] = nb; tab[f"foreign_{nm}"] = np.where(typed, foreign, np.nan)
    if side == "sim":
        og = pd.read_parquet(f"{d}/origin.parquet")
        bpos = pd.Series(np.arange(len(bins)), index=bins.bin_id.to_numpy()); ob = bpos.reindex(og.bin_id).to_numpy()
        ol = og.label.to_numpy(); on = og.n.to_numpy().astype(float)
        asg = lab[ob]; inreg = reg[ob]
        right = (asg == ol) & (asg > 0) & inreg
        res[f"truth_purity_{nm}"] = float(on[right].sum() / umi[reg & (lab > 0)].sum())            # UMI assegnati che vengono dalla cellula giusta
        inc = np.isin(ol, cell_set)
        res[f"truth_recall_{nm}"] = float(on[inc & (asg == ol)].sum() / on[inc].sum())               # UMI catturati delle cellule dell'area finiti nella loro cellula
        ty_of = cells.set_index("label").cell_type
        oty = ty_of.reindex(ol).to_numpy(); aty = ty_of.reindex(np.where(asg > 0, asg, -1)).to_numpy()
        m = (asg > 0) & inreg & np.isin(asg, cell_set)
        res[f"truth_foreign_type_{nm}"] = float(on[m & (oty != aty)].sum() / umi[reg & np.isin(lab, cell_set)].sum())
if side == "sim" and lcol == "label" and os.path.exists(f"{d}/profiles.parquet"):   # misura esplicativa: marcatori estranei nei profili veri generati
    pf = pd.read_parquet(f"{d}/profiles.parquet"); pf = pf[pf.label.isin(cell_set)]
    pos = pd.Series(np.arange(len(cell_set)), index=cell_set)
    P = sp.csr_matrix((pf.n.to_numpy().astype(float), (pos.reindex(pf.label).to_numpy(), pf.gene.to_numpy())), shape=(len(cell_set), len(genes)))
    MU = np.asarray(P @ Mmat); tot = MU.sum(axis=1); own = np.full(len(cell_set), np.nan)
    for k, ty in enumerate(types): own[ctype == ty] = MU[ctype == ty, k]
    typed = np.isin(ctype, types) & (tot >= 5); fp = 1 - own / np.where(tot > 0, tot, np.nan)
    res["foreign_profile"] = float(np.nanmean(fp[typed])); res["n_typed_profile"] = int(typed.sum()); tab["foreign_profile"] = np.where(typed, fp, np.nan)
json.dump(res, open(f"{d}/b2c_{lcol}.json", "w")); tab.to_parquet(f"{d}/cells_{lcol}.parquet", index=False)
# --- blocchi di 8 um per RCTD (una volta, con la label principale) ---
if lcol == "label":
    r8, c8 = bins.row.to_numpy() // 4, bins.col.to_numpy() // 4
    key = r8.astype(np.int64) * 100000 + c8
    bk, binv = np.unique(key, return_inverse=True)
    nreg = np.bincount(binv, weights=reg.astype(float), minlength=len(bk))
    keepb = nreg >= 8
    Bm = sp.csr_matrix((np.ones(len(bins)), (binv, np.arange(len(bins)))), shape=(len(bk), len(bins)))
    XB = (Bm @ X)[keepb]
    blk = pd.DataFrame({"block": [f"b{k}" for k in bk[keepb]], "r8": (bk[keepb] // 100000), "c8": (bk[keepb] % 100000), "n_reg": nreg[keepb]})
    xs = np.bincount(binv, weights=bins.x_um.to_numpy(), minlength=len(bk)) / np.bincount(binv, minlength=len(bk))
    ys = np.bincount(binv, weights=bins.y_um.to_numpy(), minlength=len(bk)) / np.bincount(binv, minlength=len(bk))
    blk["x"] = xs[keepb]; blk["y"] = ys[keepb]; blk["umi"] = np.asarray(XB.sum(axis=1)).ravel()
    # nucleo dominante del blocco (quadratini del solo nucleo): almeno 4 dei 16 quadratini
    ty_of = cells.set_index("label").cell_type
    dl = pd.DataFrame({"b": binv, "l": nuc})
    dl = dl[dl.l > 0].groupby(["b", "l"]).size().reset_index(name="n").sort_values(["b", "n"], ascending=[True, False]).drop_duplicates("b")
    dom = pd.Series(dl.l.to_numpy(), index=dl.b.to_numpy()); domn = pd.Series(dl.n.to_numpy(), index=dl.b.to_numpy())
    idx = np.flatnonzero(keepb)
    blk["nuc_label"] = dom.reindex(idx).fillna(0).astype(np.int64).to_numpy(); blk["nuc_bins"] = domn.reindex(idx).fillna(0).to_numpy()
    blk["nuc_type"] = np.where(blk.nuc_bins >= 4, ty_of.reindex(blk.nuc_label).to_numpy(), None)
    if side == "sim":
        og = pd.read_parquet(f"{d}/origin.parquet"); bpos = pd.Series(np.arange(len(bins)), index=bins.bin_id.to_numpy())
        og["b"] = binv[bpos.reindex(og.bin_id).to_numpy()]; og["ty"] = ty_of.reindex(og.label).to_numpy()
        g = og.groupby(["b", "ty"]).n.sum().reset_index().sort_values(["b", "n"], ascending=[True, False])
        tot = g.groupby("b").n.sum(); top = g.drop_duplicates("b").set_index("b")
        blk["true_type"] = top.ty.reindex(idx).to_numpy(); blk["true_share"] = (top.n.reindex(idx) / tot.reindex(idx)).to_numpy()
        blk["cell_umi"] = tot.reindex(idx).fillna(0).to_numpy()
    blk.to_csv(f"{d}/blocks.csv", index=False)
    sio.mmwrite(f"{d}/blocks_counts.mtx", XB.T.tocsc().astype(np.int64))
res["sec"] = round(time.time() - t0, 1); json.dump(res, open(f"{d}/b2c_{lcol}.json", "w"))
print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items() if k != "types"}), flush=True)
