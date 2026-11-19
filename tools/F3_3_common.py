"""tools/F3_3_common.py - F3.3: etichette trasferite e misure M-1...M-7 dalle uscite di una catena (formato di F3.1), per un'area
intera o per un riquadro (quarti di 500 um del linfonodo, DEC-093 R-9). Definizioni di M-1...M-6 identiche a tools/F3.1_b2c.py e
tools/F3.1_rctd.R (controllo G-1 d: sulle aree di F3.1 con le etichette originali i valori tornano entro 1e-9); M-7 = anisotropia del
campo log(1 + UMI) dei blocchi di 8 um come tools/F3.1_PH5_aniso.py.
Etichettatore trasferito (scelta tecnica di F3.3 dopo l'esplorazione E-2): profilo multinomiale per tipo dalle cellule del
riferimento di F3.1 (somma dei conteggi + 0.5 per gene), tipo di massima verosimiglianza sui conteggi del nucleo allargato di
bin2cell (4 um), nuclei con almeno 20 UMI; per le aree vere di taratura i profili escludono le cellule del riferimento della stessa
area (fuori area), per le aree di verifica e per il simulato si usano tutte. Stessa funzione sui due lati."""
import os, json, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
REF = "/mnt/micron/geo_spatialtrans/F3.1/main/ref"; UMI_MIN_TL = 20; PSEUDO = 0.5
TH8 = np.arange(8) * np.pi / 8

def load_area(d):
    bins = pd.read_parquet(f"{d}/bins.parquet"); genes = pd.read_csv(f"{d}/genes.tsv", header=None)[0].astype(str).to_numpy()
    T = pd.read_parquet(f"{d}/X.parquet")
    X = sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=(len(bins), len(genes)))
    return bins, genes, X

def expand(d, bins, genes, X, cache):
    """bin2cell destripe() + expand_labels(max_bin_distance = 2), come tools/F3.1_b2c.py; risultato in cache (npy) fuori da <d>."""
    if cache and os.path.exists(cache): return np.load(cache)
    import anndata as ad, bin2cell as b2c
    nuc = bins["label"].to_numpy().astype(np.int64); umi = np.asarray(X.sum(axis=1)).ravel()
    a = ad.AnnData(X.astype(np.float32)); a.obs_names = bins.bin_id.astype(str).to_numpy(); a.var_names = pd.Index(genes).astype(str)
    a.var_names_make_unique(); a.obs["array_row"] = bins.row.to_numpy(); a.obs["array_col"] = bins.col.to_numpy()
    a.obs["labels"] = nuc; a.obs["n_counts"] = umi
    b2c.destripe(a); b2c.expand_labels(a, labels_key="labels", expanded_labels_key="labels_expanded", max_bin_distance=2)
    exp = a.obs["labels_expanded"].to_numpy().astype(np.int64)
    if cache: os.makedirs(os.path.dirname(cache), exist_ok=True); np.save(cache, exp)
    return exp

def ref_profiles(A, genes, exclude_roi=None):
    """log-profili (tipi x geni dell'area) dal riferimento di F3.1; geni del riferimento assenti nell'area esclusi."""
    Xr = sio.mmread(f"{REF}/ref_{A}_counts.mtx").tocsr(); rg = pd.read_csv(f"{REF}/ref_{A}_genes.tsv", header=None)[0].astype(str).to_numpy()
    rc = pd.read_csv(f"{REF}/ref_{A}_cells.csv"); types = sorted(pd.read_csv(f"{REF}/markers_{A}.csv").cell_type.unique())
    use = np.ones(len(rc), bool) if exclude_roi is None else (rc.roi_id.astype(str) != str(exclude_roi)).to_numpy()
    gidx = pd.Series(np.arange(len(genes)), index=genes); gi = gidx.reindex(rg).to_numpy(); okg = ~np.isnan(gi)
    P = np.zeros((len(types), int(okg.sum())))
    for k, ty in enumerate(types):
        cols = np.flatnonzero(use & (rc.cell_type == ty).to_numpy())
        P[k] = np.asarray(Xr[np.flatnonzero(okg)][:, cols].sum(axis=1)).ravel() + PSEUDO
    return types, gi[okg].astype(int), np.log(P / P.sum(axis=1, keepdims=True)), int(use.sum())

def transfer_labels(X, exp, labels, A, genes, exclude_roi=None):
    """tipo trasferito per ogni etichetta di nucleo in <labels> (conteggi del nucleo allargato); None sotto 20 UMI."""
    types, gi, L, nref = ref_profiles(A, genes, exclude_roi)
    labels = np.asarray(labels, np.int64); pos = pd.Series(np.arange(len(labels)), index=labels); sel = np.isin(exp, labels)
    I = sp.csr_matrix((np.ones(sel.sum()), (pos.reindex(exp[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(labels), X.shape[0]))
    C = (I @ X).tocsr(); u = np.asarray(C.sum(axis=1)).ravel()
    ll = np.asarray(C[:, gi] @ L.T); tl = np.array(types, dtype=object)[ll.argmax(axis=1)]; tl[u < UMI_MIN_TL] = None
    return pd.DataFrame({"label": labels, "cell_type": tl, "umi_exp_all": u}), nref

def aniso_umi(m):
    """M-7: anisotropia del campo log(1 + UMI) dei blocchi (max - min su 8 direzioni della correlazione media a distanza 1-3)."""
    key = pd.Series(np.arange(len(m)), index=m.r8.astype(np.int64).to_numpy() * 100000 + m.c8.to_numpy())
    lu = np.log1p(m.umi.to_numpy()); corr = []
    for th in TH8:
        offs = sorted({(int(round(d * np.sin(th))), int(round(d * np.cos(th)))) for d in (1, 2, 3)} - {(0, 0)}); c_ = []
        for dr, dc in offs:
            j = key.reindex((m.r8.astype(np.int64).to_numpy() + dr) * 100000 + m.c8.to_numpy() + dc).to_numpy(); has = ~np.isnan(j)
            i = np.flatnonzero(has); j = j[has].astype(int); c_.append(np.corrcoef(lu[i], lu[j])[0, 1])
        corr.append(np.mean(c_))
    corr = np.array(corr); return float(corr.max() - corr.min())

def measures(d, A, cell_types, exp, bins, genes, X, box=None, with_aniso=False):
    """M-1...M-6 (e M-7) come F3.1. cell_types: Series label -> tipo (etichettatura scelta). box = (x0, y0, x1, y1) in um dell'area."""
    mk = pd.read_csv(f"{REF}/markers_{A}.csv"); types = sorted(mk.cell_type.unique())
    reg = bins.in_reg.to_numpy(); nuc = bins["label"].to_numpy().astype(np.int64); umi = np.asarray(X.sum(axis=1)).ravel()
    xb, yb = bins.x_um.to_numpy(), bins.y_um.to_numpy()
    inb = np.ones(len(bins), bool) if box is None else (xb >= box[0]) & (xb < box[2]) & (yb >= box[1]) & (yb < box[3])
    lab_ids, inv = np.unique(nuc[nuc > 0], return_inverse=True)
    nb_all = np.bincount(inv, minlength=len(lab_ids)); nb_reg = np.bincount(inv, weights=reg[nuc > 0].astype(float), minlength=len(lab_ids))
    cs = nb_reg >= 0.5 * nb_all
    if box is not None:
        cx = np.bincount(inv, weights=xb[nuc > 0], minlength=len(lab_ids)) / nb_all; cy = np.bincount(inv, weights=yb[nuc > 0], minlength=len(lab_ids)) / nb_all
        cs &= (cx >= box[0]) & (cx < box[2]) & (cy >= box[1]) & (cy < box[3])
    cell_set = lab_ids[cs]; ctype = cell_types.reindex(cell_set).to_numpy()
    gidx = pd.Series(np.arange(len(genes)), index=genes); mk = mk[mk.gene.isin(gidx.index)]
    Mmat = np.zeros((len(genes), len(types)))
    for k, ty in enumerate(types): Mmat[gidx.reindex(mk.gene[mk.cell_type == ty]).dropna().astype(int).to_numpy(), k] = 1
    r = {"n_cells": int(len(cell_set))}
    rb = reg & inb
    r["M1"] = float(umi[rb & (exp > 0)].sum() / umi[rb].sum())
    sel = np.isin(exp, cell_set); pos = pd.Series(np.arange(len(cell_set)), index=cell_set)
    I = sp.csr_matrix((np.ones(sel.sum()), (pos.reindex(exp[sel]).to_numpy(), np.flatnonzero(sel))), shape=(len(cell_set), len(bins)))
    C = I @ X; u = np.asarray(C.sum(axis=1)).ravel(); MU = np.asarray(C @ Mmat); tot = MU.sum(axis=1)
    r["M2"] = float(u.mean()) if len(u) else np.nan
    own = np.full(len(cell_set), np.nan)
    for k, ty in enumerate(types): own[ctype == ty] = MU[ctype == ty, k]
    typed = np.isin(ctype, types) & (tot >= 5); fo = 1 - own / np.where(tot > 0, tot, np.nan)
    r["M3"] = float(np.nanmean(fo[typed])) if typed.any() else np.nan; r["n_typed"] = int(typed.sum())
    m = pd.read_csv(f"{d}/rctd_blocks.csv")
    if box is not None: m = m[(m.x >= box[0]) & (m.x < box[2]) & (m.y >= box[1]) & (m.y < box[3])].reset_index(drop=True)
    run = m.spot_class.notna().to_numpy(); ok = run & (m.spot_class != "reject").to_numpy()
    r["M4"] = float((m.spot_class[ok] == "singlet").mean())
    nt = pd.Series(np.where(m.nuc_bins.to_numpy() >= 4, cell_types.reindex(m.nuc_label).to_numpy(), None), index=m.index)
    hn = ok & nt.notna().to_numpy() & nt.isin(types).to_numpy()
    r["M5"] = float((m["first"].to_numpy()[hn] == nt.to_numpy()[hn]).mean()); r["n_nuc_blocks"] = int(hn.sum())
    key = (m.r8.astype(str) + " " + m.c8.astype(str)).to_numpy(); ft = pd.Series(np.where(ok, m["first"].to_numpy(), None), index=key)
    nb = np.concatenate([ft.reindex((m.r8.astype(str) + " " + (m.c8 + 1).astype(str)).to_numpy()).to_numpy(),
                         ft.reindex(((m.r8 + 1).astype(str) + " " + m.c8.astype(str)).to_numpy()).to_numpy()])
    sf = np.concatenate([ft.to_numpy(), ft.to_numpy()]); pp = pd.notna(nb) & pd.notna(sf)
    r["M6"] = float((nb[pp] == sf[pp]).mean()); r["n_pairs"] = int(pp.sum())
    if with_aniso: r["M7"] = aniso_umi(m)
    return r
