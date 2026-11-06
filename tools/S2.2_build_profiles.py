"""S2.2: costruisce l'insieme dei profili per simulate_cell_profiles() dai vincitori della gara (results/S2.2/S2.2_race_winners.csv),
stimati su tutte e 5 le aree (R6) o su tutto l'atlante. Per ogni tessuto scrive in /mnt/micron/geo_spatialtrans/S2.2/profiles/:
<A>_sums.tsv.gz (sottotipo, gene, somma dei conteggi: interi, i profili sono somme / totale), <A>_subtypes.tsv (sottotipo, tipo, cellule, quota nel tipo),
<A>_genes.tsv (geni di Visium HD del tessuto, nell'ordine della matrice), <A>_phi.tsv (phi per gene), <A>_totals.tsv (quantili 0, 0.001, ..., 1 dei
totali per tipo dei nuclei veri classificabili di R6, impronta intera), <A>_meta.json (fonte, modello, fattore di profondità).
Uso: .venv/bin/python tools/S2.2_build_profiles.py   Pre-registrazione: docs/steps/S2.2.md."""
import os, sys, json, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "tools"); import S22_common as S
M = "/mnt/micron/geo_spatialtrans"; OUT = f"{M}/S2.2/profiles"; os.makedirs(OUT, exist_ok=True)
KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
ATLAS = {"A2": "crc", "A3": "crc", "A4": "ln", "A5": "brain"}
MIN_SUB = 30
win = pd.read_csv("results/S2.2/S2.2_race_winners.csv").set_index("archetype").winner
fac = pd.read_csv("results/S2.2/S2.2_depth_factor.csv").set_index("archetype")
rois = pd.read_csv("results/R2/R2_rois_checked.csv")


def merge_small(ty, sub):
    s = pd.Series(sub); n = s.map(s.value_counts())
    return np.where(n.to_numpy() < MIN_SUB, np.char.add(np.asarray(ty, dtype=str), "|rest"), sub)


def r6_all(A, fp):
    """Nuclei classificabili con tipo del preset delle 5 aree; un nucleo che cade nel riquadro di un'area precedente (le aree di R2 si
    sovrappongono, revisione G-5) si conta una volta sola, nella prima area."""
    Xs, tys, subs = [], [], []
    ra = rois[rois.archetype == A].reset_index(drop=True)
    for i, roi in enumerate(ra.roi_id):
        tab = pd.read_parquet(f"{M}/R6/counts/{A}_{roi}_nuclei.parquet"); C = sp.load_npz(f"{M}/R6/counts/{A}_{roi}_{fp}.npz").tocsr()
        lab = pd.read_parquet(f"{M}/R6/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet"); lab = lab[lab.roi_id == roi].set_index("label")
        ty = lab.cell_type.reindex(tab.label).replace(S.MERGE).to_numpy(); cl = lab.cluster.reindex(tab.label).to_numpy()
        ok = np.isin(ty, S.PRESET_TYPES[A]) & pd.notna(cl)
        for j in range(i):
            ok &= ~S.inside_box(ra.iloc[i], tab.x_um.to_numpy(), tab.y_um.to_numpy(), ra.iloc[j])
        Xs.append(C[ok]); tys.append(ty[ok]); subs.append(np.array([f"{t}|{c}" for t, c in zip(ty[ok], cl[ok])]))
    ty = np.concatenate(tys)
    return sp.vstack(Xs).tocsr(), ty, merge_small(ty, np.concatenate(subs))


for A in S.PRESET_TYPES:
    w = win[A]; src, noise = w.split("_")[0], w.split("_")[1]
    vis = pd.read_csv(f"{M}/R6/counts/{KEY[A]}_genes.tsv", sep="\t")
    if src in ("exp", "nuc"):
        X, ty, sub = r6_all(A, src); Xv = X                                   # colonne = geni di Visium
    else:
        ac = pd.read_parquet(f"{M}/S2.2/atlas/{ATLAS[A]}_cells.parquet"); AX = sp.load_npz(f"{M}/S2.2/atlas/{ATLAS[A]}_X.npz").tocsc()
        ag = pd.read_csv(f"{M}/S2.2/atlas/{ATLAS[A]}_genes.tsv", sep="\t"); pos = pd.Series(np.arange(len(ag)), index=ag.gene_id)
        ok = ac.cell_type.isin(S.PRESET_TYPES[A]).to_numpy(); ac = ac[ok].reset_index(drop=True); AX = AX[ok]
        idx = pos.reindex(vis.gene_id).to_numpy(); have = ~np.isnan(idx)
        Xv = AX[:, idx[have].astype(int)].tocoo()                            # solo i geni comuni ...
        Xv = sp.csr_matrix((Xv.data, (Xv.row, np.flatnonzero(have)[Xv.col])), shape=(AX.shape[0], len(vis)))   # ... rimessi nell'ordine di Visium, zeri altrove
        ty = ac.cell_type.to_numpy(); sub = merge_small(ty, (ac.cell_type + "|" + ac.subtype_raw.astype(str)).to_numpy())
    Xv = sp.csr_matrix(Xv, dtype=np.float64)
    us, inv = np.unique(sub, return_inverse=True)
    I = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(us), len(inv)))
    sums = (I @ Xv).tocoo()
    pd.DataFrame({"subtype": us[sums.row], "gene_index": sums.col + 1, "count": sums.data.astype(np.int64)}).to_csv(f"{OUT}/{A}_sums.tsv.gz", sep="\t", index=False)
    cnt = np.bincount(inv); sty = np.array([u.split("|")[0] for u in us])
    share = cnt / pd.Series(cnt).groupby(sty).transform("sum").to_numpy()
    pd.DataFrame({"subtype": us, "cell_type": sty, "n_cells": cnt, "share": share}).to_csv(f"{OUT}/{A}_subtypes.tsv", sep="\t", index=False)
    vis.to_csv(f"{OUT}/{A}_genes.tsv", sep="\t", index=False)
    # phi per gene (momenti dentro i sottotipi, come nella gara)
    P = np.asarray((I @ Xv).todense()); tot = P.sum(axis=1); P = P / np.where(tot > 0, tot, 1)[:, None]
    Nc = np.asarray(Xv.sum(axis=1)).ravel(); sx2 = np.asarray(Xv.multiply(Xv).sum(axis=0)).ravel()
    sxN = np.asarray((I @ sp.diags(Nc) @ Xv).todense()); sN2 = np.asarray(I @ (Nc ** 2)).ravel(); sN = np.asarray(I @ Nc).ravel()
    num = sx2 - 2 * (P * sxN).sum(axis=0) + (P ** 2 * sN2[:, None]).sum(axis=0) - (P * sN[:, None]).sum(axis=0); den = (P ** 2 * sN2[:, None]).sum(axis=0)
    phi = np.where(den > 0, np.maximum(num, 0) / np.where(den > 0, den, 1), 0.0)
    pd.DataFrame({"gene_index": np.arange(1, len(vis) + 1), "phi": phi}).to_csv(f"{OUT}/{A}_phi.tsv", sep="\t", index=False)
    # totali per tipo: nuclei veri classificabili di R6 (impronta intera), tutte le aree
    Xe, ty_e, sub_e = r6_all(A, "exp"); tot_e = np.asarray(Xe.sum(axis=1)).ravel()     # totali dei nuclei veri, impronta intera, senza doppi
    qs = np.linspace(0, 1, 1001)
    tq = pd.DataFrame({t: np.quantile(tot_e[ty_e == t], qs) for t in np.unique(ty_e)}); tq.insert(0, "q", qs)
    tq.to_csv(f"{OUT}/{A}_totals.tsv", sep="\t", index=False)
    if src in ("exp", "nuc"):                                    # totali per sottotipo (PH-3): stessi nuclei dei sottotipi, impronta intera
        tqs = pd.DataFrame({s_: np.quantile(tot_e[sub_e == s_], qs) for s_ in us}); tqs.insert(0, "q", qs)
        tqs.to_csv(f"{OUT}/{A}_totals_subtype.tsv", sep="\t", index=False)
    elif os.path.exists(f"{OUT}/{A}_totals_subtype.tsv"):
        os.remove(f"{OUT}/{A}_totals_subtype.tsv")
    meta = dict(archetype=A, winner=w, source=src, noise=noise, depth_factor=float(fac.loc[A, "f"]), depth_factor_lo=float(fac.loc[A, "lo"]),
                depth_factor_hi=float(fac.loc[A, "hi"]), n_genes=int(len(vis)), n_genes_nonzero=int((Xv.sum(axis=0) > 0).sum()), n_subtypes=int(len(us)),
                n_cells_source=int(Xv.shape[0]), atlas=ATLAS.get(A) if src == "atlas" else None,
                totals_by="subtype" if src in ("exp", "nuc") else "type")       # DEC-046 (Luca, 2026-11-06)
    json.dump(meta, open(f"{OUT}/{A}_meta.json", "w"), indent=1)
    print(meta, flush=True)
