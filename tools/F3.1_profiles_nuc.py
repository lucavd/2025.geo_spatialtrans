"""tools/F3.1_profiles_nuc.py - F3.1, controprova CP-3: insieme di profili per simulate_cell_profiles() uguale a quello di S2.2
(/mnt/micron/geo_spatialtrans/S2.2/profiles) tranne le somme per sottotipo (<A>_sums.tsv.gz), ricalcolate sui conteggi del SOLO NUCLEO
(R6, impronta 'nuc') degli stessi nuclei, con gli stessi sottotipi (stessa selezione e stessa fusione dei sottotipi piccoli di
tools/S2.2_build_profiles.py, r6_all e merge_small). Trascritti per cellula, dispersioni, quote dei sottotipi e fattore di profondita'
restano quelli di S2.2: cambia solo la composizione in geni. Solo i tessuti con profili dai nuclei (fonte 'exp').
Uso: .venv/bin/python tools/F3.1_profiles_nuc.py <out_dir> A1 A2 A4 A5 A6"""
import os, sys, json, shutil, glob, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "tools"); import S22_common as S
M = "/mnt/micron/geo_spatialtrans"; SRC = f"{M}/S2.2/profiles"; MIN_SUB = 30
out = sys.argv[1]; archs = sys.argv[2:]; os.makedirs(out, exist_ok=True)
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
def merge_small(ty, sub):   # copia di tools/S2.2_build_profiles.py
    s = pd.Series(sub); n = s.map(s.value_counts())
    return np.where(n.to_numpy() < MIN_SUB, np.char.add(np.asarray(ty, dtype=str), "|rest"), sub)
def r6_all(A, fp):          # copia di tools/S2.2_build_profiles.py
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
for A in archs:
    meta = json.load(open(f"{SRC}/{A}_meta.json")); assert meta["source"] == "exp", A
    for f in glob.glob(f"{SRC}/{A}_*"): shutil.copy(f, out)
    res = {}
    for fp in ["exp", "nuc"]:
        X, ty, sub = r6_all(A, fp)
        us, inv = np.unique(sub, return_inverse=True)
        I = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(us), len(inv)))
        res[fp] = ((I @ sp.csr_matrix(X, dtype=np.float64)).tocoo(), us, np.bincount(inv))
    (Se, ue, ne), (Sn, un, nn) = res["exp"], res["nuc"]
    assert np.array_equal(ue, un) and np.array_equal(ne, nn), f"{A}: sottotipi diversi fra exp e nuc"
    # controllo: le somme exp ricostruite = file di S2.2
    old = pd.read_csv(f"{SRC}/{A}_sums.tsv.gz", sep="\t")
    rec = pd.DataFrame({"subtype": ue[Se.row], "gene_index": Se.col + 1, "count": Se.data.astype(np.int64)})
    m = old.merge(rec, on=["subtype", "gene_index"], how="outer", suffixes=("_old", "_rec")).fillna(0)
    same = bool((m.count_old == m.count_rec).all())
    pd.DataFrame({"subtype": un[Sn.row], "gene_index": Sn.col + 1, "count": Sn.data.astype(np.int64)}).to_csv(f"{out}/{A}_sums.tsv.gz", sep="\t", index=False)
    meta["F3.1"] = {"sums": "nucleus only (R6 'nuc')", "exp_sums_reproduced": same}
    json.dump(meta, open(f"{out}/{A}_meta.json", "w"), indent=1)
    print(A, "subtypes", len(un), "exp sums reproduced:", same, "total nuc / exp:", round(Sn.data.sum() / Se.data.sum(), 3), flush=True)
