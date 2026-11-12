"""tools/F3.1_reference.py - F3.1: riferimento a cellula singola per RCTD e marcatori per tipo, uguali per i due lati.
Tipi: quelli del preset del tessuto (S22_common.MERGE e PRESET_TYPES, come i profili di S2.2).
Cellule: nuclei veri delle 5 aree del tessuto con la corona di 8 um (R6, footprint exp), tipo dalle etichette adottate di R6
(exp, umi20, seme 20261022, ph1; le stesse da cui S2.2 ha costruito i profili del simulatore); classificabili, tipo diverso da
'non risolto', UMI >= 20; al piu' N_MAX cellule per tipo (estrazione con seme); tipi con meno di 25 cellule esclusi (minimo di RCTD).
Marcatori: per tipo, i TOP geni per log2 del rapporto fra la media normalizzata (CP10k) nel tipo e la media negli altri tipi,
con media nel tipo >= 1 CP10k. Uso: .venv/bin/python tools/F3.1_reference.py <out_dir> <seed> A1 A2 ...
Scrive ref_<A>_counts.mtx (geni x cellule), ref_<A>_cells.csv, ref_<A>_genes.tsv, markers_<A>.csv, ref_<A>.json."""
import sys, os, json, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
sys.path.insert(0, "tools"); import S22_common as S   # MERGE e PRESET_TYPES: il vocabolario dei preset, come S2.2 (DEC-014)
M = "/mnt/micron/geo_spatialtrans"; C = f"{M}/R6/counts"
KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
N_MAX = 3000; MIN_CELLS = 25; TOP = 50; UMI_MIN = 20
out = sys.argv[1]; seed = int(sys.argv[2]); archs = sys.argv[3:]; os.makedirs(out, exist_ok=True)
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
for A in archs:
    rng = np.random.default_rng([seed, int(A[1])])
    genes = pd.read_csv(f"{C}/{KEY[A]}_genes.tsv", sep="\t")
    lab = pd.read_parquet(f"{M}/R6/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet")
    mats, meta = [], []
    for roi in rois.loc[rois.archetype == A, "roi_id"]:
        nuc = pd.read_parquet(f"{C}/{A}_{roi}_nuclei.parquet")
        E = sp.load_npz(f"{C}/{A}_{roi}_exp.npz").tocsr()
        assert E.shape[0] == len(nuc)
        L = lab[lab.roi_id == roi].set_index("label")
        t = nuc[["label"]].copy(); t["cell_type"] = L.cell_type.reindex(t.label).replace(S.MERGE).to_numpy(); t["classifiable"] = L.classifiable.reindex(t.label).to_numpy()
        t["umi"] = np.asarray(E.sum(axis=1)).ravel(); t["roi_id"] = roi
        ok = (t.classifiable == True) & t.cell_type.isin(S.PRESET_TYPES[A]) & (t.umi >= UMI_MIN)
        mats.append(E[np.flatnonzero(ok.to_numpy())]); meta.append(t[ok.to_numpy()])
    X = sp.vstack(mats).tocsr(); meta = pd.concat(meta, ignore_index=True)
    vc = meta.cell_type.value_counts(); types = sorted(vc[vc >= MIN_CELLS].index)
    pick = []
    for ty in types:
        ix = np.flatnonzero(meta.cell_type.to_numpy() == ty)
        pick.append(np.sort(rng.choice(ix, size=min(N_MAX, len(ix)), replace=False)))
    pick = np.concatenate(pick); Xr = X[pick]; mr = meta.iloc[pick].reset_index(drop=True)
    mr["cell_id"] = [f"{A}_{r}_{l}" for r, l in zip(mr.roi_id, mr.label)]
    sio.mmwrite(f"{out}/ref_{A}_counts.mtx", Xr.T.tocsc().astype(np.int64))
    mr[["cell_id", "cell_type", "umi", "roi_id", "label"]].to_csv(f"{out}/ref_{A}_cells.csv", index=False)
    genes.gene_id.to_csv(f"{out}/ref_{A}_genes.tsv", index=False, header=False)   # identificativi Ensembl (unici)
    # marcatori (su tutte le cellule tenute, non solo il campione)
    keepm = meta.cell_type.isin(types).to_numpy(); Xm = X[keepm]; tm = meta.cell_type.to_numpy()[keepm]
    lib = np.asarray(Xm.sum(axis=1)).ravel(); N = sp.diags(1e4 / lib) @ Xm
    means = np.vstack([np.asarray(N[tm == ty].mean(axis=0)).ravel() for ty in types])
    rows = []
    for k, ty in enumerate(types):
        other = np.delete(means, k, axis=0).mean(axis=0)
        lfc = np.log2((means[k] + 0.1) / (other + 0.1)); lfc[means[k] < 1] = -np.inf
        top = np.argsort(-lfc)[:TOP]
        rows += [{"cell_type": ty, "gene": genes.gene_id.iloc[g], "gene_name": genes.gene_name.iloc[g], "rank": r + 1, "lfc": float(lfc[g]), "mean_type": float(means[k, g])}
                 for r, g in enumerate(top) if np.isfinite(lfc[g])]
    pd.DataFrame(rows).to_csv(f"{out}/markers_{A}.csv", index=False)
    info = {"archetype": A, "n_cells_pool": int(len(meta)), "types_pool": vc.to_dict(), "types_ref": types,
            "n_ref": int(len(mr)), "n_genes": int(X.shape[1]), "seed": seed}
    json.dump(info, open(f"{out}/ref_{A}.json", "w")); print(json.dumps(info), flush=True)
