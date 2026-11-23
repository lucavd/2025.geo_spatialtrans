"""tools/sonda/S3_reference.py - S3: riferimento RCTD e marcatori dal HLCA core (docs/sonda/REGOLA_S3.md, «Riferimento RCTD»).
HLCA core (CELLxGENE 066943a2-fdac-4b29-b348-40cede398e4e), conteggi grezzi raw/X; etichetta ann_level_2; cellule con tissue != 'nose';
geni = identificativi Ensembl comuni con le feature della matrice filtrata di H2 (uguali a quelle di H1, verificato); UMI >= 20 sui geni comuni; tipi con
almeno 25 cellule; al piu' 3000 cellule per tipo estratte con numpy default_rng(20261123), tipi in ordine alfabetico, come
tools/F3.1_reference.py. Marcatori come F3.1_reference.py su tutte le cellule tenute: CP10k medio per tipo, lfc = log2((m_k + 0.1) /
(media degli altri tipi + 0.1)), m_k >= 1, primi 50. Sigla del riferimento: L. Due passate a blocchi di righe su raw/X (CSR).
Uso: .venv/bin/python -B tools/sonda/S3_reference.py  ->  /mnt/micron/geo_spatialtrans/S3_sonda/ref/"""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio, h5py, anndata as ad
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S3_common as C
F = f"{C.M}/ref_lung/HLCA_core_066943a2.h5ad"; A = "L"; N_MAX, MIN_CELLS, TOP, UMI_MIN, CH = 3000, 25, 50, 20, 20000
os.makedirs(C.REF, exist_ok=True); t0 = time.time()
ids2, nm2 = C.features(f"{C.bdir('H2', '002')}/filtered_feature_bc_matrix.h5"); ids1, _ = C.features(f"{C.bdir('H1', '002')}/filtered_feature_bc_matrix.h5")
assert np.array_equal(ids1, ids2), "H1 and H2 feature lists differ"
obs = ad.read_h5ad(F, backed="r").obs[["ann_level_2", "tissue", "dataset", "assay"]].copy()
with h5py.File(F, "r") as h:
    rv = h["raw/var"]; rid = rv[rv.attrs["_index"]][:].astype(str)
    ip = h["raw/X/indptr"][:]; ng = int(h["raw/X"].attrs["shape"][1])
gpos = pd.Series(np.arange(len(ids2)), index=ids2); gi = gpos.reindex(rid).to_numpy()      # colonna HLCA -> posizione del gene comune
common_ids = [g for g in ids2 if g in set(rid)]; cpos = pd.Series(np.arange(len(common_ids)), index=common_ids)
cmap = np.full(ng, -1, np.int64); okc = ~np.isnan(gi); cmap[np.flatnonzero(okc)] = cpos.reindex(rid[okc]).to_numpy()
G = len(common_ids); n = len(obs); ty_all = obs.ann_level_2.astype(str).to_numpy(); nose = (obs.tissue.astype(str) == "nose").to_numpy()
def rows(h, r0, r1):
    d = h["raw/X/data"][ip[r0]:ip[r1]]; ix = h["raw/X/indices"][ip[r0]:ip[r1]]; c = cmap[ix]; k = c >= 0
    rr = np.repeat(np.arange(r1 - r0), np.diff(ip[r0:r1 + 1]))
    return sp.csr_matrix((d[k].astype(np.float64), (rr[k], c[k])), shape=(r1 - r0, G))
# passata 1: UMI sui geni comuni
umi = np.zeros(n)
with h5py.File(F, "r") as h:
    for r0 in range(0, n, CH):
        r1 = min(n, r0 + CH); umi[r0:r1] = np.asarray(rows(h, r0, r1).sum(axis=1)).ravel()
keep = (~nose) & (umi >= UMI_MIN)
vc = pd.Series(ty_all[keep]).value_counts(); types = sorted(vc[vc >= MIN_CELLS].index); keep &= np.isin(ty_all, types)
# passata 2: medie CP10k per tipo (marcatori, tutte le cellule tenute)
S = np.zeros((len(types), G)); cnt = np.zeros(len(types)); tix = pd.Series(np.arange(len(types)), index=types)
with h5py.File(F, "r") as h:
    for r0 in range(0, n, CH):
        r1 = min(n, r0 + CH); kk = keep[r0:r1]
        if not kk.any(): continue
        X = rows(h, r0, r1)[kk]; N = sp.diags(1e4 / umi[r0:r1][kk]) @ X; t = tix.reindex(ty_all[r0:r1][kk]).to_numpy()
        Ind = sp.csr_matrix((np.ones(len(t)), (t, np.arange(len(t)))), shape=(len(types), len(t)))
        S += np.asarray((Ind @ N).todense()); cnt += np.bincount(t, minlength=len(types))
means = S / cnt[:, None]
# campione per RCTD
rng = np.random.default_rng(C.SEED_REF); pick = []
for ty in types:
    ix = np.flatnonzero(keep & (ty_all == ty)); pick.append(np.sort(rng.choice(ix, size=min(N_MAX, len(ix)), replace=False)))
pick = np.sort(np.concatenate(pick))
with h5py.File(F, "r") as h:
    Xr = sp.vstack([rows(h, int(i), int(i) + 1) for i in pick]).tocsr()
assert np.allclose(np.asarray(Xr.sum(axis=1)).ravel(), umi[pick])
cells = pd.DataFrame({"cell_id": [f"HLCA_{i}" for i in pick], "cell_type": ty_all[pick], "umi": umi[pick], "roi_id": "HLCA",
                      "label": pick, "obs_name": obs.index.to_numpy()[pick], "dataset": obs.dataset.astype(str).to_numpy()[pick]})
sio.mmwrite(f"{C.REF}/ref_{A}_counts.mtx", Xr.T.tocsc().astype(np.int64))
cells.to_csv(f"{C.REF}/ref_{A}_cells.csv", index=False)
pd.Series(common_ids).to_csv(f"{C.REF}/ref_{A}_genes.tsv", index=False, header=False)
name_of = pd.Series(nm2, index=ids2); rows_m = []
for k, ty in enumerate(types):
    other = np.delete(means, k, axis=0).mean(axis=0)
    lfc = np.log2((means[k] + 0.1) / (other + 0.1)); lfc[means[k] < 1] = -np.inf
    for r, g in enumerate(np.argsort(-lfc)[:TOP]):
        if np.isfinite(lfc[g]): rows_m.append({"cell_type": ty, "gene": common_ids[g], "gene_name": name_of[common_ids[g]], "rank": r + 1, "lfc": float(lfc[g]), "mean_type": float(means[k, g])})
pd.DataFrame(rows_m).to_csv(f"{C.REF}/markers_{A}.csv", index=False)
info = {"source": F, "md5": C.md5(F), "label": "ann_level_2", "n_cells_total": int(n), "n_nose_excluded": int(nose.sum()), "n_pool": int(keep.sum()),
        "types_pool": {t: int(v) for t, v in vc.items()}, "types_ref": types, "n_ref": int(len(pick)), "n_ref_by_type": cells.cell_type.value_counts().to_dict(),
        "n_genes_common": G, "n_genes_visium": int(len(ids2)), "n_genes_hlca": int(ng), "seed": C.SEED_REF, "sec": round(time.time() - t0, 1)}
json.dump(info, open(f"{C.REF}/ref_{A}.json", "w"), indent=1); print(json.dumps(info), flush=True)
