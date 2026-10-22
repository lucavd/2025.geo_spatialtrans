"""tools/S1.5b_a5_missed.py - S1.5b, decisione DEC-017: tipo dei nuclei del cervello contati a mano e persi da Cellpose.
Per ogni punto manuale delle finestre di R2b nel cervello si sommano i trascritti dei bin di 2 um con centro entro 8 um dal punto;
il tipo e' quello della maggioranza dei 15 nuclei di R6 (etichette del criterio corretto, DEC-013) piu' vicini nello spazio delle
30 componenti principali di R6 (normalizzazione per profondita' e log1p, 2 000 geni variabili, seme 20261022, come tools/R6_cluster.py).
Controllo: sui punti trovati da Cellpose il tipo trasferito deve coincidere con quello di R6 in almeno 3 casi su 4.
Uso: .venv/bin/python tools/S1.5b_a5_missed.py   Uscite: results/S1.5b/S1.5b_a5_missed_points.csv, S1.5b_a5_composition.csv"""
import sys, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, pyarrow.parquet as pq
from scipy.spatial import cKDTree
sys.path.insert(0, "tools"); from R6_common import read_h5_csc
M = "/mnt/micron/geo_spatialtrans"; R6 = f"{M}/R6"; DS = "Visium_HD_6p5mm_Mouse_Brain"; RES = "results/S1.5b"
SEED, K, RADIUS, UMI_MIN = 20261022, 15, 8.0, 20
UNTYPED = {"non risolto", "non classificabile"}
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype == "A5"].set_index("roi_id")
# ---- riferimento: nuclei di R6 (stessa preparazione di R6_cluster.py) ----
genes = pd.read_csv(f"{R6}/counts/A5_genes.tsv", sep="\t")
mats, obs = [], []
for roi in rois.index:
    obs.append(pd.read_parquet(f"{R6}/counts/A5_{roi}_nuclei.parquet")); mats.append(sp.load_npz(f"{R6}/counts/A5_{roi}_exp.npz"))
obs = pd.concat(obs, ignore_index=True); X = sp.vstack(mats).tocsr()
lab = pd.read_parquet(f"{R6}/labels/A5_exp_umi20_seed20261022_nuclei_ph1.parquet")
obs = obs.merge(lab[["roi_id", "label", "cell_type"]], on=["roi_id", "label"], how="left", validate="one_to_one")
obs["umi"] = np.asarray(X.sum(axis=1)).ravel()
ref = ad.AnnData(X=X.astype(np.float32), obs=obs.reset_index(drop=True), var=pd.DataFrame(index=genes.gene_name.values)); ref.var_names_make_unique()
ref = ref[(ref.obs.umi >= UMI_MIN).values].copy()
target = float(np.median(ref.obs.umi))
sc.pp.normalize_total(ref, target_sum=target); sc.pp.log1p(ref)
sc.pp.highly_variable_genes(ref, n_top_genes=2000, flavor="seurat")
sc.pp.pca(ref, n_comps=30, mask_var="highly_variable", svd_solver="arpack", random_state=SEED)
hv = ref.var.highly_variable.values; PCs = ref.varm["PCs"][hv]; mu = np.asarray(ref.X[:, hv].mean(axis=0)).ravel()
typed = ~ref.obs.cell_type.isin(UNTYPED).values
tree = cKDTree(ref.obsm["X_pca"][typed]); ref_types = ref.obs.cell_type.values[typed]
# ---- punti manuali e bin ----
X2, bc, _, _ = read_h5_csc(f"{M}/data_real/{DS}/binned_outputs_x/binned_outputs/square_002um/filtered_feature_bc_matrix.h5")
tp = pq.read_table(f"{M}/data_real/{DS}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet").to_pandas()
col_of = pd.Series(np.arange(len(bc)), index=bc)
P = pd.read_csv("results/R2b/R2b_points.csv"); W = pd.read_csv("results/R2b/R2b_windows.csv")
P = P.merge(W[["win_id", "archetype", "roi_id"]], on="win_id"); P = P[(P.archetype == "A5") & P.rater.isin(["Luca", "Federica"])].copy()
P["matched"] = P.match_cellpose_rgb.astype(bool)
P = P.merge(lab[["roi_id", "label", "cell_type"]].rename(columns={"label": "label_cellpose_rgb", "cell_type": "type_r6"}),
            on=["roi_id", "label_cellpose_rgb"], how="left")
out = []
for roi, Q in P.groupby("roi_id"):
    R = rois.loc[roi]; um = R.um_per_px
    box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)]
    box = box[box.barcode.isin(col_of.index)]
    bx = (box.pxl_col_in_fullres.to_numpy() - R.c0) * um; by = (box.pxl_row_in_fullres.to_numpy() - R.r0) * um
    bt = cKDTree(np.c_[bx, by]); cols = col_of.reindex(box.barcode).to_numpy()
    for _, p in Q.iterrows():
        idx = bt.query_ball_point([p.x_um, p.y_um], RADIUS)
        v = np.asarray(X2[:, cols[idx]].sum(axis=1)).ravel() if idx else np.zeros(X2.shape[0])
        umi = v.sum(); rec = dict(win_id=p.win_id, rater=p.rater, idx=p.idx, roi_id=roi, matched=p.matched, type_r6=p.type_r6, n_bins=len(idx), umi=umi)
        if umi >= UMI_MIN:
            g = np.log1p(v / umi * target)
            vv = pd.Series(g, index=genes.gene_name.values); vv = vv[~vv.index.duplicated()]
            gx = vv.reindex(ref.var_names[hv]).fillna(0).to_numpy()
            z = (gx - mu) @ PCs
            d, j = tree.query(z, k=K); vc = pd.Series(ref_types[j]).value_counts()
            top = vc[vc == vc.max()].index
            # pareggio: il tipo con il vicino piu' vicino
            rec["type_transfer"] = min(top, key=lambda t: d[np.where(ref_types[j] == t)[0][0]])
            rec["vote_share"] = vc.max() / K
        else:
            rec["type_transfer"] = "non classificabile"; rec["vote_share"] = np.nan
        out.append(rec)
D = pd.DataFrame(out); D.to_csv(f"{RES}/S1.5b_a5_missed_points.csv", index=False)
comp_r6 = pd.read_csv(f"{RES}/S1.5b_type_calibration.csv").query("archetype == 'A5'").set_index("cell_type").fraction_r6
rows = []
for rater, G in D.groupby("rater"):
    m = G[G.matched & G.type_r6.notna() & ~G.type_r6.isin(UNTYPED) & (G.type_transfer != "non classificabile")]
    acc = float((m.type_r6 == m.type_transfer).mean()); n_all = len(G); n_miss = int((~G.matched).sum())
    miss = G[~G.matched & (G.type_transfer != "non classificabile")]
    gcomp = miss.type_transfer.value_counts(normalize=True).reindex(comp_r6.index).fillna(0)
    q = n_miss / n_all; corr = (1 - q) * comp_r6 + q * gcomp
    for t in comp_r6.index:
        rows.append(dict(rater=rater, cell_type=t, n_points=n_all, n_missed=n_miss, n_missed_typed=len(miss), n_matched_eval=len(m), accuracy_matched=acc,
                         q=q, fraction_r6=comp_r6[t], fraction_missed=gcomp[t], fraction_corrected=corr[t], passes=acc >= 0.75))
C = pd.DataFrame(rows); C.to_csv(f"{RES}/S1.5b_a5_composition.csv", index=False)
pd.set_option("display.width", 220); print(C.round(4).to_string(index=False))
print(pd.crosstab(D[D.matched].type_r6, D[D.matched].type_transfer))
