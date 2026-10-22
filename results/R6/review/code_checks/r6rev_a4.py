
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, os, time
from scipy.sparse.csgraph import connected_components
M="/mnt/micron/geo_spatialtrans/R6"; os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
A="A4"; SEED=20261022; t0=time.time()
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype==A]
genes = pd.read_csv(f"{M}/counts/A4_genes.tsv", sep="\t")
lab = pd.read_parquet(f"{M}/labels/A4_exp_umi20_seed20261022_nuclei.parquet")
X = sp.vstack([sp.load_npz(f"{M}/counts/A4_{r}_exp.npz") for r in rois.roi_id]).tocsr()
adata = ad.AnnData(X=X.astype(np.float32), obs=lab.reset_index(drop=True), var=pd.DataFrame(index=genes.gene_name.values)); adata.var_names_make_unique()
sub = adata[adata.obs.classifiable.values].copy()
sc.pp.normalize_total(sub); sc.pp.log1p(sub); sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
sc.pp.pca(sub, n_comps=30, mask_var="highly_variable", svd_solver="arpack", random_state=SEED)
sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30, random_state=SEED)
ncc, cc = connected_components(sub.obsp["connectivities"], directed=False)
sizes = np.bincount(cc)
print(f"A4: classifiable {sub.n_obs}; kNN graph connected components: {ncc}; component sizes: largest {sizes.max()}, n components <50 nuclei {int((sizes<50).sum())} holding {int(sizes[sizes<50].sum())} nuclei; {time.time()-t0:.0f}s")
cl = sub.obs.cluster.astype(str).values
tab = pd.crosstab(cc, cl)
# how many leiden clusters are exactly one small component
one = sum((tab.loc[i] > 0).sum() == 1 for i in range(ncc) if sizes[i] < 50)
print(f"small components (<50) that map to exactly one leiden cluster: {one}; leiden clusters with <50 nuclei: {int((pd.Series(cl).value_counts()<50).sum())}")
# umi of nuclei in small components vs giant
umi = sub.obs.umi.values; small = sizes[cc] < 50
print(f"median umi: in small components {np.median(umi[small]):.0f} vs giant component {np.median(umi[~small]):.0f}; sr? no. n small {small.sum()}")
print(f"ROI of small-component nuclei: {pd.Series(sub.obs.roi_id.values[small]).value_counts().to_dict()}")
# does leiden with the giant component only give fewer clusters? just report distribution
print("leiden cluster sizes quantiles:", pd.Series(cl).value_counts().quantile([0,.25,.5,.75,1]).astype(int).tolist())
