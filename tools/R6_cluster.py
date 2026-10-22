"""R6: cluster nuclei per archetype from aggregated 2 um counts and label clusters with the coarse type catalogue
(PanglaoDB markers via OmniPath). Usage: python tools/R6_cluster.py <A> [fp=exp] [umi_min=20] [seed=20261022]
Outputs (/mnt/micron/geo_spatialtrans/R6/labels/): <A>_<fp>_umi<min>_seed<seed>_nuclei.parquet, _clusters.csv, _scores.csv"""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
A = sys.argv[1]; FP = sys.argv[2] if len(sys.argv) > 2 else "exp"; UMI_MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 20
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261022
M = "/mnt/micron/geo_spatialtrans/R6"; OUT = f"{M}/labels"; os.makedirs(OUT, exist_ok=True)
KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}[A]
SPECIES = {"A1": "mouse", "A2": "human", "A3": "human", "A4": "human", "A5": "mouse", "A6": "human"}[A]
# coarse catalogue (approved by Luca 2026-10-23) -> PanglaoDB cell types
CAT = {
 "A1": {"epiteliale": ["Enterocytes", "Goblet cells", "Paneth cells", "Enteroendocrine cells", "Tuft cells", "Crypt cells", "Epithelial cells"],
        "stromale e muscolo liscio": ["Fibroblasts", "Smooth muscle cells", "Stromal cells", "Myofibroblasts", "Pericytes"],
        "immunitario": ["B cells", "T cells", "Plasma cells", "Macrophages", "Dendritic cells", "Mast cells", "Monocytes", "NK cells", "Neutrophils", "T memory cells", "T helper cells"],
        "endoteliale": ["Endothelial cells"]},
 "A2": {"epiteliale": ["Enterocytes", "Goblet cells", "Crypt cells", "Epithelial cells", "Enteroendocrine cells", "Tuft cells"],
        "fibroblasti": ["Fibroblasts", "Myofibroblasts", "Stromal cells"],
        "immunitario": ["B cells", "T cells", "Plasma cells", "Macrophages", "Dendritic cells", "Mast cells", "Monocytes", "NK cells", "Neutrophils", "T memory cells", "T helper cells"],
        "endoteliale": ["Endothelial cells"],
        "muscolo liscio": ["Smooth muscle cells", "Pericytes"]},
 "A4": {"linfociti B": ["B cells", "B cells naive", "B cells memory"],
        "linfociti T": ["T cells", "T helper cells", "T memory cells", "T cytotoxic cells", "T regulatory cells", "T follicular helper cells"],
        "plasmacellule": ["Plasma cells"],
        "mieloidi": ["Macrophages", "Dendritic cells", "Monocytes", "Plasmacytoid dendritic cells", "Mast cells", "Neutrophils"],
        "stromali": ["Fibroblasts", "Stromal cells", "Smooth muscle cells", "Pericytes"],
        "endoteliali": ["Endothelial cells"]},
 "A5": {"neuroni": ["Neurons", "Interneurons", "Pyramidal cells", "GABAergic neurons", "Glutaminergic neurons"],
        "astrociti": ["Astrocytes"], "oligodendrociti": ["Oligodendrocytes", "Oligodendrocyte progenitor cells"], "microglia": ["Microglia"],
        "vascolari": ["Endothelial cells", "Endothelial cells (blood brain barrier)", "Pericytes", "Smooth muscle cells"]},
 "A6": {"cardiomiociti": ["Cardiomyocytes"], "fibroblasti": ["Fibroblasts"], "endoteliali": ["Endothelial cells", "Endothelial cells (aorta)"],
        "immunitari": ["Macrophages", "T cells", "B cells", "Mast cells", "Monocytes", "Dendritic cells", "Neutrophils"],
        "muscolo liscio e periciti": ["Smooth muscle cells", "Vascular smooth muscle cells", "Pericytes"]}}
CAT["A3"] = CAT["A2"]
# tumour vs normal epithelium (A2/A3 only): crypt/stem-like score vs differentiated (enterocyte + goblet) score; technical rule, declared
EPI_SPLIT = {"tumorale": ["Crypt cells"], "normale": ["Enterocytes", "Goblet cells"]}
t0 = time.time()
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype == A]
genes = pd.read_csv(f"{M}/counts/{KEY}_genes.tsv", sep="\t")
mats, obs = [], []
for roi in rois.roi_id:
    tab = pd.read_parquet(f"{M}/counts/{A}_{roi}_nuclei.parquet"); C = sp.load_npz(f"{M}/counts/{A}_{roi}_{FP}.npz")
    mats.append(C); obs.append(tab)
obs = pd.concat(obs, ignore_index=True); X = sp.vstack(mats).tocsr()
obs["umi"] = obs[f"umi_{FP}"]; obs["classifiable"] = obs.umi >= UMI_MIN
obs.index = obs.archetype + "_" + obs.roi_id + "_" + obs.label.astype(str)
adata = ad.AnnData(X=X.astype(np.float32), obs=obs, var=pd.DataFrame(index=genes.gene_name.values))
adata.var_names_make_unique()
sub = adata[adata.obs.classifiable].copy()
print(f"[{A} {FP} umi>={UMI_MIN}] nuclei {adata.n_obs}, classifiable {sub.n_obs}, genes {adata.n_vars}; umi median {np.median(sub.obs.umi):.0f}", flush=True)
sc.pp.normalize_total(sub); sc.pp.log1p(sub)
sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
sc.pp.pca(sub, n_comps=30, mask_var="highly_variable", svd_solver="arpack", random_state=SEED)
sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30, random_state=SEED)
sc.tl.leiden(sub, resolution=1.0, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
print(f"[{A}] leiden: {sub.obs.leiden.nunique()} clusters in {time.time()-t0:.0f}s", flush=True)
# marker scores: cluster-mean log-normalised expression, z-scored across clusters per gene, averaged over the type's genes
pan = pd.read_csv(f"{M}/PanglaoDB_omnipath_{SPECIES}.tsv", sep="\t")
pan = pan[pan[SPECIES] == True]
sym = pan.genesymbol.astype(str)
if SPECIES == "mouse": sym = sym.str.capitalize()
pan = pan.assign(gene=sym.values); present = set(sub.var_names)
cl = sub.obs.leiden.astype(str).values; clusters = sorted(set(cl), key=int)
Xl = sub.X.tocsr()
means = np.vstack([np.asarray(Xl[cl == k].mean(axis=0)).ravel() for k in clusters])   # clusters x genes
sd = means.std(axis=0); sd[sd == 0] = 1.0; Z = (means - means.mean(axis=0)) / sd
gi = pd.Series(np.arange(sub.n_vars), index=sub.var_names)
def type_score(ptypes):
    g = sorted(set(pan.gene[pan.cell_type.isin(ptypes)]) & present)
    return (Z[:, gi[g].values].mean(axis=1) if g else np.full(len(clusters), np.nan)), len(g)
scores = {}; ngen = {}
for ty, pts in CAT[A].items(): scores[ty], ngen[ty] = type_score(pts)
S = pd.DataFrame(scores, index=clusters)
best = S.idxmax(axis=1); bestv = S.max(axis=1); second = S.apply(lambda r: r.drop(r.idxmax()).max(), axis=1)
ctype = best.where(bestv >= 0, "non risolto")
if A in ("A2", "A3"):
    st, nt = type_score(EPI_SPLIT["tumorale"]); sn, nn = type_score(EPI_SPLIT["normale"])
    S["epi_tumorale_score"] = st; S["epi_normale_score"] = sn
    ctype = ctype.where(ctype != "epiteliale", pd.Series(np.where(st > sn, "epiteliale tumorale", "epiteliale normale"), index=clusters))
S["n_nuclei"] = pd.Series(cl).value_counts().reindex(clusters).values; S["type"] = ctype; S["best_score"] = bestv; S["margin"] = bestv - second
# top differential genes per cluster for human control
sc.tl.rank_genes_groups(sub, "leiden", method="t-test", n_genes=10)
top = {k: ", ".join(sub.uns["rank_genes_groups"]["names"][k][:10]) for k in clusters}
S["top10_genes"] = pd.Series(top)
S.index.name = "cluster"; S.to_csv(f"{OUT}/{A}_{FP}_umi{UMI_MIN}_seed{SEED}_clusters.csv")
pd.Series(ngen, name="n_marker_genes_present").to_csv(f"{OUT}/{A}_{FP}_umi{UMI_MIN}_seed{SEED}_markers.csv")
lab = adata.obs[["archetype", "roi_id", "label", "x_um", "y_um", "area_um2", "eccentricity", "umi", "classifiable"] + [c for c in adata.obs.columns if c.startswith("sr_")]].copy()
lab["cluster"] = pd.Series(cl, index=sub.obs_names).reindex(lab.index).values
lab["cell_type"] = pd.Series(ctype.reindex(lab.cluster).values, index=lab.index).fillna("non classificabile")
lab.loc[lab.cell_type.isna() | (lab.cluster.isna()), "cell_type"] = "non classificabile"
lab.to_parquet(f"{OUT}/{A}_{FP}_umi{UMI_MIN}_seed{SEED}_nuclei.parquet", index=False)
print(f"[{A}] types: {lab.cell_type.value_counts().to_dict()}; {time.time()-t0:.0f}s", flush=True)
print(S[["n_nuclei", "type", "best_score", "margin", "top10_genes"]].to_string())
