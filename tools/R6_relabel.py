"""R6 post hoc PH-1: relabel the pre-registered clusters by hypergeometric enrichment of each cluster's top-50 differential
genes (t-test vs all other clusters) in the PanglaoDB marker sets of the catalogue; type = set with the smallest p if
p < 0.01 / n_types (Bonferroni), else 'non risolto'. Usage: python tools/R6_relabel.py <A> [fp=exp] [umi_min=20] [seed=20261022]"""
import sys, os, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
from scipy.stats import hypergeom
sys.path.insert(0, "tools"); from R6_catalogue import CAT, EPI_SPLIT, SPECIES
A = sys.argv[1]; FP = sys.argv[2] if len(sys.argv) > 2 else "exp"; UMI_MIN = int(sys.argv[3]) if len(sys.argv) > 3 else 20
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261022; TOPN = 50; ALPHA = 0.01
M = "/mnt/micron/geo_spatialtrans/R6"; KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}[A]; SP = SPECIES[A]
base = f"{M}/labels/{A}_{FP}_umi{UMI_MIN}_seed{SEED}"
lab = pd.read_parquet(base + "_nuclei.parquet")
genes = pd.read_csv(f"{M}/counts/{KEY}_genes.tsv", sep="\t")
X = sp.vstack([sp.load_npz(f"{M}/counts/{A}_{roi}_{FP}.npz") for roi in lab.roi_id.unique()]).tocsr(); assert X.shape[0] == len(lab)
adata = ad.AnnData(X=X.astype(np.float32), obs=lab.reset_index(drop=True), var=pd.DataFrame(index=genes.gene_name.values)); adata.var_names_make_unique()
sub = adata[adata.obs.classifiable.values].copy(); sub.obs["cluster"] = sub.obs.cluster.astype(str)
sc.pp.normalize_total(sub); sc.pp.log1p(sub)
sc.pp.filter_genes(sub, min_cells=max(3, int(0.001 * sub.n_obs)))      # universe: genes detected in >= 0.1 % of nuclei
N = sub.n_vars; present = set(sub.var_names)
pan = pd.read_csv(f"{M}/PanglaoDB_omnipath_{SP}.tsv", sep="\t"); pan = pan[pan[SP] == True]
sym = pan.genesymbol.astype(str); sym = sym.str.capitalize() if SP == "mouse" else sym; pan = pan.assign(gene=sym.values)
def gset(ptypes): return set(pan.gene[pan.cell_type.isin(ptypes)]) & present
sets = {ty: gset(pts) for ty, pts in CAT[A].items()}; nty = len(sets)
sc.tl.rank_genes_groups(sub, "cluster", method="t-test", n_genes=TOPN)
clusters = list(sub.obs.cluster.unique())
rows = {}
for k in clusters:
    top = set(sub.uns["rank_genes_groups"]["names"][k][:TOPN])
    r = {}
    for ty, g in sets.items():
        o = len(top & g); r[f"p_{ty}"] = hypergeom.sf(o - 1, N, len(g), len(top)); r[f"overlap_{ty}"] = o
    if A in ("A2", "A3"):
        for e, pts in EPI_SPLIT.items():
            g = gset(pts); o = len(top & g); r[f"p_epi_{e}"] = hypergeom.sf(o - 1, N, len(g), len(top)); r[f"overlap_epi_{e}"] = o
    r["top50"] = ", ".join(sub.uns["rank_genes_groups"]["names"][k][:15]) + ", ..."
    rows[k] = r
S = pd.DataFrame(rows).T
pcols = [f"p_{ty}" for ty in sets]; P = S[pcols].astype(float)
best = P.idxmin(axis=1).str[2:]; pbest = P.min(axis=1)
ctype = best.where(pbest < ALPHA / nty, "non risolto")
if A in ("A2", "A3"):
    ctype = ctype.where(ctype != "epiteliale", pd.Series(np.where(S.p_epi_tumorale.astype(float) < S.p_epi_normale.astype(float), "epiteliale tumorale", "epiteliale normale"), index=S.index))
S["n_nuclei"] = sub.obs.cluster.value_counts().reindex(S.index); S["type_ph1"] = ctype; S["p_best"] = pbest
S.index.name = "cluster"; S = S.sort_index(key=lambda s: s.astype(int)); S.to_csv(base + "_clusters_ph1.csv")
out = lab.copy(); out["cell_type_prereg"] = out.cell_type
cl = out.cluster.astype("float").astype("Int64").astype(str)
out["cell_type"] = cl.map(ctype).where(out.cluster.notna(), "non classificabile").fillna("non classificabile")
out.to_parquet(base + "_nuclei_ph1.parquet", index=False)
pd.Series({ty: len(g) for ty, g in sets.items()}, name="n_marker_genes_in_universe").to_csv(base + "_markers_ph1.csv")
print(f"[{A} {FP} umi{UMI_MIN} seed{SEED}] PH-1 types: {out.cell_type.value_counts().to_dict()}; universe {N} genes")
print(f"  changed vs pre-registered: {(out.cell_type != out.cell_type_prereg).mean():.3f}")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(S[["n_nuclei", "type_ph1", "p_best", "top50"]].to_string())
