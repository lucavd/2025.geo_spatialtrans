
import pandas as pd, numpy as np, os, sys, re, ast
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0,"tools")
from R6_catalogue import CAT as CAT_cat, EPI_SPLIT as EPI_cat, SPECIES as SP_cat
src = open("tools/R6_cluster.py").read()
m = re.search(r"CAT = (\{.*?\})\nCAT\[\"A3\"\] = CAT\[\"A2\"\]", src, re.S); CAT_cl = ast.literal_eval(m.group(1)); CAT_cl["A3"]=CAT_cl["A2"]
m2 = re.search(r"EPI_SPLIT = (\{.*?\})\n", src); EPI_cl = ast.literal_eval(m2.group(1))
print("CAT identical cluster.py vs catalogue.py:", CAT_cl == CAT_cat, "| EPI_SPLIT identical:", EPI_cl == EPI_cat)
for A in CAT_cl:
    if CAT_cl[A] != CAT_cat.get(A): print("  diff", A, {k:(CAT_cl[A].get(k), CAT_cat[A].get(k)) for k in set(CAT_cl[A])|set(CAT_cat[A]) if CAT_cl[A].get(k)!=CAT_cat[A].get(k)})
M="/mnt/micron/geo_spatialtrans/R6"
for A in ["A1","A4","A5"]:
    cl = pd.read_csv(f"{M}/labels/{A}_exp_umi20_seed20261022_clusters.csv")
    print(A, "clusters", len(cl), "n_nuclei quantiles", cl.n_nuclei.quantile([0,.1,.5,.9,1]).astype(int).tolist(), "clusters <50 nuclei:", int((cl.n_nuclei<50).sum()), "nuclei in them:", int(cl.n_nuclei[cl.n_nuclei<50].sum()))
    print("  markers present:", pd.read_csv(f"{M}/labels/{A}_exp_umi20_seed20261022_markers.csv", index_col=0).iloc[:,0].to_dict())
    print("  markers in universe (ph1):", pd.read_csv(f"{M}/labels/{A}_exp_umi20_seed20261022_markers_ph1.csv", index_col=0).iloc[:,0].to_dict())
cl4 = pd.read_csv(f"{M}/labels/A4_exp_umi20_seed20261022_clusters_ph1.csv")
print("A4 ph1: non risolto clusters", int((cl4.type_ph1=="non risolto").sum()), "nuclei", int(cl4.n_nuclei[cl4.type_ph1=="non risolto"].sum()), "; median n in non risolto clusters", cl4.n_nuclei[cl4.type_ph1=="non risolto"].median(), "vs resolved", cl4.n_nuclei[cl4.type_ph1!="non risolto"].median())
pan = pd.read_csv(f"{M}/PanglaoDB_omnipath_mouse.tsv", sep="\t"); print(pan.columns.tolist()); print(pan.head(3).to_string())
print("mouse table genesymbol examples:", pan.genesymbol.head(10).tolist(), "| frac all-uppercase symbols:", pan.genesymbol.astype(str).str.isupper().mean().round(3))
