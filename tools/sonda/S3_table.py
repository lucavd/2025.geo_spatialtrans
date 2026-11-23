"""tools/sonda/S3_table.py - S3: tabella per area di tutti i bracci (V, T, P, B) dai JSON di tools/sonda/S3_measures.py (etichette trasferite,
labels = tl). Colonne nel formato letto da tools/riesame_D.py: arm, genes (= all), A (tumour|stroma), unit (area: per B senza il seme, per P
l'area di V), run, set, k, M1...M6 e conteggi. Uso: python3 tools/sonda/S3_table.py [measures_root] [out_csv]"""
import os, sys, glob, json, pandas as pd
R = sys.argv[1] if len(sys.argv) > 1 else "/mnt/micron/geo_spatialtrans/S3_sonda/measures"
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser("~/2025.geo_spatialtrans/results/sonda/S3/area_values.csv")
rows = []
for f in sorted(glob.glob(f"{R}/*/*.json")):
    for r in json.load(open(f)):
        if r["labels"] != "tl": continue
        r = dict(r); r["genes"] = "all"; r["run"] = r["unit"]
        if r["arm"] == "B": r["unit"] = r["unit"].rsplit("_k", 1)[0]
        rows.append(r)
U = pd.DataFrame(rows)
cols = ["arm", "genes", "A", "unit", "run", "set", "k", "M1", "M2", "M3", "M4", "M5", "M6", "n_cells", "n_typed", "n_nuc_blocks", "n_pairs", "umi_reg", "rctd_n_blocks", "rctd_n_run", "rctd_n_genes", "rctd_sec"]
U = U[cols].sort_values(["arm", "A", "run"]).reset_index(drop=True); U.to_csv(OUT, index=False)
print(U.groupby(["arm", "A"]).size().to_dict(), "missing M:", int(U[["M1", "M2", "M3", "M4", "M5", "M6"]].isna().sum().sum()))
