"""S2.2, esplorazione E-1: tabella per ROI delle misure vere del banco e intervalli di rilevanza per tessuto
(differenza massima fra due aree vere dello stesso tessuto, relativa alla media: (max - min) / media).
m1 e m2 normalizzati per la profondità (geni per UMI), m3 e m4 grezzi. Uso: python3 tools/S2.2_E1_summary.py"""
import json, glob, pandas as pd
rows = [json.load(open(f)) for f in sorted(glob.glob("results/S2.2/E1/*_real.json"))]
d = pd.DataFrame(rows).sort_values(["archetype", "roi_id"])
d["m1n"] = d.m1 / d.umi_per_bin
d["m2n"] = d.m2 / d.umi_per_block
keep = ["archetype", "roi_id", "n_cells", "n_bins", "frac_bins_eval", "umi_per_bin", "m1", "m2", "m3", "m4_within", "m4_between", "m4_ratio",
        "m1n", "m2n", "umi_per_block", "n_blocks", "n_pairs", "n_cells_m4", "n_genes", "types"]
d[keep].to_csv("results/S2.2/S2.2_E1_real.csv", index=False)
g = d.groupby("archetype")
iv = pd.DataFrame({m: (g[m].max() - g[m].min()) / g[m].mean() for m in ["m1n", "m2n", "m3", "m4_ratio"]})
iv.columns = ["rel_m1", "rel_m2", "rel_m3", "rel_m4"]
iv.round(4).to_csv("results/S2.2/S2.2_E1_relevance.csv")
print(iv.round(3).to_string())
