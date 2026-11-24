"""tools/sonda/S5_table.py - S5 (docs/sonda/REGOLA_S5.md, «Metro sul pannello»): tabella per area dei bracci ristretti a G_p (V, T, B, P, XB, XV)
con tools/sonda/S3_table.py invariato sulle misure di /mnt/micron/geo_spatialtrans/S5_sonda/measures; la sola colonna genes passa da "all" a "panel".
Uscita: results/sonda/S5/area_values.csv. Uso: python3 tools/sonda/S5_table.py"""
import os, subprocess, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); tmp = "results/sonda/S5/area_values_tmp.csv"
subprocess.run(["python3", "tools/sonda/S3_table.py", "/mnt/micron/geo_spatialtrans/S5_sonda/measures", tmp], check=True)
U = pd.read_csv(tmp); os.remove(tmp); assert set(U.genes) == {"all"}; U["genes"] = "panel"
U.to_csv("results/sonda/S5/area_values.csv", index=False); print(U.groupby(["arm", "A"]).size().to_dict(), "rctd_n_genes", sorted(U.rctd_n_genes.unique().tolist()))
