"""tools/sonda/S3b_table.py - S3b: tabella per area = righe V, B, P di results/sonda/S3/area_values.csv (invariate) + T omologo
(tools/sonda/S3_table.py sulle misure di /mnt/micron/geo_spatialtrans/S3b_sonda/measures). Uscita: results/sonda/S3b/area_values.csv"""
import os, subprocess, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
subprocess.run(["python3", "tools/sonda/S3_table.py", "/mnt/micron/geo_spatialtrans/S3b_sonda/measures", "results/sonda/S3b/T_hom_values.csv"], check=True)
old = pd.read_csv("results/sonda/S3/area_values.csv"); new = pd.read_csv("results/sonda/S3b/T_hom_values.csv")
assert set(new.arm) == {"T"} and len(new) == 10, (set(new.arm), len(new))
U = pd.concat([old[old.arm != "T"], new], ignore_index=True).sort_values(["arm", "A", "run"]).reset_index(drop=True); U.to_csv("results/sonda/S3b/area_values.csv", index=False)
print(U.groupby(["arm", "A"]).size().to_dict())
