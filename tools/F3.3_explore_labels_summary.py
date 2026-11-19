"""tools/F3.3_explore_labels_summary.py - F3.3 E-2: riepilogo per tessuto e lato delle etichette trasferite (media delle aree)."""
import glob, json, pandas as pd
O = "/mnt/micron/geo_spatialtrans/F3.3/explore/E2"
d = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{O}/*.json"))])
d.to_csv(f"{O}/E2_areas.csv", index=False)
cols = ["agree_with_orig", "M3_f31", "M3_orig", "M3_transfer", "M5_f31", "M5_orig", "M5_transfer"]
s = d.groupby(["A", "side"])[cols].mean().round(4)
s["b2c_identical"] = d.groupby(["A", "side"]).b2c_identical.all(); s["n"] = d.groupby(["A", "side"]).size()
s.to_csv(f"{O}/E2_summary.csv"); print(s.to_string())
