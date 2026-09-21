#!/usr/bin/env python
"""tools/R2_consensus.py — stima di consenso della densita' nucleare per ROI: nuclei appaiati Cellpose-RGB/StarDist (IoU>=0.5)
+ esclusivi di ciascun metodo pesati per la frazione con evidenza di ematossilina (R2_exclusive_check.csv). Output: results/R2/R2_consensus_density.csv"""
import sys, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); import R2_common as C
R = C.ROOT / "results/R2"
ex = pd.read_csv(R / "R2_exclusive_check.csv"); rs = pd.read_csv(R / "R2_roi_summary.csv")
nat = rs[rs.scale == 1].set_index(["archetype", "roi_id", "method"])
rows = []
for (A, roi), g in ex.groupby(["archetype", "roi_id"]):
    cp = g[(g.method == "cellpose_rgb") & (g.other == "stardist_he")].set_index("cls"); sd = g[(g.method == "stardist_he") & (g.other == "cellpose_rgb")].set_index("cls")
    n_m = cp.loc["matched", "n"]; n_cpx = cp.loc["exclusive", "n"] * cp.loc["exclusive", "frac_nuclear_evidence"]; n_sdx = sd.loc["exclusive", "n"] * sd.loc["exclusive", "frac_nuclear_evidence"]
    area = nat.loc[(A, roi, "cellpose_rgb"), "tissue_mm2"]
    rows.append(dict(archetype=A, roi_id=roi, tissue_mm2=area, n_matched=int(n_m), n_cp_excl_with_evidence=round(n_cpx), n_sd_excl_with_evidence=round(n_sdx),
                     density_consensus=(n_m + n_cpx + n_sdx) / area, density_cellpose=nat.loc[(A, roi, "cellpose_rgb"), "density_per_mm2"], density_stardist=nat.loc[(A, roi, "stardist_he"), "density_per_mm2"],
                     density_spaceranger=nat.loc[(A, roi, "spaceranger"), "density_per_mm2"] if (A, roi, "spaceranger") in nat.index else float("nan")))
d = pd.DataFrame(rows); d.to_csv(R / "R2_consensus_density.csv", index=False)
print(d.groupby("archetype").agg(consensus_median=("density_consensus", "median"), consensus_min=("density_consensus", "min"), consensus_max=("density_consensus", "max")).round(0))
