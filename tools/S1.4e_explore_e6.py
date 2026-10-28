#!/usr/bin/env python
"""tools/S1.4e_explore_e6.py - S1.4e, esplorazione dichiarata (docs/steps/S1.4e.md), E6: scarto fra i punti manuali (R2b) e il centroide
del nucleo del segmentatore agganciato (colonne match_/label_ di R2b_points.csv; nuclei di R2 a scala nativa), per tessuto, segmentatore
e annotatore. Limite superiore del rumore dei punti: comprende l'errore del centroide del segmentatore.
Uscita: results/S1.4e/explore/S1.4e_E6_click_offset.csv"""
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/S1.4e/explore"
SHAPE_SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}

def main():
    p = pd.read_csv(ROOT / "results/R2b/R2b_points.csv"); w = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv").set_index("win_id")
    out = []
    for A in SHAPE_SEG:
        for seg in ["cellpose_rgb", "stardist_he"]:
            q = p[(p.win_id.str[:2] == A) & p[f"match_{seg}"].fillna(False).astype(bool) & p.rater.isin(["Luca", "Federica"])]
            for (wid, rater), g in q.groupby(["win_id", "rater"]):
                nuc = pd.read_parquet(ROOT / f"results/R2/nuclei/{A}_{w.loc[wid].roi_id}_{seg}_native.parquet", columns=["label", "x_um", "y_um", "scale"])
                nuc = nuc[nuc.scale == 1].set_index("label")
                xy = nuc.loc[g[f"label_{seg}"].astype(int), ["x_um", "y_um"]].to_numpy(); dx = g[["x_um", "y_um"]].to_numpy() - xy
                out += [dict(archetype=A, seg=seg, rater=rater, dx=a, dy=b) for a, b in dx]
    d = pd.DataFrame(out); d["r"] = np.hypot(d.dx, d.dy)
    s = d.groupby(["archetype", "seg", "rater"]).agg(n=("r", "size"), r_med=("r", "median"), sd_x=("dx", "std"), sd_y=("dy", "std"),
                                                       mean_x=("dx", "mean"), mean_y=("dy", "mean")).reset_index()
    s["sd_axis"] = (s.sd_x + s.sd_y) / 2; s["shape_seg"] = s.seg == s.archetype.map(SHAPE_SEG)
    s.to_csv(OUT / "S1.4e_E6_click_offset.csv", index=False)
    print(s[(s.rater == "Luca") & s.shape_seg][["archetype", "seg", "n", "r_med", "sd_axis"]].round(3).to_string())

if __name__ == "__main__":
    main()
