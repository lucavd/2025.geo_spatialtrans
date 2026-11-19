#!/usr/bin/env python3
"""tools/F3.2_select.py - F3.2: aree di verifica finali (DEC-095 con la modifica del 2026-11-19) dalle candidate e dalla tabella delle tre letture.
Legge results/F3.2/F3.2_area_agreement.csv (colonne sheet, n, agent, blind, rule3, level, selected_class, kept) e le candidate
(F3.2_pool_<label>.csv, F3.2_pool500_A4.csv). Scrive results/F3.2/F3.2_selected_<label>.csv (aree tenute e tessera di prova) e F3.2_selected_all.csv.
Uso: python3 tools/F3.2_select.py
"""
import pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results/F3.2"
ag = pd.read_csv(R / "F3.2_area_agreement.csv")
allrows = []
for sheet, g in ag.groupby("sheet"):
    label = "A4" if sheet == "A4_500um" else sheet
    pool = pd.read_csv(R / ("F3.2_pool500_A4.csv" if sheet == "A4_500um" else f"F3.2_pool_{label}.csv"))
    test = pd.read_csv(R / f"F3.2_pool_{label}.csv"); test = test[test.level == "test_tile"]
    k = g[g.kept].merge(pool[pool.level != "test_tile"], left_on="n", right_on="order", how="left", suffixes=("", "_pool"))
    assert k.x0.notna().all(), sheet
    k = k.assign(label=label, side_um=(k.x1 - k.x0).round(3), klass=k.selected_class, area_id=[f"{label}_{c}_{i+1}" for i, c in enumerate(k.selected_class)])
    cols = ["label", "level", "klass", "area_id", "n", "x0", "y0", "x1", "y1", "side_um", "agent", "blind", "rule3", "blind_uncertain"]
    out = pd.concat([k[cols], test.assign(klass="stardist_test", n=0, side_um=250.0, agent="", blind="", rule3="", blind_uncertain=False)[cols]])
    fo = R / f"F3.2_selected_{sheet}.csv"; out.to_csv(fo, index=False); allrows.append(out[out.level != "test_tile"])
    print(sheet, len(k), "aree ->", fo.name)
pd.concat(allrows).to_csv(R / "F3.2_selected_all.csv", index=False)
