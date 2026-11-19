#!/usr/bin/env python3
"""tools/F3.2_sexcheck.py - F3.2: il tumore del colon da 6,5 mm e' di un paziente diverso da P2 (taratura)? Geni dei cromosomi sessuali.
Per sezione: UMI di XIST e dei geni del cromosoma Y (RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY, USP9Y, ZFY) per 10 000 UMI, quadratini di 8 um in tessuto.
Non e' una grandezza confrontata con il simulatore (identita' del campione). Se il sesso differisce da P2, e' un altro paziente; se coincide resta incerto.
Uso: .venv/bin/python tools/F3.2_sexcheck.py -> results/F3.2/F3.2_sexcheck.csv
"""
import h5py, numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; MIC = Path("/mnt/micron/geo_spatialtrans")
S = {"P2 (taratura, 4.0.1)": MIC / "data_real/Visium_HD_Human_Colon_Cancer", "P1 (3.0.0)": MIC / "data_holdout/CRC_P1_SR300",
     "P5 (3.0.0)": MIC / "data_holdout/CRC_P5_SR300", "6,5 mm (4.1.0)": MIC / "data_holdout/CRC_6p5mm_SR410"}
Y = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY"]
rows = []
for k, d in S.items():
    f = h5py.File(d / "binned_outputs_x/binned_outputs/square_008um/filtered_feature_bc_matrix.h5", "r")["matrix"]
    names = [n.decode() for n in f["features/name"][:]]; ind = f["indices"][:]; dat = f["data"][:]
    tot = float(dat.sum()); gsum = np.bincount(ind, weights=dat, minlength=len(names))
    rec = dict(section=k, total_umi=tot, XIST_per_1e4=gsum[names.index("XIST")] / tot * 1e4 if "XIST" in names else np.nan)
    ys = [g for g in Y if g in names]; rec["Y_genes_found"] = len(ys); rec["Y_per_1e4"] = sum(gsum[names.index(g)] for g in ys) / tot * 1e4
    rec["sex_call"] = "M" if rec["Y_per_1e4"] > 10 * max(rec["XIST_per_1e4"], 1e-3) else ("F" if rec["XIST_per_1e4"] > 10 * max(rec["Y_per_1e4"], 1e-3) else "incerto")
    rows.append(rec); print(rec, flush=True)
pd.DataFrame(rows).to_csv(ROOT / "results/F3.2/F3.2_sexcheck.csv", index=False)
