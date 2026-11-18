"""tools/F3.1_PH6_variability.py - F3.1, indagine PH-6 (dichiarata il 2026-11-18): variabilita' fra le 5 aree, simulato contro vero.
CV fra aree (sd con ddof = 1 / media) delle misure M-1...M-6 e delle quantita' di base: nuclei o cellule con meta' del nucleo in tessuto per mm^2 di tessuto
(quadratini in tessuto x 4 um^2), UMI per quadratino in tessuto, quota per tipo (tipi del riferimento). Simulato: media dei semi per area (U6) e Q6 (k = 0).
Uso: .venv/bin/python tools/F3.1_PH6_variability.py <main_dir> <results_dir>"""
import sys, os, glob, json, numpy as np, pandas as pd
main, res = sys.argv[1:3]; rows = []
MEAS = {"M-1": "share_assigned_exp", "M-2": "umi_per_cell_exp", "M-3": "foreign_exp", "M-4": "rctd_frac_singlet", "M-5": "rctd_agree_nucleus", "M-6": "rctd_neighbour_same"}
for side, pat in [("real", f"{main}/real/*"), ("sim", f"{main}/sim/*_k*"), ("simq6", f"{main}/simq6/*_k0")]:
    for d in sorted(glob.glob(pat)):
        b = os.path.basename(d); A = b[:2]; j = json.load(open(f"{d}/b2c_label.json")); rj = json.load(open(f"{d}/rctd.json"))
        r = {"side": side, "tag": b.split("_k")[0], "A": A, "density_per_mm2": j["n_cells"] / (j["n_reg"] * 4e-6), "umi_per_bin": j["umi_reg"] / j["n_reg"]}
        for k, v in MEAS.items(): r[k] = j[v] if v in j else rj[v.replace("rctd_", "")]
        types = set(pd.read_csv(f"{main}/ref/markers_{A}.csv").cell_type)
        ct = pd.read_parquet(f"{d}/cells.parquet").cell_type; ct = ct[ct.isin(types)]
        for t, v in ct.value_counts(normalize=True).items(): r[f"share::{t}"] = v
        rows.append(r)
T = pd.DataFrame(rows).groupby(["side", "A", "tag"]).mean(numeric_only=True).reset_index(); T.to_csv(f"{res}/F3.1_PH6_areas.csv", index=False)
cols = [c for c in T.columns if c not in ("side", "A", "tag")]
cv = T.groupby(["A", "side"])[cols].agg(lambda x: x.std(ddof=1) / x.mean() if x.notna().sum() > 1 and x.mean() != 0 else np.nan)
cv.to_csv(f"{res}/F3.1_PH6_cv.csv")
base = ["density_per_mm2", "umi_per_bin"] + [c for c in cols if c.startswith("share::")]
summ = []
for A in sorted(T.A.unique()):
    for side in ["sim", "simq6"]:
        if (A, side) not in cv.index: continue
        rb = cv.loc[(A, side), base] / cv.loc[(A, "real"), base]; rm = cv.loc[(A, side), list(MEAS)] / cv.loc[(A, "real"), list(MEAS)]
        summ.append({"A": A, "side": side, "median_cv_ratio_base": float(np.nanmedian(rb)), "cv_ratio_density": float(rb["density_per_mm2"]),
                     "median_cv_ratio_types": float(np.nanmedian(rb[[c for c in base if c.startswith("share::")]])), "median_cv_ratio_methods": float(np.nanmedian(rm))})
S = pd.DataFrame(summ); S.to_csv(f"{res}/F3.1_PH6_summary.csv", index=False); pd.set_option("display.width", 220)
print(cv.round(3).T.to_string()); print(S.round(3).to_string(index=False))
