"""tools/S2.4_E12_summary.py — riepilogo delle esplorazioni E-1 (aree vere) ed E-2 (pannello di S2.3) per tessuto.
Uso: .venv/bin/python tools/S2.4_E12_summary.py <explore_dir> -> results/S2.4/S2.4_E1_real.csv, S2.4_E2_sim.csv, S2.4_E12_summary.csv"""
import sys, glob, json, os
import numpy as np, pandas as pd
O = sys.argv[1]; os.makedirs("results/S2.4", exist_ok=True)
R = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{O}/real/*.json"))])
S = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{O}/sim/*.json"))])
R.to_csv("results/S2.4/S2.4_E1_real.csv", index=False); S.to_csv("results/S2.4/S2.4_E2_sim.csv", index=False)
BANDS = ["2_4", "4_8", "8_16", "16_32", "32_64"]
def band_ratio(d, b):
    n = d[f"lumen_band_{b}_n"].sum(); return (d[f"lumen_band_{b}_sum"].sum() / n) / (d["lumen_ref_sum"].sum() / d["lumen_ref_n"].sum()) if n > 0 else np.nan
rows = []
qs = [f"stripes_{a}_{k}" for a in ("rows", "cols") for k in ("odd_over_even", "sd_repro_hp", "lag1_repro_hp", "corr_halves_hp", "n_hp")]
for A, d in R.groupby("archetype"):
    s = S[S.archetype == A]; sm = s.groupby("tag").mean(numeric_only=True)   # media dei due semi per area
    row = {"A": A, "n_real": len(d)}
    for q in qs + ["mean_umi_tissue", "frac_zero_tissue"]:
        v = d[q].to_numpy(float); row[f"real_{q}_mean"] = np.nanmean(v); row[f"real_{q}_maxdiff"] = np.nanmax(v) - np.nanmin(v)
        row[f"real_{q}_2se"] = 2 * np.nanstd(v, ddof=1) / np.sqrt(np.isfinite(v).sum())
        if q in sm: w = sm[q].to_numpy(float); row[f"simS23_{q}_mean"] = np.nanmean(w); row[f"simS23_{q}_maxdiff"] = np.nanmax(w) - np.nanmin(w)
    row["n_lumens_real"] = int(d.lumen_n_lumens.sum())
    for b in BANDS:
        row[f"real_band_{b}"] = band_ratio(d, b); row[f"real_band_{b}_n"] = int(d[f"lumen_band_{b}_n"].sum())
        row[f"simS23_band_{b}"] = band_ratio(s, b) if len(s) else np.nan
    rows.append(row)
T = pd.DataFrame(rows); T.to_csv("results/S2.4/S2.4_E12_summary.csv", index=False)
pd.set_option("display.width", 250)
show = ["A"] + [c for c in T.columns if ("odd_over_even" in c or "sd_repro_hp" in c or "lag1" in c) and ("_mean" in c or "maxdiff" in c or "2se" in c) and "simS23" not in c]
print(T[show].round(3).T.to_string())
print(T[["A"] + [c for c in T.columns if "simS23_stripes" in c and "_mean" in c and ("odd" in c or "sd_repro_hp" in c or "lag1" in c)]].round(3).T.to_string())
print(T[["A", "n_lumens_real"] + [c for c in T.columns if "band" in c]].round(3).T.to_string())
print(T[["A"] + [c for c in T.columns if "frac_zero" in c or "mean_umi" in c]].round(4).T.to_string())
