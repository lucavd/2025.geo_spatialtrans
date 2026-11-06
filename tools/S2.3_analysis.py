"""tools/S2.3_analysis.py — S2.3: giudizio delle previsioni del pannello (P-5...P-9) e delle controprove CP-1...CP-3.
Simulato (media dei semi) / vero per area; stima = media delle 5 aree, +- 2 errori standard fra le aree; esiti del protocollo.
Uso: .venv/bin/python tools/S2.3_analysis.py <panel_dir> <e2_dir> -> results/S2.3/S2.3_panel_roi.csv, S2.3_P5_P9.csv, S2.3_CP.csv"""
import sys, os, json, glob
import numpy as np, pandas as pd
PD, E2 = sys.argv[1], sys.argv[2]
REAL = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{E2}/real/*.json"))]).set_index("tag")
STATS = ["panel_star", "panel_star_far", "panel_radial_R", "panel_nucleus_psi", "panel_corona_psi", "lumen_band_2_4", "lumen_band_4_8"]
LUMEN_TISSUES = ["A1", "A2", "A5", "A6"]
def lum(d, lo, hi):
    return (d[f"lumen_band_{lo}_{hi}_sum"] / d[f"lumen_band_{lo}_{hi}_n"]) / (d["lumen_ref_sum"] / d["lumen_ref_n"]) if d["lumen_ref_n"] > 0 and d[f"lumen_band_{lo}_{hi}_n"] > 0 else np.nan
def load(var):
    fs = sorted(glob.glob(f"{PD}/{var}/stats/*.json"))
    if not fs: return None
    T = pd.DataFrame([json.load(open(f)) for f in fs])
    if "error" in T.columns and T["error"].notna().any(): raise RuntimeError(T.loc[T["error"].notna(), "file"].tolist())
    T["lumen_band_2_4"] = [lum(r, 2, 4) for _, r in T.iterrows()]; T["lumen_band_4_8"] = [lum(r, 4, 8) for _, r in T.iterrows()]
    return T.groupby(["tag", "archetype"])[STATS + ["panel_nucleus_mbar", "panel_mean_umi", "lumen_n_lumens"]].mean().reset_index()
def ratios(var):
    T = load(var)
    if T is None: return None
    R = REAL.copy(); R["lumen_band_2_4"] = [lum(r, 2, 4) for _, r in R.iterrows()]; R["lumen_band_4_8"] = [lum(r, 4, 8) for _, r in R.iterrows()]
    out = T[["tag", "archetype"]].copy()
    for s in STATS:
        out[s + "_sim"] = T[s].to_numpy(); out[s + "_real"] = R.loc[T.tag, s].to_numpy(); out[s] = out[s + "_sim"] / out[s + "_real"]
    out["mean_umi_sim"] = T.panel_mean_umi.to_numpy(); out["mean_umi_real"] = R.loc[T.tag, "panel_mean_umi"].to_numpy()
    out["nucleus_mbar_sim"] = T.panel_nucleus_mbar.to_numpy(); out["nucleus_mbar_real"] = R.loc[T.tag, "panel_nucleus_mbar"].to_numpy()
    out["n_lumens"] = R.loc[T.tag, "lumen_n_lumens"].to_numpy(); out["variant"] = var
    return out
def judge(df, s, lo, hi, tissues=None):
    rows = []
    for A, d in df.groupby("archetype"):
        if tissues is not None and A not in tissues: continue
        x = d[s].to_numpy(float); x = x[np.isfinite(x)]; n = len(x)
        if n == 0: rows.append({"archetype": A, "stat": s, "n_roi": 0, "verdict": "non calcolabile"}); continue
        est = x.mean(); se = x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
        a, b = est - 2 * se, est + 2 * se
        v = "confermata" if (a >= lo and b <= hi) else ("smentita" if (b < lo or a > hi) else "indeterminata")
        rows.append({"archetype": A, "stat": s, "n_roi": n, "estimate": est, "lo2se": a, "hi2se": b, "int_lo": lo, "int_hi": hi, "verdict": v})
    return rows
os.makedirs("results/S2.3", exist_ok=True)
allr = [x for x in (ratios(v) for v in ("main", "field_off", "narrow", "r1")) if x is not None]
pd.concat(allr).to_csv("results/S2.3/S2.3_panel_roi.csv", index=False)
M = [x for x in allr if x.variant.iloc[0] == "main"]
if M:
    M = M[0]; J = []
    for s in ("panel_star", "panel_star_far", "panel_radial_R", "panel_nucleus_psi", "panel_corona_psi"): J += judge(M, s, 0.909, 1.1)
    for s in ("lumen_band_2_4", "lumen_band_4_8"): J += judge(M, s, 0.8, 1.25, LUMEN_TISSUES)
    J = pd.DataFrame(J); J.to_csv("results/S2.3/S2.3_P5_P9.csv", index=False); pd.set_option("display.width", 200); print(J.round(4).to_string())
CP = []
for x in allr:
    v = x.variant.iloc[0]
    if v == "field_off": CP += [dict(r, cp="CP-1") for r in judge(x, "panel_nucleus_psi", 0.909, 1.1)]
    if v == "narrow": CP += [dict(r, cp="CP-2") for s in ("lumen_band_4_8",) for r in judge(x, s, 0.8, 1.25, LUMEN_TISSUES)]
    if v == "r1": CP += [dict(r, cp="CP-3") for s in ("panel_star", "panel_radial_R") for r in judge(x, s, 0.909, 1.1)]
if CP: C = pd.DataFrame(CP); C.to_csv("results/S2.3/S2.3_CP.csv", index=False); print(C.round(4).to_string())
