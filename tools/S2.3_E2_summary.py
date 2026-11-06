"""tools/S2.3_E2_summary.py — S2.3, esplorazione E-2: tabella delle statistiche vere per area e intervalli di P-2, P-3 (differenza
massima fra due aree vere dello stesso tessuto) e 2 errori standard tipici (2 x deviazione standard fra le aree / sqrt(n), relativa alla
media per i rapporti). Uso: .venv/bin/python tools/S2.3_E2_summary.py <e2dir> -> results/S2.3/S2.3_E2_real.csv, S2.3_E2_relevance.csv"""
import sys, json, glob, os
import numpy as np, pandas as pd
D = sys.argv[1]; os.makedirs("results/S2.3", exist_ok=True)
T = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{D}/real/*.json"))])
assert "error" not in T.columns or T["error"].isna().all(), T.loc[T["error"].notna(), "tag"].tolist()
for lo, hi in ((2, 4), (4, 8)):
    T[f"lumen_band_{lo}_{hi}_ratio"] = (T[f"lumen_band_{lo}_{hi}_sum"] / T[f"lumen_band_{lo}_{hi}_n"]) / (T["lumen_ref_sum"] / T["lumen_ref_n"])
T.to_csv("results/S2.3/S2.3_E2_real.csv", index=False)
rows = []
def add(q, kind, col):
    for A, g in T.groupby("archetype"):
        x = g[col].to_numpy(float); x = x[np.isfinite(x)]
        if len(x) < 2: rows.append({"quantity": q, "archetype": A, "n": len(x)}); continue
        mean = x.mean(); sd = x.std(ddof=1)
        r = {"quantity": q, "archetype": A, "n": len(x), "mean": mean, "min": x.min(), "max": x.max(), "sd": sd}
        if kind == "diff":
            r["interval_halfwidth"] = x.max() - x.min(); r["two_se_typical"] = 2 * sd / np.sqrt(len(x))
        else:
            r["two_se_typical_rel"] = 2 * sd / np.sqrt(len(x)) / mean
        rows.append(r)
for grp in ("nucleus", "corona", "acellular"):
    add(f"bench_{grp}_psi", "ratio", f"bench_{grp}_psi"); add(f"bench_{grp}_rho1", "diff", f"bench_{grp}_rho1"); add(f"bench_{grp}_rho2", "diff", f"bench_{grp}_rho2")
    add(f"bench_{grp}_n_groups", "ratio", f"bench_{grp}_n_groups")
for q in ("panel_star", "panel_star_far", "panel_radial_R", "panel_nucleus_psi", "panel_corona_psi", "lumen_band_2_4_ratio", "lumen_band_4_8_ratio"):
    add(q, "ratio", q)
R = pd.DataFrame(rows); R.to_csv("results/S2.3/S2.3_E2_relevance.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(R.round(4).to_string())
print(T.groupby("archetype")[["lumen_n_lumens"]].sum())
