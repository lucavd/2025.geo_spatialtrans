#!/usr/bin/env python3
"""R5 PH-2 summary: intra-tissue Gaussian blur sigma_in at which the two-density cell model
(r fitted to the R4 PH-1 statistic) also reproduces the far statistic (bins >2 um from nucleus edge).
Usage: R5_ph2_summary.py OUT [SUBDIR=ph2]  -> OUT/tables/R5_PH2_inner[_SUBDIR].csv, R5_PH2_inner_summary[_SUBDIR].csv"""
import json, glob, sys, numpy as np, pandas as pd
out = sys.argv[1]; sub = sys.argv[2] if len(sys.argv) > 2 else "ph2"
sfx = "" if sub == "ph2" else "_" + sub
D = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{out}/{sub}/*.json"))])
D["d"] = D.pred_far - D.obs_far
D = D.sort_values(["tag", "sigma_in"])
D.to_csv(f"{out}/tables/R5_PH2_inner{sfx}.csv", index=False)
rows = []
for tag, g in D.groupby("tag"):
    g = g[np.isfinite(g.d)]
    s, d = g.sigma_in.values, g.d.values
    star, kind = np.nan, ""
    for i in range(len(s) - 1):
        if d[i] == 0: star, kind = s[i], "crossing"; break
        if d[i] * d[i + 1] < 0:
            star, kind = s[i] + (s[i + 1] - s[i]) * (-d[i]) / (d[i + 1] - d[i]), "crossing"; break
    if not kind:
        j = int(np.argmin(np.abs(d)))
        kind = f"no crossing (closest at {s[j]:g} um, d={d[j]:+.3f})"
        if j == 0: kind = "no crossing, below grid? " + kind
    r_star = np.interp(star, s, g.r.values) if np.isfinite(star) else np.nan
    rows.append(dict(tag=tag, key=tag.split("_")[0], obs_star=g.obs_star.iloc[0] if "obs_star" in g else np.nan,
                     obs_far=g.obs_far.iloc[0], sigma_in_star=star, r_at_star=r_star, n_sigma=len(s),
                     sigma_max_ok=s.max(), result=kind))
S = pd.DataFrame(rows)
S.to_csv(f"{out}/tables/R5_PH2_inner_summary{sfx}.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(S.round(3).to_string(index=False))
