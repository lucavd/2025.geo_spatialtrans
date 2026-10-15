
import pandas as pd, numpy as np
me = pd.read_csv("results/S1.4b/review/code/cp1_indep_raw.csv"); me = me[me.cond != "_diag"]
p = me.pivot_table(index=["archetype", "roi"], columns="cond", values=["m1", "m2"])
perm = [c for c in p["m1"].columns if c.startswith("perm")]
g = pd.DataFrame({"m1_vor": p["m1"]["vor"], "m1_lag": p["m1"]["lag"], "m1_perm_med": p["m1"][perm].median(1), "m2_vor": p["m2"]["vor"], "m2_lag": p["m2"]["lag"], "m2_perm_med": p["m2"][perm].median(1)}).reset_index()
a = g.groupby("archetype").mean(numeric_only=True)
a["red_m2_ratio_of_means"] = 1 - a.m2_lag / a.m2_vor
a["red_m2_mean_of_ratios"] = g.assign(r=1 - g.m2_lag / g.m2_vor).groupby("archetype").r.mean()
print(a.round(4).to_string())
