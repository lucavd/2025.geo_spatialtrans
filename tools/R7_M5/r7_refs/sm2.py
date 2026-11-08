import pandas as pd, numpy as np
W = "/tmp/r7_refs/"
d = pd.read_csv(W + "R7_Mc_refs.csv"); o = pd.read_csv(W + "R7_Mc_refs_v1.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print("components", len(d), "ids", d.id.nunique(), "same keys as v1:", bool((d[["id", "metric"]].values == o[["id", "metric"]].values).all()))
num = ["lo_new", "hi_new", "lo_cv", "hi_cv", "kappa"]
chg = np.zeros(len(d), bool)
for c in num:
    chg |= ~np.isclose(d[c].fillna(-9e9), o[c].fillna(-9e9), rtol=0, atol=1e-12)
print("changed components by tissue:", d[chg].archetype.value_counts().to_dict(), "by kind:", d[chg].interval_kind.value_counts().to_dict())
print("unchanged outside A3/A4:", bool(not d[chg & ~d.archetype.isin(["A3", "A4"])].shape[0]))
c = np.where(d.interval_kind == "range", (d.lo + d.hi) / 2, np.where(d.est.notna(), d.est, (d.lo + d.hi) / 2))
k = d.kappa.fillna(1.0); print("recheck max err", float(np.nanmax(abs(c - k * (c - d.lo) - d.lo_new))), float(np.nanmax(abs(c + k * (d.hi - c) - d.hi_new))))
print(d.groupby(["interval_kind", "archetype"]).width_ratio.agg(["min", "max", "count"]).round(4).to_string())
h = d[d.cv_between_areas > 0.45]
print("HIGH CV", len(h)); print(h[["id", "archetype", "interval_kind", "metric", "cv_between_areas", "lo_new", "hi_new", "lo_cv", "hi_cv"]].round(4).to_string(max_colwidth=50))
s = d[d.is_S25_target & d.archetype.isin(["A3", "A4"])].merge(o[["id", "metric", "lo_new", "hi_new"]], on=["id", "metric"], suffixes=("", "_v1"))
print("S25 A3/A4"); print(s[["id", "archetype", "metric", "est", "lo", "hi", "lo_new_v1", "hi_new_v1", "kappa", "lo_new", "hi_new"]].round(4).to_string(max_colwidth=40))
a5 = d[(d.archetype == "A5") & (d.qtype == "nuc")]
a5 = a5.assign(w_new=a5.hi_new - a5.lo_new, w5=a5.hi_A5nuc - a5.lo_A5nuc)
print("A5 nuc", len(a5), a5.interval_kind.value_counts().to_dict(), "A5nuc filled", int(d.lo_A5nuc.notna().sum()), "non-A5nuc filled", int(d[~((d.archetype == "A5") & (d.qtype == "nuc"))].lo_A5nuc.notna().sum()))
print(a5[["id", "interval_kind", "metric", "lo", "hi", "lo_new", "hi_new", "lo_A5nuc", "hi_A5nuc"]].assign(r=(a5.w5 / a5.w_new)).round(4).to_string(max_colwidth=50))
print("lo_A5nuc<0:", a5[a5.lo_A5nuc < 0].id.tolist())
print("negative lo_new:", d[d.notes.str.contains("lo_new < 0", na=False)].id.tolist())
print("lo_cv<0:", d[d.notes.str.contains("lo_cv < 0", na=False)].id.tolist())
