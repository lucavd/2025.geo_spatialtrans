import pandas as pd, numpy as np
W = "/tmp/r7_refs/"
d = pd.read_csv(W + "R7_Mc_refs.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print("components", len(d), "ids", d.id.nunique())
print(pd.crosstab(d.interval_kind, d.archetype, margins=True))
print("ids by primary-kind set:"); print(d.groupby("id").interval_kind.agg(lambda s: "+".join(sorted(set(s)))).value_counts())
print(pd.crosstab(d.session, d.interval_kind))
print("qtype", d.qtype.value_counts().to_dict())
# independent recheck of the scaling rule
c = np.where(d.interval_kind == "range", (d.lo + d.hi) / 2, np.where(d.est.notna(), d.est, (d.lo + d.hi) / 2))
k = d.kappa.fillna(1.0)
lo2 = c - k * (c - d.lo); hi2 = c + k * (d.hi - c)
print("max |lo_new - recheck|", float(np.nanmax(abs(lo2 - d.lo_new))), "max |hi_new - recheck|", float(np.nanmax(abs(hi2 - d.hi_new))))
print("width_ratio by kind:"); print(d.groupby(["interval_kind", "archetype"]).width_ratio.agg(["min", "max", "count"]).round(4).to_string())
print("CV available:", int(d.cv_between_areas.notna().sum()))
h = d[d.cv_between_areas > 0.45]
print("HIGH CV", len(h))
print(h[["id", "archetype", "interval_kind", "metric", "cv_between_areas", "lo", "hi", "lo_new", "hi_new", "lo_cv", "hi_cv"]].round(4).to_string(max_colwidth=60))
print("capped:", h[h.notes.str.contains("capped", na=False)].id.tolist())
s = d[d.is_S25_target]
print("S25", len(s))
print(s[["id", "archetype", "metric", "est", "lo", "hi", "kappa", "lo_new", "hi_new", "cv_between_areas"]].round(4).to_string(max_colwidth=45))
print("negative lo_new notes:", d[d.notes.str.contains("lo_new < 0", na=False)][["id", "metric", "lo_new"]].round(4).to_string(max_colwidth=50))
print("lo_cv<0:", d[d.notes.str.contains("lo_cv < 0", na=False)].id.tolist())
print("zero width:", d[d.notes.str.contains("zero-width", na=False)][["id", "metric"]].to_string(max_colwidth=50))
print("S1.4b matched per-ROI:", d[d.notes.str.contains("b8_circle", na=False)].id.tolist())
