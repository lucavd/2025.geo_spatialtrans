import pandas as pd, numpy as np
W = "/tmp/r7_refs/"
d = pd.read_csv(W + "R7_Mc_refs.csv"); o = pd.read_csv(W + "R7_Mc_refs_v2.csv"); u = pd.read_csv(W + "use.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print("components", len(d), "ids", d.id.nunique())
print(pd.crosstab(d.interval_construction, d.archetype, margins=True))
print(pd.crosstab(d.session, d.interval_construction))
print("side", d.side.value_counts().to_dict(), "m_se", d.m_se.dropna().round(4).unique().tolist())
# recheck rule
cen = np.where(d.interval_construction == "range", (d.lo + d.hi) / 2, np.where(d.est.notna(), d.est, (d.lo + d.hi) / 2))
k = d.kappa.fillna(1.0)
print("recheck max err", float(np.nanmax(abs(cen - k * (cen - d.lo) - d.lo_new))), float(np.nanmax(abs(cen + k * (d.hi - cen) - d.hi_new))))
print("width_ratio: n_all = all components, n_finite = components with non-zero old width (zero-width ones have no ratio)")
print(d.groupby(["interval_construction", "archetype"]).width_ratio.agg(["min", "max", "size", "count"]).rename(columns={"size": "n_all", "count": "n_finite"}).round(4).to_string())
# v2 -> v4 by matching (id, metric); R2c changed composition
key = ["id", "metric"]; mm = d.merge(o[key + ["lo_new", "hi_new", "kappa"]], on=key, how="outer", suffixes=("", "_v2"), indicator=True)
print("merge", mm._merge.value_counts().to_dict())
print("only in v2:", mm[mm._merge == "right_only"][key].values.tolist()); print("only in v4:", mm[mm._merge == "left_only"][key].values.tolist())
b = mm[mm._merge == "both"]
ch = ~(np.isclose(b.lo_new, b.lo_new_v2, atol=1e-12, rtol=0) & np.isclose(b.hi_new, b.hi_new_v2, atol=1e-12, rtol=0))
print("changed v2->v4:", int(ch.sum()), "by tissue", b[ch].archetype.value_counts().to_dict(), "by construction", b[ch].interval_construction.value_counts().to_dict())
b = b.assign(w4=b.hi_new - b.lo_new, w2=b.hi_new_v2 - b.lo_new_v2)
print("width v4/v2 by construction x tissue (changed only):"); print(b[ch].assign(r=b.w4 / b.w2).groupby(["interval_construction", "archetype"]).r.agg(["min", "max", "size", "count"]).rename(columns={"size": "n_all", "count": "n_finite"}).round(4).to_string())
# S25
s = d[d.is_S25_target].merge(o[key + ["lo_new", "hi_new"]], on=key, suffixes=("", "_v2"))
print("S25", len(s)); print(s[["id", "archetype", "metric", "est", "lo", "hi", "kappa", "lo_new_v2", "hi_new_v2", "lo_new", "hi_new", "lo_stress", "hi_stress", "lo_geomA5", "hi_geomA5"]].round(4).to_string(max_colwidth=34))
# high CV
h = d[d.cv_between_areas > 0.45]
print("HIGH CV", len(h)); print(h[["id", "archetype", "interval_construction", "metric", "cv_between_areas", "lo_new", "hi_new", "lo_cv", "hi_cv"]].round(4).to_string(max_colwidth=45))
# windows
print("WINDOWS"); print(d[d.interval_kind == "window"][["id", "archetype", "interval_construction", "est", "lo", "hi", "kappa", "lo_new", "hi_new", "lo_stress", "hi_stress"]].round(4).to_string())
# geometric A5
a5 = d[d.archetype == "A5"]; corr = a5.interval_construction.isin(["boot", "se_m", "range"])
print("A5 comps", len(a5), "geomA5 filled", int(d.lo_geomA5.notna().sum()), "non-A5 filled", int(d[d.archetype != "A5"].lo_geomA5.notna().sum()))
g = a5[corr].assign(r=(a5.hi_geomA5 - a5.lo_geomA5) / (a5.hi_new - a5.lo_new))
print("A5 geom/new width ratio by qtype x construction"); print(g.groupby(["qtype", "interval_construction"]).r.agg(["min", "max", "size", "count"]).rename(columns={"size": "n_all", "count": "n_finite"}).round(4).to_string())
print("A5 lo_new<0:", a5[(a5.lo_new < 0) & (a5.lo >= 0)][["id", "metric"]].values.tolist())
print("A5 lo_geom<0:", a5[(a5.lo_geomA5 < 0) & (a5.lo >= 0)][["id", "metric"]].values.tolist())
print("A5 S25 geom:"); print(a5[a5.is_S25_target][["id", "metric", "lo_new", "hi_new", "lo_geomA5", "hi_geomA5"]].round(4).to_string(max_colwidth=40))
# stress
st = d.interval_construction.isin(["boot", "se_m", "window_one_stage"])
print("stress-affected", int(st.sum()), "unaffected", int((~st).sum()))
print("lo_new<0 (non-neg lo):", d[(d.lo_new < 0) & (d.lo >= 0)].id.tolist())
print("lo_stress<0 newly:", d[(d.lo_stress < 0) & (d.lo_new >= 0) & (d.lo >= 0)][["id", "archetype", "metric", "lo_new", "lo_stress"]].round(4).to_string(max_colwidth=50))
print("lo_cv<0:", d[(d.lo_cv < 0) & (d.lo >= 0)].id.tolist())
# usage of changed references
du = d.merge(u, on="id", how="left")
mt = du[du.margin_use.astype(str).str.upper() == "TRUE"]
print("margin_use TRUE components:"); print(mt[["id", "archetype", "interval_construction", "metric", "lo", "hi", "lo_new", "hi_new"]].round(4).to_string(max_colwidth=50))
print("margin_use by construction", du.groupby("interval_construction").margin_use.apply(lambda s: s.astype(str).value_counts().to_dict()).to_dict())
print("used_downstream of S25 ids:", du[du.is_S25_target][["id", "used_downstream"]].drop_duplicates().values.tolist())
