
import json, pandas as pd
s=json.load(open("/tmp/r7_S13_S14/summary.json"))
print("REPRO", s["repro"], "NOTREP", s["log"]["not_reproduced"], "NA5", s["n_a5rows"])
print("AFF", s["aff_by_step"])
for x in s["trans_by_step"]: print("TBS", x)
for r in s["changed"]: print("CH", r["step"],"|",r["prediction_id"],"|",r["archetype"],"|",str(r["verdict_recorded"])[:55],"->",str(r["verdict_new"])[:55],"| st:",str(r["verdict_new_stress"])[:55])
for r in s["changed_stress_only"]: print("ST", r["step"],"|",r["prediction_id"],"|",r["archetype"],"|",str(r["verdict_recorded"])[:50],"->",str(r["verdict_new"])[:50],"| st:",str(r["verdict_new_stress"])[:50])
for r in s["a5"]: print("G5", r)
f=pd.read_csv("/tmp/r7_S13_S14/folds_S14_S14b.csv"); print("fold changes new", (f.v_old!=f.v_new).sum(), "geomA5", (f.v_old!=f.v_a5).sum())
d=pd.read_csv("/tmp/r7_S13_S14/R7_M5_S13_S14.csv", keep_default_na=False)
print(d.transition.isin(["unaffected","unchanged"]).sum(), d.transition.value_counts().head(3).to_dict(), d.reproduces.value_counts().to_dict())
print(d.interval_construction.value_counts().head(20).to_dict())
