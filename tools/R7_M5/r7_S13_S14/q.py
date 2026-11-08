
import pandas as pd, json
d=pd.read_csv("/tmp/r7_S13_S14/R7_M5_S13_S14.csv", keep_default_na=False)
aff=d[d.transition!="unaffected"]
# verdict-level transitions over affected rows (short labels)
def short(x):
    x=str(x); return x.split(" (")[0] if x.startswith(("K2","V2","choice","superata","smentita (","confermata (")) else x[:40]
t=(aff.verdict_recorded.map(short)+"->"+aff.verdict_new.map(short)).value_counts()
print("TRANS", "; ".join("%s: %d"%(k,v) for k,v in t.items()))
ts=(aff.verdict_recorded.map(short)+"->"+aff.verdict_new_stress.map(short)).value_counts()
ch=aff[aff.verdict_recorded.astype(str)!=aff.verdict_new.astype(str)]
print("NCH", len(ch), "NCH_stress", (aff.verdict_recorded.astype(str)!=aff.verdict_new_stress.astype(str)).sum())
print("CHSTEP", ch.step.value_counts().to_dict())
# real verdict changes (exclude aggregate rows whose winner unchanged)
for _,r in ch[ch.step.isin(["S1.4e","S1.4f"])].iterrows():
    print("E", r.step,"|",r.prediction_id,"|",r.archetype,"|",r.verdict_recorded,"->",r.verdict_new,"|",r.verdict_new_stress,"|",r.touches_conclusion[:60])
print("UNAFF by step", d[d.transition=="unaffected"].step.value_counts().to_dict())
print("est_unit", aff.est_unit.value_counts().to_dict(), "margin_unit", aff.margin_unit.value_counts().to_dict())
print("loao rows", (d.est_unit=="loao").sum())
print(d[d.prediction_id.str.startswith("Regola (decisione)")].notes.tolist())
