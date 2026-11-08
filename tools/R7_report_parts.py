"""tools/R7_report_parts.py — R7: tabelle riassuntive per il report e per il riassunto «In breve». -> results/R7/report/*.csv"""
import os, numpy as np, pandas as pd
OUT = "results/R7/report"; os.makedirs(OUT, exist_ok=True)
T = pd.concat([pd.read_csv(f"results/R7/M5/{f}.csv", low_memory=False).assign(table=f) for f in ["R7_M5_S13_S14", "R7_M5_S15_S16", "R7_M5_R_S2"]], ignore_index=True)
THREE = {"confermata", "indeterminata", "smentita"}
def kind(t):
    if t in ("unaffected", "unchanged"): return t
    a, _, b = str(t).partition("->")
    return "three" if a in THREE and b in THREE else "rule/other"
T["kind"] = T.transition.map(kind)
rows = []
for s, g in T.groupby("step"):
    th = g[g.kind == "three"]
    rows.append(dict(step=s, rows=len(g), affected=int((g.kind != "unaffected").sum()), changed_three=len(th), rule_or_other_changes=int((g.kind == "rule/other").sum()),
                     **{k: int((th.transition == k).sum()) for k in ["confermata->indeterminata", "smentita->indeterminata", "smentita->confermata", "indeterminata->confermata", "indeterminata->smentita", "confermata->smentita"]}))
S = pd.DataFrame(rows); tot = S.drop(columns="step").sum(); tot["step"] = "TOTAL"; S = pd.concat([S, tot.to_frame().T], ignore_index=True)
S.to_csv(f"{OUT}/R7_transitions_by_step.csv", index=False); print(S.to_string(index=False))
R = T[T.kind == "rule/other"][["table", "step", "prediction_id", "archetype", "quantity", "verdict_recorded", "verdict_new", "verdict_new_stress", "verdict_geom_A5", "transition", "touches_conclusion"]]
R.to_csv(f"{OUT}/R7_rule_changes.csv", index=False); print(R[["step", "prediction_id", "transition", "verdict_new_stress", "verdict_geom_A5"]].to_string(index=False))
print("reproduces:", T.reproduces.astype(str).str.upper().value_counts().to_dict())
M = pd.read_csv("results/R7/Mc/R7_Mc_refs.csv")
M[M.is_S25_target == True][["id", "archetype", "metric", "est", "lo", "hi", "kappa", "lo_new", "hi_new", "width_ratio"]].to_csv(f"{OUT}/R7_S25_targets.csv", index=False)
print(M.groupby(["interval_construction"]).size().to_dict(), "high CV:", int((M.cv_between_areas > 0.45).sum()))
