#!/usr/bin/env python
"""results/S1.4c/review/recalc/r07_margins.py — revisione avversariale S1.4c (recalc): fragilita' della decisione.
Dalle tabelle di r02 (recalc_fold_verdicts_<tag>.csv): celle archetipo x metrica in cui K-2 e K-1 differiscono; verdetti
d'archetipo PASS che dipendono da un solo fold con d entro il 10% della mediana del metro; scelta se le metriche S
non fossero declassate. Uso: .venv/bin/python r07_margins.py [tag]
"""
import sys, numpy as np, pandas as pd
O = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/"; tag = sys.argv[1] if len(sys.argv) > 1 else "obs"
fv = pd.read_csv(O + f"recalc_fold_verdicts_{tag}.csv"); vt = pd.read_csv(O + f"recalc_verdicts_{tag}.csv")
panel = ["B1", "B2", "B7", "B8"]
w = vt[vt.metric.isin(panel)].pivot_table(index=["archetype", "metric"], columns="variant", values="verdict", aggfunc="first")
print("Celle del pannello (archetipo x metrica):"); print(w.to_string())
print("\nK-2 != K-1:"); print(w[w["K-2"] != w["K-1"]].to_string())
fv["ratio"] = fv.dist_median / fv.yard_median
fr = []
for (A, V, m), z in fv[fv.metric.isin(panel)].groupby(["archetype", "variant", "metric"]):
    npass = (z.verdict == "PASS").sum()
    if npass == 4:
        near = z[(z.verdict == "PASS") & (z.ratio > 0.9)]
        if len(near): fr.append(dict(archetype=A, variant=V, metric=m, n_pass=npass, fold=",".join(map(str, near.fold)), ratio=",".join(f"{r:.3f}" for r in near.ratio)))
    if npass == 3:
        near = z[(z.verdict != "PASS") & (z.ratio < 1.1)]
        if len(near): fr.append(dict(archetype=A, variant=V, metric=m, n_pass=npass, fold=",".join(map(str, near.fold)), ratio=",".join(f"{r:.3f}" for r in near.ratio)))
print("\nVerdetti d'archetipo a un fold dal cambio (|d/mediana - 1| <= 0.1):"); print(pd.DataFrame(fr).to_string())
print("\nPASS per variante (pannello):", vt[vt.metric.isin(panel) & (vt.verdict == "PASS")].groupby("variant").size().to_dict())
print("PASS per variante (7 metriche):", vt[vt.verdict == "PASS"].groupby("variant").size().to_dict())
