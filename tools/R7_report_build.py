"""tools/R7_report_build.py — R7: parti markdown del report (dal file dello step) e figure. -> reports/R7_parts/*.md, reports/figures/R7/*.png"""
import os, re, shutil, glob
import pandas as pd
def md(df):
    df = df.astype(str)
    return "\n".join(["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)] + ["| " + " | ".join(r) + " |" for r in df.values.tolist()])
PD = "reports/R7_parts"; FD = "reports/figures/R7"; os.makedirs(PD, exist_ok=True); os.makedirs(FD, exist_ok=True)
step = open("docs/steps/R7.md").read()
def section(title):
    m = re.search(r"^# %s\n(.*?)(?=^# )" % re.escape(title), step, re.S | re.M); return m.group(1).strip() if m else ""
for t, f in (("In breve", "in_breve.md"), ("Dubbi rinviati", "dubbi.md"), ("Esiti", "esiti.md"), ("Analisi post hoc (dichiarate dopo i risultati; non cambiano gli esiti)", "posthoc.md"),
             ("Revisione avversariale", "revisione.md"), ("Deviazioni dichiarate", "deviazioni.md"), ("Correzioni dopo la revisione avversariale", "correzioni.md"),
             ("Esplorazioni dichiarate (prima della pre-registrazione)", "esplorazioni.md"), ("Attese pre-registrate", "attese.md"), ("Perimetro", "perimetro.md")):
    txt = section(t); assert txt, t
    txt = re.sub(r"\]\(/(decisions|steps)/", r"](../docs/\1/", txt); txt = re.sub(r"^## ", "### ", txt, flags=re.M)
    txt = re.sub(r"\[\^[a-z0-9]+\]", "", txt)
    open(f"{PD}/{f}", "w").write(txt + "\n")
for f in glob.glob("results/R7/figures/*.png"): shutil.copy(f, FD)
T = pd.read_csv("results/R7/report/R7_transitions_by_step.csv")
open(f"{PD}/transitions.md", "w").write(md(T) + "\n")
S = pd.read_csv("results/R7/report/R7_S25_targets.csv")
S["vecchio"] = [f"{e:.4g} [{l:.4g}, {h:.4g}]" for e, l, h in zip(S.est, S.lo, S.hi)]; S["corretto"] = [f"[{l:.4g}, {h:.4g}]" for l, h in zip(S.lo_new, S.hi_new)]
open(f"{PD}/s25.md", "w").write(md(S[["id", "archetype", "metric", "vecchio", "corretto", "kappa"]].round(3)) + "\n")
R = pd.read_csv("results/R7/report/R7_rule_changes.csv")
open(f"{PD}/rules.md", "w").write(md(R[["step", "prediction_id", "transition", "verdict_new_stress", "verdict_geom_A5"]].fillna("")) + "\n")
print("ok")
