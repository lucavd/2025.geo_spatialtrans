#!/usr/bin/env python3
"""tools/S2.8_report_parts.py - S2.8: estrae le sezioni di docs/steps/S2.8.md in reports/S2.8_parts/ (il testo esiste una volta sola,
nel file dello step) e disegna le figure del report in reports/figures/S2.8/ dalle tabelle di results/S2.8/. Uso (root del repo):
.venv/bin/python tools/S2.8_report_parts.py"""
import os, re
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.makedirs("reports/S2.8_parts", exist_ok=True); os.makedirs("reports/figures/S2.8", exist_ok=True)
s = open("docs/steps/S2.8.md").read(); body = s.split("\n---\n", 1)[1]
secs = {}; cur = None
for l in body.splitlines(keepends=True):
    if l.startswith("# "): cur = l[2:].strip(); secs[cur] = ""; continue
    if cur: secs[cur] += l
def fix_links(t): return re.sub(r"\]\(/(decisions|steps)/", r"](../docs/\1/", t)
def demote(t): return re.sub(r"^(#+) ", lambda m: "#" + m.group(1) + " ", t, flags=re.M)
want = {"In breve": "in_breve", "Dubbi rinviati": "dubbi", "Perimetro": "perimetro", "Regola dell'ingresso (R-1, scritta prima della misura)": "regola",
        "Attese pre-registrate": "attese", "Correzioni e aggiunte dopo la pre-registrazione (dichiarate, prima dei verdetti)": "correzioni",
        "Esiti": "esiti", "Revisione avversariale": "revisione", "Deviazioni dichiarate": "deviazioni", "Verifiche di Luca": "verifiche",
        "Stime di preparazione (lette dai file di S2.7 prima della pre-registrazione; non sono risultati dello step)": "stime"}
for k, f in want.items():
    if k not in secs: raise SystemExit(f"sezione mancante: {k}")
    open(f"reports/S2.8_parts/{f}.md", "w").write(demote(fix_links(secs[k])))
S = pd.read_csv("results/S2.8/S2.8_summary.csv")
TIS = {"A1": "intestino", "A2": "tumore del colon", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
CF = [("ref", "motore di S2.5", "#7f7f7f"), ("e504", "a pezzi, 504 µm", "#1f77b4"), ("e1008", "a pezzi, 1 008 µm", "#9ecae1"), ("einf", "a pezzi, pezzo solo", "#2ca02c")]
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.size": 9})
for col, err, lim, yl, fn in [("time_mean", "time_es2", 300, "tempo di orologio (s)", "fig1_time.png"), ("peak_max", "peak_es2", 8, "picco di memoria (GB)", "fig2_memory.png")]:
    fig, axs = plt.subplots(2, 3, figsize=(9, 5.2), sharey=True)
    for ax, (A, name) in zip(axs.ravel(), TIS.items()):
        for i, (cf, lab, colr) in enumerate(CF):
            r = S[(S.conf == cf) & (S.A == A)].iloc[0]
            ax.bar(i, r[col], yerr=r[err], color=colr, label=lab, capsize=2)
        ax.axhline(lim, ls="--", color="k", lw=0.8); ax.set_title(name); ax.set_xticks([])
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axs[:, 0]: ax.set_ylabel(yl)
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(f"reports/figures/S2.8/{fn}", dpi=150); plt.close(fig)
print("ok", sorted(os.listdir("reports/S2.8_parts")))
