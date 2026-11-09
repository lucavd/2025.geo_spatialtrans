"""tools/S2.5_report_parts.py - S2.5: figure e parti markdown del report (reports/S2.5.qmd). Legge solo tabelle e file dello step.
Uso: .venv/bin/python tools/S2.5_report_parts.py -> reports/figures/S2.5/*.png, reports/S2.5_parts/*.md"""
import os, re, json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 9})
FD = "reports/figures/S2.5"; PD = "reports/S2.5_parts"; os.makedirs(FD, exist_ok=True); os.makedirs(PD, exist_ok=True)
TIS = {"A1": "intestine", "A2": "colon tumour", "A3": "tumour stroma", "A4": "lymph node", "A5": "brain", "A6": "heart"}
TIS_IT = {"A1": "intestino", "A2": "tumore del colon", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
PANEL = "/mnt/micron/geo_spatialtrans/S2.5/panel"
step = open("docs/steps/S2.5.md").read()
def section(title):
    m = re.search(r"^# %s\n(.*?)(?=^# )" % re.escape(title), step, re.S | re.M); return m.group(1).strip() if m else ""
for t, f in (("In breve", "in_breve.md"), ("Dubbi rinviati", "dubbi.md"), ("Esiti", "esiti.md"),
             ("Analisi post hoc (dichiarate dopo i risultati; non cambiano gli esiti)", "posthoc.md"), ("Revisione avversariale", "revisione.md"),
             ("Deviazioni dichiarate", "deviazioni.md"), ("Correzioni (dopo il calcolo, prima del verdetto)", "correzioni.md"),
             ("Esplorazioni dichiarate (prima della pre-registrazione)", "esplorazioni.md"), ("Attese pre-registrate", "attese.md"),
             ("Decisioni chieste a Luca prima del calcolo", "decisioni.md")):
    txt = section(t); assert txt, t
    txt = re.sub(r"\]\(/(decisions|steps)/", r"](../docs/\1/", txt)
    txt = re.sub(r"^(#+) ", lambda m: "#" * (len(m.group(1)) + 1) + " ", txt, flags=re.M)
    open(f"{PD}/{f}", "w").write(txt + "\n")
# tabella di tutti i giudizi: main e main_r2 affiancati
VM = pd.read_csv("results/S2.5/S2.5_verdicts_main.csv"); VR = pd.read_csv("results/S2.5/S2.5_verdicts_main_r2.csv")
key = ["id", "archetype", "metric"]; V = VM.merge(VR, on=key, suffixes=("_m", "_r"))
def f(x):
    if pd.isna(x): return ""
    x = float(x); return f"{x:.4g}" if abs(x) < 1000 else f"{x:.0f}"
lines = ["| id | tessuto | misura | vero (stima) | intervallo del vero (R7) | simulato, fattore di S2.2 | esito | simulato, fattore di R-2 | esito |", "|---|---|---|---|---|---|---|---|---|"]
for _, r in V.iterrows():
    lines.append(f"| {r.id} | {TIS_IT[r.archetype]} | {r.metric} | {f(r.real_est_m)} | [{f(r.rel_lo_m)}, {f(r.rel_hi_m)}] | {f(r.sim_est_m)} [{f(r.sim_lo_m)}, {f(r.sim_hi_m)}] | {r.esito_m} | {f(r.sim_est_r)} [{f(r.sim_lo_r)}, {f(r.sim_hi_r)}] | {r.esito_r} |")
open(f"{PD}/verdicts.md", "w").write("\n".join(lines) + "\n")
VD = pd.read_csv("results/S2.5/S2.5_verdicts_depth_off.csv")
lines = ["| tessuto | intervallo del vero | simulato senza fattore di profondità | esito |", "|---|---|---|---|"]
for _, r in VD[VD.id == "P-1"].iterrows():
    lines.append(f"| {TIS_IT[r.archetype]} | [{f(r.rel_lo)}, {f(r.rel_hi)}] | {f(r.sim_est)} [{f(r.sim_lo)}, {f(r.sim_hi)}] | {r.esito} |")
open(f"{PD}/cp1.md", "w").write("\n".join(lines) + "\n")
# figura 1: simulato / vero per misura e tessuto, prima e dopo R-2
MET = [("P-1", "mean_umi2", "UMI per 2 um bin"), ("P-2", "frac_empty2", "empty 2 um bins"), ("P-3", "mean_genes2", "genes per 2 um bin"),
       ("P-4", "mean_umi8", "UMI per 8 um block"), ("P-6", "mean_genes8", "genes per 8 um block"), ("P-7", "star", "nucleus / ring ratio")]
fig, axs = plt.subplots(2, 3, figsize=(11, 6), sharex=True)
for ax, (pid, m, lab) in zip(axs.flat, MET):
    for k, A in enumerate(TIS):
        a = VM[(VM.id == pid) & (VM.archetype == A)].iloc[0]; b = VR[(VR.id == pid) & (VR.archetype == A)].iloc[0]; e = a.real_est
        lo = max(a.rel_lo / e, 1e-3); hi = a.rel_hi / e
        ax.add_patch(plt.Rectangle((k - 0.35, lo), 0.7, hi - lo, color="0.85", lw=0))
        for dx, row, c, mk, lb in ((-0.12, a, "tab:red", "x", "simulated, S2.2 depth factor"), (0.12, b, "tab:blue", "o", "simulated, R-2 depth factor")):
            ax.errorbar(k + dx, row.sim_est / e, yerr=[[(row.sim_est - row.sim_lo) / e], [(row.sim_hi - row.sim_est) / e]], fmt=mk, c=c, ms=4, label=lb if k == 0 else None)
    ax.axhline(1, c="0.4", lw=0.7); ax.set_yscale("log"); ax.set_title(f"{pid}: {lab}")
    ax.set_xticks(range(6)); ax.set_xticklabels(list(TIS.values()), rotation=40, ha="right")
axs[0, 0].set_ylabel("simulated / real (log)"); axs[1, 0].set_ylabel("simulated / real (log)")
axs[0, 0].legend(fontsize=7, frameon=False, loc="lower left")
fig.text(0.5, 0.005, "grey bar: corrected interval of the real value (R7), divided by the real estimate", ha="center", fontsize=8)
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(f"{FD}/fig1_ratios.png", dpi=150); plt.close(fig)
# figura 2: orientamento, aree vere e simulate
RA = pd.read_csv("results/S2.5/S2.5_real_areas.csv"); SA = pd.read_csv("results/S2.5/S2.5_sim_main_r2.csv")
fig, ax = plt.subplots(figsize=(7, 3.4))
for k, A in enumerate(TIS):
    r = RA[RA.archetype == A].reg_aniso8; s = SA[SA.archetype == A].reg_aniso8
    b = VR[(VR.id == "P-9") & (VR.archetype == A)].iloc[0]
    ax.add_patch(plt.Rectangle((k - 0.35, b.rel_lo), 0.7, b.rel_hi - b.rel_lo, color="0.85", lw=0))
    ax.scatter(np.full(len(r), k - 0.12), r, s=16, c="k", label="real area" if k == 0 else None)
    ax.scatter(np.full(len(s), k + 0.12), s, s=16, c="tab:blue", marker="^", label="simulated area (R-2 panel, 2 seeds)" if k == 0 else None)
ax.set_xticks(range(6)); ax.set_xticklabels(list(TIS.values()), rotation=30, ha="right"); ax.set_ylabel("orientation index of 8 um counts")
ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig(f"{FD}/fig2_orientation.png", dpi=150); plt.close(fig)
# figura 3: tempo e memoria contro i trascritti generati
rows = []
for v in ("main", "main_r2", "depth_off"):
    for fj in glob.glob(f"{PANEL}/{v}/*.json"):
        d = json.load(open(fj)); rows.append(dict(src=f"panel {v} (4e6, up to 6 processes)", A=d["archetype"], n=d["n_generated"], t=d["elapsed_sim_s"], mb=d["peak_mb"]))
for fj in sorted(glob.glob("results/S2.5/PH1/*.json")):
    d = json.load(open(fj)); rows.append(dict(src=f"PH-1, A2_r4, {d['chunk_transcripts']:.0e} per block", A="A2", n=d["n_generated"], t=d["elapsed_s"], mb=d["peak_mb"]))
for fj in sorted(glob.glob("results/S2.5/PH2/*.json")):
    d = json.load(open(fj)); rows.append(dict(src="PH-2, colon tumour, 2e6, idle machine", A="A2", n=d["n_generated"], t=d["elapsed_s"], mb=d["peak_mb"]))
C = pd.DataFrame(rows); C.to_csv(f"{PD}/cost_points.csv", index=False)
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
sty = {"panel main (4e6, up to 6 processes)": ("0.55", "o", 10), "panel main_r2 (4e6, up to 6 processes)": ("0.75", "s", 10), "panel depth_off (4e6, up to 6 processes)": ("0.35", "d", 10)}
others = [s for s in C.src.unique() if s not in sty]; pal = ["tab:red", "tab:orange", "tab:purple", "tab:blue"]
for i, s in enumerate(others): sty[s] = (pal[i % 4], "*" if s.startswith("PH-1") else "^", 40)
for ax, q, lim, lab in ((axs[0], "t", 300, "seconds per 1 mm2 area"), (axs[1], "mb", 8192, "peak resident memory (MB)")):
    for s, g in C.groupby("src"):
        c, mk, sz = sty[s]; ax.scatter(g.n / 1e6, g[q], c=c, marker=mk, s=sz, label=s)
    ax.axhline(lim, c="tab:red", ls="--", lw=1); ax.set_xlabel("transcripts generated (millions)"); ax.set_ylabel(lab)
axs[1].legend(fontsize=6.5, frameon=False, loc="lower right"); fig.tight_layout(); fig.savefig(f"{FD}/fig3_cost.png", dpi=150); plt.close(fig)
print("ok", len(C), "punti di costo")
