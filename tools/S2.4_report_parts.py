"""tools/S2.4_report_parts.py — S2.4: figure e parti markdown del report (reports/S2.4.qmd). Legge solo tabelle e file dello step.
Uso: .venv/bin/python tools/S2.4_report_parts.py -> reports/figures/S2.4/*.png, reports/S2.4_parts/*.md"""
import os, re, json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 9})
FD = "reports/figures/S2.4"; PD = "reports/S2.4_parts"; os.makedirs(FD, exist_ok=True); os.makedirs(PD, exist_ok=True)
TIS = {"A1": "intestine", "A2": "colon tumour", "A3": "tumour stroma", "A4": "lymph node", "A5": "brain", "A6": "heart"}
SEC = {"A1": "intestine", "A23": "colon tumour + stroma", "A4": "lymph node", "A5": "brain", "A6": "heart"}
step = open("docs/steps/S2.4.md").read()
def section(title):
    m = re.search(r"^# %s\n(.*?)(?=^# )" % re.escape(title), step, re.S | re.M); return m.group(1).strip() if m else ""
for t, f in (("In breve", "in_breve.md"), ("Dubbi rinviati", "dubbi.md"), ("Esiti", "esiti.md"), ("Analisi post hoc (dichiarate dopo i risultati; non cambiano gli esiti)", "posthoc.md"),
             ("Revisione avversariale", "revisione.md"), ("Deviazioni dichiarate", "deviazioni.md"), ("Correzioni prima del calcolo", "correzioni.md"), ("Esplorazioni dichiarate (prima della pre-registrazione)", "esplorazioni.md")):
    txt = section(t)
    txt = re.sub(r"\]\(/(decisions|steps)/", r"](../docs/\1/", txt)
    txt = re.sub(r"^## ", "### ", txt, flags=re.M)
    open(f"{PD}/{f}", "w").write(txt + "\n")
# tabella di tutte le previsioni
V = pd.read_csv("results/S2.4/S2.4_verdicts.csv")
def fmt(x, d=3): return "" if pd.isna(x) else f"{x:.{d}f}"
lines = ["| id | tessuto | quantità | stima | ± 2 SE | intervallo | n | esito |", "|---|---|---|---|---|---|---|---|"]
for _, r in V.iterrows():
    two = fmt(2 * r.se) if pd.isna(r.get("ci_lo", np.nan)) else f"[{r.ci_lo:.3f}, {r.ci_hi:.3f}]"
    lines.append(f"| {r.id} | {r.tessuto} | {r.quantita} | {fmt(r.stima, 4 if abs(r.stima) < 1 else 2) if not pd.isna(r.stima) else ''} | {two} | [{r.lo:g}, {r.hi:g}] | {int(r.n)} | {r.esito} |")
open(f"{PD}/verdicts.md", "w").write("\n".join(lines) + "\n")
# figura 1: strisce
R = pd.read_csv("results/S2.4/S2.4_E1_real.csv"); S = pd.DataFrame([json.load(open(f)) for f in glob.glob("/mnt/micron/geo_spatialtrans/S2.4/panel/stats_main/sim/*.json")])
S = S.groupby(["tag", "archetype"]).mean(numeric_only=True).reset_index()
fig, axs = plt.subplots(2, 3, figsize=(10, 5.5), sharex=True)
for i, ax_ in enumerate(("rows", "cols")):
    for j, (q, lab) in enumerate((("odd_over_even", "odd/even ratio"), ("sd_repro_hp", "SD of fine part (log)"), ("lag1_repro_hp", "lag-1 of fine part"))):
        ax = axs[i, j]; col = f"stripes_{ax_}_{q}"
        for k, A in enumerate(TIS):
            r = R[R.archetype == A][col]; s = S[S.archetype == A][col]
            ax.scatter(np.full(len(r), k - 0.12), r, s=14, c="k", marker="o", label="real area" if k == 0 else None)
            ax.scatter(np.full(len(s), k + 0.12), s, s=14, c="tab:red", marker="^", label="simulated area (mean of 2 seeds)" if k == 0 else None)
        ax.set_title(f"{lab}, {ax_}"); ax.set_xticks(range(6)); ax.set_xticklabels(list(TIS.values()), rotation=40, ha="right")
axs[0, 0].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{FD}/fig1_stripes.png", dpi=150); plt.close(fig)
# figura 2: vetro, prima e dopo R-3
G = pd.read_csv("results/S2.4/S2.4_P5_glass.csv"); bands = ["16-32", "32-64", "64-128", "128-256", "256-512", "512-100000"]
fig, axs = plt.subplots(1, 5, figsize=(13, 3), sharey=True)
for ax, k in zip(axs, SEC):
    g = G[G.key == k].set_index("band").reindex(bands)
    e = json.load(open(f"/mnt/micron/geo_spatialtrans/S2.4/explore/glass/{k}_pad1500.json"))
    eb = {b["band"]: b["ratio_with_tail"] for b in e["bands"]}
    x = np.arange(len(bands))
    ax.axhspan(0.8, 1.25, color="0.9"); ax.axhline(1, c="0.5", lw=0.7)
    ax.errorbar(x, g.ratio, yerr=[g.ratio - g.lo, g.hi - g.ratio], fmt="o", c="tab:blue", ms=4, label="background + tail fitted together (R-3)")
    ax.plot(x, [eb.get(b.replace("512-100000", "512-10000"), np.nan) for b in bands], "x", c="tab:red", label="far-glass background + S2.3 tail (E-3)")
    ax.set_xticks(x); ax.set_xticklabels(["16-32", "32-64", "64-128", "128-256", "256-512", ">=512"], rotation=45, fontsize=7); ax.set_title(SEC[k])
axs[0].set_ylabel("glass observed / predicted"); axs[0].legend(fontsize=6.5, frameon=False, loc="lower left")
fig.tight_layout(); fig.savefig(f"{FD}/fig2_glass.png", dpi=150); plt.close(fig)
# figura 3: lumi profondi, simulato / vero, S2.3 (senza coda) e S2.4
E = pd.read_csv("results/S2.4/S2.4_E12_summary.csv").set_index("A"); P4 = V[V.id == "P-4"]
fig, axs = plt.subplots(1, 3, figsize=(10, 3), sharey=True)
for ax, A in zip(axs, ("A1", "A2", "A6")):
    bl = ["8_16", "16_32", "32_64"]; x = np.arange(3)
    s23 = [E.loc[A, f"simS23_band_{b}"] / E.loc[A, f"real_band_{b}"] for b in bl]
    p = P4[P4.A == A].set_index("quantita").reindex([f"fascia {b.replace('_', '-')} um" for b in bl])
    ax.axhspan(0.8, 1.25, color="0.9"); ax.axhline(1, c="0.5", lw=0.7)
    ax.plot(x - 0.1, s23, "x", c="tab:red", label="S2.3 panel (no second tail)")
    ax.errorbar(x + 0.1, p.stima, yerr=2 * p.se, fmt="o", c="tab:blue", ms=4, label="S2.4 panel (± 2 SE)")
    ax.set_xticks(x); ax.set_xticklabels(["8-16 um", "16-32 um", "32-64 um"]); ax.set_title(TIS[A]); ax.set_yscale("log")
axs[0].set_ylabel("lumen band, simulated / real"); axs[0].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{FD}/fig3_lumens.png", dpi=150); plt.close(fig)
# figura 4: tempo e memoria
J = pd.DataFrame([json.load(open(f)) for f in glob.glob("/mnt/micron/geo_spatialtrans/S2.4/panel/main/*_k?.json")])
fig, axs = plt.subplots(1, 2, figsize=(9, 3))
for ax, q, lim, lab in ((axs[0], "elapsed_total_s", 300, "seconds per 1 mm2 area"), (axs[1], "peak_mb", 8192, "peak resident memory (MB)")):
    for k, A in enumerate(TIS): v = J[J.archetype == A][q]; ax.scatter(np.full(len(v), k), v, s=12, c="k")
    ax.axhline(lim, c="tab:red", ls="--", lw=1, label="design budget"); ax.set_xticks(range(6)); ax.set_xticklabels(list(TIS.values()), rotation=40, ha="right"); ax.set_ylabel(lab)
axs[0].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{FD}/fig4_cost.png", dpi=150); plt.close(fig)
print("ok")
