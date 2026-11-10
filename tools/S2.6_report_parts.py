"""tools/S2.6_report_parts.py - S2.6: figure e parti markdown del report (reports/S2.6.qmd). Legge solo tabelle di results/S2.6 e docs/steps/S2.6.md.
Uso (root del repo): .venv/bin/python tools/S2.6_report_parts.py -> reports/figures/S2.6/*.png, reports/S2.6_parts/*.md"""
import os, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 8})
R = "results/S2.6"; FD = "reports/figures/S2.6"; PD = "reports/S2.6_parts"; os.makedirs(FD, exist_ok=True); os.makedirs(PD, exist_ok=True)
NAME = {"A1": "intestino", "A2": "tumore del colon", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
A6 = list(NAME)
# --- fig 1: crescita con l'area
q1 = pd.read_csv(f"{R}/S2.6_q1_runs.csv")
fig, ax = plt.subplots(2, 6, figsize=(13, 4.4))
for i, A in enumerate(A6):
    g = q1[q1.archetype == A]
    ax[0, i].scatter(g.mm2, g.sim_s / g.mm2, s=12, c="tab:blue", label="simulazione")
    ax[0, i].scatter(g.mm2, g.step1_s / g.mm2, s=12, c="tab:orange", marker="^", label="livello cellulare")
    ax[0, i].axhline(300, ls="--", c="grey", lw=0.8); ax[0, i].set_title(NAME[A]); ax[0, i].set_ylim(0, 360)
    ax[1, i].scatter(g.mm2, g.sim_peak_gb, s=12, c="tab:blue"); ax[1, i].set_xlabel("area (mm²)")
    xx = np.array([0, 9]); b, a = np.polyfit(g.mm2, g.sim_peak_gb, 1); ax[1, i].plot(xx, a + b * xx, c="tab:blue", lw=0.7)
    ax[1, i].text(0.3, 0.9 * ax[1, i].get_ylim()[1], f"{b:.1f} GB/mm²", fontsize=7)
ax[0, 0].set_ylabel("s per mm² (un processo)"); ax[1, 0].set_ylabel("picco di memoria (GB)"); ax[0, 0].legend(fontsize=6, loc="center right")
fig.tight_layout(); fig.savefig(f"{FD}/fig1_growth.png", dpi=150); plt.close(fig)
# --- fig 2: J-1 per fascia
j1 = pd.read_csv(f"{R}/S2.6_q2_J1_J2.csv"); bands = ["0-8", "8-16", "16-40", "40-100", "100-250"]
fig, ax = plt.subplots(2, 3, figsize=(10, 5.2), sharey=True)
for i, A in enumerate(A6):
    a = ax.flat[i]; a.axhspan(0.98, 1.02, color="0.88")
    for k, (arm, col, dx) in enumerate([("P1000", "tab:blue", -0.15), ("P504", "tab:green", 0), ("P1000H0", "tab:red", 0.15)]):
        g = j1[(j1.archetype == A) & (j1.arm == arm) & (j1.band.isin(bands))].set_index("band").reindex(bands)
        a.errorbar(np.arange(5) + dx, g.ratio, yerr=g.two_se, fmt="o", ms=3, c=col, label={"P1000": "1 000 µm con bordo", "P504": "504 µm con bordo", "P1000H0": "1 000 µm senza bordo"}[arm])
    a.set_xticks(range(5)); a.set_xticklabels(bands); a.set_title(NAME[A]); a.axhline(1, c="k", lw=0.5)
    if i % 3 == 0: a.set_ylabel("UMI per blocco, a pezzi / pezzo solo")
    if i >= 3: a.set_xlabel("distanza dalla giuntura (µm)")
ax.flat[0].legend(fontsize=6, loc="lower right"); ax.flat[0].set_ylim(0.8, 1.15)
fig.tight_layout(); fig.savefig(f"{FD}/fig2_junctions.png", dpi=150); plt.close(fig)
# --- fig 3: J-3
j3 = pd.read_csv(f"{R}/S2.6_q2_J3.csv")
fig, ax = plt.subplots(figsize=(8, 3))
arms = [("ONE", 1000, "pezzo solo", "0.4"), ("P1000", 1000, "1 000 µm con bordo", "tab:blue"), ("P504", 504, "504 µm con bordo", "tab:green"), ("P1000H0", 1000, "senza bordo", "tab:red")]
for k, (arm, T, lab, col) in enumerate(arms):
    v = [j3[(j3.archetype == A) & (j3.arm == arm) & (j3.T_um == T)].pearson.values[0] for A in A6]
    ax.bar(np.arange(6) + (k - 1.5) * 0.2, v, 0.2, color=col, label=lab)
ax.set_xticks(range(6)); ax.set_xticklabels([NAME[A] for A in A6]); ax.axhline(0, c="k", lw=0.5)
ax.set_ylabel("correlazione degli UMI\nfra blocchi adiacenti\nattraverso la giuntura"); ax.legend(fontsize=6, ncol=4)
fig.tight_layout(); fig.savefig(f"{FD}/fig3_cross_seam.png", dpi=150); plt.close(fig)
# --- fig 4: sezioni
fig, ax = plt.subplots(1, 3, figsize=(10, 3.2), sharey=True)
order = ["tissue", "glass_16_64", "glass_64_256", "glass_256_512", "glass_512plus"]; lab = ["tessuto", "16-64", "64-256", "256-512", "≥ 512"]
for i, A in enumerate(["A1", "A4", "A6"]):
    s = pd.read_csv(f"{R}/S2.6_section_{A}.csv").set_index("band").reindex(order); a = ax[i]
    a.axhspan(0.8, 1.25, xmin=0.2, color="0.9"); a.axhspan(0.98, 1.02, xmax=0.2, color="0.8")
    y = s.ratio_pieces_one.values; lo = s.lo_pieces_one.values; hi = s.hi_pieces_one.values
    a.errorbar(range(5), y, yerr=[y - lo, hi - y], fmt="o", c="tab:purple", ms=4); a.axhline(1, c="k", lw=0.5)
    a.set_xticks(range(5)); a.set_xticklabels(lab); a.set_title(NAME[A]); a.set_xlabel("tessuto e vetro per distanza dal tessuto (µm)")
ax[0].set_ylabel("UMI per blocco, a pezzi / pezzo solo"); ax[0].set_ylim(0, 2.5)
fig.tight_layout(); fig.savefig(f"{FD}/fig4_sections.png", dpi=150); plt.close(fig)
# --- parti dal file dello step
st = open("docs/steps/S2.6.md").read()
def section(h1, nxt):
    i = st.index(h1); j = st.index(nxt, i + len(h1)); return st[i:j].strip() + "\n"
ib = section("## In breve:", "# Dubbi rinviati"); open(f"{PD}/in_breve.md", "w").write(ib.replace("## In breve:", "## In breve:", 1))
open(f"{PD}/dubbi.md", "w").write(section("# Dubbi rinviati", "# Verifiche di Luca").replace("# Dubbi rinviati\n", ""))
for name, a, b in [("esiti", "# Esiti", "# Correzioni (dopo il calcolo"), ("correzioni", "# Correzioni (dopo il calcolo", "# Analisi post hoc"),
                   ("posthoc", "# Analisi post hoc", "# Revisione avversariale"), ("revisione", "# Revisione avversariale", "# Deviazioni dichiarate"),
                   ("deviazioni", "# Deviazioni dichiarate", "# In breve"), ("attese", "# Attese pre-registrate", "# Decisioni chieste a Luca"),
                   ("esplorazioni", "# Esplorazioni dichiarate", "# Attese pre-registrate")]:
    t = section(a, b); t = re.sub(r"^# ", "## ", t, flags=re.M) if name in ("esiti",) else t
    t = "\n".join(l[1:] if l.startswith("#") else l for l in t.splitlines()[1:]) + "\n"
    t = re.sub(r"\]\(/decisions/", "](../docs/decisions/", t); t = re.sub(r"\]\(/steps/", "](../docs/steps/", t)
    open(f"{PD}/{name}.md", "w").write(t)
open(f"{PD}/in_breve.md", "w").write(re.sub(r"^## In breve:", "## In breve:", ib))
v = pd.read_csv(f"{R}/S2.6_verdicts.csv")
v["tessuto"] = v.archetype.map(NAME)
rows = ["| id | tessuto | quantità | stima | intervallo (± 2 SE o percentile) | intervallo di rilevanza | esito |", "|---|---|---|---|---|---|---|"]
for _, r in v.iterrows():
    f = lambda x: "" if pd.isna(x) else f"{x:.4g}"
    rows.append(f"| {r.id} | {r.tessuto} | {r.quantity} | {f(r.estimate)} | [{f(r.low)}, {f(r.high)}] | [{f(r.rel_lo)}, {f(r.rel_hi)}] | {r.outcome} |")
open(f"{PD}/verdicts.md", "w").write("\n".join(rows) + "\n")
print("parti e figure scritte")
