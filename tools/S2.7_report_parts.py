"""tools/S2.7_report_parts.py - S2.7: figure e parti markdown del report (reports/S2.7.qmd). Legge solo tabelle di results/S2.7 e docs/steps/S2.7.md.
Uso (root del repo): .venv/bin/python tools/S2.7_report_parts.py -> reports/figures/S2.7/*.png, reports/S2.7_parts/*.md"""
import os, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 8})
R = "results/S2.7"; FD = "reports/figures/S2.7"; PD = "reports/S2.7_parts"; os.makedirs(FD, exist_ok=True); os.makedirs(PD, exist_ok=True)
NAME = {"A1": "intestino", "A2": "tumore del colon", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
# fig 1: P-1 e CP-1
P = pd.read_csv(f"{R}/S2.7_P1.csv"); PH = pd.read_csv(f"{R}/PH1/S2.7_P1.csv"); PH = PH[PH.n_new > 0]
mets = [("mean_umi2", "UMI per quadratino 2 µm"), ("mean_genes2", "geni per quadratino 2 µm"), ("mean_genes8", "geni per blocco 8 µm"), ("mean_umi2_thinned_0.95", "controprova: cattura -5 %")]
fig, ax = plt.subplots(1, 4, figsize=(12, 3.2), sharey=True)
for k, (m, lab) in enumerate(mets):
    a = ax[k]; a.axvspan(0.98, 1.02, color="0.88"); a.axvline(1, c="k", lw=0.5)
    for i, A in enumerate(NAME):
        r = P[(P.archetype == A) & (P.metric == m)].iloc[0]
        a.errorbar(r.ratio, i, xerr=[[r.ratio - r.lo], [r.hi - r.ratio]], fmt="o", ms=3, c="tab:blue" if m != "mean_umi2_thinned_0.95" else "tab:red")
        if A == "A5" and m != "mean_umi2_thinned_0.95":
            q = PH[PH.metric == m].iloc[0]; a.errorbar(q.ratio, i + 0.3, xerr=[[q.ratio - q.lo], [q.hi - q.ratio]], fmt="s", ms=3, c="tab:green")
    a.set_title(lab); a.set_xlim(0.92, 1.05); a.set_xlabel("a pezzi / motore di S2.5")
ax[0].set_yticks(range(6)); ax[0].set_yticklabels(NAME.values()); ax[0].invert_yaxis()
fig.tight_layout(); fig.savefig(f"{FD}/fig1_law.png", dpi=150); plt.close(fig)
# fig 2: sezioni
order = ["tissue", "glass_16_64", "glass_64_256", "glass_256_512", "glass_512plus"]; lab = ["tessuto", "16-64", "64-256", "256-512", "≥ 512"]
fig, ax = plt.subplots(1, 3, figsize=(10, 3.2), sharey=True)
for i, A in enumerate(["A1", "A4", "A6"]):
    s = pd.read_csv(f"{R}/S2.7_section_{A}.csv").set_index("band").reindex(order); a = ax[i]
    a.axhspan(0.8, 1.25, xmin=0.2, color="0.9"); a.axhspan(0.98, 1.02, xmax=0.2, color="0.8")
    y = s.ratio_pieces_one.values; lo = s.lo_pieces_one.values; hi = s.hi_pieces_one.values
    a.errorbar(range(5), y, yerr=[y - lo, np.minimum(hi, 2.0) - y], fmt="o", c="tab:purple", ms=4); a.axhline(1, c="k", lw=0.5)
    a.set_xticks(range(5)); a.set_xticklabels(lab); a.set_title(NAME[A]); a.set_xlabel("tessuto e vetro per distanza dal tessuto (µm)")
ax[0].set_ylabel("UMI per blocco, pezzi (S2.7) / pezzo solo (S2.6)"); ax[0].set_ylim(0.6, 2.0)
fig.tight_layout(); fig.savefig(f"{FD}/fig2_sections.png", dpi=150); plt.close(fig)
# fig 3: risorse
T = pd.read_csv(f"{R}/S2.7_P3_P4.csv"); T["A"] = T.section.str.extract(r"sec_(A\d)")
fig, ax = plt.subplots(1, 2, figsize=(9, 3))
for j, (col, lim, yl) in enumerate([("piece_peak_gb", 8, "picco di un pezzo (GB)"), ("s_per_mm2_tissue", 300, "s di processo per mm² di tessuto")]):
    a = ax[j]; x = np.arange(5); As = ["A1", "A2", "A4", "A5", "A6"]
    v1 = [T[(T.A == A) & (T.piece_um == 1008)][col].values[0] for A in As]
    v2 = [T[(T.A == A) & (T.piece_um == 504)][col].values[0] if ((T.A == A) & (T.piece_um == 504)).any() else np.nan for A in As]
    a.bar(x - 0.2, v1, 0.4, label="pezzi di 1 008 µm"); a.bar(x + 0.2, v2, 0.4, label="pezzi di 504 µm (solo tumore e cervello)")
    a.axhline(lim, ls="--", c="grey"); a.set_xticks(x); a.set_xticklabels([NAME[A] for A in As], rotation=20); a.set_ylabel(yl)
ax[0].legend(fontsize=6)
fig.tight_layout(); fig.savefig(f"{FD}/fig3_resources.png", dpi=150); plt.close(fig)
# parti dal file dello step
st = open("docs/steps/S2.7.md", encoding="utf-8").read()
def section(a, b):
    i = st.index(a); j = st.index(b, i + len(a)); return st[i:j].strip() + "\n"
def clean(t, drop_head=True):
    L = t.splitlines()[1:] if drop_head else t.splitlines()
    t = "\n".join(L) + "\n"
    t = re.sub(r"\]\(/decisions/", "](../docs/decisions/", t); t = re.sub(r"\]\(/steps/", "](../docs/steps/", t); return t
open(f"{PD}/in_breve.md", "w").write(clean(section("## In breve:", "# Dubbi rinviati"), False))
for name, a, b in [("dubbi", "# Dubbi rinviati", "# Verifiche di Luca"), ("esiti", "# Esiti", "# Analisi post hoc"), ("posthoc", "# Analisi post hoc", "# Revisione avversariale"),
                   ("revisione", "# Revisione avversariale", "# Deviazioni dichiarate"), ("deviazioni", "# Deviazioni dichiarate", "# In breve"),
                   ("attese", "# Attese pre-registrate", "# Correzioni prima dei calcoli"), ("correzioni", "# Correzioni prima dei calcoli", "# Esiti"), ("perimetro", "# Perimetro", "# Stime di preparazione")]:
    open(f"{PD}/{name}.md", "w").write(clean(section(a, b)).replace("\n## ", "\n### "))
rows = ["| id | tessuto | quantità | rapporto | ± 2 SE | rilevanza | esito |", "|---|---|---|---|---|---|---|"]
for _, r in P.iterrows():
    rows.append(f"| {r.id} | {NAME[r.archetype]} | {r.metric} | {r.ratio:.4f} | [{r.lo:.4f}, {r.hi:.4f}] | [{r.rel_lo:.3g}, {r.rel_hi:.3g}] | {r.esito} |")
open(f"{PD}/p1.md", "w").write("\n".join(rows) + "\n")
rows = ["| tessuto | fascia | blocchi | rapporto | intervallo | esito |", "|---|---|---|---|---|---|"]
for A in ["A1", "A4", "A6"]:
    s = pd.read_csv(f"{R}/S2.7_section_{A}.csv")
    for _, r in s.iterrows():
        if not r.n_blocks: continue
        lo_, hi_ = (0.98, 1.02) if r.band == "tissue" else (0.8, 1.25)
        es = "confermata" if (r.lo_pieces_one >= lo_ and r.hi_pieces_one <= hi_) else ("smentita" if (r.hi_pieces_one < lo_ or r.lo_pieces_one > hi_) else "indeterminata")
        rows.append(f"| {NAME[A]} | {r.band} | {int(r.n_blocks)} | {r.ratio_pieces_one:.4f} | [{r.lo_pieces_one:.4f}, {r.hi_pieces_one:.4f}] | {es} |")
open(f"{PD}/p2.md", "w").write("\n".join(rows) + "\n")
print("parti e figure scritte")
