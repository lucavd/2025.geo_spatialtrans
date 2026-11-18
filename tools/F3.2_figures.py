#!/usr/bin/env python3
"""tools/F3.2_figures.py - F3.2: figure del report (sfondo bianco, un pannello per tessuto dove pertinente).
fig1 censimento: valori per famiglia e insieme di aree; fig2 prova a secco di R-4 (F3.1, dati di taratura); fig3 mappa delle aree per sezione.
Uso: .venv/bin/python tools/F3.2_figures.py [fig1|fig2|fig3|all] -> reports/F3.2_parts/F3.2_fig*.png
"""
import sys, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; R = ROOT / "results/F3.2"; OUT = ROOT / "reports/F3.2_parts"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white"})
SETS = [("R2_ROI", "le 30 aree di R2 (5 per tessuto)", "#1f4e79"), ("R2b_WINDOW", "finestre contate a mano, dentro le aree di R2", "#6fa8dc"),
        ("R5_EDGE_TILE", "314 tessere di bordo di R5", "#e69f00"), ("SECTION", "sezione intera", "#999999"), ("DESIGN", "impostazione, non stima", "#dddddd")]
FAM = {"cell_layer.type_catalogue": "catalogo dei tipi (densità, taglia del nucleo)", "cell_layer.tissue": "preset del tessuto",
       "cell_layer.compartments(option)": "compartimenti (opzione)", "cell_layer.composition": "composizione dei tipi",
       "cell_layer.heterogeneity": "vuoti e densità locale", "cell_layer.size_resid(option)": "residui del modello di taglia (opzione)",
       "cell_layer.defaults": "impostazioni predefinite", "step2.capture": "cattura: strisce, fondo, seconda coda",
       "step2.transcripts": "trascritti: spostamento, nucleo/citoplasma, campo fine", "step2.depth": "fattore di profondità",
       "step2.profiles": "profili di espressione (tabelle)", "step2.grid": "griglia", "step2.pieces": "simulazione a pezzi"}

def fig1():
    c = pd.read_csv(R / "F3.2_census.csv")
    t = c.groupby(["family", "area_set"]).size().unstack(fill_value=0)
    t = t.loc[t.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    left = np.zeros(len(t))
    for k, lab, col in SETS:
        v = t[k].values if k in t else np.zeros(len(t))
        ax.barh(range(len(t)), v, left=left, color=col, edgecolor="white", linewidth=0.5, label=lab); left += v
    ax.set_yticks(range(len(t))); ax.set_yticklabels([FAM.get(f, f) for f in t.index])
    ax.set_xlabel("valori usati per difetto dal simulatore (numero)")
    n_same = int(c.area_set.isin(["R2_ROI", "R2b_WINDOW"]).sum())
    ax.set_title(f"{n_same} dei {len(c)} valori vengono dalle stesse 30 aree di 1 mm²\n(una sola sezione per tessuto)", loc="left")
    ax.legend(frameon=False, loc="lower right"); ax.margins(x=0.02)
    fig.tight_layout(); fig.savefig(OUT / "F3.2_fig1_census.png", dpi=150); plt.close(fig)

QLAB = {"density_per_mm2": "densità dei nuclei", "umi_per_bin": "UMI per quadratino", "M-1": "UMI assegnati", "M-2": "UMI per cellula",
        "M-3": "marcatori estranei", "M-4": "singoletti di RCTD", "M-5": "concordanza con il nucleo", "M-6": "coerenza spaziale"}
TIS = {"A1": "intestino di topo", "A2": "tumore del colon", "A3": "stroma del tumore", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}

def fig2():
    d = pd.read_csv(ROOT / "results/F3.1/F3.1_PH6_areas.csv"); dr = pd.read_csv(R / "F3.2_R4_dryrun.csv")
    Q = list(QLAB)
    fig, axs = plt.subplots(2, 3, figsize=(10.5, 6.4), sharex=True)
    for ax, A in zip(axs.ravel(), TIS):
        g = d[d.A == A]
        for qi, q in enumerate(Q):
            m = g[g.side == "real"][q].mean()
            for side, col, off, mk in [("real", "#555555", -0.22, "o"), ("sim", "#1f4e79", 0.0, "s"), ("simq6", "#e69f00", 0.22, "^")]:
                v = g[g.side == side][q].dropna().values / m
                ax.scatter(np.full(len(v), qi + off), v, s=14, color=col, marker=mk, alpha=0.85, linewidths=0, zorder=3)
        ax.axhline(1, color="#bbbbbb", lw=0.8, zorder=1)
        k = dr[(dr.arm == "sim") & (dr.archetype == A)].real_in_sim_range.sum()
        ax.set_title(f"{TIS[A]}\n{k} {'area vera' if k == 1 else 'aree vere'} su 40 dentro l'intervallo simulato", loc="left")
        ax.set_xticks(range(len(Q))); ax.set_xticklabels([QLAB[q] for q in Q], rotation=55, ha="right")
        ax.margins(x=0.04, y=0.08)
    for ax in axs[:, 0]: ax.set_ylabel("valore / media delle aree vere")
    h = [plt.Line2D([], [], color=c_, marker=m_, ls="", label=l_) for c_, m_, l_ in [("#555555", "o", "aree vere (5)"), ("#1f4e79", "s", "simulate, semina uniforme (5)"),
                                                                                    ("#e69f00", "^", "simulate, semina a compartimenti (5)")]]
    fig.legend(handles=h, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("Le aree simulate variano molto meno delle vere: 13 aree vere su 240 cadono dentro l'intervallo delle simulate", x=0.01, ha="left", y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.91)); fig.savefig(OUT / "F3.2_fig2_R4_dryrun.png", dpi=150); plt.close(fig)

def fig3():
    import importlib.util
    s = importlib.util.spec_from_file_location("a", str(ROOT / "tools/F3.2_areas.py")); A = importlib.util.module_from_spec(s); s.loader.exec_module(A)
    u = pd.read_csv(R / "F3.2_used_areas.csv"); secs = ["A1", "A23", "A4", "A5", "A6"]
    NAME = {"A1": "intestino di topo", "A23": "tumore e stroma del colon (P2)", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}
    fig, axs = plt.subplots(2, 3, figsize=(11, 7.4)); axs = axs.ravel()
    for ax, sec in zip(axs, secs):
        B, _ = A.bin_scores(sec); b = B.sample(min(len(B), 60000), random_state=1)
        ax.scatter(b.x / 1000, b.y / 1000, s=0.2, color="#d9d9d9", rasterized=True)
        for _, r in u[(u.section == sec) & u.x0.notna()].iterrows():
            col = {"R2_ROI": "#1f4e79", "R2b_WINDOW": "#6fa8dc", "R5_EDGE_TILE": "#e69f00"}[r.set]
            ax.add_patch(Rectangle((r.x0 / 1000, r.y0 / 1000), (r.x1 - r.x0) / 1000, (r.y1 - r.y0) / 1000, fill=False, ec=col, lw=0.8))
        f = R / f"F3.2_pool_{sec}.csv"
        if f.exists():
            p = pd.read_csv(f)
            for _, r in p.iterrows():
                col = "#cc0000" if r.level == "test_tile" else "#2e7d32"
                ax.add_patch(Rectangle((r.x0 / 1000, r.y0 / 1000), (r.x1 - r.x0) / 1000, (r.y1 - r.y0) / 1000, fill=False, ec=col, lw=1.0, ls="--"))
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(NAME[sec], loc="left"); ax.set_xlabel("mm"); ax.margins(0.03)
    axs[-1].set_axis_off()
    h = [Rectangle((0, 0), 1, 1, fill=False, ec=c_, ls=l_, label=t_) for c_, l_, t_ in [("#1f4e79", "-", "aree di R2"), ("#6fa8dc", "-", "finestre contate a mano"),
         ("#e69f00", "-", "tessere di bordo di R5 (con margine)"), ("#2e7d32", "--", "candidate per la verifica"), ("#cc0000", "--", "tessera di prova di StarDist")]]
    axs[-1].legend(handles=h, frameon=False, loc="center")
    fig.suptitle("Le aree candidate per la verifica stanno ad almeno 280 µm da ogni area già usata", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(OUT / "F3.2_fig3_areas_map.png", dpi=150); plt.close(fig)

if __name__ == "__main__":
    w = sys.argv[1] if len(sys.argv) > 1 else "all"
    for k, f in [("fig1", fig1), ("fig2", fig2), ("fig3", fig3)]:
        if w in ("all", k): f(); print(k, "ok")
