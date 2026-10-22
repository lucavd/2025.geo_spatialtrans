"""R6 figures (white background, one panel per archetype). Usage: python tools/R6_figures.py [variant=ph1]"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
V = sys.argv[1] if len(sys.argv) > 1 else "ph1"; SUF = "" if V == "prereg" else "_" + V
M = "/mnt/micron/geo_spatialtrans/R6"; OUT = "results/R6/figures"; os.makedirs(OUT, exist_ok=True)
TIS = {"A1": "A1 intestino di topo", "A2": "A2 tumore del colon", "A3": "A3 stroma del tumore", "A4": "A4 linfonodo", "A5": "A5 cervello di topo", "A6": "A6 cuore"}
ARCH = list(TIS)
labs = {A: pd.read_parquet(f"{M}/labels/{A}_exp_umi20_seed20261022_nuclei{SUF}.parquet") for A in ARCH}
pal = plt.get_cmap("tab10")
def colors(types):
    types = sorted(types, key=lambda t: (t in ("non risolto", "non classificabile"), t))
    return {t: ("lightgrey" if t == "non classificabile" else "darkgrey" if t == "non risolto" else pal(i % 10)) for i, t in enumerate(types)}
# Fig 1: fractions per ROI (stacked)
fig, axes = plt.subplots(2, 3, figsize=(14, 8), facecolor="white")
for ax, A in zip(axes.ravel(), ARCH):
    d = labs[A]; ct = pd.crosstab(d.roi_id, d.cell_type, normalize="index"); col = colors(ct.columns)
    bottom = np.zeros(len(ct))
    for t in col:
        ax.bar(ct.index, ct[t].values, bottom=bottom, color=col[t], label=t, edgecolor="white", linewidth=0.5); bottom += ct[t].values
    ax.set_title(TIS[A]); ax.set_ylim(0, 1); ax.set_ylabel("frazione dei nuclei"); ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
fig.suptitle(f"Composizione per tipo nelle 5 ROI di ogni tessuto (etichette {V})"); fig.tight_layout(); fig.savefig(f"{OUT}/R6_fig1_fractions{SUF}.png", dpi=150); plt.close(fig)
# Fig 2: nuclear area by type
fig, axes = plt.subplots(2, 3, figsize=(14, 8), facecolor="white")
for ax, A in zip(axes.ravel(), ARCH):
    d = labs[A]; d = d[~d.cell_type.isin(["non classificabile"])]; fr = d.cell_type.value_counts(normalize=True); ty = [t for t in fr.index if fr[t] >= 0.01]
    data = [d.area_um2[d.cell_type == t].values for t in ty]
    bp = ax.boxplot(data, tick_labels=[f"{t}\n(n={len(x)})" for t, x in zip(ty, data)], showfliers=False, patch_artist=True)
    col = colors(ty)
    for patch, t in zip(bp["boxes"], ty): patch.set_facecolor(col[t]); patch.set_alpha(0.7)
    ax.set_yscale("log"); ax.set_ylabel("area del nucleo (µm²)"); ax.set_title(TIS[A]); ax.tick_params(axis="x", labelsize=7)
fig.suptitle(f"Area del nucleo per tipo (tipi con almeno l'1 % dei nuclei; etichette {V})"); fig.tight_layout(); fig.savefig(f"{OUT}/R6_fig2_area_by_type{SUF}.png", dpi=150); plt.close(fig)
# Fig 3: spatial map of one ROI per archetype (first ROI), types as colours
fig, axes = plt.subplots(2, 3, figsize=(15, 10), facecolor="white")
for ax, A in zip(axes.ravel(), ARCH):
    d = labs[A]; roi = sorted(d.roi_id.unique())[0]; d = d[d.roi_id == roi]; col = colors(d.cell_type.unique())
    for t, c in col.items():
        g = d[d.cell_type == t]; ax.scatter(g.x_um, g.y_um, s=1.5 if A == "A4" else 4, c=[c], label=f"{t} ({len(g)})", linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(f"{TIS[A]}, ROI {roi}"); ax.legend(fontsize=6, markerscale=4, loc="upper right", framealpha=0.9); ax.set_xlabel("µm")
fig.suptitle(f"Mappa dei tipi in una ROI per tessuto (etichette {V})"); fig.tight_layout(); fig.savefig(f"{OUT}/R6_fig3_maps{SUF}.png", dpi=150); plt.close(fig)
# Fig 4: confusion with Space Ranger (A5, A6)
fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor="white")
for ax, A in zip(axes, ["A5", "A6"]):
    f = f"results/R6/R6_sr_confusion_{A}{SUF}.csv"
    if not os.path.exists(f): ax.set_visible(False); continue
    ct = pd.read_csv(f, index_col=0); im = ax.imshow(ct.values, cmap="Blues")
    ax.set_xticks(range(ct.shape[1])); ax.set_xticklabels(ct.columns, rotation=45, ha="right", fontsize=8); ax.set_yticks(range(ct.shape[0])); ax.set_yticklabels(ct.index, fontsize=8)
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]): ax.text(j, i, int(ct.values[i, j]), ha="center", va="center", fontsize=7, color="black" if ct.values[i, j] < ct.values.max() / 2 else "white")
    ax.set_xlabel("tipo Space Ranger (10x Cloud, mappato al catalogo)"); ax.set_ylabel("nostro tipo"); ax.set_title(TIS[A])
fig.suptitle(f"Confronto con l'annotazione di Space Ranger sui nuclei appaiati (etichette {V})"); fig.tight_layout(); fig.savefig(f"{OUT}/R6_fig4_sr_confusion{SUF}.png", dpi=150); plt.close(fig)
print("figures written", OUT)
