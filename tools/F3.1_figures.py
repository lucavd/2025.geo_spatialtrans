"""tools/F3.1_figures.py - F3.1: figure del report (sfondo bianco; vero in grigio, simulato in blu, profili dal solo nucleo in arancio).
Uso: .venv/bin/python tools/F3.1_figures.py <results_dir> <out_dir>"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
res, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.dpi": 200})
GREY, BLUE, ORANGE = "#7f7f7f", "#1f5fa8", "#d9822b"
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; TN = ["intestino", "tumore\ndel colon", "stroma", "linfonodo", "cervello", "cuore"]
MEAS = [("M-1", "share_assigned_exp", "UMI assegnati alle cellule (quota)"), ("M-2", "umi_per_cell_exp", "UMI per cellula"),
        ("M-3", "foreign_exp", "UMI su marcatori di altri tipi (quota)"), ("M-4", "rctd_frac_singlet", "blocchi chiamati singoletto (quota)"),
        ("M-5", "rctd_agree_nucleus", "tipo di RCTD = tipo del nucleo (quota)"), ("M-6", "rctd_neighbour_same", "blocchi vicini con lo stesso tipo (quota)")]
V = pd.read_csv(f"{res}/F3.1_verdicts.csv"); real = pd.read_csv(f"{res}/F3.1_real_runs.csv"); sim = pd.read_csv(f"{res}/F3.1_sim_areas.csv")
nuc = pd.read_csv(f"{res}/F3.1_simnuc_runs.csv")
for t in (real, sim, nuc): t["A"] = t.tag.str[:2]
MARK = {"confermata": "=", "smentita": "×", "indeterminata": "?", "non valutabile": "–"}
fig, axes = plt.subplots(2, 3, figsize=(10, 6.2)); rng = np.random.default_rng(0)
for ax, (mid, col, lab) in zip(axes.ravel(), MEAS):
    for i, A in enumerate(ARCH):
        v = V[(V.id == mid) & (V.archetype == A)].iloc[0]
        ax.add_patch(plt.Rectangle((i - 0.38, v.rel_lo), 0.34, v.rel_hi - v.rel_lo, color=GREY, alpha=0.25, lw=0))
        r = real.loc[real.A == A, col]; s = sim.loc[sim.A == A, col]
        ax.scatter(i - 0.21 + rng.uniform(-0.06, 0.06, len(r)), r, s=9, color=GREY, zorder=3)
        if np.isfinite(v.get("sim_lo", np.nan)):
            ax.plot([i + 0.21] * 2, [v.sim_lo, v.sim_hi], color=BLUE, lw=2.2, alpha=0.55, solid_capstyle="butt")
        ax.scatter(i + 0.21 + rng.uniform(-0.06, 0.06, len(s)), s, s=9, color=BLUE, zorder=3)
        ax.annotate(MARK.get(v.verdict, "?"), (i, 1.0), xycoords=("data", "axes fraction"), ha="center", va="top", fontsize=8,
                    color=("black" if v.verdict == "smentita" else GREY))
    ax.set_xticks(range(6)); ax.set_xticklabels(TN); ax.set_title(f"{lab}", loc="left"); ax.margins(y=0.12)
    if mid == "M-2": ax.set_yscale("log")
axes[0, 0].scatter([], [], s=9, color=GREY, label="aree vere (banda: intervallo di rilevanza)")
axes[0, 0].scatter([], [], s=9, color=BLUE, label="aree simulate, media di 2 semi (barra: intervallo)")
fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Bin2cell e RCTD sul simulato e sul vero, 5 aree per tessuto (× smentita, ? indeterminata, = confermata)", x=0.01, ha="left", fontsize=8)
fig.tight_layout(rect=(0, 0.04, 1, 0.97)); fig.savefig(f"{out}/F3.1_fig1_verdicts.png"); plt.close(fig)
# CP-3: marcatori estranei, concordanza, coerenza: vero, simulato (k = 0), profili dal solo nucleo (stessi semi)
simr = pd.read_csv(f"{res}/F3.1_sim_runs.csv"); simr["A"] = simr.tag.str[:2]; s0 = simr[simr.k == 0]
fig, axes = plt.subplots(1, 3, figsize=(10, 3.1))
for ax, (mid, col, lab) in zip(axes, [MEAS[2], MEAS[4], MEAS[5]]):
    for i, A in enumerate(["A1", "A2", "A4", "A5", "A6"]):
        r = real.loc[real.A == A, col]; a = s0[s0.A == A].set_index("tag")[col]; b = nuc[nuc.A == A].set_index("tag")[col]; tg = a.index.intersection(b.index)
        ax.scatter(np.full(len(r), i - 0.25), r, s=9, color=GREY, zorder=3)
        for t in tg: ax.plot([i, i + 0.25], [a[t], b[t]], color="#bbbbbb", lw=0.6, zorder=1)
        ax.scatter(np.full(len(tg), i), a[tg], s=9, color=BLUE, zorder=3); ax.scatter(np.full(len(tg), i + 0.25), b[tg], s=9, color=ORANGE, zorder=3)
    ax.set_xticks(range(5)); ax.set_xticklabels(["intestino", "tumore\ndel colon", "linfonodo", "cervello", "cuore"]); ax.set_title(lab, loc="left"); ax.margins(y=0.1)
for c_, l_ in [(GREY, "vero"), (BLUE, "simulato, profili dai nuclei con la corona"), (ORANGE, "simulato, profili dal solo nucleo (stessi semi)")]:
    axes[0].scatter([], [], s=9, color=c_, label=l_)
fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Con i profili dal solo nucleo le cellule simulate restano più contaminate del vero; cambia solo il cervello (linee: stessa area, stessi semi)", x=0.01, ha="left", fontsize=8)
fig.tight_layout(rect=(0, 0.08, 1, 0.95)); fig.savefig(f"{out}/F3.1_fig2_CP3.png"); plt.close(fig)
# CP-2: misure senza verita' contro verita' (10 esecuzioni simulate per tessuto)
fig, axes = plt.subplots(1, 2, figsize=(8, 3.3)); cols = plt.get_cmap("tab10")
for k, A in enumerate(ARCH):
    S = simr[simr.A == A]
    axes[0].scatter(S.truth_foreign_type_exp, S.foreign_exp, s=10, color=cols(k), label=TN[k].replace("\n", " "))
    axes[1].scatter(S.rctd_truth_acc_first, S.rctd_agree_nucleus, s=10, color=cols(k))
axes[0].set_xlabel("vero: quota di UMI da cellule di altro tipo"); axes[0].set_ylabel("misurato: UMI su marcatori di altri tipi")
axes[1].set_xlabel("vero: accuratezza del primo tipo di RCTD"); axes[1].set_ylabel("misurato: concordanza con il nucleo")
axes[0].legend(frameon=False, fontsize=6.5); fig.suptitle("Controprova: le misure calcolabili sul vero seguono la verità del simulato?", x=0.01, ha="left", fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(f"{out}/F3.1_fig3_CP2.png"); plt.close(fig)
# PH-1: mescolanza dei vicini e riponderazione
mx = pd.read_csv(f"{res}/F3.1_PH1_mixing.csv"); mx["tag"] = mx.run.str.replace(r"_k\d$", "", regex=True); rw = pd.read_csv(f"{res}/F3.1_PH1_mixing_reweight.csv")
fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
for i, A in enumerate(ARCH):
    for side, col_, dx in [("real", GREY, -0.15), ("sim", BLUE, 0.15)]:
        v = mx[(mx.side == side) & (mx.A == A)].groupby("tag").mix_mean.mean()
        axes[0].scatter(np.full(len(v), i + dx), v, s=9, color=col_, zorder=3)
    r = rw[(rw.A == A) & (rw.side == "sim")].iloc[0]
    axes[1].plot([i, i], [r.M3_sim, r.M3_sim_reweighted_to_real_mixing], color="#bbbbbb", lw=1, zorder=1)
    axes[1].scatter([i], [r.M3_real], s=14, color=GREY, zorder=3); axes[1].scatter([i], [r.M3_sim], s=14, color=BLUE, zorder=3)
    axes[1].scatter([i], [r.M3_sim_reweighted_to_real_mixing], s=14, facecolor="white", edgecolor=BLUE, zorder=3)
for ax in axes: ax.set_xticks(range(6)); ax.set_xticklabels(TN); ax.margins(y=0.1)
axes[0].set_title("Vicini di altro tipo (quota dei 6 nuclei più vicini)", loc="left")
axes[1].set_title("UMI su marcatori di altri tipi: vero, simulato, simulato con i vicini del vero", loc="left")
for c_, l_, f_ in [(GREY, "vero", GREY), (BLUE, "simulato", BLUE), (BLUE, "simulato ripesato con la mescolanza del vero", "white")]:
    axes[1].scatter([], [], s=14, color=f_, edgecolor=c_, label=l_)
axes[1].legend(frameon=False, fontsize=6.5, loc="upper right")
fig.suptitle("Nel simulato i tipi sono mescolati più che nel vero, e questo spiega gran parte dei marcatori estranei in intestino e tumore (analisi a posteriori)", x=0.01, ha="left", fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(f"{out}/F3.1_fig4_mixing.png"); plt.close(fig)
print("figures written")
