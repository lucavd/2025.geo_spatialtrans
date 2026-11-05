"""tools/S2.1_analysis.py - S2.1: tabelle per tessuto dei descrittivi D-1...D-4 (media delle 5 aree, intervallo bootstrap delle aree,
B = 2000, seme 20261105) e figure del report. Uso (root del repo): .venv/bin/python tools/S2.1_analysis.py"""
import glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES = "results/S2.1"; FIG = f"{RES}/figures"; import os; os.makedirs(FIG, exist_ok=True)
SEED, B = 20261105, 2000
NAME = {"A1": "intestine (mouse)", "A2": "colon tumour", "A3": "tumour stroma", "A4": "lymph node", "A5": "brain (mouse)", "A6": "heart"}
sim = pd.read_csv(f"{RES}/S2.1_desc_sim_roi.csv"); real = pd.read_csv(f"{RES}/S2.1_desc_real_roi.csv")
# regola del centro: nuclei o cellule che contengono il centro di almeno uno dei 16 figli (come Space Ranger); tocco: area > 0 (revisione F01)
Q = {"d1": ("d1_centre_in_nucleus", "seg_d1_centre_in_nucleus", "sr_d1_centre_in_nucleus"),
     "d2c": ("d2c_nuclei_mean", "seg_d2c_nuclei_mean", "sr_d2_nuclei_mean"),
     "d2t": ("d2_nuclei_mean", "seg_d2_nuclei_mean", None),
     "d3c": ("d3c_cells_mean", None, "sr_d3_cells_mean"), "d3_pure": ("d3_pure_bins", None, None)}
LAB = {"d1": "2 um bins in tissue with centre in a nucleus", "d2c": "nuclei per 8 um bin, centre rule (mean)",
       "d2t": "nuclei per 8 um bin, touch rule (mean)", "d3c": "cells per 8 um bin, centre rule (mean)",
       "d3_pure": "2 um bins inside one territory (>= 99 %)"}
rng = np.random.default_rng(SEED); rows = []
def boot(v):
    v = np.asarray(v, float); idx = rng.integers(0, len(v), size=(B, len(v))); m = v[idx].mean(1)
    return v.mean(), np.percentile(m, 2.5), np.percentile(m, 97.5)
for A in NAME:
    for q, (cs, cg, cr) in Q.items():
        for src, df, col in (("simulated", sim, cs), ("segmenter (Step 1)", real, cg), ("Space Ranger 4", real, cr)):
            if col is None or col not in df: continue
            v = df.loc[df.archetype == A, col].dropna()
            if len(v) == 0: continue
            m, lo, hi = boot(v); rows.append(dict(archetype=A, tissue=NAME[A], quantity=q, label=LAB[q], source=src, n_roi=len(v), mean=m, lo=lo, hi=hi))
tab = pd.DataFrame(rows); tab.to_csv(f"{RES}/S2.1_desc_tissue.csv", index=False)
# D-4: tempi e memoria
RUN = os.environ.get("S21_RUN", "S2.1/run1")
tt = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"/mnt/micron/geo_spatialtrans/{RUN}/grid/*_time.csv"))])
tt["archetype"] = tt.tag.str[:2]; tt.to_csv(f"{RES}/S2.1_D4_time.csv", index=False)
# figura descrittivi: un pannello per quantita', tessuti in ascissa, una serie per fonte
COL = {"simulated": "#1f77b4", "segmenter (Step 1)": "#7f7f7f", "Space Ranger 4": "#d62728"}
fig, axs = plt.subplots(1, 5, figsize=(19, 4.2), facecolor="white")
for ax, q in zip(axs, Q):
    t = tab[tab.quantity == q]
    for k, (src, c) in enumerate(COL.items()):
        u = t[t.source == src]; x = np.array([list(NAME).index(a) for a in u.archetype]) + (k - 1) * 0.22
        ax.errorbar(x, u["mean"], yerr=[u["mean"] - u.lo, u.hi - u["mean"]], fmt="o", color=c, label=src, ms=5, capsize=2)
        roi = (sim if src == "simulated" else real)
        col = {"simulated": Q[q][0], "segmenter (Step 1)": Q[q][1], "Space Ranger 4": Q[q][2]}[src]
        if col and col in roi:
            for a in NAME:
                v = roi.loc[roi.archetype == a, col].dropna(); ax.scatter(np.full(len(v), list(NAME).index(a) + (k - 1) * 0.22), v, s=8, color=c, alpha=0.35)
    ax.set_xticks(range(6)); ax.set_xticklabels([NAME[a] for a in NAME], rotation=40, ha="right", fontsize=8); ax.set_title(LAB[q], fontsize=9)
    ax.grid(axis="y", lw=0.3)
axs[0].legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig(f"{FIG}/S2.1_desc.png", dpi=150)
# figura D-2: distribuzione dei nuclei per blocco di 8 um, un pannello per tessuto
fig, axs = plt.subplots(2, 3, figsize=(12, 6.5), facecolor="white", sharey=True)
for ax, A in zip(axs.ravel(), NAME):
    for k, (src, pre, df) in enumerate((("simulated", "d2c_nuclei", sim), ("segmenter (Step 1)", "seg_d2c_nuclei", real), ("Space Ranger 4", "sr_d2_nuclei", real))):
        cols = [f"{pre}_p{i}" for i in range(5)]
        if cols[0] not in df or df.loc[df.archetype == A, cols[0]].isna().all(): continue
        v = df.loc[df.archetype == A, cols].mean().to_numpy()
        ax.bar(np.arange(5) + (k - 1) * 0.27, v, width=0.27, color=list(COL.values())[k], label=src)
    ax.set_xticks(range(5)); ax.set_xticklabels(["0", "1", "2", "3", ">=4"]); ax.set_title(NAME[A], fontsize=10); ax.grid(axis="y", lw=0.3)
axs[0, 0].legend(fontsize=7, frameon=False); fig.supxlabel("nuclei per 8 um bin in tissue (centre rule: nucleus contains the centre of at least one 2 um child)"); fig.supylabel("fraction of 8 um bins (mean of 5 areas)")
fig.tight_layout(); fig.savefig(f"{FIG}/S2.1_D2_dist.png", dpi=150)
# figura CP-1: accordo con Space Ranger per area e per griglia, un pannello per tessuto
cp = pd.read_csv(f"{RES}/S2.1_CP1_roi.csv"); GR = {"none": "our grid", "transpose": "rows <-> columns", "half": "half-bin shift", "reverse": "rows reversed"}
fig, axs = plt.subplots(1, 5, figsize=(16, 3.8), facecolor="white", sharey=True)
for ax, A in zip(axs, ["A2", "A3", "A4", "A5", "A6"]):
    for k, gname in enumerate(GR):
        u = cp[(cp.archetype == A) & (cp.grid == gname)]
        ax.scatter(np.full(len(u), k), u.jaccard, s=18, color="#d62728" if gname == "none" else "#7f7f7f")
    ax.axhline(0.999, color="#d62728", lw=0.6, ls="--"); ax.axhline(0.9, color="#7f7f7f", lw=0.6, ls=":")
    ax.set_xticks(range(4)); ax.set_xticklabels(list(GR.values()), rotation=35, ha="right", fontsize=8); ax.set_title(NAME[A], fontsize=10); ax.grid(axis="y", lw=0.3)
axs[0].set_ylabel("agreement with Space Ranger (Jaccard)"); fig.tight_layout(); fig.savefig(f"{FIG}/S2.1_CP1.png", dpi=150)
print(tab.round(4).to_string()); print(tt.groupby("archetype")[["n_cells", "elapsed_s", "r_max_used_mb", "vm_hwm_mb"]].agg(["median", "max"]))
