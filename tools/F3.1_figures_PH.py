"""tools/F3.1_figures_PH.py - F3.1: figure delle indagini PH-4...PH-6 (stile di tools/F3.1_figures.py: vero grigio, simulato blu, profili dal nucleo arancio,
compartimenti Q6 verde). Uso: .venv/bin/python tools/F3.1_figures_PH.py <results_dir> <out_dir>"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
res, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.dpi": 200})
GREY, BLUE, ORANGE, GREEN = "#7f7f7f", "#1f5fa8", "#d9822b", "#3a9a5b"
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; TN = ["intestino", "tumore\ndel colon", "stroma", "linfonodo", "cervello", "cuore"]
rng = np.random.default_rng(0)
def strip(ax, df, col, sides, dx):
    for (side, color, lab), off in zip(sides, dx):
        for i, A in enumerate(ARCH):
            v = df.loc[(df.side == side) & (df.A == A), col].dropna()
            if len(v): ax.scatter(i + off + rng.uniform(-0.04, 0.04, len(v)), v, s=8, color=color, zorder=3, label=lab if i == 0 else None)
    ax.set_xticks(range(6)); ax.set_xticklabels(TN)
# ---- figura 5: marcatori estranei nel solo nucleo; linfonodo, immunoglobuline e origine dei trascritti
sp_ = pd.read_csv(f"{res}/F3.1_PH4d_spill_areas.csv"); sn = pd.read_csv(f"{res}/F3.1_PH4d_spill_simnuc_areas.csv"); sn["side"] = "simnuc"
D = pd.concat([sp_, sn], ignore_index=True)
fig, axes = plt.subplots(1, 3, figsize=(11, 3.9), gridspec_kw={"width_ratios": [1.5, 1, 1.05]})
ax = axes[0]
strip(ax, D, "foreign_nuc", [("real", GREY, "vero"), ("sim", BLUE, "simulato (2 semi)"), ("simq6", GREEN, "simulato, compartimenti"), ("simnuc", ORANGE, "simulato, profili dal nucleo")], [-0.3, -0.1, 0.1, 0.3])
for i, A in enumerate(ARCH):
    fp = sp_.loc[(sp_.side == "sim") & (sp_.A == A), "foreign_profile"].mean(); ax.plot([i - 0.4, i + 0.4], [fp, fp], color="black", lw=0.8, ls="--", label="profili generati, prima dello spostamento" if i == 0 else None)
ax.set_ylabel("UMI su marcatori di altri tipi (quota)"); ax.set_title("(a) marcatori estranei nel solo nucleo (senza allargamento)", loc="left"); ax.legend(frameon=False, loc="upper left", fontsize=6)
ig = pd.read_csv(f"{res}/F3.1_PH4c_ig_share.csv").set_index("type")
types = ["linfociti T", "endoteliali", "mieloidi", "stromali", "linfociti B", "plasmacellule"]; ax = axes[1]; x = np.arange(len(types)); w = 0.19
for k, (c, color, lab) in enumerate([("reference", "black", "riferimento (nuclei veri + corona)"), ("real_nucleus_bins", GREY, "vero, solo nucleo"),
                                      ("sim_profiles", "#9ecae1", "simulato, profili generati"), ("sim_nucleus_bins", BLUE, "simulato, solo nucleo"), ("simq6_nucleus_bins", GREEN, "compartimenti, solo nucleo")]):
    ax.bar(x + (k - 2) * w * 0.85, ig.loc[types, c], width=w * 0.85, color=color, label=lab)
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([t.replace("linfociti ", "") for t in types], rotation=30, ha="right")
ax.set_ylabel("UMI delle immunoglobuline (quota)"); ax.set_title("(b) linfonodo: immunoglobuline per tipo", loc="left"); ax.legend(frameon=False, fontsize=5.5, loc="upper left")
og = pd.read_csv(f"{res}/F3.1_PH4c_origin.csv"); ax = axes[2]
tt = ["linfociti T", "endoteliali", "mieloidi", "stromali", "linfociti B", "plasmacellule"]; xs = []
for j, side in enumerate(["sim", "simq6"]):
    g = og[og.side == side].set_index("type").loc[tt]; own = g.own / g.umi; bpl = g.from_B_or_plasma / g.umi; oth = 1 - own - bpl
    xx = np.arange(len(tt)) * 2.6 + j * 1.1; xs.append(xx)
    ax.bar(xx, own, color=BLUE if side == "sim" else GREEN, width=1.0, label=f"dalla cellula stessa ({'uniforme' if side == 'sim' else 'compartimenti'})")
    ax.bar(xx, bpl, bottom=own, color="#c0392b", width=1.0, alpha=0.85 if side == "sim" else 0.55, label="da B o plasmacellule" if j == 0 else None)
    ax.bar(xx, oth, bottom=own + bpl, color="#dddddd", width=1.0, label="da altre cellule" if j == 0 else None)
ax.set_xticks(np.arange(len(tt)) * 2.6 + 0.55); ax.set_xticklabels([t.replace("linfociti ", "") for t in tt], rotation=30, ha="right")
ax.set_ylabel("quota degli UMI del nucleo"); ax.set_ylim(0, 1); ax.set_title("(c) linfonodo simulato: da chi vengono\ngli UMI del nucleo", loc="left"); ax.legend(frameon=False, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
fig.tight_layout(); fig.savefig(f"{out}/F3.1_fig5_spill.png"); plt.close(fig)
# ---- figura 6: anisotropia (PH-5) e variabilita' fra aree (PH-6)
an = pd.read_csv(f"{res}/F3.1_PH5_aniso.csv"); cv = pd.read_csv(f"{res}/F3.1_PH6_summary.csv")
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
for ax, col, ttl in [(axes[0], "aniso_umi", "(a) anisotropia del campo degli UMI"), (axes[1], "aniso_type", "(b) anisotropia del campo dei tipi di RCTD")]:
    strip(ax, an, col, [("real", GREY, "vero"), ("sim", BLUE, "simulato"), ("simq6", GREEN, "compartimenti")], [-0.22, 0, 0.22])
    ax.set_title(ttl, loc="left"); ax.set_ylabel("massimo meno minimo sulle direzioni")
axes[0].legend(frameon=False, loc="upper left")
ax = axes[2]; x = np.arange(6)
for k, (side, color, lab) in enumerate([("sim", BLUE, "uniforme"), ("simq6", GREEN, "compartimenti")]):
    c = cv[cv.side == side].set_index("A").loc[ARCH]
    ax.scatter(x + (k - 0.5) * 0.3, c.median_cv_ratio_types, marker="s", color=color, s=14, label=f"composizione ({lab})")
    ax.scatter(x + (k - 0.5) * 0.3, c.cv_ratio_density, marker="^", color=color, s=14, label=f"densità ({lab})")
    ax.scatter(x + (k - 0.5) * 0.3, c.median_cv_ratio_methods, marker="o", facecolor="white", edgecolor=color, s=16, label=f"misure M-1...M-6 ({lab})")
ax.axhline(1, color=GREY, lw=0.7); ax.axhline(0.5, color=GREY, lw=0.7, ls="--"); ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(TN)
ax.set_ylabel("variabilità fra aree, simulato / vero (CV)"); ax.set_title("(c) variabilità fra le 5 aree (mediane)", loc="left"); ax.legend(frameon=False, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
fig.tight_layout(); fig.savefig(f"{out}/F3.1_fig6_aniso_var.png"); plt.close(fig)
print("ok")
