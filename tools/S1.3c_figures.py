"""tools/S1.3c_figures.py — figure del report S1.3c. Uso: .venv/bin/python tools/S1.3c_figures.py  (scrive results/S1.3c/figures/)"""
import shutil, glob
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/S1.3c"; FIGD = RES / "figures"; FIGD.mkdir(exist_ok=True)
OUT = Path("/mnt/micron/geo_spatialtrans/S1.3c")
for f in glob.glob(str(OUT / "posthoc" / "PH*.csv")): shutil.copy(f, RES / ("S1.3c_" + Path(f).name))
for f in glob.glob(str(OUT / "posthoc" / "PH4_*.png")): shutil.copy(f, FIGD / Path(f).name)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 7, "figure.facecolor": "white", "axes.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False})
ARCH = [f"A{i}" for i in range(1, 7)]
NAME = {"A1": "A1 epitelio", "A2": "A2 tumore", "A3": "A3 stroma", "A4": "A4 linfoide", "A5": "A5 cervello", "A6": "A6 cuore"}
FOC, GREY, ORA = "#1f4e79", "#9a9a9a", "#d17a00"
R = pd.read_csv(RES / "S1.3c_roi.csv")

# F1: E per ROI
fig, ax = plt.subplots(figsize=(7.0, 3.0))
ser = [("primary", 3, "primaria, s = 6 µm", FOC, "o", 0.0), ("primary", 5, "primaria, s = 10 µm", GREY, "s", 0.18),
       ("sens_seg", 3, "altro segmentatore, s = 6", ORA, "^", -0.18), ("sens_allgenes", 3, "tutti i geni, s = 6", "#5aa0d0", "v", 0.09)]
for v, k, lab, c, m, off in ser:
    g = R[(R.variant == v) & (R.k == k)]
    x = np.array([ARCH.index(a) for a in g.archetype]) + off
    ax.scatter(x, g.E, s=14 if v == "primary" and k == 3 else 9, color=c, marker=m, label=lab, zorder=3 if k == 3 and v == "primary" else 2, alpha=0.9)
ax.axhline(0.05, color="black", lw=0.8, ls="--"); ax.text(5.45, 0.051, "soglia di rilevanza 0.05", ha="right", va="bottom", fontsize=7)
ax.set_xticks(range(6)); ax.set_xticklabels([NAME[a] for a in ARCH]); ax.set_ylabel("E = φ_obs − φ_null")
ax.set_ylim(-0.003, 0.06); ax.set_title("L'eccesso di finestre vuote resta sotto la soglia in tutti i 30 ROI e in tutte le varianti", loc="left")
ax.legend(frameon=False, loc="upper left", ncol=2); fig.tight_layout(); fig.savefig(FIGD / "F1_E_per_roi.png", dpi=200); plt.close(fig)

# F2: calibrazione PH-3
C = pd.read_csv(RES / "S1.3c_PH3_calibration.csv")
fig, axes = plt.subplots(1, 6, figsize=(10, 2.5), sharey=True)
for ax, A in zip(axes, ARCH):
    for typ, c, m in (("A", FOC, "o"), ("B", ORA, "s")):
        g = C[(C.archetype == A) & (C.type == typ)]
        ax.scatter(g.f_target + (0.006 if typ == "B" else -0.006), g.E, s=6, color=c, marker=m, alpha=0.5)
        md = g.groupby("f_target").E.median(); ax.plot(md.index, md.values, color=c, lw=1.2, label=f"vuoto tipo {typ}")
    ax.plot([0, 0.32], [0, 0.32], color=GREY, lw=0.7, ls=":", label="recupero perfetto"); ax.axhline(0.05, color="black", lw=0.7, ls="--")
    ax.set_title(NAME[A], loc="left"); ax.set_xlim(0, 0.33); ax.set_ylim(0, 0.33); ax.set_xticks([0.05, 0.1, 0.2, 0.3])
    ax.set_xticklabels(["5", "10", "20", "30"]); ax.set_xlabel("vuoto inserito (% dei bin)")
axes[0].set_ylabel("E misurata"); axes[0].legend(frameon=False, loc="upper left")
fig.suptitle("Anche con vuoto sintetico perfetto (bin azzerati), E recupera dal 10 % (A4) al 97 % (A6) della frazione inserita (n = 5 ROI; linea = mediana)", x=0.01, ha="left", fontsize=8)
fig.tight_layout(); fig.savefig(FIGD / "F2_calibrazione_E.png", dpi=200); plt.close(fig)

# F3: profilo dei lumi PH-1
P = pd.read_csv(RES / "S1.3c_PH1_profile.csv"); L = pd.read_csv(RES / "S1.3c_PH1_levels.csv")
bands = ["0-2", "2-4", "4-8", "8-16", "16-32", ">32"]
fig, axes = plt.subplots(1, 6, figsize=(10, 2.6), sharey=True)
for ax, A in zip(axes, ARCH):
    g = P[(P.archetype == A) & (P.n >= 20)]
    for roi, gg in g.groupby("roi_id"):
        gg = gg.set_index("band").reindex(bands); ax.plot(range(6), gg.umi_mean, color=FOC, lw=0.8, alpha=0.6, marker="o", ms=2)
    l = L[L.archetype == A]
    ax.axhline(l.umi_tissue_mean.median(), color=ORA, lw=1, label="tessuto (mediana ROI)")
    ax.axhline(l.umi_bg_far_per_bin.median(), color="black", lw=0.8, ls="--", label="sfondo ≥ 50 µm fuori tessuto")
    ax.set_yscale("log"); ax.set_xticks(range(6)); ax.set_xticklabels(bands, rotation=60); ax.set_title(NAME[A], loc="left")
    ax.set_xlabel("distanza dal tessuto (µm)")
    if A in ("A5", "A6"): ax.set_title(NAME[A] + " †", loc="left")
axes[0].set_ylabel("UMI medi per bin «lume»"); h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, loc="upper right", ncol=2, fontsize=7)
fig.suptitle("Nei lumi veri (A1–A4) i bin non scendono al livello dello sfondo lontano (una linea per ROI, fasce con ≥ 20 bin)\n† in A5 e A6 il «lume» è soprattutto tessuto pallido escluso dalla maschera R2 (PH-4): non è vuoto vero", x=0.01, ha="left", fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(FIGD / "F3_profilo_lumi.png", dpi=200); plt.close(fig)

# F4: UMI vs distanza dal nucleo
import json
fig, axes = plt.subplots(1, 6, figsize=(10, 2.5), sharey=True)
for ax, A in zip(axes, ARCH):
    Rs, Rp = [], []
    for f in sorted(glob.glob(str(OUT / "roi" / f"{A}_*_profile.parquet"))):
        p = pd.read_parquet(f); base = p[p.d_bin < 2].eval("mean*size").sum() / p[p.d_bin < 2]["size"].sum()
        ax.plot(p.d_bin + 0.5, p["mean"] / base, color=FOC, lw=0.8, alpha=0.7)
        roi = Path(f).name.split("_")[1]; d = json.load(open(OUT / "roi" / f"{A}_{roi}_primary.json"))
        Rs.append(d["by_k"]["3"]["R_profile"]); Rp.append(d["CP1_R_profile"])
    ax.axhline(1, color=GREY, lw=0.7, ls=":")
    ax.text(0.97, 0.95, f"R reale {np.median(Rs):.2f}\nR placebo {np.median(Rp):.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=6)
    ax.set_title(NAME[A], loc="left"); ax.set_xlabel("distanza dal nucleo (µm)"); ax.set_xlim(0, 20)
axes[0].set_ylabel("UMI per bin / UMI a 0–2 µm")
fig.suptitle("Gli UMI calano con la distanza dal nucleo in A1, A5, A6 e il placebo li appiattisce; in A3 calano anche col placebo (simulatore attuale: 1; n = 5 ROI)", x=0.01, ha="left", fontsize=8)
fig.tight_layout(); fig.savefig(FIGD / "F4_profilo_nucleo.png", dpi=200); plt.close(fig)
print("figure ok", sorted(p.name for p in FIGD.iterdir()))
