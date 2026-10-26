"""tools/S1.5d_figures.py - S1.5d: figure del report (sfondo bianco, un pannello per tessuto). Ingressi: results/S1.5d/*.csv, celle del banco
(seme 1), nuclei veri (R6), ROI H&E (R2). Uso: .venv/bin/python -B tools/S1.5d_figures.py   Uscite: results/S1.5d/figures/*.png"""
import sys, os, importlib.util, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "tools"); import R2_common as C
s = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(s); s.loader.exec_module(cal)
RES = "results/S1.5d"; FIG = f"{RES}/figures"; os.makedirs(FIG, exist_ok=True); OUT = "/mnt/micron/geo_spatialtrans/S1.5d"
TIS = {"A1": "intestino di topo", "A2": "tumore del colon", "A3": "stroma del tumore", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False})
COL = plt.rcParams["axes.prop_cycle"].by_key()["color"]
ST = {"U": ("semina uniforme (predefinito di oggi)", "0.45", "x"), "N": ("compartimenti di S1.5c", COL[0], "s"),
      "RV": ("compartimenti con i vuoti (vincente)", COL[1], "o"), "W3": ("vuoti, residui empirici", COL[4], "D"),
      "W3p": ("vuoti, residui senza decimi (controprova)", COL[8], "v"), "U3": ("uniforme, residui empirici", COL[5], "^"),
      "U3c": ("uniforme, residui empirici, decimo ricentrato", COL[6], "P")}
V = pd.read_csv(f"{RES}/S1.5d_verdicts.csv"); ROIS = pd.read_csv("results/R2/R2_rois_checked.csv")
def band(ax, y, lo, hi): ax.plot([lo, hi], [y, y], color="0.82", lw=6, solid_capstyle="butt", zorder=0)
# ---- figura 1: misure di tessuto
arms = ["U", "N", "RV"]; off = dict(zip(arms, [-0.22, 0, 0.22]))
fig, axs = plt.subplots(1, 4, figsize=(11, 3.4), sharey=True, gridspec_kw={"wspace": 0.08})
for ax, pid, lab in [(axs[0], "P-11", "spazio vuoto (frazione a > 3 distanze tipiche)"), (axs[1], "P-13", "territori > 3 volte la mediana"),
                     (axs[2], "P-2", "CV dei territori, simulato / vero - 1"), (axs[3], "P-1", "nuclei a contatto (finestre contate a mano)")]:
    v = V[V.id == pid]
    for y, A in enumerate(TIS):
        r = v[v.archetype == A]; band(ax, y, r.lo.iloc[0], r.hi.iloc[0])
        for arm in arms:
            q = r[r.arm == arm].iloc[0]; lb, c_, mk = ST[arm]
            ax.errorbar(q.est, y + off[arm], xerr=2 * q.se, fmt=mk, color=c_, ms=4 if arm == "RV" else 3.3, lw=0.9, label=lb if y == 0 else None, zorder=3)
    ax.set_xlabel(lab, fontsize=7.5); ax.margins(x=0.06)
axs[0].set_yticks(range(6)); axs[0].set_yticklabels(list(TIS.values())); axs[0].invert_yaxis()
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, fontsize=7, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
fig.subplots_adjust(left=0.12, right=0.99, top=0.88, bottom=0.15, wspace=0.08); fig.savefig(f"{FIG}/S1.5d_fig1_tissue.png", dpi=200); plt.close(fig)
# ---- figura 2: mappe, un ROI per tessuto (H&E, vero, compartimenti di S1.5c, vincente), riquadro di 300 um
fig, axs = plt.subplots(6, 4, figsize=(10, 15.5))
for i, A in enumerate(TIS):
    _, d = cal.load(A); r = sorted(d.roi_id.unique())[0]; e = d[d.roi_id == r]; x0, y0 = e.x.min() + 100, e.y.min() + 100
    tys = MAIN[A]; cmap = {t: COL[[2, 9, 7, 3][k]] for k, t in enumerate(tys)}
    rgb = C.load_roi(A, r); upp = ROIS[(ROIS.archetype == A) & (ROIS.roi_id == r)].um_per_px.iloc[0]
    ax = axs[i, 0]; ax.imshow(rgb, extent=(0, rgb.shape[1] * upp, rgb.shape[0] * upp, 0), interpolation="lanczos")
    srcs = [("vero (nuclei segmentati)", e.x.values, e.y.values, e.type.values)]
    for arm in ["N", "RV"]:
        c_ = pd.read_parquet(f"{OUT}/cells/{arm}/{A}_{r}_s01.parquet"); srcs.append((ST[arm][0], c_.x.values, c_.y.values, c_.cell_type.values))
    axs[i, 0].set_title(f"{TIS[A]} ({r}): H&E", fontsize=7, loc="left")
    for j, (nm, x, y, t) in enumerate(srcs, start=1):
        ax = axs[i, j]; k = (x > x0) & (x < x0 + 300) & (y > y0) & (y < y0 + 300)
        ax.scatter(x[k], y[k], s=3 if A in ("A1", "A4") else 8, c=[cmap.get(tt, "0.75") for tt in t[k]], lw=0)
        ax.set_title(nm, fontsize=7, loc="left")
    for ax in axs[i]:
        ax.set_xlim(x0, x0 + 300); ax.set_ylim(y0 + 300, y0); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    axs[i, 3].legend(handles=[plt.Line2D([], [], marker="o", ls="", color=cmap[t], ms=4, label=t) for t in tys] +
                     [plt.Line2D([], [], marker="o", ls="", color="0.75", ms=4, label="altri / senza tipo")],
                     frameon=False, fontsize=6, loc="center left", bbox_to_anchor=(1.0, 0.5))
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5d_fig2_maps.png", dpi=150); plt.close(fig)
# ---- figura 3: taglie. (a) cellule quasi senza citoplasma nel tumore del colon; (b) 90 percentile dell'area per tipo
fig, axs = plt.subplots(1, 2, figsize=(11, 5.6), gridspec_kw={"width_ratios": [1, 1.5]})
P = pd.read_csv(f"{RES}/S1.5d_posthoc_p8c_real.csv"); real = P[P.source == "vero"]; ax = axs[0]
aa = ["U", "U3", "U3c", "N", "RV", "W3", "W3p"]
ax.axvspan(0, 0.0287, color="0.9", lw=0, zorder=0); ax.axvspan(real.value.min(), real.value.max(), color="0.75", lw=0, zorder=0)
for y, arm in enumerate(aa):
    q = V[(V.id == "P-8c") & (V.arm == arm)].iloc[0]; lb, c_, mk = ST[arm]
    ax.errorbar(q.est, y, xerr=2 * q.se, fmt=mk, color=c_, ms=4.5, lw=0.9)
ax.scatter(real.value, [len(aa)] * len(real), marker="|", s=120, color="k", lw=1.2)
ax.set_yticks(range(len(aa) + 1)); ax.set_yticklabels([ST[a][0] for a in aa] + ["5 ROI vere"], fontsize=7); ax.invert_yaxis()
ax.set_xlabel("frazione con nucleo ≥ 0.95 del territorio (tumore del colon)"); ax.set_xlim(0, 0.056)
ax.set_title("Intervallo pre-registrato (chiaro) e gamma delle ROI vere (scuro)", fontsize=7.5, loc="left")
ax = axs[1]; v = V[(V.id == "P-14")]; keys = v[v.arm == "RV"][["archetype", "unit"]].values.tolist()
for y, (A, t) in enumerate(keys):
    r = v[(v.archetype == A) & (v.unit == t)]; band(ax, y, r.lo.iloc[0], r.hi.iloc[0])
    for arm, o in [("RV", -0.17), ("W3", 0.17)]:
        q = r[r.arm == arm].iloc[0]; lb, c_, mk = ST[arm]
        ax.errorbar(q.est, y + o, xerr=2 * q.se, fmt=mk, color=c_, ms=3.5, lw=0.8, label=lb if y == 0 else None)
ax.set_yticks(range(len(keys))); ax.set_yticklabels([f"{TIS[A]}: {t}" for A, t in keys], fontsize=6.5); ax.invert_yaxis()
ax.set_xscale("log"); ax.set_xticks([10, 20, 50, 100, 200]); ax.set_xticklabels(["10", "20", "50", "100", "200"])
ax.set_xlabel("90° percentile dell'area del nucleo (µm²), celle interne"); ax.legend(frameon=False, fontsize=7, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2)
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5d_fig3_size.png", dpi=200); plt.close(fig)
print(sorted(os.listdir(FIG)))
