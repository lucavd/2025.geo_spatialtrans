"""tools/S1.5e_figures.py - S1.5e: figure del report (sfondo bianco, un pannello per tessuto dove pertinente). Ingressi: results/S1.5e/*.csv, celle del
banco (seme 1), nuclei veri (R6), ROI H&E (R2), esplorazione E3. Uso: .venv/bin/python -B tools/S1.5e_figures.py   Uscite: results/S1.5e/figures/*.png"""
import sys, os, importlib.util, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "tools"); import R2_common as C
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py"); m5e = mod("m5e", "tools/S1.5e_metrics.py")
RES = "results/S1.5e"; FIG = f"{RES}/figures"; os.makedirs(FIG, exist_ok=True); OUT = "/mnt/micron/geo_spatialtrans/S1.5e"
TIS = {"A1": "intestino di topo", "A2": "tumore del colon", "A3": "stroma del tumore", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False})
COL = plt.rcParams["axes.prop_cycle"].by_key()["color"]
ST = {"U": ("semina uniforme (predefinito)", "0.45", "x"), "RV": ("vuoti che spostano le cellule (S1.5d)", COL[0], "s"),
      "X": ("vuoti tondi che tolgono cellule", COL[2], "^"), "XA": ("vuoti allungati che tolgono cellule", COL[1], "o"),
      "XAD": ("vuoti allungati + densità con tetto", COL[4], "D"), "Q": ("vincente senza quote", COL[8], "v")}
V = pd.read_csv(f"{RES}/S1.5e_verdicts.csv"); ROIS = pd.read_csv("results/R2/R2_rois_checked.csv"); WIN = pd.read_csv(f"{RES}/S1.5e_ruleA.csv").winner.iloc[0]
ARMS = [a for a in ["U", "RV", "X", "XA", "XAD"] if a in set(V.arm)]
def band(ax, y, lo, hi): ax.plot([lo, hi], [y, y], color="0.82", lw=6, solid_capstyle="butt", zorder=0)
def panel(fig_name, items, arms):
    off = dict(zip(arms, np.linspace(-0.3, 0.3, len(arms))))
    fig, axs = plt.subplots(1, len(items), figsize=(2.9 * len(items), 3.5), sharey=True)
    for ax, (pid, lab) in zip(axs, items):
        v = V[(V.id == pid) & (V.unit == "tessuto")]
        for y, A in enumerate(TIS):
            r = v[v.archetype == A]
            if not len(r): continue
            band(ax, y, r.lo.iloc[0], r.hi.iloc[0])
            for arm in arms:
                q = r[r.arm == arm]
                if not len(q): continue
                q = q.iloc[0]; lb, c_, mk = ST[arm]
                ax.errorbar(q.est, y + off[arm], xerr=2 * q.se, fmt=mk, color=c_, ms=4.2 if arm == WIN else 3.2, lw=0.9, label=lb if y == 0 else None, zorder=3)
        ax.set_xlabel(lab, fontsize=7.5); ax.margins(x=0.06)
    axs[0].set_yticks(range(6)); axs[0].set_yticklabels(list(TIS.values())); axs[0].invert_yaxis()
    h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, fontsize=7, loc="upper center", ncol=len(arms), bbox_to_anchor=(0.5, 1.0))
    fig.subplots_adjust(left=0.1, right=0.99, top=0.86, bottom=0.16, wspace=0.08); fig.savefig(f"{FIG}/{fig_name}", dpi=200); plt.close(fig)
panel("S1.5e_fig1_contacts_density.png", [("P-1", "nuclei a contatto (finestre contate a mano)"), ("P-15", "densità del ROI (nuclei / mm² di tessuto)"),
      ("P-2", "CV dei territori, simulato / vero - 1"), ("P-7", "CV del numero di vicini")], ARMS)
panel("S1.5e_fig2_voids.png", [("P-11", "spazio vuoto (bersaglio)"), ("P-13", "territori > 3 volte la mediana"),
      ("P-16", "quota dell'area vuota in vuoti allungati (> 3)"), ("P-17", "numero di vuoti per ROI")], ARMS)
# ---- figura 3: contatto per classe di distanza dal vicino: manuale (E3), RV, vincente
E3 = pd.read_csv(f"{RES}/S1.5e_explore_contacts_summary.csv"); WN = pd.read_csv(f"{RES}/S1.5e_windows_nn.csv")
lab = ["< 5", "5-6", "6-7", "7-8", "8-10", "> 10"]
fig, axs = plt.subplots(2, 3, figsize=(10, 5.6), sharey=True); axs = axs.ravel()
for ax, A in zip(axs, TIS):
    m = E3[(E3.archetype == A) & (E3.source == "manual")].iloc[0]
    ax.plot(range(6), [m[f"s1_b{k}"] for k in range(6)], "k-o", ms=4, lw=1.2, label="vero (risposte di Luca, 60 focali)")
    for arm in dict.fromkeys(["RV", WIN]):
        z = WN[(WN.arm == arm) & (WN.archetype == A)].groupby("bin")[["n", "n_contact"]].sum()
        lb, c_, mk = ST[arm]; ax.plot(z.index, z.n_contact / z.n.where(z.n > 0), marker=mk, color=c_, ms=3.5, lw=1, label=lb)
    ax2 = ax.twinx(); ax2.bar(np.arange(6) - 0.2, [m[f"frac_b{k}"] for k in range(6)], width=0.35, color="0.85", zorder=0)
    z = WN[(WN.arm == WIN) & (WN.archetype == A)].groupby("bin").n.sum(); ax2.bar(z.index + 0.2, z / z.sum(), width=0.35, color=ST[WIN][1], alpha=0.25, zorder=0)
    ax2.set_ylim(0, 1.6); ax2.set_yticks([]); ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    ax.set_title(TIS[A], fontsize=8, loc="left"); ax.set_xticks(range(6)); ax.set_xticklabels(lab, fontsize=7); ax.set_ylim(-0.02, 1.02)
for ax in axs[3:]: ax.set_xlabel("distanza dal vicino più vicino (µm)")
axs[0].set_ylabel("frazione a contatto"); axs[3].set_ylabel("frazione a contatto")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, fontsize=7, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(f"{FIG}/S1.5e_fig3_contact_by_distance.png", dpi=200); plt.close(fig)
# ---- figura 4: mappe con i vuoti, un ROI per tessuto (H&E, vero, RV, vincente), riquadro di 400 um; pixel vuoti in grigio
fig, axs = plt.subplots(6, 4, figsize=(10, 15.5))
for i, A in enumerate(TIS):
    _, d = cal.load(A); r = sorted(d.roi_id.unique())[0]; e = d[d.roi_id == r]; x0, y0 = e.x.min() + 100, e.y.min() + 100; g = m5e.RoiGeomE(A, r)
    rgb = C.load_roi(A, r); upp = ROIS[(ROIS.archetype == A) & (ROIS.roi_id == r)].um_per_px.iloc[0]
    axs[i, 0].imshow(rgb, extent=(0, rgb.shape[1] * upp, rgb.shape[0] * upp, 0), interpolation="lanczos"); axs[i, 0].set_title(f"{TIS[A]} ({r}): H&E", fontsize=7, loc="left")
    srcs = [("vero (nuclei segmentati)", e.x.values, e.y.values, e.interior.values)]
    for arm in dict.fromkeys(["RV", WIN]):
        c_ = pd.read_parquet(f"{OUT}/cells/{arm}/{A}_{r}_s01.parquet"); srcs.append((ST[arm][0], c_.x.values, c_.y.values, c_.interior.values))
    for j, (nm, x, y, it) in enumerate(srcs, start=1):
        ax = axs[i, j]; v = m5e.void_image(x, y, it.astype(bool), g)
        ax.imshow(np.where(v, 0.55, np.where(g.m, 1.0, 0.85)), cmap="gray", vmin=0, vmax=1, extent=(0, g.m.shape[1] * g.px, g.m.shape[0] * g.px, 0), interpolation="nearest")
        k = (x > x0) & (x < x0 + 400) & (y > y0) & (y < y0 + 400); ax.scatter(x[k], y[k], s=2 if A in ("A1", "A4") else 5, c="k", lw=0)
        ax.set_title(nm, fontsize=7, loc="left")
    for ax in axs[i]: ax.set_xlim(x0, x0 + 400); ax.set_ylim(y0 + 400, y0); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5e_fig4_maps.png", dpi=150); plt.close(fig)
print(sorted(os.listdir(FIG)))
