"""tools/S1.5c_figures.py - S1.5c: figure del report (sfondo bianco, un pannello per tessuto). Ingressi: results/S1.5c/*.csv, celle del banco (seme 1),
nuclei veri (etichette R6 PH-1). Uso: .venv/bin/python tools/S1.5c_figures.py <braccio vincente>   Uscite: results/S1.5c/figures/*.png"""
import sys, os, importlib.util, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.ticker
s = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(s); s.loader.exec_module(cal)
WIN = sys.argv[1] if len(sys.argv) > 1 else "N-a"; RES = "results/S1.5c"; FIG = f"{RES}/figures"; os.makedirs(FIG, exist_ok=True)
OUT = "/mnt/micron/geo_spatialtrans/S1.5c"
TIS = {"A1": "intestino di topo", "A2": "tumore del colon", "A3": "stroma del tumore", "A4": "linfonodo", "A5": "cervello di topo", "A6": "cuore"}
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False})
COL = plt.rcParams["axes.prop_cycle"].by_key()["color"]
UM = pd.read_csv(f"{RES}/S1.5c_unit_metrics.csv", keep_default_na=False, na_values=[""]); UM[["cell_type", "bin"]] = UM[["cell_type", "bin"]].fillna("")
RR = pd.read_csv(f"{RES}/S1.5c_real_refs.csv", keep_default_na=False, na_values=[""]); RR[["cell_type", "bin"]] = RR[["cell_type", "bin"]].fillna("")
V = pd.read_csv(f"{RES}/S1.5c_verdicts.csv")
mid = {"0-10": 5, "10-20": 15, "20-40": 30, "40-80": 60, "80-160": 120}
ARML = {"U": "semina uniforme (oggi)", WIN: "compartimenti per tipo", "G": "stesse posizioni, tipi a caso"}
# ---- figura 1: rapporto di connessione dei tipi, reale (banda fra ROI) contro bracci
fig, axs = plt.subplots(2, 3, figsize=(10, 6.2), sharex=True)
for ax, A in zip(axs.ravel(), TIS):
    for k, t in enumerate(MAIN[A]):
        r = RR[(RR.archetype == A) & (RR.metric == "markconn") & (RR.cell_type == t)].copy(); r["m"] = r.bin.map(mid)
        g = r.groupby("m").value; ax.fill_between(g.min().index, g.min(), g.max(), color=COL[k], alpha=0.18, lw=0)
        for arm, ls, mk in [(WIN, "-", "o"), ("U", ":", None)]:
            x = UM[(UM.arm == arm) & (UM.archetype == A) & (UM.metric == "markconn") & (UM.cell_type == t)].copy(); x["m"] = x.bin.map(mid)
            y = x.groupby(["m", "roi_id"]).value.mean().groupby("m").mean()
            ax.plot(y.index, y.values, ls=ls, color=COL[k], marker=mk, ms=2.5, lw=1.2, label=t if arm == WIN else None)
    ax.axhline(1, color="grey", lw=0.5); ax.set_xscale("log"); ax.set_xticks(list(mid.values())); ax.set_xticklabels([str(v) for v in mid.values()])
    ax.set_title(TIS[A], loc="left"); ax.legend(frameon=False, fontsize=6.5, loc="upper right")
for ax in axs[1]: ax.set_xlabel("distanza fra le cellule (µm)")
for ax in axs[:, 0]: ax.set_ylabel("vicini dello stesso tipo / attesi a caso")
fig.text(0.01, 0.005, "Bande: minimo-massimo delle 5 aree vere. Linee piene: compartimenti per tipo (media delle 5 aree, 10 semi). Puntinate: semina uniforme di oggi.", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f"{FIG}/S1.5c_fig1_markconn.png", dpi=180); plt.close(fig)
# ---- figura 2: misure per tipo (distanza relativa dal vicino, rapporto nucleo/territorio): intervallo reale contro bracci
fig, axs = plt.subplots(1, 2, figsize=(10, 6.5))
for ax, pid, lab in [(axs[0], "P-5", "distanza dal vicino del tipo / distanza tipica del tessuto"), (axs[1], "P-4", "rapporto nucleo / territorio del tipo")]:
    v = V[V.id == pid]; keys = v[v.arm == "U"][["archetype", "unit"]].values.tolist()
    for y, (A, t) in enumerate(keys):
        r = v[(v.archetype == A) & (v.unit == t)]; lo, hi = r.lo.iloc[0], r.hi.iloc[0]
        ax.plot([lo, hi], [y, y], color="0.75", lw=5, solid_capstyle="butt")
        for arm, c_, mk in [("U", "0.35", "x"), ("G", COL[2], "s"), (WIN, COL[3], "o")]:
            q = r[r.arm == arm].iloc[0]; ax.errorbar(q.est, y + {"U": -0.22, "G": 0, WIN: 0.22}[arm], xerr=2 * q.se, fmt=mk, color=c_, ms=3.5, lw=0.8,
                                                      label=ARML[arm] if y == 0 else None)
    ax.set_yticks(range(len(keys))); ax.set_yticklabels([f"{TIS[A]}: {t}" for A, t in keys], fontsize=6.5); ax.invert_yaxis(); ax.set_xlabel(lab)
    if pid == "P-4":
        ax.set_xscale("log"); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5]); ax.set_xticklabels(["0.02", "0.05", "0.1", "0.2", "0.5"])
axs[0].legend(frameon=False, fontsize=7, loc="lower right"); axs[1].set_yticklabels([])
fig.text(0.01, 0.005, "Barre grigie: minimo-massimo delle 5 aree vere. Punti: stima ± 2 errori standard (5 aree, 10 semi ciascuna).", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f"{FIG}/S1.5c_fig2_by_type.png", dpi=180); plt.close(fig)
# ---- figura 3: misure di tessuto: variabilita' dei territori, numero di vicini, contatti
fig, axs = plt.subplots(1, 3, figsize=(10, 3.2))
for ax, pid, lab in [(axs[0], "P-2", "CV dei territori, simulato / vero - 1"), (axs[1], "P-7", "CV del numero di vicini (4 distanze tipiche)"), (axs[2], "P-1", "nuclei a contatto (finestre contate a mano)")]:
    v = V[V.id == pid]
    for y, A in enumerate(TIS):
        r = v[v.archetype == A]; lo, hi = r.lo.iloc[0], r.hi.iloc[0]; ax.plot([lo, hi], [y, y], color="0.75", lw=5, solid_capstyle="butt")
        for arm, c_, mk in [("U", "0.35", "x"), (WIN, COL[3], "o")]:
            q = r[r.arm == arm].iloc[0]; ax.errorbar(q.est, y + (-0.18 if arm == "U" else 0.18), xerr=2 * q.se, fmt=mk, color=c_, ms=3.5, lw=0.8, label=ARML[arm] if y == 0 else None)
    ax.set_yticks(range(6)); ax.set_yticklabels(list(TIS.values()) if pid == "P-2" else []); ax.invert_yaxis(); ax.set_xlabel(lab, fontsize=7.5)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=7, loc="lower center", ncol=len(l))
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(f"{FIG}/S1.5c_fig3_tissue.png", dpi=180); plt.close(fig)
# ---- figura 4: mappe dei tipi, un ROI per tessuto (vero, semina di oggi, compartimenti), riquadro di 300 um
fig, axs = plt.subplots(6, 3, figsize=(8.5, 17))
for i, A in enumerate(TIS):
    lab, d = cal.load(A); r = sorted(d.roi_id.unique())[0]; e = d[d.roi_id == r]
    x0, y0 = e.x.min() + 100, e.y.min() + 100
    tys = MAIN[A]; cmap = {t: COL[k] for k, t in enumerate(tys)}
    srcs = [("vero", e.x.values, e.y.values, e.type.values)]
    for arm in ["U", WIN]:
        c_ = pd.read_parquet(f"{OUT}/cells/{arm}/{A}_{r}_s01.parquet"); srcs.append((ARML[arm], c_.x.values, c_.y.values, c_.cell_type.values))
    for j, (nm, x, y, t) in enumerate(srcs):
        ax = axs[i, j]; k = (x > x0) & (x < x0 + 300) & (y > y0) & (y < y0 + 300)
        col = [cmap.get(tt, "0.8") for tt in t[k]]; ax.scatter(x[k], y[k], s=3 if A in ("A1", "A4") else 8, c=col, lw=0)
        ax.set_xlim(x0, x0 + 300); ax.set_ylim(y0, y0 + 300); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{TIS[A]} ({r}): {nm}", fontsize=7, loc="left")
    axs[i, 2].legend(handles=[plt.Line2D([], [], marker="o", ls="", color=cmap[t], ms=4, label=t) for t in tys] + [plt.Line2D([], [], marker="o", ls="", color="0.8", ms=4, label="altri / senza tipo")],
                     frameon=False, fontsize=6, loc="center left", bbox_to_anchor=(1.0, 0.5))
fig.tight_layout(); fig.savefig(f"{FIG}/S1.5c_fig4_maps.png", dpi=150); plt.close(fig)
print("figure scritte")
