"""tools/F3.3_figures.py - F3.3: figure del report (sfondo bianco; taratura in grigio, aree nuove della stessa sezione in verde, altre sezioni o
pazienti in arancio e viola, simulato in blu). Uso: .venv/bin/python tools/F3.3_figures.py"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "reports/F3.3_parts"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.dpi": 200})
GREY, BLUE, GREEN = "#7f7f7f", "#1f5fa8", "#2e8b57"
L2C = {"CRC_P1_SR300": ("#d9822b", "o", "paziente 1"), "CRC_P5_SR300": ("#b2451f", "s", "paziente 5"), "CRC_6p5mm_SR410": ("#8a5bb5", "^", "tumore 6,5 mm"),
       "MouseBrain_FFPE_SR401": ("#8a5bb5", "D", "cervello FFPE")}
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; TN = ["intestino", "tumore\ndel colon", "stroma", "linfonodo\n(500 µm)", "cervello", "cuore"]
MEAS = [("M1", "UMI assegnati alle cellule (quota)"), ("M2", "UMI per cellula"), ("M3", "UMI su marcatori di altri tipi (quota)"),
        ("M4", "blocchi chiamati singoletto (quota)"), ("M5", "tipo di RCTD = tipo del nucleo (quota)"), ("M6", "blocchi vicini con lo stesso tipo (quota)")]
U = pd.read_csv("results/F3.3/F3.3_area_values.csv"); U = U[U.src != "IN_panel"]
SIM = U[U.src == "SIM"].groupby(["A", "u2"])[[m for m, _ in MEAS]].mean().reset_index()
rng = np.random.default_rng(0)
fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.4))
for ax, (m, lab) in zip(axes.ravel(), MEAS):
    for i, A in enumerate(ARCH):
        r = U[(U.A == A) & (U.src == "IN")][m]; s = SIM[SIM.A == A][m]
        ax.scatter(i - 0.3 + rng.uniform(-0.05, 0.05, len(r)), r, s=8, color=GREY, zorder=3, label="taratura" if i == 0 else None)
        o1 = U[(U.A == A) & (U.src == "OUT") & (U.level == "L1")][m]
        ax.scatter(i - 0.1 + rng.uniform(-0.05, 0.05, len(o1)), o1, s=9, color=GREEN, zorder=3, label="aree nuove, stessa sezione" if i == 0 else None)
        for k, (sec, (col, mk, nm)) in enumerate(L2C.items()):
            o2 = U[(U.A == A) & (U.src == "OUT") & (U.set.str.startswith(sec))][m]
            if len(o2): ax.scatter(i + 0.1 + rng.uniform(-0.05, 0.05, len(o2)), o2, s=10, color=col, marker=mk, zorder=3, label=nm if A in ("A2", "A5") and m == "M1" else None)
        ax.scatter(i + 0.3 + rng.uniform(-0.05, 0.05, len(s)), s, s=8, color=BLUE, zorder=3, label="simulato" if i == 0 else None)
    ax.set_xticks(range(6)); ax.set_xticklabels(TN); ax.set_title(lab, loc="left"); ax.margins(y=0.1)
    if m == "M2": ax.set_yscale("log")
h, l = [], []
for ax in axes.ravel():
    for hh, ll in zip(*ax.get_legend_handles_labels()):
        if ll not in l: h.append(hh); l.append(ll)
fig.legend(h, l, loc="lower center", ncol=7, frameon=False); fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(f"{OUT}/F3.3_fig1_values.png"); plt.close(fig)
# figura 2: spostamento del vero, fuori meno dentro, in unita' di h
C = pd.read_csv("results/F3.3/F3.3_cells.csv")
NAME = {"A1_L1_A1": "intestino, aree nuove", "A23_L1_A2": "tumore, aree nuove", "A4_L1_A4": "linfonodo, aree nuove", "A5_L1_A5": "cervello, aree nuove",
        "A6_L1_A6": "cuore, aree nuove", "CRC_P1_SR300_L2_A2": "tumore, paziente 1", "CRC_P5_SR300_L2_A2": "tumore, paziente 5", "CRC_6p5mm_SR410_L2_A2": "tumore, 6,5 mm",
        "CRC_P1_SR300_L2_A3": "stroma, paziente 1", "CRC_P5_SR300_L2_A3": "stroma, paziente 5", "CRC_6p5mm_SR410_L2_A3": "stroma, 6,5 mm",
        "MouseBrain_FFPE_SR401_L2_A5": "cervello, sezione FFPE"}
order = [k for k in NAME if k in set(C.set)]
P = C.pivot_table(index="set", columns="measure", values="shift_over_h").loc[order]
Vd = C.pivot_table(index="set", columns="measure", values="verdict_shift", aggfunc="first").loc[order]
fig, ax = plt.subplots(figsize=(6.6, 5.0)); lim = 5
im = ax.imshow(P.clip(-lim, lim).to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        v = P.iat[i, j]; ax.text(j, i, f"{v:+.1f}" + ("*" if Vd.iat[i, j] == "smentita" else ""), ha="center", va="center", fontsize=6.5, color="black")
ax.set_xticks(range(6)); ax.set_xticklabels(["UMI assegnati", "UMI per cellula", "marcatori estranei", "singoletti", "tipo = nucleo", "vicini stesso tipo"], rotation=30, ha="right")
ax.set_yticks(range(len(order))); ax.set_yticklabels([NAME[k] for k in order]); ax.axhline(4.5, color="black", lw=0.8)
cb = fig.colorbar(im, ax=ax, shrink=0.8); cb.set_label("vero fuori meno vero dentro, in margini h")
fig.tight_layout(); fig.savefig(f"{OUT}/F3.3_fig2_shift.png"); plt.close(fig)
# figura 3: stime aggregate
S = pd.read_csv("results/F3.3/F3.3_summary.csv")
rows = [("P-1", "L1", "spostamento del vero\naree nuove"), ("P-3", "L1", "degrado del simulato\naree nuove"), ("P-2", "L2", "spostamento del vero\naltre sezioni o pazienti"),
        ("P-3", "L2", "degrado del simulato\naltre sezioni o pazienti")]
fig, ax = plt.subplots(figsize=(6.0, 3.0))
for k, (pid, lev, lab) in enumerate(rows):
    r = S[(S.id == pid) & (S.level == lev)].iloc[0]
    ax.plot([r.lo, r.hi], [k, k], color=BLUE, lw=2.5, solid_capstyle="butt"); ax.scatter([r.est], [k], color=BLUE, zorder=3)
    if np.isfinite(r.null_lo): ax.plot([r.null_lo, r.null_hi], [k - 0.22, k - 0.22], color=GREY, lw=1.5)
ax.axvspan(-0.5, 0.5, color=GREY, alpha=0.15, lw=0); ax.axvline(0, color=GREY, lw=0.6)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[2] for r in rows]); ax.invert_yaxis()
ax.set_xlabel("margini h (radice con segno del quadrato medio); grigio: senza alcuno spostamento")
fig.tight_layout(); fig.savefig(f"{OUT}/F3.3_fig3_aggregate.png"); plt.close(fig)
print("ok")
