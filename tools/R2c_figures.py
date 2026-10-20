#!/usr/bin/env python
"""tools/R2c_figures.py — R2c: figure del report (sfondo bianco, un pannello per archetipo).
Fig. 1 results/R2c/fig/R2c_s1.png: S-1 manuale (Luca; IC a grappolo; intervallo delle varianti unsure; valori per finestra),
  cellpose, StarDist, Space Ranger (A4/A6), K-2 (mediana e intervallo dei 10 semi), Federica (sue finestre).
Fig. 2 results/R2c/fig/R2c_distance.png: frazione 'touch' (Luca, coppie 'pair', unsure esclusi) per classi di distanza, con
  frazione di contatto di cellpose e StarDist sulle coppie appaiate in quel segmentatore; linea d_max."""
import importlib.util, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/R2c"; (RES / "fig").mkdir(exist_ok=True)
spec = importlib.util.spec_from_file_location("ana", ROOT / "tools/R2c_analysis.py"); A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
s1 = pd.read_csv(RES / "R2c_s1.csv").set_index("archetype"); sF = pd.read_csv(RES / "R2c_s1_federica.csv")
sF = sF[sF.set == "Federica"].set_index("archetype")
k2 = pd.read_csv(RES / "R2c_k2_windows.csv"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype")
foc = pd.read_csv(RES / "R2c_focal.csv"); foc = foc[foc.focal]
L = A.load_answers("Luca"); P = L[L.kind == "pair"][["win_id", "i", "j", "answer", "d_um", "pair_id"]]
# ---------------------------------------------------------------- Fig. 1
fig, axs = plt.subplots(1, 6, figsize=(17, 4.6), facecolor="white", sharey=False)
for ax, a in zip(axs, ARCH):
    r = s1.loc[a]; x = {"manuale": 0, "Federica": 1, "cellpose": 2, "StarDist": 3, "Space Ranger": 4, "K-2": 5}
    ax.axhspan(r.lo_primary, r.hi_primary, color="#90caf9", alpha=0.35, lw=0)
    ax.errorbar(0, r.s1_primary, yerr=[[r.s1_primary - r.lo_primary], [r.hi_primary - r.s1_primary]], fmt="o", color="#0d47a1", capsize=4, ms=7, zorder=4)
    ax.plot([0.18, 0.18], [r.s1_unsure_sep, r.s1_unsure_touch], color="#0d47a1", lw=3, alpha=0.5)
    f = foc[foc.archetype == a]
    for w, g in f.groupby("win_id"):
        y = A.focal_s1(P[P.win_id == w], g, "primary"); ax.plot(-0.18, np.nanmean(y), "_", color="#1565c0", ms=10, mew=1.5)
    if a in sF.index:
        q = sF.loc[a]; ax.errorbar(1, q.s1_primary, yerr=[[q.s1_primary - q.lo_primary], [q.hi_primary - q.s1_primary]], fmt="s", color="#6a1b9a", capsize=3)
    ax.plot(2, r.seg_cellpose_rgb, "D", color="#e65100", ms=7); ax.plot(3, r.seg_stardist_he, "D", color="#2e7d32", ms=7)
    if np.isfinite(r.get("seg_spaceranger", np.nan)):
        ax.plot(4, r.seg_spaceranger, "D", color="#616161", ms=7)
    kk = k2[k2.archetype == a].groupby("seed")[["n", "n_contact"]].sum(); v = kk.n_contact / kk.n
    ax.errorbar(5, np.median(v), yerr=[[np.median(v) - v.min()], [v.max() - np.median(v)]], fmt="^", color="#c62828", capsize=3, ms=8)
    ax.set_xticks(list(x.values())); ax.set_xticklabels(list(x.keys()), rotation=45, ha="right", fontsize=9)
    ax.set_title(f"{a}  B-dec: {r.b_dec} · K-2 {r.k2_verdict}", fontsize=10); ax.set_xlim(-0.6, 5.6)
    ax.set_ylim(-0.03, max(1.0, ax.get_ylim()[1])); ax.grid(axis="y", alpha=0.3)
axs[0].set_ylabel("S-1: frazione di nuclei a contatto")
fig.suptitle("R2c — S-1 manuale (Luca: punto e IC 95 % a grappolo, banda azzurra; barra = varianti 'non giudicabile'; trattini = finestre) "
             "contro Federica, segmentatori sugli stessi focali appaiati e K-2 (mediana, min–max dei 10 semi)", fontsize=10)
fig.tight_layout(); fig.savefig(RES / "fig/R2c_s1.png", dpi=140, facecolor="white"); fig.savefig(RES / "fig/R2c_s1.pdf", facecolor="white")
# ---------------------------------------------------------------- Fig. 2
pr = pd.read_csv(RES / "R2c_pairs.csv")[["pair_id", "both_cellpose_rgb", "both_stardist_he", "contact_cellpose_rgb", "contact_stardist_he"]]
Q = P[P.answer != "unsure"].merge(pr, on="pair_id", how="left")
fig, axs = plt.subplots(1, 6, figsize=(17, 3.9), facecolor="white")
for ax, a in zip(axs, ARCH):
    g = Q[Q.win_id.str.startswith(a + "_")]; d = float(dm.d_max_um[a])
    bins = np.linspace(0, d, 7); mid = (bins[:-1] + bins[1:]) / 2
    def frac(m, y):
        out, n = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            s = (g.d_um > lo) & (g.d_um <= hi) & m; out.append(y[s].mean() if s.sum() else np.nan); n.append(int(s.sum()))
        return np.array(out), n
    t = (g.answer == "touch").astype(float)
    fm, nm = frac(np.ones(len(g), bool), t)
    ax.plot(mid, fm, "o-", color="#0d47a1", label="manuale (Luca)")
    for seg, col, lab in [("cellpose_rgb", "#e65100", "cellpose"), ("stardist_he", "#2e7d32", "StarDist")]:
        m = (g[f"both_{seg}"] == True).to_numpy(); fs, _ = frac(m, (g[f"contact_{seg}"] == True).astype(float))
        ax.plot(mid, fs, "D--", color=col, ms=4, label=lab + " (appaiate)")
    for xx, nn in zip(mid, nm):
        ax.text(xx, 1.02, str(nn), ha="center", fontsize=7, color="#555")
    ax.axvline(d, color="k", ls=":", lw=1); ax.set_ylim(-0.03, 1.08); ax.set_xlim(0, d * 1.03)
    ax.set_title(f"{a} (d_max {d:.1f} µm)", fontsize=10); ax.set_xlabel("distanza fra i punti (µm)"); ax.grid(alpha=0.3)
axs[0].set_ylabel("frazione a contatto"); axs[0].legend(fontsize=7, loc="center right")
fig.suptitle("R2c — contatto in funzione della distanza fra i centri (coppie 'pair'; numeri = coppie manuali per classe)", fontsize=10)
fig.tight_layout(); fig.savefig(RES / "fig/R2c_distance.png", dpi=140, facecolor="white"); fig.savefig(RES / "fig/R2c_distance.pdf", facecolor="white")
print("ok")
