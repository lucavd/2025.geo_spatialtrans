"""tools/R5_figures.py — R5, figure del report (sfondo bianco, un pannello per sezione o tessuto). Uso: .venv/bin/python tools/R5_figures.py <outdir>
Scrive results/R5/figures/R5_F*.png."""
import sys, os, json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
O = sys.argv[1]; FD = "results/R5/figures"; os.makedirs(FD, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False})
KEYS = ["A1", "A23", "A4", "A5", "A6"]
NAME = {"A1": "mouse intestine", "A23": "colon cancer (tumour and stroma)", "A4": "lymph node", "A5": "mouse brain", "A6": "heart"}
TIS = {"A1": "mouse intestine", "A2": "colon tumour", "A3": "tumour stroma", "A4": "lymph node", "A5": "mouse brain", "A6": "heart"}
COL = {"obs": "black", "K1": "#999999", "K2": "#0072B2", "K3": "#E69F00"}
T = pd.read_csv(f"{O}/tables/R5_tiles.csv")

# F1: profilo pooled per sezione, osservato e previsto
fig, axs = plt.subplots(1, 5, figsize=(12, 2.8), sharey=True)
for ax, key in zip(axs, KEYS):
    acc = None
    for f in sorted(glob.glob(f"{O}/tiles/{key}/t*.json")):
        d = json.load(open(f))
        if "error" in d: continue
        p = d["profile"]; e = np.array(p["edges"])
        arr = {k: np.array(p[k], float) for k in ("n", "y", "pred_K1", "pred_K2", "pred_K3")}
        acc = arr if acc is None else {k: acc[k] + arr[k] for k in acc}
    x = np.r_[e[0] - 2, e] + 1.0                       # centro della fascia (digitize: fascia i = [e[i-1], e[i]))
    n = np.maximum(acc["n"], 1)
    sel_ref = (x >= -10) & (x < -6)          # tessuto entro 10 µm dal bordo: e' la sola parte di tessuto nella verosimiglianza
    ref = acc["y"][sel_ref].sum() / n[sel_ref].sum()
    ok = acc["n"] > 50
    ax.semilogy(x[ok], acc["y"][ok] / n[ok] / ref, "o", ms=2, color=COL["obs"], label="observed")
    for m in ("K1", "K2", "K3"):
        ax.semilogy(x[ok], acc[f"pred_{m}"][ok] / n[ok] / ref, "-", lw=1.0 if m != "K2" else 1.6, color=COL[m],
                    label={"K1": "Gaussian", "K2": "Gaussian + tail", "K3": "anisotropic + tail"}[m])
    ax.axvline(0, color="k", lw=0.5); ax.set_title(f"{NAME[key].replace(' (tumour and stroma)', '')}, {int((T.key == key).sum())} tiles", loc="left")
axs[0].set_ylabel("UMI per 2 µm bin / tissue level"); axs[0].legend(frameon=False, loc="upper right")
fig.supxlabel("signed distance from the H&E tissue edge (µm; negative = inside tissue)", fontsize=8)
fig.tight_layout(); fig.savefig(f"{FD}/R5_F1_edge_profiles.png", dpi=200)

# F2: stime per tessera
P = pd.read_csv(f"{O}/tables/R5_predictions_edges.csv")
fig, axs = plt.subplots(1, 3, figsize=(11, 3.0))
for i, key in enumerate(KEYS):
    t = T[T.key == key]; jit = np.random.default_rng(i).uniform(-0.15, 0.15, len(t))
    axs[0].plot(i + jit, t.K2_sigma, "o", ms=2, color="#7f7f7f", alpha=0.6)
    axs[1].plot(i + jit, t.K2_w, "o", ms=2, color="#7f7f7f", alpha=0.6)
    for ax, pid in ((axs[0], "P-1"), (axs[1], "P-2")):
        r = P[(P.key == key) & (P.prediction == pid)].iloc[0]
        ax.errorbar(i + 0.3, r.estimate, yerr=[[r.estimate - r.lo], [r.hi - r.estimate]], fmt="s", ms=4, color="#0072B2", capsize=2)
    axs[2].plot(t.K2_dc, t.K2_dr, "o", ms=2, alpha=0.5, label=NAME[key])
axs[0].axhspan(1, 4, color="#0072B2", alpha=0.08); axs[0].set_ylabel("Gaussian width σ (µm)"); axs[0].set_title("Diffusion width per tile, median and 95 % interval", loc="left")
axs[1].axhspan(0, 0.01, color="#999999", alpha=0.15); axs[1].set_ylabel("tail mass w"); axs[1].set_title("Tail mass per tile (grey band: negligible)", loc="left")
for ax in axs[:2]:
    ax.set_xticks(range(5)); ax.set_xticklabels([NAME[k].replace(" (tumour and stroma)", "") for k in KEYS], rotation=30, ha="right")
axs[2].axhline(0, color="k", lw=0.5); axs[2].axvline(0, color="k", lw=0.5)
axs[2].set_xlabel("shift along columns (µm)"); axs[2].set_ylabel("shift along rows (µm)"); axs[2].set_title("Image-to-grid shift per tile", loc="left")
axs[2].legend(frameon=False, fontsize=6, markerscale=2)
fig.tight_layout(); fig.savefig(f"{FD}/R5_F2_tile_estimates.png", dpi=200)

# F3: rapporto nucleo/corona: osservato (R4), prima della diffusione con il nucleo dei bordi (sensibilita' sulle ROI riuscite; come scritto
# non calcolabile fuori dal cervello) e con la sfocatura interna di PH-2 (per ROI)
if os.path.exists(f"{O}/tables/R5_predictions_nc.csv"):
    Q = pd.read_csv(f"{O}/tables/R5_predictions_nc.csv")
    A = ["A1", "A2", "A3", "A4", "A5", "A6"]; x = np.arange(6)
    sv = Q[(Q.prediction == "P-5") & Q.variant.str.startswith("sensitivity")].set_index("archetype").reindex(A)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
    axs[0].plot(x - 0.18, sv.r_obs_mean, "o", color="black", label="observed (R4, after diffusion)")
    ok = sv.estimate.notna()
    axs[0].errorbar((x + 0.0)[ok], sv.estimate[ok], yerr=[(sv.estimate - sv.lo)[ok].clip(lower=0), (sv.hi - sv.estimate)[ok].clip(lower=0)], fmt="s", color="#0072B2",
                    capsize=2, label="before diffusion, edge kernel (successful ROIs only)")
    ph = f"{O}/tables/R5_PH2_inner_summary.csv"
    if os.path.exists(ph):
        H = pd.read_csv(ph); H["a"] = H.tag.str[:2]
        for i, a in enumerate(A):
            h = H[(H.a == a) & H.r_at_star.notna()]
            axs[0].plot(np.full(len(h), i + 0.18), h.r_at_star, "D", ms=3, color="#D55E00", label="before diffusion, inner blur (PH-2, per ROI)" if i == 0 else None)
    axs[0].axhspan(0.909, 1.1, color="#999999", alpha=0.15); axs[0].set_yscale("log")
    axs[0].set_xticks(x); axs[0].set_xticklabels([TIS[a] for a in A], rotation=30, ha="right"); axs[0].set_ylabel("nucleus / ring density ratio")
    axs[0].set_title("Nucleus vs cytoplasm density (grey: uniform)", loc="left"); axs[0].legend(frameon=False, fontsize=6)
    cp = Q[Q.prediction.str.startswith("CP-1") & Q.variant.str.startswith("sensitivity")].set_index("archetype").reindex(A)
    ok = cp.estimate.notna()
    axs[1].errorbar(x[ok], cp.estimate[ok], yerr=[(cp.estimate - cp.lo)[ok].clip(lower=0), (cp.hi - cp.estimate)[ok].clip(lower=0)], fmt="s", color="#0072B2", capsize=2)
    axs[1].axhspan(0.909, 1.1, color="#0072B2", alpha=0.08); axs[1].axhline(1, color="k", lw=0.5)
    axs[1].set_xticks(x); axs[1].set_xticklabels([TIS[a] for a in A], rotation=30, ha="right")
    axs[1].set_ylabel("predicted / observed ratio\nwithout bins within 2 µm of the nucleus edge"); axs[1].set_title("CP-1 with the edge kernel (band: within a tenth)", loc="left")
    fig.tight_layout(); fig.savefig(f"{FD}/R5_F3_nucleus_ratio.png", dpi=200)

# F5: sfocatura interna (PH-2, per ROI) contro la larghezza ai bordi (P-1, per sezione)
if os.path.exists(f"{O}/tables/R5_PH2_inner_summary.csv"):
    H = pd.read_csv(f"{O}/tables/R5_PH2_inner_summary.csv"); H["a"] = H.tag.str[:2]
    P = pd.read_csv(f"{O}/tables/R5_predictions_edges.csv"); sec = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
    A = ["A1", "A2", "A3", "A4", "A5", "A6"]
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    for i, a in enumerate(A):
        r = P[(P.key == sec[a]) & (P.prediction == "P-1")].iloc[0]
        ax.errorbar(i - 0.15, r.estimate, yerr=[[r.estimate - r.lo], [r.hi - r.estimate]], fmt="s", color="#0072B2", capsize=2, label="edge kernel σ (tissue → glass)" if i == 0 else None)
        h = H[H.a == a]; hs = h[h.sigma_in_star.notna()]; hn = h[h.sigma_in_star.isna()]
        ax.plot(np.full(len(hs), i + 0.15), hs.sigma_in_star, "D", ms=3, color="#D55E00", label="inner blur σ_in (PH-2, per ROI)" if i == 0 else None)
        below = hn[hn.result.str.contains("below grid")]
        ax.plot(np.full(len(below), i + 0.15), np.full(len(below), 0.18), "v", ms=4, mfc="none", color="#D55E00", label="no solution ≥ 0.25 µm (drawn at 0.18)" if i == 0 else None)
        if len(below): ax.text(i + 0.27, 0.18, f"×{len(below)}", va="center", fontsize=6, color="#D55E00")
        other = len(hn) - len(below)
        if other: ax.text(i + 0.15, 6.5, f"{other} ROI\nno solution", ha="center", fontsize=5)
    ax.set_yscale("log"); ax.set_ylim(0.15, 9); ax.set_yticks([0.25, 0.5, 1, 2, 4, 8]); ax.set_yticklabels(["0.25", "0.5", "1", "2", "4", "8"]); ax.minorticks_off(); ax.set_xticks(range(6)); ax.set_xticklabels([TIS[a] for a in A], rotation=30, ha="right")
    ax.set_ylabel("Gaussian width (µm)"); ax.set_title("Blur at tissue edges vs inside tissue", loc="left"); ax.legend(frameon=False, fontsize=6)
    fig.tight_layout(); fig.savefig(f"{FD}/R5_F5_inner_vs_edge.png", dpi=200)

# F4: lumi
if os.path.exists(f"{O}/tables/R5_CP2_lumens.csv"):
    L = pd.read_csv(f"{O}/tables/R5_CP2_lumens.csv")
    fig, ax = plt.subplots(figsize=(6, 3.0))
    bands = ["2-4", "4-8", "8-16", "16-32"]
    for i, key in enumerate(KEYS):
        l = L[L.key == key].set_index("band_um").reindex(bands)
        xx = np.arange(4) + (i - 2) * 0.12
        ax.errorbar(xx, l.ratio, yerr=[l.ratio - l.lo, l.hi - l.ratio], fmt="o", ms=3, capsize=2, label=f"{NAME[key]} ({int(l.n_lumens.max()) if l.n_lumens.notna().any() else 0} lumens)")
    ax.axhspan(0.8, 1.25, color="#999999", alpha=0.15); ax.axhline(1, color="k", lw=0.5); ax.set_yscale("log")
    ax.set_xticks(range(4)); ax.set_xticklabels([f"{b} µm" for b in bands]); ax.set_xlabel("distance from lumen edge")
    ax.set_ylabel("observed / predicted UMI"); ax.set_title("Empty lumens inside tissue vs edge kernel (grey: within a quarter)", loc="left")
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout(); fig.savefig(f"{FD}/R5_F4_lumens.png", dpi=200)
print("figures written")
