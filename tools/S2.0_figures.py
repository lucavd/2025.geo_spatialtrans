"""tools/S2.0_figures.py — S2.0: figure del report (sfondo bianco, un pannello per sezione).
Uso: .venv/bin/python tools/S2.0_figures.py <outdir> <part>   part: stripes | edges | inner | disp | lumens"""
import sys, os, json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "tools")
OUT, PART = sys.argv[1], sys.argv[2]
FIGD = "reports/figures/S2.0"; os.makedirs(FIGD, exist_ok=True)
SEC = ["A1", "A23", "A4", "A5", "A6"]
NAME = {"A1": "mouse intestine", "A23": "colon cancer (tumour, stroma)", "A4": "lymph node", "A5": "mouse brain", "A6": "heart"}
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

def stripes():
    R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
    fig, ax = plt.subplots(2, 5, figsize=(16, 6.4))
    for j, k in enumerate(SEC):
        F = np.load(f"{OUT}/stripes/{k}_factors.npz"); ph = json.load(open(f"{OUT}/stripes/PH1_{k}.json")); js = json.load(open(f"{OUT}/stripes/{k}.json"))
        a, b = F["rows_half1"], F["rows_half2"]; ok = np.isfinite(a) & np.isfinite(b)
        ax[0, j].scatter(a[ok], b[ok], s=2, alpha=0.3, color="#1f77b4", rasterized=True)
        lim = np.nanpercentile(np.abs(np.r_[a[ok], b[ok]]), 99.5); ax[0, j].plot([-lim, lim], [-lim, lim], color="grey", lw=0.8)
        ax[0, j].set_xlim(-lim, lim); ax[0, j].set_ylim(-lim, lim)
        ax[0, j].set_title(f"{NAME[k]}\nrows: SD {js['rows']['sd_repro']:.3f}, r = {js['rows']['corr_halves']:.2f}")
        ax[0, j].set_xlabel("row factor, left half (log)"); ax[0, j].set_ylabel("row factor, right half (log)" if j == 0 else "")
        z = np.load(f"{R5RUN}/section/{k}.npz"); U = z["U"].astype(float); T = z["in_tissue"].astype(bool)
        from scipy import ndimage
        import S2_0_common as S
        dist = ndimage.distance_transform_edt(~T) * 2.0; G = (~T) & (dist >= 100); G[:50] = G[-50:] = False; G[:, :50] = G[:, -50:] = False
        n = G.sum(1); s = (U * G).sum(1)
        with np.errstate(divide="ignore", invalid="ignore"): g = np.log(s / n / U[G].mean())
        _, gnp = S.periodic_split(np.where((n >= 300) & (s > 0), g, np.nan), np.arange(len(g)))
        t = 0.5 * (a + b); okg = np.isfinite(t) & np.isfinite(gnp)
        ax[1, j].scatter(t[okg], gnp[okg], s=2, alpha=0.3, color="#d62728", rasterized=True)
        ax[1, j].plot([-lim, lim], [-lim, lim], color="grey", lw=0.8); ax[1, j].set_xlim(-lim, lim)
        ax[1, j].set_title(f"far glass vs tissue: r = {ph['rows']['corr_tissue_glass']:.2f} ({ph['rows']['n_rows_both']} rows)")
        ax[1, j].set_xlabel("row factor in tissue (log)"); ax[1, j].set_ylabel("row factor on glass >= 100 µm (log)" if j == 0 else "")
    fig.tight_layout(); fig.savefig(f"{FIGD}/S2.0_F1_stripes_rows.png", dpi=130)
    # spettro e mappa della parita' per tessera
    fig, ax = plt.subplots(1, 5, figsize=(16, 3.2))
    for j, k in enumerate(SEC):
        F = np.load(f"{OUT}/stripes/{k}_factors.npz")
        for nm, col in (("rows", "#1f77b4"), ("cols", "#ff7f0e")):
            a, b = F[f"{nm}_half1"], F[f"{nm}_half2"]; ok = np.isfinite(a) & np.isfinite(b)
            a2 = np.where(ok, a - np.nanmean(a[ok]), 0); b2 = np.where(ok, b - np.nanmean(b[ok]), 0)
            X = np.real(np.fft.rfft(a2) * np.conj(np.fft.rfft(b2))); fr = np.fft.rfftfreq(len(a2))
            m = fr > 0; Xs = pd.Series(X[m]).rolling(25, center=True, min_periods=1).mean().to_numpy()
            ax[j].plot(1 / fr[m], Xs / Xs.max(), color=col, lw=0.8, label=nm)
        ax[j].set_xscale("log"); ax[j].set_xlabel("period (bins of 2 µm)"); ax[j].set_title(NAME[k]); ax[j].legend(frameon=False)
        ax[j].set_ylabel("cross-power between halves (smoothed, max = 1)" if j == 0 else "")
    fig.tight_layout(); fig.savefig(f"{FIGD}/S2.0_F2_stripes_spectrum.png", dpi=130)
    print("stripes figures done")

RES = os.environ.get("S20_RES", "results/S2.0")
COL = {"C0": "#7f7f7f", "C1": "#1f77b4", "C2": "#d62728", "C3": "#2ca02c"}
RLAB = {"C0": "C0 R5 core", "C1": "C1 two populations", "C2": "C2 edge zone", "C3": "C3 glass transport"}

def edges():
    T = pd.read_csv(f"{RES}/S2.0_edges_tiles.csv"); T = T[~T.error.astype(bool)]
    fig, ax = plt.subplots(1, 5, figsize=(16, 3.4), sharey=True)
    for j, k in enumerate(SEC):
        t = T[T.key == k]
        data = [t[f"dAIC_{r}"].dropna().to_numpy() for r in ("C1", "C2", "C3")]
        bp = ax[j].boxplot(data, widths=0.6, showfliers=False, patch_artist=True)
        for p_, r in zip(bp["boxes"], ("C1", "C2", "C3")): p_.set_facecolor(COL[r]); p_.set_alpha(0.4)
        for i, d in enumerate(data):
            ax[j].scatter(np.random.default_rng(1).normal(i + 1, 0.06, len(d)), d, s=3, color=COL[("C1", "C2", "C3")[i]], alpha=0.5, rasterized=True)
        ax[j].axhline(0, color="black", lw=0.6); ax[j].set_yscale("symlog", linthresh=10)
        ax[j].set_xticks([1, 2, 3]); ax[j].set_xticklabels(["C1", "C2", "C3"]); ax[j].set_title(f"{NAME[k]} ({len(t)} windows)")
    ax[0].set_ylabel("ΔAIC vs C0 per window (negative = better)")
    fig.tight_layout(); fig.savefig(f"{FIGD}/S2.0_F3_edges_dAIC.png", dpi=130)
    E = pd.read_csv(f"{RES}/S2.0_edges_sections.csv")
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
    for a_, (rl, pm, lab) in zip(ax, (("C0", "sigma", "C0 blur σ (µm)"), ("C2", "sigma_e", "C2 edge-zone blur σ_e (µm)"), ("C2", "D", "C2 edge-zone depth D (µm)"))):
        e = E[(E.rule == rl) & (E.param == pm)].set_index("key").reindex(SEC)
        a_.errorbar(range(5), e["median"], yerr=[e["median"] - e.lo, e.hi - e["median"]], fmt="o", color=COL[rl], capsize=3)
        a_.set_xticks(range(5)); a_.set_xticklabels([NAME[k].split(" (")[0] for k in SEC], rotation=25, ha="right"); a_.set_ylabel(lab)
        a_.set_title("section median of windows [95% CI]")
        if pm == "D": a_.set_yscale("log"); a_.axhline(1e4, color="grey", ls=":", lw=0.8); a_.text(0, 6e3, "upper bound of the fit", fontsize=7, color="grey")
    fig.tight_layout(); fig.savefig(f"{FIGD}/S2.0_F4_edges_params.png", dpi=130)
    print("edges figures done")

def inner():
    I = pd.read_csv(f"{RES}/S2.0_inner_summary.csv")
    AR = ["A1", "A2", "A3", "A4", "A5", "A6"]; AN = ["intestine", "tumour", "stroma", "lymph node", "brain", "heart"]
    fig, ax = plt.subplots(figsize=(9, 3.8))
    for i, rl in enumerate(("C0", "C1", "C2", "C3")):
        for j, A_ in enumerate(AR):
            x = I[(I.archetype == A_) & (I.rule == rl)]
            if x.empty: continue
            x = x.iloc[0]; xp = j + (i - 1.5) * 0.18
            if np.isfinite(x.ratio_mean):
                ax.errorbar(xp, x.ratio_mean, yerr=[[x.ratio_mean - x.lo], [x.hi - x.ratio_mean]], fmt="o", color=COL[rl], capsize=2, label=RLAB[rl] if j == 0 else None)
            elif np.isfinite(x.ratio_mean_ok_only):
                ax.errorbar(xp, x.ratio_mean_ok_only, yerr=[[x.ratio_mean_ok_only - x.lo_ok_only], [x.hi_ok_only - x.ratio_mean_ok_only]], fmt="o", mfc="white", color=COL[rl], capsize=2)
    from matplotlib.lines import Line2D
    hd = [Line2D([], [], marker="o", ls="", color=COL[r_], label=RLAB[r_]) for r_ in ("C0", "C1", "C2", "C3")]
    ax.axhspan(0.909, 1.1, color="#cccccc", alpha=0.4, lw=0); ax.set_xticks(range(6)); ax.set_xticklabels(AN)
    ax.set_ylabel("predicted / observed, deep nucleus bins"); ax.set_ylim(0.4, 1.4); ax.legend(handles=hd, frameon=False, fontsize=7, ncol=4, loc="lower left")
    ax.set_title("inner check: grey = [0.909, 1.1]; open symbol = only ROIs where r could be inverted")
    fig.tight_layout(); fig.savefig(f"{FIGD}/S2.0_F5_inner.png", dpi=130); print("inner figure done")

def disp():
    def load(d): return pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/{d}/*.json"))])
    a = load("disp"); b = load("disp_stripes")
    m = a.merge(b[["tag", "ratio", "v"]], on="tag", suffixes=("", "_st"))
    AR = ["A1", "A2", "A3", "A4", "A5", "A6"]; AN = ["intestine", "tumour", "stroma", "lymph node", "brain", "heart"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    for j, A_ in enumerate(AR):
        t = m[m.archetype == A_]
        ax[0].scatter(np.full(len(t), j - 0.12), t.ratio, color="#7f7f7f", s=14, label="pre-registered prototype (parity only)" if j == 0 else None)
        ax[0].scatter(np.full(len(t), j + 0.12), t.ratio_st, color="#d62728", s=14, label="post hoc: + measured row/column factors" if j == 0 else None)
        TC = plt.cm.tab10(j)
        ax[1].scatter(t.mbar_obs, t.ratio - 1, facecolors="none", edgecolors=TC, s=16); ax[1].scatter(t.mbar_obs, t.ratio_st - 1, color=TC, s=16, label=AN[j])
    ax[0].axhspan(0.909, 1.1, color="#cccccc", alpha=0.5, lw=0); ax[0].set_xticks(range(6)); ax[0].set_xticklabels(AN); ax[0].set_ylabel("dispersion index, observed / simulated")
    ax[0].legend(frameon=False, fontsize=7); ax[0].set_title("M3: one point per ROI; grey = [0.909, 1.1]")
    ax[1].set_xlabel("mean UMI per deep nucleus bin (2 µm)"); ax[1].set_ylabel("observed / simulated - 1")
    ax[1].set_title("open = parity only, filled = + row/column factors"); ax[1].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIGD}/S2.0_F6_disp.png", dpi=130); print("disp figure done")

def lumens():
    L = pd.read_csv(f"{RES}/S2.0_CP1_lumens.csv")
    fig, ax = plt.subplots(1, 5, figsize=(16, 3.2), sharey=True)
    for j, k in enumerate(SEC):
        for i, rl in enumerate(("C0", "C2")):
            t = L[(L.key == k) & (L.rule == rl)]
            x = np.arange(len(t)) + (i - 0.5) * 0.2
            ax[j].errorbar(x, t.obs_over_pred, yerr=[t.obs_over_pred - t.lo, t.hi - t.obs_over_pred], fmt="o", color=COL[rl], capsize=2, label=RLAB[rl])
        ax[j].axhspan(0.8, 1.25, color="#cccccc", alpha=0.5, lw=0); ax[j].axvspan(-0.5, 1.5, color="#fff2cc", alpha=0.5, lw=0)
        ax[j].set_xticks(range(4)); ax[j].set_xticklabels(["2-4", "4-8", "8-16", "16-32"]); ax[j].set_xlabel("distance from lumen edge (µm)")
        nl = int(L[L.key == k].n_lumens.iloc[0]); ax[j].set_title(f"{NAME[k]} ({nl} lumens)" + (" descriptive" if k == "A5" else ""))
    ax[0].set_ylabel("observed / predicted UMI"); ax[0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIGD}/S2.0_F7_lumens.png", dpi=130); print("lumens figure done")

if __name__ == "__main__":
    {"stripes": stripes, "edges": edges, "inner": inner, "disp": disp, "lumens": lumens}[PART]()
