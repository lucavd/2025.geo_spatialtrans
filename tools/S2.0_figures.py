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

if __name__ == "__main__":
    {"stripes": stripes}[PART]()
