"""tools/S2.0_ph1_stripes_glass.py — S2.0, post hoc PH-1 (dichiarata dopo M1): i fattori non periodici di riga e colonna stimati nel tessuto
si ritrovano sul vetro lontano dal tessuto (>= 100 µm, a >= 100 µm dal bordo dell'area)? E spettro della parte riproducibile.
Uso: .venv/bin/python tools/S2.0_ph1_stripes_glass.py <outdir> -> <outdir>/stripes/PH1_<KEY>.json"""
import sys, json
import numpy as np
from scipy import ndimage
sys.path.insert(0, "tools")
import S2_0_common as S
R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
OUT = sys.argv[1]
for KEY in ("A1", "A23", "A4", "A5", "A6"):
    z = np.load(f"{R5RUN}/section/{KEY}.npz"); U = z["U"].astype(float); T = z["in_tissue"].astype(bool)
    F = np.load(f"{OUT}/stripes/{KEY}_factors.npz")
    dist = ndimage.distance_transform_edt(~T) * S.BIN_UM
    G = (~T) & (dist >= 100); G[:50] = G[-50:] = False; G[:, :50] = G[:, -50:] = False
    res = {"key": KEY, "n_glass_bins": int(G.sum()), "glass_mean": float(U[G].mean())}
    for ax, nm in ((0, "rows"), (1, "cols")):
        idx = np.arange(U.shape[ax])
        tis = 0.5 * (F[f"{nm}_half1"] + F[f"{nm}_half2"])          # parte non periodica nel tessuto (media delle meta')
        n = G.sum(1 - ax); s = (U * G).sum(1 - ax)
        with np.errstate(divide="ignore", invalid="ignore"):
            g = np.log(s / n / U[G].mean())
        okg = (n >= 300) & (s > 0)
        _, gnp = S.periodic_split(np.where(okg, g, np.nan), idx)
        ok = np.isfinite(gnp) & np.isfinite(tis)
        # rumore di Poisson atteso sul vetro: var(log) ~ 1 / conteggi
        pois_sd = float(np.sqrt(np.mean(1.0 / s[ok]))) if ok.any() else float("nan")
        res[nm] = {"n_rows_both": int(ok.sum()), "corr_tissue_glass": float(np.corrcoef(tis[ok], gnp[ok])[0, 1]) if ok.sum() > 10 else float("nan"),
                   "slope_glass_on_tissue": float(np.polyfit(tis[ok], gnp[ok], 1)[0]) if ok.sum() > 10 else float("nan"),
                   "sd_glass_np": float(np.nanstd(gnp[ok])), "poisson_sd_glass": pois_sd, "mean_glass_counts_per_row": float(np.mean(s[ok])) if ok.any() else float("nan")}
        # spettro incrociato delle due meta' (parte riproducibile)
        a, b = F[f"{nm}_half1"], F[f"{nm}_half2"]; okab = np.isfinite(a) & np.isfinite(b)
        a2 = np.where(okab, a - np.nanmean(a[okab]), 0.0); b2 = np.where(okab, b - np.nanmean(b[okab]), 0.0)
        X = np.real(np.fft.rfft(a2) * np.conj(np.fft.rfft(b2))); fr = np.fft.rfftfreq(len(a2))
        top = np.argsort(X[1:])[::-1][:6] + 1
        res[nm]["cross_spectrum_top_periods_bins"] = [round(1 / fr[i], 2) for i in top]
        res[nm]["cross_power_frac_period_le_4"] = float(X[fr >= 0.25].sum() / X[1:].sum())
        res[nm]["cross_power_frac_period_gt_20"] = float(X[(fr > 0) & (fr < 0.05)].sum() / X[1:].sum())
    json.dump(res, open(f"{OUT}/stripes/PH1_{KEY}.json", "w"), indent=1)
    print(json.dumps(res), flush=True)
