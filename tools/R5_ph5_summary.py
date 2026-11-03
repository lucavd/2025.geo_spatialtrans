"""tools/R5_ph5_summary.py — riassunto di PH-5 (bordi con la parte gaussiana fissata a sigma = SIG): per sezione mediana di Delta NLL
(fissato - libero) per tessera e per 1000 bin, quota di tessere con Delta NLL > 2 (una sola differenza di parametri: soglia di Wilks 1.92 al 5 %),
mediane di w, ell, scarto con il nucleo stretto. Uso: .venv/bin/python tools/R5_ph5_summary.py <outdir> [SIG=0.7] -> <outdir>/tables/R5_PH5_narrowcore.csv"""
import sys, json, glob, numpy as np, pandas as pd
OUT = sys.argv[1]; SIG = float(sys.argv[2]) if len(sys.argv) > 2 else 0.7
rows = []
for f in sorted(glob.glob(f"{OUT}/ph5_s{SIG:g}/*.json")):
    d = json.load(open(f)); ft = d["fit"]
    rows.append({"key": d["key"], "tile": d["tile"], "dnll": d["delta_nll"], "n_bins": d["n_bins"], "w": ft["w"], "ell": ft["ell"], "dr": ft["dr"], "dc": ft["dc"],
                 "sigma_free": d["sigma_K2_free"], "w_free": d["w_K2_free"]})
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/tables/R5_PH5_narrowcore_tiles.csv", index=False)
S = T.groupby("key").apply(lambda g: pd.Series({"n_tiles": len(g), "sigma_fixed": SIG, "median_dnll": g.dnll.median(), "median_dnll_per_1000_bins": (1000 * g.dnll / g.n_bins).median(),
     "frac_dnll_gt_2": float((g.dnll > 2).mean()), "frac_dnll_lt_0": float((g.dnll < -0.01).mean()), "median_w_narrow": g.w.median(), "median_ell_narrow": g.ell.median(),
     "median_sigma_free": g.sigma_free.median(), "median_w_free": g.w_free.median()}), include_groups=False).reset_index()
S.to_csv(f"{OUT}/tables/R5_PH5_narrowcore.csv", index=False); pd.set_option("display.width", 250); print(S.round(3).to_string(index=False))
