#!/usr/bin/env python
"""tools/S1.4b_posthoc_arbiter.py — S1.4b, analisi POST HOC di C-L2 (dopo il verdetto; non lo cambia).
PH-A: arbitro ufficiale con 20x punti (4e6) sugli insiemi falliti per il sotto-criterio delle aree.
PH-B: calibrazione esatta della statistica «celle con |z| > 4»: conteggi multinomiali simulati con le stesse aree e lo
stesso numero di punti (500 repliche); p = frazione delle repliche con statistica >= osservata.
Uso: .venv/bin/python tools/S1.4b_posthoc_arbiter.py <dir arb_in> <arbiter.csv ufficiale> <out_prefix>
"""
import sys, os
import numpy as np, pandas as pd, shapely
sys.path.insert(0, "tools"); import importlib.util
spec = importlib.util.spec_from_file_location("arb", "tools/S1.4b_arbiter.py"); arb = importlib.util.module_from_spec(spec); spec.loader.exec_module(arb)
d_in, off, pref = sys.argv[1:4]
o = pd.read_csv(off); fail = o.loc[~o["pass_"], "set"].tolist()
rows_a = [arb.run_set(os.path.join(d_in, s + ".parquet"), 4_000_000) for s in fail]
pd.DataFrame(rows_a).to_csv(pref + "_PHA.csv", index=False)
rng = np.random.default_rng(20261017); rows_b = []
for _, r in o.iterrows():
    d = pd.read_parquet(os.path.join(d_in, r["set"] + ".parquet"))
    reg = shapely.from_wkb(open(os.path.join(d_in, r["set"] + "_region.wkb"), "rb").read())
    cells = shapely.from_wkb(d["wkb"].to_numpy()); xmin, ymin, xmax, ymax = reg.bounds
    A_box = (xmax - xmin) * (ymax - ymin); p = shapely.area(cells) / A_box; ok = p > 0
    N = 200000; pv = np.append(p, max(0.0, 1 - p.sum()))
    sims = []
    for b in range(500):
        n = rng.multinomial(N, pv)[:-1]
        z = (n[ok] - N * p[ok]) / np.sqrt(N * p[ok] * (1 - p[ok])); sims.append(int((np.abs(z) > 4).sum()))
    sims = np.array(sims)
    rows_b.append(dict(set=r["set"], n_cells=int(ok.sum()), mean_pts_per_cell=N * p[ok].mean(), observed=int(r["n_absz_gt4"]),
                       null_mean=sims.mean(), null_q95=np.quantile(sims, 0.95), null_q99=np.quantile(sims, 0.99), p_value=float(np.mean(sims >= r["n_absz_gt4"])),
                       gauss_E=float(r["exp_absz_gt4"])))
pd.DataFrame(rows_b).to_csv(pref + "_PHB.csv", index=False)
print(pd.DataFrame(rows_a)[["set", "n_pts_region", "misassigned", "uncovered", "double", "n_absz_gt4", "exp_absz_gt4", "pass_"]].to_string())
print(pd.DataFrame(rows_b).to_string())
