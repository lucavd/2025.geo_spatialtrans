"""tools/S1.5d_explore_size.py - S1.5d, esplorazione dichiarata del modello di taglia (prima della pre-registrazione dell'aggiunta sul modello di taglia).
Tumore del colon (A2) e, per confronto, gli altri tessuti. Celle interne. Reale: nuclei veri con territorio (insieme di S1.5b/S1.5c).
Simulato: layer del seme 1 dei bracci U e N-a di S1.5c (/mnt/micron/geo_spatialtrans/S1.5c/layers), con l'area bersaglio del modello di taglia.
Domande: (S1) quanto pesa la coda N/C >= 0.95 per tipo, vero contro simulato; (S2) il residuo del modello log-log per tipo dipende da A_T
(eteroschedasticita'), e la coda del vero e' piu' corta di una lognormale con la stessa sigma; (S3) quanto troncherebbero tetti phi diversi.
Uscita: results/S1.5d/S1.5d_explore_size.csv. Uso: .venv/bin/python tools/S1.5d_explore_size.py"""
import importlib.util, subprocess, json, numpy as np, pandas as pd
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py")
rows = []
for A in cal.ARCH:
    _, d = cal.load(A); e = d[d.interior.astype(bool)]
    for t, g in [("_all", e)] + list(e.groupby("type")):
        if len(g) < 200: continue
        nc = g.a_n / g.a_t; lx = np.log(g.a_t); ly = np.log(g.a_n)
        b = np.polyfit(lx, ly, 1); res = ly - np.polyval(b, lx)
        q = pd.qcut(g.a_t, 4, labels=False)
        o = dict(archetype=A, source="real", cell_type=t, n=len(g), frac_nc095=np.mean(nc >= 0.95), q99_nc=np.quantile(nc, 0.99), med_nc=np.median(nc),
                 beta=b[0], sigma=res.std(ddof=2))
        for k in range(4): o[f"sigma_q{k+1}"] = res[q.values == k].std(ddof=1)
        # coda attesa da una lognormale con la stessa sigma, condizionata a A_T: P(A_N >= 0.95 A_T)
        from scipy.stats import norm
        mu = np.polyval(b, lx); o["frac_nc095_lognormal"] = float(np.mean(norm.sf((np.log(0.95) + lx - mu) / res.std(ddof=2))))
        rows.append(o)
pd.DataFrame(rows).to_csv("results/S1.5d/S1.5d_explore_size_real.csv", index=False)
print("ok", len(rows))
