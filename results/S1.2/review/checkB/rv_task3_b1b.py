# rv_task3_b1b.py — revisione avversariale S1.2 / checkB, compito 3.
# Frazione esclusa e densita' attesa = totale x (1 - esclusa) sui 30 ROI graphclust, con codice indipendente:
# componenti 4-connesse per cluster (scipy.ndimage.label, struttura a croce), escluse se area < 100 µm²
# (pixel 8 µm -> 64 µm²/pixel), denominatore = tutti i pixel della mappa. Confronto con results/S1.2/S1.2_B1b.csv.
# Controprova: con 8-connessione la frazione esclusa deve cambiare (diagonali unite).
# Uso: cd ~/2025.geo_spatialtrans && .venv/bin/python results/S1.2/review/checkB/rv_task3_b1b.py
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from math import gcd
from functools import reduce
RV = "results/S1.2/review/checkB"
b1 = pd.read_csv("results/S1.2/S1.2_B1b.csv"); meta = pd.read_csv(f"{RV}/data/clust_meta.csv")
TOT = dict(A1=12419, A2=8924, A3=3100, A4=28096, A5=1185, A6=988); EV = dict(A1=10291, A2=7376, A3=2311, A4=26844, A5=953, A6=906)
rows = []
for _, mr in meta.iterrows():
    A, roi, px = mr.archetype, mr.roi_id, mr.pixel_size_um
    cl = pd.read_csv(f"{RV}/data/clust_{A}_{roi}.csv", dtype={"cluster": str})
    assert not cl.duplicated(["x", "y"]).any()
    stride = reduce(gcd, np.unique(np.diff(np.unique(cl.x))).tolist() + np.unique(np.diff(np.unique(cl.y))).tolist())
    x = cl.x.to_numpy() - cl.x.min(); y = cl.y.to_numpy() - cl.y.min()
    res = {}
    for name, st in (("4conn", ndi.generate_binary_structure(2, 1)), ("8conn", ndi.generate_binary_structure(2, 2))):
        exc = 0; ncomp = 0; nexc = 0
        for k in cl.cluster.unique():
            g = np.zeros((y.max() + 1, x.max() + 1), bool); s = (cl.cluster == k).to_numpy(); g[y[s], x[s]] = True
            lab, nl = ndi.label(g, structure=st)
            sz = np.bincount(lab.ravel())[1:]
            small = sz * px * px < 100
            exc += sz[small].sum(); ncomp += nl; nexc += small.sum()
        res[name] = (exc / len(cl), ncomp, nexc)
    p = b1[(b1.archetype == A) & (b1.roi_id == roi)].iloc[0]
    fe = res["4conn"][0]
    rows.append(dict(archetype=A, roi_id=roi, stride=stride, n_px=len(cl), n_clusters=cl.cluster.nunique(),
                     area_total_mm2_rv=len(cl) * px * px / 1e6, area_total_mm2_ref=p.area_total_mm2,
                     frac_excl_rv=fe, frac_excl_ref=p.frac_excluded, diff_frac=fe - p.frac_excluded,
                     n_comp_4=res["4conn"][1], n_excl_4=res["4conn"][2], frac_excl_8conn=res["8conn"][0],
                     expected_rv=TOT[A] * (1 - fe), expected_ref=p.expected, dens_seed42=p.dens_seed42,
                     seed42_minus_expected_cells=(p.dens_seed42 - TOT[A] * (1 - fe)) * len(cl) * px * px / 1e6,
                     expected_in_interval=EV[A] <= TOT[A] * (1 - fe) <= TOT[A], seed42_in_interval=EV[A] <= p.dens_seed42 <= TOT[A],
                     frac_excl_max_for_PASS=1 - EV[A] / TOT[A]))
out = pd.DataFrame(rows); out.to_csv(f"{RV}/rv_task3_b1b.csv", index=False)
print(out.round(5).to_string())
print("max |diff frac|", out.diff_frac.abs().max(), " max |diff area|", (out.area_total_mm2_rv - out.area_total_mm2_ref).abs().max())
print(out.groupby("archetype").agg(exp_in=("expected_in_interval", "sum"), s42_in=("seed42_in_interval", "sum"),
      med_fe=("frac_excl_rv", "median"), med_fe8=("frac_excl_8conn", "median"), excess_mean=("seed42_minus_expected_cells", "mean")))
