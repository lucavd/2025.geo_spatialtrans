"""tools/S1.5b_posthoc.py - S1.5b, analisi post hoc (dopo i verdetti; non cambiano verdetti ne' preset).
PH-1: composizione dei nuclei del cervello persi da Cellpose corretta per l'errore del trasferimento (matrice di confusione sui punti trovati,
inversione con minimi quadrati non negativi), valutatori Luca e Federica. Uscita: results/S1.5b/S1.5b_posthoc_a5.csv"""
import numpy as np, pandas as pd
from scipy.optimize import nnls
D = pd.read_csv("results/S1.5b/S1.5b_a5_missed_points.csv"); C = pd.read_csv("results/S1.5b/S1.5b_a5_composition.csv")
types = ["neuroni", "astrociti", "oligodendrociti", "microglia", "vascolari"]; rows = []
for rater in ["Luca", "Federica", "both"]:
    G = D if rater == "both" else D[D.rater == rater]
    m = G[G.matched & G.type_r6.isin(types) & G.type_transfer.isin(types)]
    Mx = pd.crosstab(m.type_r6, m.type_transfer).reindex(index=types, columns=types, fill_value=0).to_numpy(float)
    P = Mx / Mx.sum(axis=1, keepdims=True)                       # P[vero, trasferito]
    obs = G[~G.matched].type_transfer.value_counts().reindex(types, fill_value=0).to_numpy(float)
    x, _ = nnls(P.T, obs); xs = x / x.sum()
    recall = np.diag(P)
    for k, t in enumerate(types):
        rows.append(dict(rater=rater, cell_type=t, n_missed=int(obs.sum()), observed_share=obs[k] / obs.sum(), corrected_share=xs[k], recall_matched=recall[k], n_matched_true=int(Mx[k].sum())))
R = pd.DataFrame(rows)
r6 = C[C.rater == "Luca"].set_index("cell_type").fraction_r6; q = C[C.rater == "Luca"].q.iloc[0]
R["composition_if_corrected"] = [ (1 - q) * r6[t] + q * s if r == "Luca" else np.nan for r, t, s in zip(R.rater, R.cell_type, R.corrected_share)]
R.to_csv("results/S1.5b/S1.5b_posthoc_a5.csv", index=False); pd.set_option("display.width", 200); print(R.round(3).to_string(index=False))

# PH-2: oggetti di area < 8 um2 nel cervello (vicini al filtro minimo di 4 um2 delle maschere): peso sul modello della glia e sulla banda
import importlib.util
spec = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(spec); spec.loader.exec_module(cal)
lab, d = cal.load("A5"); ev = d[d.interior & d.matched & (d.a_n > 0)]
glia = ["astrociti", "microglia", "oligodendrociti", "vascolari"]; aref = 1e6 / 1039.7; thr = 59.1
out = []
for name, e in [("tutti", ev), ("senza a_n < 8", ev[ev.a_n >= 8])]:
    g = e[e.type.isin(glia)]; rho, b, s, n = cal.fit(g.a_n.values, g.a_t.values, aref)
    band = e.groupby("roi_id").apply(lambda x: ((x.a_n > thr / np.sqrt(2)) & (x.a_n < thr * np.sqrt(2))).mean())
    out.append(dict(set=name, n_glia=n, rho_glia=rho, beta_glia=b, sigma_glia=s, band_min=band.min(), band_max=band.max(),
                    frac_small_r1=float((ev[ev.roi_id == "r1"].a_n < 8).mean()), frac_small_oligo=float((ev[ev.type == "oligodendrociti"].a_n < 8).mean())))
O = pd.DataFrame(out); O.to_csv("results/S1.5b/S1.5b_posthoc_a5_small.csv", index=False); print(O.round(4).to_string(index=False))
