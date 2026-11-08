"""tools/R7_kappa_implied.py — R7, regola R-1 dopo P-2: kappa e kappa_range calibrati sulla correlazione implicita
della struttura del tessuto (M-4, stima puntuale) dove P-2 è smentita (stroma, linfonodo); calcolati anche negli altri
tessuti come sensibilità. Media di 5 semi. -> results/R7/coverage/R7_kappa_implied.csv, R7_kappa_final_v2.csv"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R7_common as C
pairs = pd.read_csv("results/R7/cp1/R7_cp1_pairs.csv"); P2 = pd.read_csv("results/R7/cp1/R7_P2.csv")
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); K = pd.read_csv("results/R7/coverage/R7_kappa_final.csv")
rows = []
for (A, f), g in pairs.groupby(["archetype", "field"]):
    ids = list(rois[rois.archetype == A].roi_id); R = np.eye(5)
    for _, r in g.iterrows():
        i, j = ids.index(r.i), ids.index(r.j); R[i, j] = R[j, i] = r.R_implied
    ks = [C.calibrate_kappa(R, np.random.default_rng([C.SEED, 1000 + int(A[1]), int(f == "nuclei"), s])) for s in range(5)]
    rows.append(dict(archetype=A, field=f, rho_bar=C.rho_bar(R), kappa=float(np.mean(ks)), kappa_sd5=float(np.std(ks, ddof=1)),
                     kappa_range=C.kappa_range(R, np.random.default_rng([C.SEED, 1100 + int(A[1]), int(f == "nuclei")])),
                     P2_verdict=P2[(P2.archetype == A) & (P2.field == f)].verdict.iloc[0]))
I = pd.DataFrame(rows); I.to_csv("results/R7/coverage/R7_kappa_implied.csv", index=False); print(I.round(4).to_string(index=False))
K2 = K.copy(); K2["source"] = "F (geometria, R-2)"
for A in ("A3", "A4"):                       # P-2 smentita per entrambi i campi
    u = I[(I.archetype == A) & (I.field == "umi")].iloc[0]; n = I[(I.archetype == A) & (I.field == "nuclei")].iloc[0]
    m = K2.archetype == A
    K2.loc[m, "kappa_bin"] = u.kappa; K2.loc[m, "kappa_nuc"] = n.kappa; K2.loc[m, "kappa_box"] = max(u.kappa, n.kappa)
    K2.loc[m, "kappa_range_bin"] = u.kappa_range; K2.loc[m, "kappa_range_nuc"] = n.kappa_range; K2.loc[m, "kappa_range_box"] = max(u.kappa_range, n.kappa_range)
    K2.loc[m, "kappa_se"] = max(u.kappa_sd5, n.kappa_sd5) / np.sqrt(5)
    K2.loc[m, "source"] = "struttura del tessuto misurata dentro le aree (R-1 dopo P-2 smentita); box = massimo dei due campi"
K2.to_csv("results/R7/coverage/R7_kappa_final_v2.csv", index=False); print(K2.round(4).to_string(index=False))
