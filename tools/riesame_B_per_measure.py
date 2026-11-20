"""tools/riesame_B_per_measure.py - Riesame Parte B: distanza con segno per tessuto x misura (in unita' di h) di T (P1, P5), S e B da V
sui geni comuni, stesso stimatore di tools/riesame_D.py (termini per cella, prima della media). Descrittivo: nessun w, nessun esito.
Ingressi: results/riesame/B_area_values.csv, results/riesame/B_h_panel.csv. Uscita: results/riesame/B_per_measure_common.csv."""
import numpy as np, pandas as pd
U = pd.read_csv("results/riesame/B_area_values.csv"); H = pd.read_csv("results/riesame/B_h_panel.csv")
U = U[U.genes == "common"]; H = H[H.genes == "common"]; rows = []
for A in ["A2", "A3"]:
    for m in ["M1", "M2", "M3", "M4", "M5", "M6"]:
        h = float(H[(H.A == A) & (H.measure == f"M-{m[1]}")].h.iloc[0])
        V = U[(U.arm == "V") & (U.A == A)].groupby("unit")[m].mean()
        Tg = {s: g.groupby("unit")[m].mean() for s, g in U[(U.arm == "T") & (U.A == A)].groupby("set")}
        allr = [V] + list(Tg.values()); sp = np.sqrt(sum(((g - g.mean()) ** 2).sum() for g in allr if len(g) > 1) / sum(len(g) - 1 for g in allr if len(g) > 1))
        r = {"A": A, "m": m, "h": h, "V": V.mean()}
        r["T_q"] = np.mean([((o.mean() - V.mean()) ** 2 - sp ** 2 / len(o) - sp ** 2 / len(V)) / h ** 2 for s, o in Tg.items() if "6p5" not in s])
        for arm in ["S", "B"]:
            a = U[(U.arm == arm) & (U.A == A)].groupby("unit")[m].mean()
            r[arm] = a.mean(); r[arm + "_q"] = ((a.mean() - V.mean()) ** 2 - sp ** 2 / len(V) - a.var(ddof=1) / len(a)) / h ** 2
        rows.append(r)
R = pd.DataFrame(rows)
for k in "TSB": R[k + "_d"] = np.sign(R[k + "_q"]) * np.sqrt(abs(R[k + "_q"]))
R.to_csv("results/riesame/B_per_measure_common.csv", index=False)
print({k: round(float(np.sign(R[k + "_q"].mean()) * np.sqrt(abs(R[k + "_q"].mean()))), 6) for k in "TSB"})
