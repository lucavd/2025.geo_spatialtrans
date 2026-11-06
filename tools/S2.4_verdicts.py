"""tools/S2.4_verdicts.py — S2.4: tabella unica degli esiti (P-1...P-7) e conteggi per previsione ed esito.
P-5 dalle fasce di R-3 (results/S2.4/S2.4_P5_glass.csv): una fascia con una sola tessera non ha incertezza ed e' «non valutabile».
Uscite: results/S2.4/S2.4_verdicts.csv, S2.4_verdict_counts.csv. Uso: .venv/bin/python tools/S2.4_verdicts.py"""
import numpy as np, pandas as pd
P = pd.read_csv("results/S2.4/S2.4_predictions.csv")
G = pd.read_csv("results/S2.4/S2.4_P5_glass.csv")
SEC = {"A1": "intestino", "A23": "tumore e stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
rows = []
for _, g in G.iterrows():
    se = (g.hi - g.lo) / 3.92
    es = "non valutabile" if g.n_tiles < 2 else g.verdict
    rows.append({"id": "P-5", "A": g.key, "tessuto": SEC[g.key], "quantita": f"vetro {g.band} um", "sim": np.nan, "real": np.nan, "stima": g.ratio, "se": se,
                 "lo": 0.8, "hi": 1.25, "n": int(g.n_tiles), "esito": es, "ci_lo": g.lo, "ci_hi": g.hi})
V = pd.concat([P, pd.DataFrame(rows)], ignore_index=True).sort_values(["id", "A"], kind="stable")
V.to_csv("results/S2.4/S2.4_verdicts.csv", index=False)
C = V.groupby(["id", "esito"]).size().unstack(fill_value=0); C.to_csv("results/S2.4/S2.4_verdict_counts.csv")
print(C.to_string())
