#!/usr/bin/env python
"""tools/S1.4e_posthoc_raters.py - S1.4e, post hoc PH-1 (dopo il risultato, dichiarato): contatto per coppia secondo Luca e secondo Federica
sulle stesse coppie di R2c (le finestre giudicate da entrambi), per classe di distanza fra i punti di Luca; coppie dei due tessuti con l'eccesso
residuo (tumore del colon, stroma) che Federica non ha giudicato; tempo mediano per item dell'annotatore (t_ms).
Uscita: results/S1.4e/S1.4e_PH1_raters.csv.  Uso: .venv/bin/python tools/S1.4e_posthoc_raters.py"""
import json, numpy as np, pandas as pd
def load(r):
    it = pd.read_csv(f"results/R2c/R2c_items_{r}.csv"); a = pd.DataFrame(json.load(open(f"results/R2c/annotations/R2c_answers_{r}.json"))["answers"])
    return it.merge(a, on="item", validate="1:1")
L, F = load("Luca"), load("Federica"); rows = []
for r, d in [("Luca", L), ("Federica", F)]: rows.append(dict(what="t_ms mediano per item", rater=r, archetype="", cls="", n=len(d), value=float(d.t_ms.median())))
Lp, Fp = L[L.kind == "pair"], F[F.kind == "pair"]
m = Lp.merge(Fp[["pair_id", "answer"]], on="pair_id", suffixes=("_L", "_F")); m["A"] = m.win_id.str[:2]
for A in ["A1", "A2", "A3", "A4"]:
    for lo, hi, c in [(0, 6, "lt6"), (6, 10, "6to10")]:
        g = m[(m.A == A) & (m.d_um >= lo) & (m.d_um < hi)]
        for r, col in [("Luca", "answer_L"), ("Federica", "answer_F")]:
            j = g[g[col] != "unsure"]; rows.append(dict(what="contatto sulle coppie comuni", rater=r, archetype=A, cls=c, n=len(j), value=float((j[col] == "touch").mean()) if len(j) else np.nan))
Lp = Lp.assign(A=Lp.win_id.str[:2]); fw = set(Fp.win_id)
for A in ["A2", "A3"]:
    rest = Lp[(Lp.A == A) & ~Lp.win_id.isin(fw)]
    rows.append(dict(what="coppie di Luca in finestre senza Federica", rater="Luca", archetype=A, cls=",".join(sorted(rest.win_id.unique())), n=len(rest), value=np.nan))
d = pd.DataFrame(rows); d.to_csv("results/S1.4e/S1.4e_PH1_raters.csv", index=False); print(d.round(3).to_string())
