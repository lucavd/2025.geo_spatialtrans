"""tools/F3.3_prereg_precision.py - F3.3: precisione attesa di S* e G* sotto nessuno spostamento, solo dai dati di taratura
(modello parametrico con la deviazione standard fra aree dentro campione e le numerosita' degli insiemi di verifica). Per scrivere
i «2 errori standard tipici» della pre-registrazione. Uso: .venv/bin/python tools/F3.3_prereg_precision.py"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "tools"); import F3_3_stats as F
T = pd.read_csv("results/F3.3/F3.3_prereg_numbers.csv")
SETS = {"A1": {"L1": [("L1", 4)]}, "A2": {"L1": [("L1", 3)], "L2": [("P1", 4), ("P5", 3), ("CRC_6p5mm", 5)]}, "A3": {"L2": [("P1", 3), ("P5", 1), ("CRC_6p5mm", 3)]},
        "A4": {"L1": [("L1", 2)]}, "A5": {"L1": [("L1", 1)], "L2": [("MouseBrain_FFPE", 4)]}, "A6": {"L1": [("L1", 5)]}}
rows = []
for lev in ("L1", "L2"):
    cells = []
    for _, r in T.iterrows():
        st = SETS[r.archetype].get(lev)
        if not st: continue
        m0 = r.real_in; n = int(r.n_in)
        cells.append({"in": np.full(n, m0), "out": {nm: np.full(nv, m0) for nm, nv in st}, "h": r.h, "sd": r.sd_in, "df": n - 1, "sim": r.sim,
                      "good": r.verdict_in_V != "smentita", "id": r.id, "A": r.archetype})
    for nm, fun in [("S", F.stat_S), ("G", F.stat_G)]:
        c = F.corrected(cells, fun, "param", 2000, 20261119)
        rng = np.random.default_rng(20261120); sims = np.array([fun(F.null_draw(cells, rng, "param")) for _ in range(2000)])
        rows.append({"level": lev, "stat": nm, "n_cells": len(cells), "n_good": sum(c_["good"] for c_ in cells),
                     "null_median": c["null_median"], "null_q025": np.percentile(sims, 2.5), "null_q975": np.percentile(sims, 97.5),
                     "half_width_typical": (np.percentile(sims, 97.5) - np.percentile(sims, 2.5)) / 2})
R = pd.DataFrame(rows); R.to_csv("results/F3.3/F3.3_prereg_precision.csv", index=False); print(R.round(3).to_string(index=False))
