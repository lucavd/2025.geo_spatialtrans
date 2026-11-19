"""tools/F3.3_prereg_precision_q.py - F3.3, correzione del 2026-11-19 (cancello G-1 f): precisione attesa delle statistiche quadratiche
(radice con segno del quadrato medio dello spostamento vero e dell'aumento del quadrato dello scarto) sotto nessuno spostamento e con uno
spostamento costruito di 1 h in ogni cella, solo dai dati di taratura (deviazione standard fra aree dentro campione, numerosita' degli
insiemi di verifica). Uso: .venv/bin/python tools/F3.3_prereg_precision_q.py"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "tools"); import F3_3_stats as F
T = pd.read_csv("results/F3.3/F3.3_prereg_numbers.csv")
SETS = {"A1": {"L1": [("L1", 4)]}, "A2": {"L1": [("L1", 3)], "L2": [("P1", 4), ("P5", 3), ("CRC_6p5mm", 5)]}, "A3": {"L2": [("P1", 3), ("P5", 1), ("CRC_6p5mm", 3)]},
        "A4": {"L1": [("L1", 2)]}, "A5": {"L1": [("L1", 1)], "L2": [("MouseBrain_FFPE", 4)]}, "A6": {"L1": [("L1", 5)]}}
rng = np.random.default_rng(20261121); rows = []
for lev in ("L1", "L2"):
    for shift in (0.0, 1.0):
        est = {"S": [], "G": []}
        for rep in range(300):
            cells = []
            for _, r in T.iterrows():
                st = SETS[r.archetype].get(lev)
                if not st: continue
                n = int(r.n_in); s = r.sd_in; m0 = r.real_in
                cells.append({"in": rng.normal(m0, s, n), "out": {nm: rng.normal(m0 + shift * r.h * rng.choice([-1, 1]), s, nv) for nm, nv in st},
                              "h": r.h, "sd": s, "df": n - 1, "sim": m0 + r.gap_in if r.verdict_in_V == "smentita" else m0, "good": r.verdict_in_V != "smentita", "se_sim": r.se_sim})
            est["S"].append(F.ssqrt(F.stat_QS(cells))); est["G"].append(F.ssqrt(F.stat_QG(cells)))
        for k, v in est.items():
            rows.append({"level": lev, "stat": k, "true_shift_h": shift, "median": np.median(v), "q025": np.percentile(v, 2.5), "q975": np.percentile(v, 97.5)})
R = pd.DataFrame(rows); R.to_csv("results/F3.3/F3.3_prereg_precision_q.csv", index=False); print(R.round(3).to_string(index=False))
