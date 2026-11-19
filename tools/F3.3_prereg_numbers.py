"""tools/F3.3_prereg_numbers.py - F3.3: numeri della pre-registrazione, solo da dati di taratura e simulazioni di F3.1.
Per tessuto (classe) e misura: valore vero dentro campione (etichette trasferite per M-3 e M-5), margine di rilevanza h (F3.1;
M-3 e M-5 ricostruiti con la stessa costruzione sulle etichette trasferite; linfonodo: quarti di 500 um dei follicoli, R-9; M-4 almeno
0.01, DEC-086), valore simulato, giudizio dentro campione con la regola V di F3.3 (scarto +- t * errore standard contro +-h) e con la
regola di F3.1, e la semi-ampiezza tipica attesa per ogni insieme di verifica (dalla variabilita' fra aree di taratura).
Uso: .venv/bin/python tools/F3.3_prereg_numbers.py  ->  results/F3.3/F3.3_prereg_numbers.csv"""
import sys, glob, json, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "tools"); import S2_5_common as CM
O = "/mnt/micron/geo_spatialtrans/F3.3/calib"; OUT = "results/F3.3"; import os; os.makedirs(OUT, exist_ok=True)
rows = []
for f in sorted(glob.glob(f"{O}/*.json")): rows += json.load(open(f))
D = pd.DataFrame(rows)
MEAS = [("M-1", "M1", "bin"), ("M-2", "M2", "nuc"), ("M-3", "M3", "nuc"), ("M-4", "M4", "bin"), ("M-5", "M5", "bin"), ("M-6", "M6", "bin")]
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; SEED = 20261115
CLASS = {"A1": "mucosa", "A2": "tumore", "A3": "stroma", "A4": "follicolo_500", "A5": "corteccia", "A6": "miocardio"}
SETS = {"A1": [("L1", 4)], "A2": [("L1", 3), ("P1", 4), ("P5", 3), ("CRC_6p5mm", 5)], "A3": [("P1", 3), ("P5", 1), ("CRC_6p5mm", 3)],
        "A4": [("L1", 2)], "A5": [("L1", 1), ("MouseBrain_FFPE", 4)], "A6": [("L1", 5)]}
I31 = pd.read_csv("results/F3.1/F3.1_intervals.csv").set_index(["archetype", "id"]); V31 = pd.read_csv("results/F3.1/F3.1_verdicts.csv").set_index(["archetype", "id"])
def units(side, A, col, lab):
    d = D[(D.side == side) & (D.A == A) & (D.labels == lab)]
    if A == "A4":
        d = d[d.part.str.startswith("q") & d.tag.str.contains("_f")]; return d.groupby(["tag", "part"])[col].mean().to_numpy()   # quarti, media dei semi
    d = d[d.part == "whole"]; return d.groupby("tag")[col].mean().to_numpy()                                                    # aree, media dei semi
def vrule(g, w, h):
    if abs(g) + w <= h: return "confermata"
    if abs(g) - w > h: return "smentita"
    return "indeterminata"
out = []
for i, A in enumerate(ARCH):
    for j, (mid, col, q) in enumerate(MEAS):
        lab = "tl" if mid in ("M-3", "M-5") else "f31"
        r = units("real", A, col, lab); s = units("sim", A, col, lab); n = len(r)
        if A == "A4" or mid in ("M-3", "M-5"):
            ci = CM.corrected_interval(r, A, np.random.default_rng([SEED, i, j]), side="real", q=q); h = (ci["hi"] - ci["lo"]) / 2
            if mid in ("M-1", "M-3", "M-4", "M-5", "M-6"): h = (min(1.0, ci["hi"]) - ci["lo"]) / 2
            hsrc = "quarti dei follicoli (R-9)" if A == "A4" else "costruzione di F3.1 sulle etichette trasferite"
        else:
            h = float(I31.loc[(A, mid), "half_width"]); hsrc = "F3.1"
        if mid == "M-4": h = max(h, 0.01)
        R = float(r.mean()); S = float(s.mean()); sd = float(r.std(ddof=1)); se_s = float(s.std(ddof=1) / np.sqrt(len(s)))
        t = stats.t.ppf(0.975, n - 1); w_in = t * np.sqrt(se_s ** 2 + sd ** 2 / n); g = S - R
        row = {"archetype": A, "class": CLASS[A], "id": mid, "labels": lab, "n_in": n, "real_in": R, "sd_in": sd, "sim": S, "se_sim": se_s, "h": h, "h_source": hsrc,
               "gap_in": g, "gap_in_over_h": g / h, "w_in": w_in, "verdict_in_V": vrule(g, w_in, h),
               "verdict_f31": V31.loc[(A, mid), "verdict"] if (A, mid) in V31.index else None, "real_in_f31": float(I31.loc[(A, mid), "real_est"])}
        for st, nv in SETS[A]:
            row[f"wD_{st}_n{nv}"] = stats.t.ppf(0.975, n - 1) * sd * np.sqrt(1 / nv + 1 / n)
            row[f"wout_{st}_n{nv}"] = stats.t.ppf(0.975, n - 1) * np.sqrt(se_s ** 2 + sd ** 2 / nv)
        out.append(row)
T = pd.DataFrame(out); T.to_csv(f"{OUT}/F3.3_prereg_numbers.csv", index=False)
pd.set_option("display.width", 250)
print(T[["archetype", "id", "real_in", "real_in_f31", "sim", "h", "gap_in_over_h", "w_in", "verdict_in_V", "verdict_f31"]].round(4).to_string(index=False))
