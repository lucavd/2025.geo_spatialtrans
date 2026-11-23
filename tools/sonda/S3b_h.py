"""tools/sonda/S3b_h.py - S3b (docs/sonda/REGOLA_S3b.md, «h»): margine h per tessuto e misura sulle 5 aree di V (H2) con kappa = 1 e fattore di
asimmetria = 1 (decisione di Luca, 2026-11-23: sola precisione della media di H2). Bootstrap percentile della media come tools/sonda/S3_h.py
(tools/S2_5_common.py corrected_interval, B = 2000, default_rng([20261115, i, j]), i = 1 tumore, 2 stroma, j = posizione della misura): si usano
lo_boot e hi_boot, non i limiti corretti. h = (min(1, hi) - max(0, lo)) / 2 per le quote, (hi - max(0, lo)) / 2 per M-2; M-4 almeno 0,01.
Uso: .venv/bin/python -B tools/sonda/S3b_h.py [area_values.csv] [out.csv]"""
import os, sys, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); import S2_5_common as CM
SEED = 20261115; MAP = {"tumour": ("A2", 1), "stroma": ("A3", 2)}
MEAS = [("M-1", "M1", "bin", True), ("M-2", "M2", "nuc", False), ("M-3", "M3", "nuc", True), ("M-4", "M4", "bin", True), ("M-5", "M5", "bin", True), ("M-6", "M6", "bin", True)]
inp = sys.argv[1] if len(sys.argv) > 1 else "results/sonda/S3b/area_values.csv"; out = sys.argv[2] if len(sys.argv) > 2 else "results/sonda/S3b/h.csv"
U = pd.read_csv(inp); V = U[U.arm == "V"]; rows = []
for A, (Ac, i) in MAP.items():
    R = V[V.A == A].sort_values("unit"); assert len(R) == 5, (A, len(R))
    for j, (mid, col, q, frac) in enumerate(MEAS):
        v = R[col].to_numpy(float); ci = CM.corrected_interval(v, Ac, np.random.default_rng([SEED, i, j]), side="real", q=q)
        lo = max(0.0, ci["lo_boot"]); hi = min(1.0, ci["hi_boot"]) if frac else ci["hi_boot"]; h_raw = (hi - lo) / 2; h = max(h_raw, 0.01) if mid == "M-4" else h_raw
        rows.append({"A": A, "measure": mid, "genes": "all", "n_areas": len(v), "V_mean": ci["est"], "V_sd": float(v.std(ddof=1)), "lo": lo, "hi": hi,
                     "h_raw": h_raw, "h": h, "kappa": 1.0, "kappa_from": "none (S3b)", "cv": ci["cv"], "cv_factor": 1.0})
H = pd.DataFrame(rows); H.to_csv(out, index=False); print(H[["A", "measure", "V_mean", "V_sd", "h"]].round(5).to_string(index=False))
