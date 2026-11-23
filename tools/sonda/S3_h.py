"""tools/sonda/S3_h.py - S3: margine h per tessuto e misura sulle 5 aree di V (H2), costruzione di tools/riesame_B_h.py: bootstrap della media
(tools/S2_5_common.py corrected_interval, B = 2000, default_rng([20261115, i, j])), kappa e fattore di asimmetria di R7 del colon: tumore come A2
(i = 1), stroma come A3 (i = 2), j = posizione della misura; q = bin per M-1, M-4, M-5, M-6, nuc per M-2, M-3; quote: hi <= 1; h = (hi - lo) / 2;
M-4 almeno 0,01. Uso: .venv/bin/python -B tools/sonda/S3_h.py [area_values.csv] [out.csv]"""
import os, sys, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); import S2_5_common as CM
SEED = 20261115; MAP = {"tumour": ("A2", 1), "stroma": ("A3", 2)}
MEAS = [("M-1", "M1", "bin", True), ("M-2", "M2", "nuc", False), ("M-3", "M3", "nuc", True), ("M-4", "M4", "bin", True), ("M-5", "M5", "bin", True), ("M-6", "M6", "bin", True)]
inp = sys.argv[1] if len(sys.argv) > 1 else "results/sonda/S3/area_values.csv"; out = sys.argv[2] if len(sys.argv) > 2 else "results/sonda/S3/h.csv"
U = pd.read_csv(inp); V = U[U.arm == "V"]; rows = []
for A, (Ac, i) in MAP.items():
    R = V[V.A == A].sort_values("unit"); assert len(R) == 5, (A, len(R))
    for j, (mid, col, q, frac) in enumerate(MEAS):
        v = R[col].to_numpy(float); ci = CM.corrected_interval(v, Ac, np.random.default_rng([SEED, i, j]), side="real", q=q)
        hi = min(1.0, ci["hi"]) if frac else ci["hi"]; h_raw = (hi - ci["lo"]) / 2; h = max(h_raw, 0.01) if mid == "M-4" else h_raw
        rows.append({"A": A, "measure": mid, "genes": "all", "n_areas": len(v), "V_mean": ci["est"], "V_sd": float(v.std(ddof=1)), "lo": ci["lo"], "hi": hi,
                     "h_raw": h_raw, "h": h, "kappa": ci["kappa"], "kappa_from": Ac, "cv": ci["cv"], "cv_factor": ci["cv_factor"]})
H = pd.DataFrame(rows); H.to_csv(out, index=False); print(H[["A", "measure", "V_mean", "V_sd", "h", "kappa", "cv_factor"]].round(5).to_string(index=False))
