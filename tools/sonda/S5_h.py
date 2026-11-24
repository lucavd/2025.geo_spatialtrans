"""tools/sonda/S5_h.py - S5 (docs/sonda/REGOLA_S5.md): h_p per tessuto e misura sulle 5 aree di Vp con tools/sonda/S3b_h.py invariato (kappa = 1,
B = 2000, default_rng([20261115, i, j])) sulla tabella del pannello; la sola colonna genes passa da "all" a "panel".
Uscita: results/sonda/S5/h.csv. Uso: .venv/bin/python -B tools/sonda/S5_h.py"""
import os, sys, subprocess, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); out = "results/sonda/S5/h.csv"
subprocess.run([sys.executable, "-B", "tools/sonda/S3b_h.py", "results/sonda/S5/area_values.csv", out], check=True)
H = pd.read_csv(out); assert set(H.genes) == {"all"}; H["genes"] = "panel"; H.to_csv(out, index=False)
