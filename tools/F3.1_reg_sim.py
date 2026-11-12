"""tools/F3.1_reg_sim.py - F3.1: aggiunge in_reg (definizione reg di S2.5) ai quadratini simulati esportati.
Centro del quadratino nelle regioni dello Step 1 dell'area e a >= 25 um dal bordo del riquadro (come tools/S2.5_stats.py).
Uso: .venv/bin/python tools/F3.1_reg_sim.py <area_dir> [<area_dir> ...]  (area_dir = .../<tag>_k<k>)"""
import sys, os, json, importlib.util, numpy as np, pandas as pd, shapely
def _load(name, path):
    s = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
RL = _load("s23real", "tools/S2.3_real.py"); FRAME_UM = 25.0
for d in sys.argv[1:]:
    info = json.load(open(f"{d}/export.json")); tag = info["tag"]; side = float(info["side_um"])
    b = pd.read_parquet(f"{d}/bins.parquet"); x, y = b.x_um.to_numpy(), b.y_um.to_numpy()
    fr = (y >= FRAME_UM) & (y <= side - FRAME_UM) & (x >= FRAME_UM) & (x <= side - FRAME_UM)
    P = np.zeros(len(b), bool); P[fr] = shapely.contains_xy(RL.region_geom(tag), x[fr], y[fr])
    b["in_reg"] = P; b.to_parquet(f"{d}/bins.parquet", index=False)
    info["n_reg"] = int(P.sum()); json.dump(info, open(f"{d}/export.json", "w")); print(tag, "n_reg", int(P.sum()), flush=True)
