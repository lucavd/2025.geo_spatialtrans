#!/usr/bin/env python
"""tools/R3b_posthoc_valid.py — R3b, dopo la revisione (rilievo F-1 del revisore «ricalcolo»): etichette con centroide fuori dal
tessuto valido nelle maschere grezze (Space Ranger, StarDist, CO), dopo il filtro cp2 (usato da P-2) e fra i generatori (R3b/R3/real).
Uscita: results/R3b/R3b_posthoc_valid.csv.  Uso (root del repo): .venv/bin/python tools/R3b_posthoc_valid.py"""
import importlib.util, sys
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "tools"))
def _load(n, f):
    s = importlib.util.spec_from_file_location(n, ROOT / f); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
K = _load("k", "tools/R2c_candidates.py"); B = _load("b", "tools/R3b_consensus.py")
R3B = Path("/mnt/micron/geo_spatialtrans/R3b"); rois = K.S.rois_table()
def n_out(lab, upp, valid, side_um):
    f = B.cp2_filter(lab, upp, valid, side_um); a = set(np.unique(lab[lab > 0])); b = set(np.unique(f[f > 0])); return len(a - b), len(a)
rows = []
for A, roi in B.roi_list():
    rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]; upp = float(rw.um_per_px); side_um = float(rw.side_px) * upp
    valid = K.S.read_valid(A, roi)
    for s in ("spaceranger", "stardist_he", "consensus"):
        lab = np.load(R3B / "masks" / f"{A}_{roi}_consensus_native.npz")["labels"] if s == "consensus" else B.load_native(A, roi, s)
        o, n = n_out(lab.astype(np.int64), upp, valid, side_um)
        f = B.cp2_filter(lab.astype(np.int64), upp, valid, side_um); o2, _ = n_out(f, upp, valid, side_um)
        g = pd.read_parquet(R3B / "R3/real" / f"{A}_{roi}_{s}_gen.parquet") if s != "spaceranger" else None
        gl = set(g.label.astype(int)) if g is not None else set()
        go = len(gl - set(np.unique(f[f > 0]))) if g is not None else np.nan
        rows.append(dict(archetype=A, roi=roi, mask=s, n_labels=n, n_out_raw=o, n_out_after_cp2=o2, n_gen=len(gl) if g is not None else np.nan, n_gen_out=go))
d = pd.DataFrame(rows); d.to_csv(ROOT / "results/R3b/R3b_posthoc_valid.csv", index=False)
print(d.groupby("mask")[["n_labels", "n_out_raw", "n_out_after_cp2", "n_gen", "n_gen_out"]].sum().to_string())
