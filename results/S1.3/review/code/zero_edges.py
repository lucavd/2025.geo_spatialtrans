#!/usr/bin/env python3
"""results/S1.3/review/code/zero_edges.py — lati di lunghezza ~0 nei tile GEOS (keep_tiles) del pacchetto: il conteggio
dei lati usato da S1.3 (nrow(anello) - 1, tools/S1.3_roi.R) li conta come lati veri."""
import numpy as np, pandas as pd, shapely
D = "/mnt/micron/geo_spatialtrans/S1.3/review_code_data"; rows = []
for lab in ["A1_r1", "A2_r1", "A3_r1", "A4_f1", "A5_r1", "A6_r1", "I1", "I2", "I3", "I4", "I5"]:
    tv = pd.read_csv(f"{D}/{lab}_tv.csv"); T = shapely.from_wkb(tv.tile_wkb.values)
    me = np.array([np.hypot(*np.diff(np.asarray(t.exterior.coords), axis=0).T).min() for t in T])
    rows.append(dict(case=lab, n_tiles=len(T), n_min_edge_lt_1e9=int((me < 1e-9).sum()), frac_lt_1e9=float((me < 1e-9).mean()),
                     n_min_edge_lt_1e6=int((me < 1e-6).sum())))
o = pd.DataFrame(rows); o.to_csv("/home/user/2025.geo_spatialtrans/results/S1.3/review/code/zero_edges_tiles.csv", index=False); print(o.to_string())
