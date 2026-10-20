#!/usr/bin/env python
"""tools/R3b_posthoc_srprec.py — R3b, POST HOC (dopo la figura di controllo della maschera CO): nelle 8 finestre, quanti oggetti CO
contati (regola di R2b) sono appaiati a un punto di Luca, separatamente per forma StarDist e poligono Space Ranger rimasto; e quanti
punti di Luca sono coperti solo grazie ai poligoni Space Ranger. Uscita: results/R3b/R3b_posthoc_srprec.csv.
Uso (root del repo): .venv/bin/python tools/R3b_posthoc_srprec.py"""
import importlib.util, json, sys
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location("p", ROOT / "tools/R3b_p1p2.py"); P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
MB = P.MB; R3B = P.R3B
wins = pd.read_csv(ROOT / "results/R2b/R2b_windows.csv").set_index("win_id")
rows = []
for wj in json.load(open(ROOT / "results/R2b/annotations/R2b_annotations_Luca.json"))["windows"]:
    if wj["archetype"] not in ("A4", "A6"): continue
    A, roi = wj["archetype"], wj["roi_id"]; wrow = wins.loc[wj["win_id"]]
    lab = np.load(R3B / "masks" / f"{A}_{roi}_consensus_native.npz")["labels"]
    keep = pd.read_parquet(R3B / "nuclei" / f"{A}_{roi}_consensus_native.parquet", columns=["label", "keep"])
    sets = {"consensus": (lab, MB.objects_table(lab))}; keeps = {"consensus": set(keep.label[keep.keep].astype(int))}
    _, mrows, prow, _ = MB.analyse_window(wj, wrow, "Luca", sets, keeps)
    # oggetti contati e appaiati, ricostruiti con le stesse funzioni di analyse_window
    m, side, upp = int(wj["margin_px"]), int(wj["side_px"]), float(wj["um_per_px"]); c0, r0 = int(wrow.c0_roi), int(wrow.r0_roi)
    excl = [(np.array(e["pts"], float) - m).tolist() for e in wj.get("exclusions", [])]; em = MB.poly_mask(excl, side)
    o = MB.window_objects(lab, sets["consensus"][1], c0, r0, side, upp, keeps["consensus"])
    xi = np.clip(np.floor(o.cx).astype(int), 0, side - 1); yi = np.clip(np.floor(o.cy).astype(int), 0, side - 1); o["in_excl"] = em[yi, xi]
    counted = o[o.keep & o.in_win & ~o.in_excl]
    mp = pd.read_parquet(R3B / "masks" / f"{A}_{roi}_consensus_map.parquet").set_index("label").src
    matched = set(prow.label_consensus[prow.match_consensus].astype(int))
    assert len(counted) == mrows[0]["n_obj"]
    for src in ("sd", "sr"):
        cc = counted[counted.label.map(mp) == src]
        rows.append(dict(win_id=wj["win_id"], archetype=A, src=src, n_counted=len(cc), n_matched_to_point=int(cc.label.isin(matched).sum()),
                         n_points_matched_here=int(prow.label_consensus[prow.match_consensus].astype(int).map(mp).eq(src).sum()), n_points=len(prow)))
d = pd.DataFrame(rows); d.to_csv(ROOT / "results/R3b/R3b_posthoc_srprec.csv", index=False)
g = d.groupby(["archetype", "src"])[["n_counted", "n_matched_to_point", "n_points_matched_here", "n_points"]].sum()
g["precision"] = g.n_matched_to_point / g.n_counted; print(g.to_string())
