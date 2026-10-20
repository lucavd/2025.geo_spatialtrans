#!/usr/bin/env python
"""tools/R3b_posthoc_sdarea_windows.py — R3b, POST HOC: frazione di contatti sulle finestre per la variante «CO, taglia dai soli
nuclei di forma StarDist» (tools/R3b_posthoc_sdarea.R). Ciclo copiato da tools/R3b_windows.py (stessa regola dei focali).
Uscita: results/R3b/R3b_posthoc_sdarea_windows.csv.  Uso: .venv/bin/python tools/R3b_posthoc_sdarea_windows.py"""
import importlib.util, json
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py"); K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
R2B = ROOT / "results/R2b"; RES = ROOT / "results/R3b"; D = Path("/mnt/micron/geo_spatialtrans/R3b/sim_out_consensus_sdarea")
SEEDS = [f"s{k:02d}" for k in range(1, 11)]
wins = pd.read_csv(R2B / "R2b_windows.csv"); wins = wins[wins.archetype.isin(["A4", "A6"])]
dm = pd.read_csv(ROOT / "results/R2c/R2c_dmax.csv").set_index("archetype")
wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]
rows = []
for wr in wins.itertuples():
    dmax = float(dm.d_max_um[wr.archetype]); upp = float(wr.um_per_px); x0, y0 = float(wr.x0_um), float(wr.y0_um); s = int(wr.side_px) * upp
    ex = K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
    fs = sorted(D.glob(f"{wr.archetype}_{wr.roi_id}_fold*_s??_K2.parquet")); assert [p.name.split("_")[3] for p in fs] == SEEDS, wr.win_id
    for p in fs:
        t = pd.read_parquet(p)
        keep = np.asarray(np.minimum.reduce([t.cx - x0, x0 + s - t.cx, t.cy - y0, y0 + s - t.cy]), float) >= dmax
        if u is not None:
            keep &= shapely.distance(shapely.points(t[["cx", "cy"]].to_numpy()), u) >= dmax
        rows.append(dict(source="CO_sdarea", win_id=wr.win_id, archetype=wr.archetype, roi=wr.roi_id, seed=p.name.split("_")[3],
                         n=int(keep.sum()), n_contact=int((t.n_c2[keep] > 0).sum())))
r = pd.DataFrame(rows); r["s1"] = r.n_contact / r.n; r.to_csv(RES / "R3b_posthoc_sdarea_windows.csv", index=False)
g = r.groupby(["archetype", "seed"]).s1.mean().groupby("archetype")
print(pd.DataFrame({"mean": g.mean(), "two_se": 2 * g.std(ddof=1) / np.sqrt(10)}).round(4).to_string())
