#!/usr/bin/env python
"""tools/S1.5e_windows.py - S1.5e, P-1 e CP-3: frazione di nuclei simulati a contatto sulle finestre di R2c, e contatto per classe di distanza dal
vicino (classi di tools/S1.5e_explore_contacts.py). Stessa regola dei focali di tools/S1.5d_windows.py (centroide raster a distanza >= d_max del
tessuto dai lati e dalle esclusioni annotate da Luca; contatto = n_c2 > 0 del raster S1.4c); distanza dal vicino fra tutti i centri del raster.
Bracci: variabile S15E_ARMS. Uscite: results/S1.5e/S1.5e_windows.csv (finestra x seme x braccio), S1.5e_windows_nn.csv (anche per classe).
Uso: S15E_ARMS=U,RV,X,XA,XAD .venv/bin/python tools/S1.5e_windows.py"""
import importlib.util, json, os
import numpy as np, pandas as pd, shapely
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py"); K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
R2B = ROOT / "results/R2b"; RES = ROOT / "results/R2c"; OUT = ROOT / "results/S1.5e"; RAS = Path("/mnt/micron/geo_spatialtrans/S1.5e/raster")
SEEDS = [f"s{k:02d}" for k in range(1, 11)]; BINS = [0, 5, 6, 7, 8, 10, np.inf]
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wins = pd.read_csv(R2B / "R2b_windows.csv"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]
    rows, nrows = [], []
    for arm in os.environ.get("S15E_ARMS", "U,RV,X,XA,XAD").split(","):
        for wr in wins.itertuples():
            dmax = float(dm.d_max_um[wr.archetype]); upp = float(wr.um_per_px); x0, y0 = float(wr.x0_um), float(wr.y0_um)
            s = int(wr.side_px) * upp
            ex = K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
            fs = sorted((RAS / arm).glob(f"{wr.archetype}_{wr.roi_id}_s??.parquet"))
            assert [p.stem.split("_")[2] for p in fs] == SEEDS, (arm, wr.win_id, [p.name for p in fs])
            for p in fs:
                t = pd.read_parquet(p)
                de = np.asarray(np.minimum.reduce([t.cx - x0, x0 + s - t.cx, t.cy - y0, y0 + s - t.cy]), float)
                keep = de >= dmax
                if u is not None:
                    keep &= shapely.distance(shapely.points(t[["cx", "cy"]].to_numpy()), u) >= dmax
                sd = p.stem.split("_")[2]; c = (t.n_c2.values > 0)
                rows.append(dict(arm=arm, win_id=wr.win_id, archetype=wr.archetype, roi=wr.roi_id, seed=sd,
                                 n=int(keep.sum()), n_contact=int(c[keep].sum()), n_conflict=int(t.n_conflict.iloc[0]) if len(t) else 0))
                if keep.any():
                    P = t[["cx", "cy"]].to_numpy(); dd, _ = cKDTree(P).query(P[keep], k=2); b = np.digitize(dd[:, 1], BINS[1:-1])
                    for k in range(len(BINS) - 1):
                        m = b == k; nrows.append(dict(arm=arm, win_id=wr.win_id, archetype=wr.archetype, seed=sd, bin=k, n=int(m.sum()), n_contact=int(c[keep][m].sum())))
    r = pd.DataFrame(rows); r["s1"] = r.n_contact / r.n
    r.to_csv(OUT / "S1.5e_windows.csv", index=False); pd.DataFrame(nrows).to_csv(OUT / "S1.5e_windows_nn.csv", index=False)
    print(r.groupby(["arm", "archetype"]).size().to_string())
if __name__ == "__main__":
    main()
