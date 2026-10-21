#!/usr/bin/env python
"""tools/S1.5_windows.py - S1.5, P-5: frazione di nuclei simulati a contatto sulle finestre di R2c, simulatore completo.
Stessa regola di tools/S1.4d_windows.py (focali di R2c: centroide raster nella finestra a distanza >= d_max del tessuto
dai lati e dalle esclusioni annotate da Luca; contatto = n_c2 > 0 del raster S1.4c), applicata ai nuclei di place_cells()
rasterizzati da tools/S1.4c_raster.py sim. Sorgenti: i due bracci del modello di taglia (preserve_area, none).
Uscita: results/S1.5/S1.5_windows.csv (una riga per finestra x seme x braccio). Uso: .venv/bin/python tools/S1.5_windows.py
"""
import importlib.util, json
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py"); K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
R2B = ROOT / "results/R2b"; RES = ROOT / "results/R2c"; OUT = ROOT / "results/S1.5"; RAS = Path("/mnt/micron/geo_spatialtrans/S1.5/raster")
SEEDS = [f"s{k:02d}" for k in range(1, 11)]

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wins = pd.read_csv(R2B / "R2b_windows.csv"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]
    rows = []
    import os
    arms = os.environ.get("S15_ARMS", "preserve_area,none").split(",")
    for arm in arms:
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
                rows.append(dict(arm=arm, win_id=wr.win_id, archetype=wr.archetype, roi=wr.roi_id, seed=p.stem.split("_")[2],
                                 n=int(keep.sum()), n_contact=int((t.n_c2[keep] > 0).sum()), n_conflict=int(t.n_conflict.iloc[0]) if len(t) else 0))
    r = pd.DataFrame(rows); r["s1"] = r.n_contact / r.n
    r.to_csv(OUT / os.environ.get("S15_WIN_OUT", "S1.5_windows.csv"), index=False)
    print(r.groupby(["arm", "archetype"]).size().to_string())

if __name__ == "__main__":
    main()
