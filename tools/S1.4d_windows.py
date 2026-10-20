#!/usr/bin/env python
"""tools/S1.4d_windows.py — S1.4d: frazione di nuclei simulati a contatto sulle finestre di R2c (P-1, P-4, controprova CP-1).
Regola dei focali di R2c (tools/R2c_refs.py, blocco K-2): nuclei simulati con centroide raster nella finestra a distanza >= d_max
del tessuto (results/R2c/R2c_dmax.csv) dai lati e dalle esclusioni annotate da Luca; contatto = n_c2 > 0 (raster S1.4c).
Sorgenti: K2_tissue (S1.4d, phi del tessuto, generatori primari), K2_stardist (S1.4d, generatori StarDist, linfonodo e cuore),
Kc (S1.4c, vecchio modello che scala il nucleo). Uscita: results/S1.4d/S1.4d_windows.csv (una riga per finestra x seme x sorgente).
Uso (root del repo): .venv/bin/python tools/S1.4d_windows.py
"""
import importlib.util, json
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py"); K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
R2B = ROOT / "results/R2b"; RES = ROOT / "results/R2c"; OUT = ROOT / "results/S1.4d"
SRC = {"K2_tissue": (Path("/mnt/micron/geo_spatialtrans/S1.4d/sim_out_tissue"), "K2"),
       "K2_stardist": (Path("/mnt/micron/geo_spatialtrans/S1.4d/sim_out_sd"), "K2"),
       "Kc": (Path("/mnt/micron/geo_spatialtrans/S1.4c/sim_out"), "Kc"),
       "K2_proto": (Path("/mnt/micron/geo_spatialtrans/S1.4c/sim_out"), "K2")}   # controllo: deve riprodurre R2c_k2_windows.csv
SEEDS = [f"s{k:02d}" for k in range(1, 11)]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wins = pd.read_csv(R2B / "R2b_windows.csv"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]
    rows = []
    for src, (d, tag) in SRC.items():
        for wr in wins.itertuples():
            if src == "K2_stardist" and wr.archetype not in ("A4", "A6"):
                continue
            dmax = float(dm.d_max_um[wr.archetype]); upp = float(wr.um_per_px); x0, y0 = float(wr.x0_um), float(wr.y0_um)
            s = int(wr.side_px) * upp
            ex = K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
            fs = sorted(d.glob(f"{wr.archetype}_{wr.roi_id}_fold*_s??_{tag}.parquet"))
            assert [p.name.split("_")[3] for p in fs] == SEEDS, (src, wr.win_id, [p.name for p in fs])
            for p in fs:
                t = pd.read_parquet(p)
                de = np.asarray(np.minimum.reduce([t.cx - x0, x0 + s - t.cx, t.cy - y0, y0 + s - t.cy]), float)
                keep = de >= dmax
                if u is not None:
                    keep &= shapely.distance(shapely.points(t[["cx", "cy"]].to_numpy()), u) >= dmax
                rows.append(dict(source=src, win_id=wr.win_id, archetype=wr.archetype, roi=wr.roi_id, seed=p.name.split("_")[3],
                                 n=int(keep.sum()), n_contact=int((t.n_c2[keep] > 0).sum()), n_conflict=int(t.n_conflict.iloc[0]) if len(t) else 0))
    r = pd.DataFrame(rows); r["s1"] = r.n_contact / r.n
    r.to_csv(OUT / "S1.4d_windows.csv", index=False)
    print(r.groupby(["source", "archetype"]).size().to_string())


if __name__ == "__main__":
    main()
