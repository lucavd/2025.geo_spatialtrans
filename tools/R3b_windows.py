#!/usr/bin/env python
"""tools/R3b_windows.py — R3b: frazione di nuclei simulati a contatto sulle finestre di R2c (P-3, CP-1), regola di
tools/S1.4d_windows.py invariata: nuclei simulati con centroide raster nella finestra a distanza >= d_max del tessuto dai lati e
dalle esclusioni annotate da Luca; contatto = n_c2 > 0 (raster di tools/S1.4c_raster.py sim).
Sorgenti: CO_own (maschera CO, phi della CO), SD_own (StarDist, phi di StarDist), SD_sr (StarDist, phi di Space Ranger = S1.4d:
regressione G-3 contro le righe K2_stardist di results/S1.4d/S1.4d_windows.csv), CP1 (punti di Luca come centri, CO altrove).
Uscite: results/R3b/R3b_windows.csv (finestra x seme x sorgente), results/R3b/R3b_G3_windows.csv
Uso (root del repo): .venv/bin/python tools/R3b_windows.py
"""
import importlib.util, json
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py"); K = importlib.util.module_from_spec(spec); spec.loader.exec_module(K)
R2B = ROOT / "results/R2b"; RES = ROOT / "results/R3b"; M = Path("/mnt/micron/geo_spatialtrans/R3b")
SRC = {"CO_own": M / "sim_out_consensus_own", "SD_own": M / "sim_out_stardist_he_own", "SD_sr": M / "sim_out_stardist_he_sr", "CP1": M / "sim_out_cp1"}
SEEDS = [f"s{k:02d}" for k in range(1, 11)]


def main():
    RES.mkdir(parents=True, exist_ok=True)
    wins = pd.read_csv(R2B / "R2b_windows.csv"); wins = wins[wins.archetype.isin(["A4", "A6"])]
    dm = pd.read_csv(ROOT / "results/R2c/R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]
    rows = []
    for src, d in SRC.items():
        for wr in wins.itertuples():
            dmax = float(dm.d_max_um[wr.archetype]); upp = float(wr.um_per_px); x0, y0 = float(wr.x0_um), float(wr.y0_um)
            s = int(wr.side_px) * upp
            ex = K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
            fs = sorted(d.glob(f"{wr.archetype}_{wr.roi_id}_fold*_s??_K2.parquet"))
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
    r.to_csv(RES / "R3b_windows.csv", index=False)
    ref = pd.read_csv(ROOT / "results/S1.4d/S1.4d_windows.csv"); ref = ref[ref.source == "K2_stardist"]
    a = r[r.source == "SD_sr"].set_index(["win_id", "seed"]).sort_index(); b = ref.set_index(["win_id", "seed"]).sort_index()
    g3 = pd.DataFrame([dict(check="S1.4d_windows K2_stardist (n, n_contact)", n=len(a), n_ref=len(b),
                            equal=bool(len(a) == len(b) and a.index.equals(b.index) and np.array_equal(a.n, b.n) and np.array_equal(a.n_contact, b.n_contact)))])
    g3.to_csv(RES / "R3b_G3_windows.csv", index=False)
    print(r.groupby(["source", "archetype", "seed"]).s1.mean().groupby(["source", "archetype"]).agg(["mean", "std"]).round(4).to_string())
    print(g3.to_string(index=False))


if __name__ == "__main__":
    main()
