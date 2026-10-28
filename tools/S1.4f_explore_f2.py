#!/usr/bin/env python
"""tools/S1.4f_explore_f2.py - S1.4f, esplorazione dichiarata (docs/steps/S1.4f.md), F2: solo unita' gia' simulate in S1.4e (banco B, bracci
ER ed ET, 300 unita' per braccio), nessuna simulazione nuova. Allineamento (quota dei nuclei allungati >= 1.2 con asse maggiore di traverso
al vicino, psi > 45 gradi) per classe di distanza dal vicino: < 5, 5-6, 6-8 um; stessa misura di P-3 di S1.4e (poly_shape, nn_dir, psi_deg
di tools/S1.4e_measure.py). Serve a scrivere l'attesa della controprova CP-2 (classi non usate nella taratura).
Uscita: results/S1.4f/explore/S1.4f_F2_units.csv. Uso: .venv/bin/python tools/S1.4f_explore_f2.py [ncores]"""
import importlib.util, glob, os, sys
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
_s = importlib.util.spec_from_file_location("m14e", ROOT / "tools/S1.4e_measure.py"); ME = importlib.util.module_from_spec(_s); _s.loader.exec_module(ME)
M = Path("/mnt/micron/geo_spatialtrans/S1.4e"); OUT = ROOT / "results/S1.4f/explore"
CLS = {"nn<5": (0, 5), "nn5-6": (5, 6), "nn6-8": (6, 8)}

def one(arm, key):
    nuc = pd.read_parquet(M / "B" / "nuc" / arm / f"{key}.parquet"); cells = pd.read_parquet(M / "B" / "cells" / f"{key}.parquet")
    S = np.array([ME.poly_shape(g) for g in shapely.from_wkb(nuc.wkb.to_numpy())]); d, j, u = ME.nn_dir(S[:, :2])
    psi = ME.psi_deg(S[:, 2], u[:, 0], u[:, 1]); base = cells.interior.to_numpy() & (S[:, 3] >= 1.2); A, roi, sd = key.split("_"); out = []
    for c, (lo, hi) in CLS.items():
        m = base & (d >= lo) & (d < hi)
        out.append(dict(arm=arm, archetype=A, roi_id=roi, seed=sd, cls=c, n=int(m.sum()), p3=float(np.mean(psi[m] > 45)) if m.any() else np.nan))
    return out

def main(nc):
    from multiprocessing import Pool
    keys = sorted(os.path.basename(f)[:-8] for f in glob.glob(str(M / "B" / "cells" / "*.parquet")))
    with Pool(nc) as p: r = p.starmap(one, [(a, k) for a in ["ER", "ET"] for k in keys])
    D = pd.DataFrame([x for y in r for x in y]); OUT.mkdir(parents=True, exist_ok=True); D.to_csv(OUT / "S1.4f_F2_units.csv", index=False)
    R = D.groupby(["arm", "archetype", "roi_id", "cls"]).p3.mean().groupby(["arm", "archetype", "cls"]).mean().unstack("arm")
    ref = pd.read_csv(ROOT / "results/S1.4f/S1.4f_refs.csv").set_index(["archetype", "cls"]).est
    R["real"] = ref.reindex(R.index); R["frac_real_ER_to_ET"] = (R.real - R.ER) / (R.ET - R.ER); print(R.round(3).to_string())

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
