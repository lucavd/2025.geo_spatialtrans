"""tools/S1.6_unit_metrics.py - S1.6: copia di tools/S1.5e_unit_metrics.py con le cartelle di S1.6 (stesse funzioni). Testo originale:
S1.5e: metriche per unita' simulata, come tools/S1.5d_unit_metrics.py (stesse funzioni di tools/S1.5c_metrics.py e
tools/S1.5d_metrics.py: P-5...P-7, P-11...P-14), piu' densita' del ROI (nuclei / mm^2 della maschera di tessuto valido, P-15) e forma dei vuoti
(tools/S1.5e_metrics.void_shape con le celle interne della regola di R3 scritte dal banco: P-16, P-17). Stadio «real»: le stesse funzioni sui nuclei
veri (results/S1.5e/S1.5e_real_void_shape.csv, per G-4). Ingresso: /mnt/micron/geo_spatialtrans/S1.5e/cells/<braccio>/*.parquet. Uscita: un csv per
unita' in .../S1.5e/umet/<braccio>/ (rilanciabile) e results/S1.5e/S1.5e_unit_metrics.csv.
Uso: .venv/bin/python tools/S1.5e_unit_metrics.py [worker] [bracci separati da virgole] | real"""
import sys, os, glob, importlib.util, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial import cKDTree
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
m5c = mod("m5c", "tools/S1.5c_metrics.py"); m5d = mod("m5d", "tools/S1.5d_metrics.py"); m5e = mod("m5e", "tools/S1.5e_metrics.py"); cal = mod("cal", "tools/S1.5b_calibration.py")
OUT = "/mnt/micron/geo_spatialtrans/S1.6"; RES = "results/S1.6"
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
GEOM = {}
def geom(A, roi): return GEOM.get((A, roi)) or GEOM.setdefault((A, roi), m5e.RoiGeomE(A, roi))
def one(f):
    arm = f.split("/")[-2]; A, roi, sd = os.path.basename(f)[:-8].split("_"); o = f"{OUT}/umet/{arm}/{A}_{roi}_{sd}.csv"
    if os.path.exists(o): return o
    e = pd.read_parquet(f); I = e.interior.values.astype(bool); rows = []
    for k, v in m5c.roi_metrics(e.x.values, e.y.values, e.cell_type.values, I, MAIN[A]).items():
        rows.append(dict(metric=k[0], cell_type=k[1], bin=k[2] if len(k) > 2 else "", value=v) if isinstance(k, tuple) else dict(metric=k, cell_type="", bin="", value=v))
    g = geom(A, roi)
    for k, v in m5d.position_metrics(e.x.values, e.y.values, I, e.a_t.values, g).items(): rows.append(dict(metric=k, cell_type="", bin="", value=v))
    for k, v in m5e.void_shape(e.x.values, e.y.values, I, g).items(): rows.append(dict(metric=k, cell_type="", bin="", value=v))
    rows.append(dict(metric="roi_density", cell_type="", bin="", value=m5e.roi_density(len(e), g)))
    for t in MAIN[A]:
        k = I & (e.cell_type.values == t)
        if k.sum(): rows.append(dict(metric="q90_an", cell_type=t, bin="", value=float(np.quantile(e.a_n.values[k], 0.9))))
    d = pd.DataFrame(rows).assign(arm=arm, archetype=A, roi_id=roi, seed=int(sd[1:])); os.makedirs(os.path.dirname(o), exist_ok=True)
    d.to_csv(o + ".tmp", index=False); os.replace(o + ".tmp", o); return o
def real_one(args):
    A, roi = args; _, d = cal.load(A); e = d[d.roi_id == roi]
    return dict(archetype=A, roi_id=roi, **m5e.void_shape(e.x.values, e.y.values, e.interior.values.astype(bool), geom(A, roi)))
if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "real":
        R = m5d.rois_table()
        with Pool(6) as p: r = pd.DataFrame(p.map(real_one, list(zip(R.archetype, R.roi_id))))
        r.to_csv(f"{RES}/S1.5e_real_void_shape.csv", index=False); print(r.to_string()); sys.exit(0)
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    arms = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    fs = sorted(glob.glob(f"{OUT}/cells/*/*.parquet")); fs = [f for f in fs if arms is None or f.split("/")[-2] in arms]
    with Pool(nw) as p: p.map(one, fs, chunksize=4)
    allf = sorted(glob.glob(f"{OUT}/umet/*/*.csv"))
    pd.concat([pd.read_csv(o, keep_default_na=False, na_values=[""]) for o in allf]).to_csv(f"{RES}/S1.6_unit_metrics.csv", index=False); print(len(allf), "unita'")
