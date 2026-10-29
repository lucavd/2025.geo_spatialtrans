"""tools/S1.6_explore_parts.py - S1.6, esplorazione E2 (BL-173), prima della pre-registrazione: aree delle parti dei nuclei del braccio U di S1.5e
(nuclei salvati per le 24 ROI con finestra x 10 semi, /mnt/micron/geo_spatialtrans/S1.5e/nuc/U) e del braccio K2 di S1.4e (300 unita', B/nuc/K2).
Uso: tools/S1.6_explore_parts.py <braccio> [cartella dei bracci]. Per ogni nucleo: numero di parti, area della parte
piu' grande, aree delle altre; numero di nuclei di tipo MULTIPOLYGON (BL-173 contava questi). Uscita: results/S1.6/S1.6_explore_parts.csv (per decade dell'area delle parti non principali) e stampa del riepilogo."""
import glob, os, sys, numpy as np, pandas as pd, shapely
from multiprocessing import Pool
ARM = sys.argv[1] if len(sys.argv) > 1 else "U"
def one(f):
    d = pd.read_parquet(f); g = shapely.from_wkb(d.wkb.values)
    parts, idx = shapely.get_parts(g, return_index=True); a = shapely.area(parts)
    df = pd.DataFrame(dict(i=idx, a=a)); df["rank"] = df.groupby("i").a.rank(ascending=False, method="first")
    sec = df[df["rank"] > 1].a.values; n_multi = int((df.groupby("i").size() > 1).sum())
    n_type_multi = int((shapely.get_type_id(g) == 6).sum())       # nuclei di tipo MULTIPOLYGON (anche con una sola parte)
    return len(g), n_multi, sec, float(df[df["rank"] == 1].a.min()), n_type_multi
BASE = sys.argv[2] if len(sys.argv) > 2 else "/mnt/micron/geo_spatialtrans/S1.5e/nuc"
fs = sorted(glob.glob(f"{BASE}/{ARM}/*.parquet"))
with Pool(8) as p: res = p.map(one, fs)
n = sum(r[0] for r in res); nm = sum(r[1] for r in res); sec = np.concatenate([r[2] for r in res]); amin = min(r[3] for r in res)
nt = sum(r[4] for r in res)
print(f"E2 {ARM}: unita' {len(fs)}, nuclei {n}, di tipo MULTIPOLYGON {nt}, nuclei in piu' parti {nm} ({nm / n:.6f}), parti non principali {len(sec)}")
print(f"E2: area minima della parte principale {amin:.4g} um2; parti non principali: max {sec.max():.3g}, min {sec.min():.3g} um2")
edges = [0] + [10.0 ** k for k in range(-16, 3)]
h = np.histogram(sec, bins=edges)[0]
out = pd.DataFrame(dict(arm=ARM, lo=edges[:-1], hi=edges[1:], n_parts=h, n_nuclei=n, n_type_multipolygon=nt, n_multipart=nm)); os.makedirs("results/S1.6", exist_ok=True)
out.to_csv(f"results/S1.6/S1.6_explore_parts_{BASE.split('/')[-3]}_{ARM}.csv", index=False); print(out[out.n_parts > 0].to_string(index=False))
