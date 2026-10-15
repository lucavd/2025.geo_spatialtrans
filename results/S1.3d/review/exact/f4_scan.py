"""review/exact/f4_scan.py — conteggio INDIPENDENTE delle celle F4 invalide: shapely.is_valid su ogni cella di ogni geom/*_f4cells.parquet."""
import glob, os, sys, numpy as np, pandas as pd, pyarrow.parquet as pq, shapely
from concurrent.futures import ProcessPoolExecutor
D = "/mnt/micron/geo_spatialtrans/S1.3d"
def one(f):
    t = pq.read_table(f).to_pandas().sort_values(["cell", "part", "ring", "k"])
    # costruzione vettoriale: anelli -> poligoni -> multipoligoni
    t["rid"] = t.groupby(["cell", "part", "ring"], sort=False).ngroup()
    rings = shapely.linearrings(t[["x", "y"]].values, indices=t.rid.values)
    rr = t.drop_duplicates("rid")[["rid", "cell", "part", "ring"]]
    out = []
    polys = []; pid = []
    for (c, p), d in rr.groupby(["cell", "part"], sort=True):
        ids = d.sort_values("ring").rid.values; polys.append(shapely.polygons(rings[ids[0]], holes=rings[ids[1:]] if len(ids) > 1 else None)); pid.append(c)
    pid = np.array(pid); polys = np.array(polys)
    geoms = shapely.multipolygons(polys, indices=pd.factorize(pid)[0]) if len(np.unique(pid)) < len(pid) else polys
    cells = np.unique(pid) if len(np.unique(pid)) < len(pid) else pid
    # celle a parte singola come poligono (stessa validità), multiple come multipoligono
    v = shapely.is_valid(geoms); bad = cells[~v]
    return os.path.basename(f)[:-len("_f4cells.parquet")], len(cells), [int(b) for b in bad]
fs = sorted(glob.glob(f"{D}/geom/*_f4cells.parquet"))
with ProcessPoolExecutor(12) as ex: R = list(ex.map(one, fs, chunksize=8))
df = pd.DataFrame([dict(key=k, n_cells=n, n_invalid=len(b), invalid=";".join(map(str, b))) for k, n, b in R])
df.to_csv(f"{D}/review_exact/f4/f4_scan_all.csv", index=False)
print("file", len(df), "celle", df.n_cells.sum(), "celle invalide", df.n_invalid.sum(), "casi con invalide", (df.n_invalid > 0).sum())
print(df[df.n_invalid > 0].to_string())
