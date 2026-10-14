"""tools/S1.3d_coverage.py — S1.3d (BL-084): validità di copertura (GEOS CoverageValidator via shapely, gap_width = 0) dei
territori V0 (tessellate_voronoi) e F4 (taglio unico) di ogni chiave. F4 per O1 = validatore + CoverageUnion: copertura F4
invalida = FAIL di F4 (dev_notes). Uso: .venv/bin/python tools/S1.3d_coverage.py <outdir> <file_chiavi> [--cores 40]"""
import sys, os, json, argparse
from multiprocessing import Pool
import numpy as np, pyarrow.parquet as pq, shapely

def load(path):
    t = pq.read_table(path).to_pandas().sort_values(["cell", "part", "ring", "k"])
    geoms = {}
    for c, dc in t.groupby("cell", sort=True):
        parts = []
        for p, dp in dc.groupby("part", sort=True):
            rings = [list(zip(dr.x.to_numpy(), dr.y.to_numpy())) for r, dr in dp.groupby("ring", sort=True)]
            parts.append(shapely.Polygon(rings[0], rings[1:]))
        geoms[int(c)] = parts[0] if len(parts) == 1 else shapely.MultiPolygon(parts)
    return [geoms[k] for k in sorted(geoms)]

def one(args):
    outdir, key = args; g = os.path.join(outdir, "geom"); res = {"key": key}
    for tag, fn in (("V0", "_cells.parquet"), ("F4", "_f4cells.parquet")):
        p = os.path.join(g, key + fn)
        if not os.path.exists(p): res[tag] = None; continue
        G = np.array(load(p), dtype=object)
        valid = bool(shapely.coverage_is_valid(G, gap_width=0.0))
        inv = shapely.coverage_invalid_edges(G, gap_width=0.0)
        n_bad = int(np.sum(~shapely.is_empty(inv) & ~shapely.is_missing(inv)))
        res[tag] = {"valid": valid, "n_cells_invalid_edges": n_bad, "len_invalid": float(np.nansum(shapely.length(inv)))}
    return res

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("outdir"); ap.add_argument("keys"); ap.add_argument("--cores", type=int, default=40); a = ap.parse_args()
    keys = [k.strip() for k in open(a.keys) if k.strip()]
    out = os.path.join(a.outdir, "coverage"); os.makedirs(out, exist_ok=True)
    todo = [k for k in keys if not os.path.exists(os.path.join(out, k + ".json"))]
    with Pool(a.cores) as p:
        for r in p.imap_unordered(one, [(a.outdir, k) for k in todo]):
            json.dump(r, open(os.path.join(out, r["key"] + ".json"), "w"))
    print(f"coverage: {len(keys)} chiavi ({len(todo)} calcolate), shapely {shapely.__version__} GEOS {shapely.geos_version_string}")
