# Revisione avversariale S1.1 — 02_components.py
# Ricalcolo INDIPENDENTE (scipy.ndimage.label, 4-vicinato, per cluster) di componenti, dimensioni,
# frazione esclusa, composizione del nullo, verita' di I5. Non usa R/04b1_extract_regions.R ne' test_S1.1.R.
import sys, json, numpy as np, pandas as pd
from functools import reduce
from math import gcd
from scipy import ndimage as ndi
EXP = "/mnt/micron/geo_spatialtrans/S1.1/review/export"
OUT = "/mnt/micron/geo_spatialtrans/S1.1/review/py_out"
import os; os.makedirs(OUT, exist_ok=True)
S4 = np.array([[0,1,0],[1,1,1],[0,1,0]])

def load(id_):
    df = pd.read_parquet(f"{EXP}/{id_}.parquet")
    meta = json.load(open(f"{EXP}/{id_}.json"))
    return df, meta

def axis_gcd(v):
    u = np.unique(v)
    if len(u) < 2: return None
    return int(np.gcd.reduce(np.diff(u).astype(np.int64)))

def components(df, stride):
    x = df.x.to_numpy().astype(np.int64); y = df.y.to_numpy().astype(np.int64)
    assert np.all((x - x.min()) % stride == 0) and np.all((y - y.min()) % stride == 0)
    i = (x - x.min()) // stride; j = (y - y.min()) // stride
    lv, code = np.unique(df.cluster.to_numpy(), return_inverse=True)
    G = np.zeros((j.max() + 1, i.max() + 1), dtype=np.int32)
    assert len(np.unique(j * (i.max() + 1) + i)) == len(i), "duplicati"
    G[j, i] = code + 1
    rows = []
    for k in range(1, len(lv) + 1):
        M = G == k
        L, n = ndi.label(M, structure=S4)
        if n == 0: continue
        flat = L.ravel()
        sizes = np.bincount(flat)[1:]
        nz = np.flatnonzero(flat)
        labs, first = np.unique(flat[nz], return_index=True)   # primo pixel in ordine di riga
        pos = nz[first]
        rj, ri = np.divmod(pos, G.shape[1])
        rows.append(pd.DataFrame({"cluster": lv[k - 1], "n_px": sizes[labs - 1],
                                  "rep_x": x.min() + ri * stride, "rep_y": y.min() + rj * stride}))
    comp = pd.concat(rows, ignore_index=True)
    return comp, G, lv

ids = sys.argv[1:]
summary = {}
for id_ in ids:
    df, meta = load(id_)
    px = float(meta["pixel_size_um"])
    sx, sy = axis_gcd(df.x), axis_gcd(df.y)
    stride = 6 if id_.startswith("I6_") else 1          # stride atteso dal disegno, non dai dati
    comp, G, lv = components(df, stride)
    e2 = (stride * px) ** 2
    comp["area_um2"] = comp.n_px * e2
    comp.to_parquet(f"{OUT}/{id_}_components.parquet")
    excl = comp.area_um2 < 100
    tot = comp.area_um2.sum()
    summary[id_] = dict(n_bins=int(len(df)), gcd_x=sx, gcd_y=sy, stride_used=stride, px=px,
                        grid=list(G.shape), n_clusters=int(len(lv)),
                        n_components=int(len(comp)), n_regions_100=int((~excl).sum()),
                        n_excluded_100=int(excl.sum()), frac_area_excluded_100=float(comp.area_um2[excl].sum() / tot),
                        area_total_um2=float(tot), sum_npx=int(comp.n_px.sum()),
                        n_single_px=int((comp.n_px == 1).sum()),
                        size_max=int(comp.n_px.max()), size_median=float(comp.n_px.median()),
                        cluster_counts=df.cluster.value_counts().sort_index().to_dict())
    print(id_, {k: summary[id_][k] for k in ["n_bins","gcd_x","gcd_y","n_components","n_regions_100","frac_area_excluded_100"]}, flush=True)
json.dump(summary, open(f"{OUT}/components_summary_{ids[0]}.json", "w"), indent=1, default=int)

