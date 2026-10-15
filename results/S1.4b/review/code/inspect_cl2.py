
import numpy as np, pandas as pd, shapely
from scipy.spatial import cKDTree
f = "/mnt/micron/geo_spatialtrans/S1.4b/arb_in/A1_cellpose_rgb_fold3.parquet"
d = pd.read_parquet(f); reg = shapely.from_wkb(open(f.replace(".parquet", "_region.wkb"), "rb").read()); cells = shapely.from_wkb(d.wkb.to_numpy())
x, y, w = d.x.to_numpy(), d.y.to_numpy(), d.w.to_numpy(); W = w.max()
tree = cKDTree(np.column_stack([x, y, np.sqrt(W - w)])); xmin, ymin, xmax, ymax = reg.bounds
for seed_off in [0, 1, 2, 3]:
    rng = np.random.default_rng(20261016 + int(d.set_id.iloc[0]) + seed_off * 7919); N = 200000
    P = np.column_stack([rng.uniform(xmin, xmax, N), rng.uniform(ymin, ymax, N)]); inr = shapely.contains_xy(reg, P[:, 0], P[:, 1]); Q = P[inr]
    _, arb = tree.query(np.column_stack([Q, np.zeros(len(Q))]), k=1)
    # assegnazione alternativa: cella di R che contiene il punto
    st = shapely.STRtree(cells); pi, ci = st.query(shapely.points(Q), predicate="intersects")
    nR = np.bincount(ci, minlength=len(cells))
    A = (xmax - xmin) * (ymax - ymin); p = shapely.area(cells) / A; n_i = np.bincount(arb, minlength=len(cells))
    z = (n_i - N * p) / np.sqrt(N * p * (1 - p)); big = np.where(np.abs(z) > 4)[0]
    print("seed_off", seed_off, "n|z|>4 =", len(big), "pts in region", len(Q), "arb==R-count all cells:", np.array_equal(n_i, nR) if seed_off == 0 else "")
    if seed_off == 0:
        for b in big: print("  cell", b, "z=%.2f n=%d exp=%.1f area=%.2f dist_to_region_boundary=%.2f" % (z[b], n_i[b], N * p[b], p[b] * A, shapely.distance(cells[b].centroid, reg.boundary)))
