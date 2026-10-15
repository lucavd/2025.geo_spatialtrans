
import glob, numpy as np, pandas as pd, shapely
from scipy.spatial import cKDTree
from scipy.stats import binom, norm
rows = []
for f in sorted(glob.glob("/mnt/micron/geo_spatialtrans/S1.4b/arb_in/*.parquet")):
    d = pd.read_parquet(f); xy = np.column_stack([d.x, d.y]); dd, _ = cKDTree(xy).query(xy, k=2)
    reg = shapely.from_wkb(open(f.replace(".parquet", "_region.wkb"), "rb").read()); cells = shapely.from_wkb(d.wkb.to_numpy())
    xmin, ymin, xmax, ymax = reg.bounds; Abox = (xmax - xmin) * (ymax - ymin); p = shapely.area(cells) / Abox; N = 200000
    ok = p > 0; p = p[ok]; mu = N * p; sd = np.sqrt(N * p * (1 - p))
    # P(|z| > 4) esatta (binomiale) per cella
    pe = binom.sf(np.floor(mu + 4 * sd), N, p) + binom.cdf(np.ceil(mu - 4 * sd) - 1, N, p)
    rows.append(dict(set=f.split("/")[-1][:-8], n=len(d), min_nn_um=dd[:, 1].min(), n_tiny_lt_1um2=int((p * Abox < 1).sum()), min_area=float((p * Abox).min()),
                     E_normal=2 * norm.sf(4) * ok.sum(), E_exact=pe.sum(), E_exact_from_tiny=pe[p * Abox < 5].sum()))
r = pd.DataFrame(rows); a = pd.read_csv("results/S1.4b/S1.4b_arbiter.csv"); a["set"] = a["set"].astype(str); r = r.merge(a[["set", "n_absz_gt4", "pass_"]], on="set")
pd.set_option("display.width", 220); print(r.round(4).to_string()); print("TOT obs", r.n_absz_gt4.sum(), "E_normal", r.E_normal.sum().round(2), "E_exact", r.E_exact.sum().round(2), "min nn", r.min_nn_um.min())
r.to_csv("results/S1.4b/review/code/cl2_zstat_exact.csv", index=False)
