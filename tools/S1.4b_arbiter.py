#!/usr/bin/env python
"""tools/S1.4b_arbiter.py — S1.4b, C-L2: arbitro indipendente del diagramma di potenza (prototipo tools/S1.4b_laguerre.R).

Punti uniformi nel riquadro della regione; per quelli dentro la regione la cella "vera" e' quella di potenza minima,
trovata con un KD-tree 3D sui generatori sollevati (x, y, sqrt(W - w_i)), W = max w: |p' - c_i'|^2 = pow_i(p) + W.
Confronto con le celle di R (WKB): un punto e' assegnato male se la cella dell'arbitro non lo contiene ed e' a piu'
di TOL = 1e-6 µm dal suo bordo; scoperto se nessuna cella lo contiene (oltre TOL); doppio se due celle lo contengono
oltre TOL dai bordi. Aree: conteggi per cella contro l'area del poligono (z binomiale).
Uso: .venv/bin/python tools/S1.4b_arbiter.py <dir_input> <out.csv> [--mutant MU5] [--n 200000]
"""
import sys, os, glob
import numpy as np, pandas as pd
import shapely
from scipy.spatial import cKDTree
from scipy.stats import norm

TOL = 1e-6


def run_set(f, n_pts, mutant=None, seed=20261016):
    d = pd.read_parquet(f)
    reg = shapely.from_wkb(open(f.replace(".parquet", "_region.wkb"), "rb").read())
    cells = shapely.from_wkb(d["wkb"].to_numpy())
    x = d["x"].to_numpy(); y = d["y"].to_numpy(); w = d["w"].to_numpy()
    W = w.max()
    z = np.sqrt(W + w) if mutant == "MU5" else np.sqrt(W - w)          # MU5: errore nel sollevamento
    tree = cKDTree(np.column_stack([x, y, z]))
    xmin, ymin, xmax, ymax = reg.bounds
    rng = np.random.default_rng(seed + int(d["set_id"].iloc[0]))
    P = np.column_stack([rng.uniform(xmin, xmax, n_pts), rng.uniform(ymin, ymax, n_pts)])
    shapely.prepare(reg)
    inreg = shapely.contains_xy(reg, P[:, 0], P[:, 1])
    Q = P[inreg]
    _, arb = tree.query(np.column_stack([Q, np.zeros(len(Q))]), k=1, workers=8)
    pts = shapely.points(Q)
    st = shapely.STRtree(cells)
    pi, ci = st.query(pts, predicate="intersects")
    ncont = np.bincount(pi, minlength=len(Q))
    own_hit = np.zeros(len(Q), bool); own_hit[pi[ci == arb[pi]]] = True
    empty = shapely.is_empty(cells)
    bd = shapely.boundary(cells)
    # distanza dal bordo della cella dell'arbitro (per i casi dubbi)
    bad = ~own_hit
    dist_b = np.full(len(Q), np.inf)
    if bad.any():
        idx = np.where(bad)[0]; a = arb[idx]
        dist_b[idx] = np.where(empty[a], np.inf, shapely.distance(pts[idx], bd[a]))
    mis = bad & (dist_b > TOL)
    unc = (ncont == 0)
    if unc.any():
        iu = np.where(unc)[0]; near_c = st.query_nearest(pts[iu], return_distance=True)
        dd = np.full(len(iu), np.inf); dd[near_c[0][0]] = near_c[1]
        unc_bad = int((dd > TOL).sum())
    else:
        unc_bad = 0
    dbl = np.where(ncont >= 2)[0]; dbl_bad = 0
    if len(dbl):
        for p in dbl:
            cs = ci[pi == p]
            if np.all(shapely.distance(pts[p], bd[cs]) > TOL): dbl_bad += 1
    # aree: n_i contro N * A_i / A_box (punti uniformi nel riquadro)
    A_box = (xmax - xmin) * (ymax - ymin)
    area = shapely.area(cells); p_i = area / A_box
    n_i = np.bincount(arb, minlength=len(cells))
    with np.errstate(invalid="ignore", divide="ignore"):
        zsc = (n_i - n_pts * p_i) / np.sqrt(n_pts * p_i * (1 - p_i))
    ok = area > 0
    n4 = int((np.abs(zsc[ok]) > 4).sum()); e4 = 2 * norm.sf(4) * ok.sum()
    return dict(set=os.path.basename(f).replace(".parquet", ""), mutant=mutant or "none", n_cells=len(cells),
                n_pts_region=int(len(Q)), misassigned=int(mis.sum()), uncovered=unc_bad, double=dbl_bad,
                n_absz_gt4=n4, exp_absz_gt4=e4, max_absz=float(np.nanmax(np.abs(zsc[ok]))),
                area_pass=bool(n4 <= e4 + 3 * np.sqrt(e4) + 1e-12),
                pass_=bool(mis.sum() == 0 and unc_bad == 0 and dbl_bad == 0 and n4 <= e4 + 3 * np.sqrt(e4) + 1e-12))


if __name__ == "__main__":
    a = sys.argv[1:]; d_in, out = a[0], a[1]
    mut = a[a.index("--mutant") + 1] if "--mutant" in a else None
    n = int(a[a.index("--n") + 1]) if "--n" in a else 200000
    rows = [run_set(f, n, mut) for f in sorted(glob.glob(os.path.join(d_in, "*.parquet")))]
    df = pd.DataFrame(rows); df.to_csv(out, index=False)
    print(df.to_string()); print("C-L2 [%s]: %d/%d PASS" % (mut or "none", df["pass_"].sum(), len(df)))
