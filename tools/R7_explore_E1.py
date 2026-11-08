"""tools/R7_explore_E1.py — R7, esplorazione E-1 (prima della pre-registrazione).
Correlazione fra le medie di due aree di 1 mm^2 con lo stesso scarto di una coppia di R2,
stimata spostando la coppia nel tessuto intorno alle aree (campo: UMI per bin di 8 um, Space Ranger).
Uscita: results/R7/explore/R7_E1_pairs.csv, R7_E1_tissue.csv"""
import os, sys, itertools, json
import numpy as np, pandas as pd, pyarrow.parquet as pq
sys.path.insert(0, "tools"); import R6_common as C
OUT = "results/R7/explore"; os.makedirs(OUT, exist_ok=True)
DR = "data/real/datasets"
STEP_UM = 50.0      # passo delle posizioni della coppia
CELL_UM = 8.0       # cella della griglia di accumulo (un bin di 8 um)
DIL_UM = 1000.0     # intorno: unione delle aree dilatata di 1 mm
MIN_TISSUE = 0.5    # frazione minima di celle con un bin in tessuto in ciascun quadrato
R = pd.read_csv("results/R2/R2_rois_checked.csv")
def integral(a):
    s = np.zeros((a.shape[0]+1, a.shape[1]+1)); s[1:,1:] = a.cumsum(0).cumsum(1); return s
def box_sum(S, r0, c0, h, w):
    return S[r0+h, c0+w] - S[r0, c0+w] - S[r0+h, c0] + S[r0, c0]
rows_pair, rows_t = [], []
for A, g in R.groupby("archetype"):
    g = g.reset_index(drop=True); ds = g.dataset[0]; upp = g.um_per_px[0]
    sq = f"{DR}/{ds}/binned_outputs_x/binned_outputs/square_008um"
    X, bc, _, _ = C.read_h5_csc(f"{sq}/filtered_feature_bc_matrix.h5")
    umi = np.asarray(X.sum(axis=0)).ravel() if X.shape[1] == len(bc) else np.asarray(X.sum(axis=1)).ravel()
    tp = pq.read_table(f"{sq}/spatial/tissue_positions.parquet").to_pandas()
    tp = tp[tp.in_tissue == 1].set_index("barcode").loc[[b.decode() if isinstance(b, bytes) else b for b in bc]]
    cpx = CELL_UM / upp
    rr = np.floor(tp.pxl_row_in_fullres.values / cpx).astype(int); cc = np.floor(tp.pxl_col_in_fullres.values / cpx).astype(int)
    H, W = rr.max()+2, cc.max()+2
    U = np.zeros((H, W)); N = np.zeros((H, W))
    np.add.at(U, (rr, cc), umi); np.add.at(N, (rr, cc), 1.0)
    N = np.minimum(N, 1.0)  # una cella = un bin; doppi rari per l'arrotondamento
    SU, SN = integral(U), integral(N)
    side = int(round(1000.0 / CELL_UM))
    box = [(int(round(x.r0 / cpx)), int(round(x.c0 / cpx))) for _, x in g.iterrows()]
    dil = int(round(DIL_UM / CELL_UM))
    rmin = max(0, min(b[0] for b in box) - dil); cmin = max(0, min(b[1] for b in box) - dil)
    rmax = min(H - side - 1, max(b[0] for b in box) + dil); cmax = min(W - side - 1, max(b[1] for b in box) + dil)
    step = int(round(STEP_UM / CELL_UM))
    def X_at(r0, c0):
        n = box_sum(SN, r0, c0, side, side)
        return (box_sum(SU, r0, c0, side, side) / n) if n >= MIN_TISSUE * side * side else np.nan
    x_real = np.array([X_at(*b) for b in box])
    # anchors
    grid = [(r0, c0) for r0 in range(rmin, rmax, step) for c0 in range(cmin, cmax, step)]
    Xg = {p: X_at(*p) for p in grid}
    vals = np.array([v for v in Xg.values() if np.isfinite(v)])
    rho = np.eye(len(g))
    for i, j in itertools.combinations(range(len(g)), 2):
        dr, dc = box[j][0]-box[i][0], box[j][1]-box[i][1]
        a, b = [], []
        for (r0, c0) in grid:
            r1, c1 = r0+dr, c0+dc
            if not (0 <= r1 < H-side-1 and 0 <= c1 < W-side-1): continue
            xa = Xg[(r0, c0)]; xb = X_at(r1, c1)
            if np.isfinite(xa) and np.isfinite(xb): a.append(xa); b.append(xb)
        a, b = np.array(a), np.array(b)
        rr_ = float(np.corrcoef(a, b)[0, 1]) if len(a) > 10 else np.nan
        f = max(0, side-abs(dr)) * max(0, side-abs(dc)) / side**2
        rho[i, j] = rho[j, i] = rr_
        rows_pair.append(dict(archetype=A, i=g.roi_id[i], j=g.roi_id[j], d_um=float(np.hypot(dr, dc)*CELL_UM), f_box=f, rho_slide=rr_, n_pos=len(a)))
    n = len(g); off = rho[~np.eye(n, dtype=bool)]; rb = np.nanmean(off)
    Fm = np.eye(n)
    for p in rows_pair:
        if p["archetype"] == A:
            i = list(g.roi_id).index(p["i"]); j = list(g.roi_id).index(p["j"]); Fm[i, j] = Fm[j, i] = p["f_box"]
    fb = Fm[~np.eye(n, dtype=bool)].mean()
    k = lambda r: float(np.sqrt((1+(n-1)*r)/(1-r)))
    rows_t.append(dict(archetype=A, n_anchor=len(vals), sd_anchor=float(np.std(vals, ddof=1)), sd_real5=float(np.std(x_real, ddof=1)),
                       mean_real5=float(np.mean(x_real)), rho_bar_slide=float(rb), f_bar=float(fb), k_slide=k(rb), k_box=k(fb)))
    print(A, json.dumps(rows_t[-1]), flush=True)
pd.DataFrame(rows_pair).to_csv(f"{OUT}/R7_E1_pairs.csv", index=False)
pd.DataFrame(rows_t).to_csv(f"{OUT}/R7_E1_tissue.csv", index=False)
