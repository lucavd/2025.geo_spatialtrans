
import sys, json, numpy as np, pyarrow.parquet as pq
sys.path.insert(0, "tools/sonda"); import S1_registration as m
for tag in ["X1", "X2"]:
    P = m.paths(tag); M, H, A = np.load(f"{m.OUT}/{tag}/composite_parts.npy")
    tp8 = pq.read_table(f"{P['b8']}/spatial/tissue_positions.parquet").to_pandas()
    cells = pq.read_table(P["cells"], columns=["x_centroid", "y_centroid"]).to_pandas().values
    nc = int(tp8.array_col.max()) + 1; nr = int(tp8.array_row.max()) + 1
    tis = np.zeros((nr, nc), bool); t = tp8[tp8.in_tissue == 1]; tis[t.array_row.values, t.array_col.values] = True
    print(tag, "cells", len(cells), "grid8", nr, nc, "in_tissue8", int(tis.sum()))
    for v in ["direct|px0.2125", "direct|um", "inverse|px0.2125", "inverse|um"]:
        a = m.apply_h(m.composite(M, H, A, v), cells); k = np.floor(np.round(a) / 4).astype(int)
        ins = (k[:, 0] >= 0) & (k[:, 0] < nc) & (k[:, 1] >= 0) & (k[:, 1] < nr)
        on = np.zeros(len(k), bool); on[ins] = tis[k[ins, 1], k[ins, 0]]
        occ = np.zeros((nr, nc), bool); occ[k[ins, 1], k[ins, 0]] = True
        from scipy.ndimage import binary_dilation
        rev = (binary_dilation(occ, np.ones((3, 3), bool)) & tis).sum() / tis.sum()
        print(f"  {v:18s} inside_grid={ins.mean():.3f} on_tissue={on.mean():.3f} on_tissue|inside={on[ins].mean() if ins.any() else float('nan'):.3f} visium_tissue_with_xenium_cell={rev:.3f}")
