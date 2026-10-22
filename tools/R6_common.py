"""R6 common: bin->nucleus footprints, count aggregation (pure functions, tested by tools/R6_test.py)."""
import numpy as np, pandas as pd, scipy.sparse as sp, h5py
from scipy import ndimage

def footprint_labels(mask, rows, cols, um_per_px, d_exp_um, mode="exp", shift_um=20.0):
    """Label (0 = ambient) for each bin centre at integer pixel (rows, cols) of the ROI label image `mask`.
    mode 'nuc': label of the pixel; 'exp': pixel label, else nearest nucleus if its edge is within d_exp_um;
    'shift': as 'exp' on a mask rolled by shift_um in both axes (negative control)."""
    rows = np.asarray(rows, dtype=np.int64); cols = np.asarray(cols, dtype=np.int64)
    H, W = mask.shape
    inside = (rows >= 0) & (rows < H) & (cols >= 0) & (cols < W)
    out = np.zeros(len(rows), dtype=np.int64)
    if mode == "shift":
        s = int(round(shift_um / um_per_px)); mask = np.roll(np.roll(mask, s, axis=0), s, axis=1); mode = "exp"
    r, cc = rows[inside], cols[inside]
    lab = mask[r, cc].astype(np.int64)
    if mode == "exp":
        dist, (ir, ic) = ndimage.distance_transform_edt(mask == 0, return_indices=True)
        amb = lab == 0
        d_um = dist[r[amb], cc[amb]] * um_per_px
        near = mask[ir[r[amb], cc[amb]], ic[r[amb], cc[amb]]].astype(np.int64)
        near[d_um > d_exp_um] = 0
        lab[amb] = near
    out[inside] = lab
    return out

def aggregate_counts(X_bins, labels, keep_labels):
    """X_bins: sparse (genes x bins); labels: label per bin (0 = ambient); keep_labels: array of nucleus labels (rows of output).
    Returns (counts csr nuclei x genes, n_bins per nucleus, umi per nucleus)."""
    keep_labels = np.asarray(keep_labels)
    pos = pd.Series(np.arange(len(keep_labels)), index=keep_labels)
    sel = np.isin(labels, keep_labels)
    idx_n = pos.reindex(labels[sel]).to_numpy()
    idx_b = np.flatnonzero(sel)
    I = sp.csr_matrix((np.ones(len(idx_b)), (idx_n, idx_b)), shape=(len(keep_labels), X_bins.shape[1]))
    C = (I @ X_bins.T.tocsr()).tocsr()
    n_bins = np.bincount(idx_n, minlength=len(keep_labels))
    umi = np.asarray(C.sum(axis=1)).ravel()
    return C, n_bins, umi

def read_h5_csc(path):
    """10x filtered_feature_bc_matrix.h5 -> (csc genes x barcodes, barcodes array, gene names, gene ids)."""
    with h5py.File(path, "r") as h:
        m = h["matrix"]
        shape = tuple(int(v) for v in m["shape"][:])
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape)
        bc = m["barcodes"][:].astype(str)
        names = m["features"]["name"][:].astype(str); ids = m["features"]["id"][:].astype(str)
    return X, bc, names, ids
