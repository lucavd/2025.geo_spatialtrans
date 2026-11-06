"""S2.2: funzioni comuni del banco della gara sulla fonte dei profili (metriche per quadratino, blocchi di 8 um,
somiglianza fra blocchi vicini, misura di Daniele dentro e fra cellule). Le stesse funzioni si applicano ai quadratini veri
e a quelli simulati. Pre-registrazione: docs/steps/S2.2.md."""
import numpy as np, scipy.sparse as sp

KEY_STRIDE = 1_000_000          # chiave del blocco = riga_blocco * KEY_STRIDE + colonna_blocco
MIN_BINS_BLOCK = 12             # un blocco di 8 um entra se ha almeno 12 dei 16 quadratini di valutazione
MIN_BINS_CELL = 4               # misura di Daniele: cellule con almeno 4 quadratini di valutazione
D_BETWEEN_UM = 30.0             # misura di Daniele: la cellula di confronto dello stesso tipo sta a più di 30 um
# tipi dei preset dello Step 1 (R/04b0_tissue_presets.R, righe 42-43; uguale a tools/S2.2_real.py, controllato da R/testing/test_S2.2.R)
PRESET_TYPES = {"A1": ["epiteliale", "immunitario", "stromale e muscolo liscio"], "A2": ["epiteliale", "fibroblasti"],
                "A3": ["endoteliale", "epiteliale", "fibroblasti", "immunitario"],
                "A4": ["endoteliali", "linfociti B", "linfociti T", "mieloidi", "plasmacellule", "stromali"],
                "A5": ["astrociti", "microglia", "oligodendrociti", "vascolari", "neuroni"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
MERGE = {"epiteliale normale": "epiteliale", "epiteliale tumorale": "epiteliale"}


def inside_box(roi_row, x_um, y_um, other_row):
    """Nuclei di un'area (coordinate in um dall'origine dell'area) che cadono nel riquadro di un'altra area (pixel dell'immagine).
    Le aree di R2 dello stesso tessuto si sovrappongono fino al 39 % (revisione G-5 di S2.2)."""
    col = roi_row.c0 + np.asarray(x_um) / roi_row.um_per_px; row = roi_row.r0 + np.asarray(y_um) / roi_row.um_per_px
    return (col >= other_row.c0) & (col < other_row.c1) & (row >= other_row.r0) & (row < other_row.r1)


def _rowcos(A, B):
    """Coseno riga per riga fra due csr della stessa forma."""
    num = np.asarray(A.multiply(B).sum(axis=1)).ravel()
    na = np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).ravel())
    nb = np.sqrt(np.asarray(B.multiply(B).sum(axis=1)).ravel())
    den = na * nb
    out = np.full(len(num), np.nan)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def blocks(Xb, rows, cols):
    """Somma 4 x 4 dei quadratini di valutazione. Ritorna (B csr blocchi x geni dei blocchi ammessi, chiavi ordinate)."""
    key = (np.asarray(rows) // 4).astype(np.int64) * KEY_STRIDE + (np.asarray(cols) // 4).astype(np.int64)
    uk, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    A = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(uk), len(inv)))
    B = (A @ Xb).tocsr()
    ok = cnt >= MIN_BINS_BLOCK
    return B[ok], uk[ok]


def adjacent_pairs(keys):
    """Coppie (i, j) di blocchi ammessi adiacenti a destra o in basso; keys ordinate."""
    I, J = [], []
    for step in (1, KEY_STRIDE):
        nb = keys + step
        j = np.searchsorted(keys, nb)
        j = np.minimum(j, len(keys) - 1)
        ok = keys[j] == nb
        I.append(np.flatnonzero(ok)); J.append(j[ok])
    return np.concatenate(I), np.concatenate(J)


def cell_halves(lab_bins, Xb, rng):
    """Divide a caso i quadratini di ogni cellula in due metà (seme dato). Ritorna (etichette, H1, H2) csr cellule x geni."""
    lab_bins = np.asarray(lab_bins)
    u, inv, cnt = np.unique(lab_bins, return_inverse=True, return_counts=True)
    keep = cnt >= MIN_BINS_CELL
    order = np.lexsort((rng.random(len(lab_bins)), inv))      # dentro ogni cellula, ordine casuale
    rank = np.empty(len(lab_bins), np.int64)
    starts = np.r_[0, np.cumsum(cnt)[:-1]]
    rank[order] = np.arange(len(lab_bins)) - np.repeat(starts, cnt)
    half = (rank >= (cnt[inv] // 2)).astype(np.int64)            # prima metà 0, seconda 1 (la seconda ha il quadratino dispari)
    sel = keep[inv]
    newidx = np.cumsum(keep) - 1
    rows_ = newidx[inv[sel]]
    n = int(keep.sum())
    H = []
    for h in (0, 1):
        m = sel & (half == h)
        A = sp.csr_matrix((np.ones(int(m.sum())), (newidx[inv[m]], np.flatnonzero(m))), shape=(n, len(lab_bins)))
        H.append((A @ Xb).tocsr())
    return u[keep], H[0], H[1]


def daniele(lab_bins, Xb, cell_type, cell_xy, rng, tries=10):
    """Misura di Daniele (BL-204): coseno fra le due metà della stessa cellula, contro il coseno fra la metà 1 di una cellula e la
    metà 2 di un'altra cellula dello stesso tipo a più di D_BETWEEN_UM. cell_type, cell_xy: dict etichetta -> tipo, (x, y)."""
    labs, H1, H2 = cell_halves(lab_bins, Xb, rng)
    if len(labs) < 2:
        return dict(m4_within=np.nan, m4_between=np.nan, m4_ratio=np.nan, n_cells_m4=len(labs))
    within = _rowcos(H1, H2)
    ty = np.array([cell_type[l] for l in labs]); xy = np.array([cell_xy[l] for l in labs], float)
    partner = np.full(len(labs), -1)
    for t in np.unique(ty):
        idx = np.flatnonzero(ty == t)
        if len(idx) < 2:
            continue
        todo = idx.copy()
        for _ in range(tries):
            cand = rng.choice(idx, size=len(todo))
            d = np.hypot(*(xy[todo] - xy[cand]).T)
            ok = (cand != todo) & (d > D_BETWEEN_UM)
            partner[todo[ok]] = cand[ok]
            todo = todo[~ok]
            if len(todo) == 0:
                break
    has = partner >= 0
    between = np.full(len(labs), np.nan)
    if has.any():
        between[has] = _rowcos(H1[np.flatnonzero(has)], H2[partner[has]])
    w = np.nanmean(within); b = np.nanmean(between)
    return dict(m4_within=float(w), m4_between=float(b), m4_ratio=float(w / b) if b > 0 else np.nan, n_cells_m4=int(len(labs)))


def bin_metrics(Xb, rows, cols, lab_bins, cell_type, cell_xy, seed):
    """Metriche della gara su una matrice csr quadratini x geni (solo quadratini di valutazione, geni dell'universo comune).
    m1 geni per quadratino di 2 um; m2 geni per blocco di 8 um; m3 coseno medio fra blocchi di 8 um adiacenti;
    m4 misura di Daniele. Lo stesso seme per il vero e per il simulato dello stesso ROI."""
    Xb = sp.csr_matrix(Xb)
    nz = (Xb > 0).astype(np.int8)
    m1 = float(np.asarray(nz.sum(axis=1)).mean())
    B, keys = blocks(Xb, rows, cols)
    m2 = float(np.asarray((B > 0).sum(axis=1)).mean()) if B.shape[0] else np.nan
    I, J = adjacent_pairs(keys)
    m3 = float(np.nanmean(_rowcos(B[I], B[J]))) if len(I) else np.nan
    G = Xb.shape[1]
    out = dict(m1=m1, m2=m2, m3=m3, sparsity2=1 - m1 / G, sparsity8=1 - m2 / G, n_bins=int(Xb.shape[0]),
               n_blocks=int(B.shape[0]), n_pairs=int(len(I)), umi_per_bin=float(Xb.sum() / Xb.shape[0]),
               umi_per_block=float(B.sum() / B.shape[0]) if B.shape[0] else np.nan, n_genes=int(G))
    out.update(daniele(lab_bins, Xb, cell_type, cell_xy, np.random.default_rng(seed)))
    return out
