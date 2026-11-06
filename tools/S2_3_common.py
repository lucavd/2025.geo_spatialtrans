"""tools/S2_3_common.py — S2.3: funzioni pure del banco e delle statistiche (testate in tools/S2.3_test.py).
Pre-registrazione: docs/steps/S2.3.md (cancelli 4994f45).

Campo di variabilita' fine F (design sezione 5, DEC-040): sulla sottogriglia di 0.5 um, Z = rumore bianco gaussiano convoluto con una
gaussiana di deviazione standard L (um), riportato a varianza 1 con la somma dei quadrati dei pesi discreti; F = quantile della gamma
di media 1 e varianza s2 in Phi(Z), calcolato per interpolazione lineare su una tabella di Z (passo 1e-3 in [-8.5, 8.5]); la stessa
tabella e' in R (R/05b3_place_and_move_transcripts.R). Statistiche per gruppo: psi (indice di dispersione di S2.0 M3, gruppi generici),
rho a distanza 1 e 2 quadratini (S2.0 PH-3, vettorizzata)."""
import numpy as np
from scipy import ndimage, stats as _stats

SEED = 20261112
SEED_BOOT = 20261113
B = 2000
SUB_UM = 0.5
ZTAB = np.round(np.arange(-8500, 8501) * 1e-3, 3)

def kernel_taps(L_um, step=SUB_UM):
    """Pesi discreti della gaussiana di deviazione standard L_um (in pixel L_um/step), troncata a 4 deviazioni standard."""
    sd = L_um / step
    h = int(np.ceil(4.0 * sd))
    x = np.arange(-h, h + 1, dtype=float)
    return np.exp(-0.5 * (x / sd) ** 2)

def z_corr_theory(L_um, lag_px, step=SUB_UM):
    """Correlazione teorica di Z fra due sottopixel a distanza lag_px lungo un asse (stessa riga): sum g_i g_{i+k} / sum g^2."""
    g = kernel_taps(L_um, step); k = int(lag_px)
    if k >= len(g): return 0.0
    return float(np.sum(g[:len(g) - k] * g[k:]) / np.sum(g * g))

def gamma_table(s2):
    """Trasformazione Z -> F sulla tabella ZTAB (gamma di media 1 e varianza s2)."""
    p = _stats.norm.cdf(ZTAB)
    return _stats.gamma.ppf(p, a=1.0 / s2, scale=s2)

def gaussian_z(shape, L_um, rng, step=SUB_UM):
    """Z gaussiano a varianza 1 con correlazione da convoluzione separabile (rumore bianco generato su un riquadro allargato)."""
    g = kernel_taps(L_um, step); h = (len(g) - 1) // 2
    W = rng.standard_normal((shape[0] + 2 * h, shape[1] + 2 * h))
    Z = ndimage.correlate1d(W, g, axis=0, mode="constant")
    Z = ndimage.correlate1d(Z, g, axis=1, mode="constant")
    Z = Z[h:h + shape[0], h:h + shape[1]] / float(np.sum(g * g))
    return Z

def gamma_field(shape, s2, L_um, rng, step=SUB_UM):
    """Campo F (media 1, varianza puntuale s2) sulla sottogriglia; s2 <= 0 -> F = 1 (il generatore non viene usato)."""
    if s2 <= 0:
        return np.ones(shape)
    Z = np.clip(gaussian_z(shape, L_um, rng, step), ZTAB[0], ZTAB[-1])
    return np.interp(Z, ZTAB, gamma_table(s2))

# ------------------------------------------------------------------ statistiche per gruppo
def _group_index(grp, min_bins, y):
    sel = grp >= 0
    g = grp[sel]
    uniq, inv, n = np.unique(g, return_inverse=True, return_counts=True)
    sy = np.bincount(inv, weights=y[sel])
    ok = (n >= min_bins) & (sy > 0)
    return sel, inv, n, sy, ok

def psi_groups(u, grp, c, min_bins=4):
    """Indice di dispersione (S2.0 M3, gruppi qualsiasi; grp = -1 fuori): per gruppo X2 = sum (y - m)^2 / m con
    m = media del gruppo * c / media di c nel gruppo; psi = sum X2 / sum (n - 1) sui gruppi con >= min_bins quadratini e >= 1 UMI.
    Ritorna (psi, m medio per quadratino, numero di gruppi)."""
    y = np.asarray(u, float); cc = np.asarray(c, float); grp = np.asarray(grp)
    sel, inv, n, sy, ok = _group_index(grp, min_bins, y)
    sc = np.bincount(inv, weights=cc[sel])
    okb = ok[inv]
    m = (sy / n)[inv] * cc[sel] / (sc / n)[inv]
    x2 = np.where(okb, (y[sel] - m) ** 2 / np.where(okb, m, 1.0), 0.0)
    X2 = np.bincount(inv, weights=x2)
    df = ((n - 1) * ok).sum()
    psi = float(X2[ok].sum() / df) if df > 0 else float("nan")
    return psi, (float(m[okb].mean()) if okb.any() else float("nan")), int(ok.sum())

def residuals(u, grp, c, min_bins=4):
    """Residui standardizzati (y - m)/sqrt(m) con m come in psi_groups; NaN fuori dai gruppi validi."""
    y = np.asarray(u, float); cc = np.asarray(c, float); grp = np.asarray(grp)
    sel, inv, n, sy, ok = _group_index(grp, min_bins, y)
    sc = np.bincount(inv, weights=cc[sel])
    okb = ok[inv]
    m = (sy / n)[inv] * cc[sel] / (sc / n)[inv]
    z = np.full(len(y), np.nan)
    z[np.flatnonzero(sel)[okb]] = (y[sel][okb] - m[okb]) / np.sqrt(m[okb])
    return z

def rho_groups(u, grp, c, rows, cols, min_bins=4):
    """Correlazione dei residui fra quadratini dello stesso gruppo a distanza 1 e 2 lungo le righe (h) e le colonne (v)
    (S2.0 PH-3). Ritorna dict h1, h2, v1, v2, rho1 = media(h1, v1), rho2 = media(h2, v2), e il numero di coppie."""
    z = residuals(u, grp, c, min_bins); grp = np.asarray(grp)
    rows = np.asarray(rows, np.int64); cols = np.asarray(cols, np.int64)
    r0, c0 = rows.min(), cols.min()
    G = np.full((rows.max() - r0 + 3, cols.max() - c0 + 3), -1, np.int64)
    G[rows - r0, cols - c0] = np.arange(len(rows))
    out = {}
    for nm, dr, dc in (("h1", 0, 1), ("h2", 0, 2), ("v1", 1, 0), ("v2", 2, 0)):
        j = G[rows - r0 + dr, cols - c0 + dc]
        i = np.arange(len(rows))
        ok = (j >= 0)
        i, j = i[ok], j[ok]
        ok = (grp[i] >= 0) & (grp[i] == grp[j]) & np.isfinite(z[i]) & np.isfinite(z[j])
        i, j = i[ok], j[ok]
        out[nm] = float(np.corrcoef(z[i], z[j])[0, 1]) if len(i) > 10 else float("nan")
        out[nm + "_n"] = int(len(i))
    out["rho1"] = float(np.nanmean([out["h1"], out["v1"]])); out["rho2"] = float(np.nanmean([out["h2"], out["v2"]]))
    return out

def group_stats(u, grp, c, rows, cols, min_bins=4):
    psi, mbar, ng = psi_groups(u, grp, c, min_bins)
    rh = rho_groups(u, grp, c, rows, cols, min_bins)
    return {"psi": psi, "mbar": mbar, "n_groups": ng, **rh}

# ------------------------------------------------------------------ gruppi sui dati (maschere allineate)
def aligned_px(pos_um, d_img, upp):
    """Pixel della maschera (immagine dell'area) letti al centro del quadratino spostato dello scarto della sezione:
    il quadratino in pos vede la sorgente in pos - d_img (convenzione di make_blur di S2.0)."""
    q = (np.asarray(pos_um) - np.asarray(d_img)[None, :]) / upp
    return np.floor(q[:, 0]).astype(np.int64), np.floor(q[:, 1]).astype(np.int64)

def nuclear_maps(mask, upp, d_exp_um=8.0):
    """Mappe in pixel nativi: etichetta del nucleo, distanza dentro il nucleo, distanza fuori (um), etichetta dell'impronta di 8 um
    (nucleo piu' vicino con bordo entro d_exp_um; pareggi: indice della trasformata di scipy, come R4/R6)."""
    d_in = ndimage.distance_transform_edt(mask > 0) * upp
    d_out, (ir, ic) = ndimage.distance_transform_edt(mask == 0, return_indices=True)
    d_out = d_out * upp
    near = mask[ir, ic]
    exp = np.where(mask > 0, mask, np.where(d_out <= d_exp_um, near, 0))
    return {"lab": mask, "d_in": d_in, "d_out": d_out, "exp": exp}

def lookup(arr, r, c, fill=0):
    H, W = arr.shape
    ins = (r >= 0) & (r < H) & (c >= 0) & (c < W)
    out = np.full(len(r), fill, dtype=arr.dtype)
    out[ins] = arr[r[ins], c[ins]]
    return out

def groups_from_maps(Mp, r, c, keep, tissue_at_bin, rows, cols, deep_um=1.0, ring_gap_um=2.0, acell_gap_um=2.0, d_exp_um=8.0, tile=8):
    """Gruppi del banco per i quadratini (pixel allineati r, c): nucleo (centro a >= deep_um dentro un nucleo tenuto), corona
    (impronta fuori dal nucleo, a >= ring_gap_um dalla maschera), spazio senza cellule (tessuto fuori dalle impronte e a >= acell_gap_um
    da esse, blocchi di tile x tile quadratini della griglia). -1 fuori."""
    keep = np.asarray(keep)
    lab = lookup(Mp["lab"], r, c).astype(np.int64); din = lookup(Mp["d_in"], r, c, 0.0); dout = lookup(Mp["d_out"], r, c, np.inf)
    exp = lookup(Mp["exp"], r, c).astype(np.int64)
    kin = np.isin(lab, keep); kex = np.isin(exp, keep)
    g_nuc = np.where(kin & (lab > 0) & (din >= deep_um), lab, -1)
    g_ring = np.where((lab == 0) & kex & (exp > 0) & (dout >= ring_gap_um), exp, -1)
    acell = tissue_at_bin & (exp == 0) & (dout >= d_exp_um + acell_gap_um)
    rows = np.asarray(rows, np.int64); cols = np.asarray(cols, np.int64)
    tid = (rows // tile) * 100000 + (cols // tile)
    g_acell = np.where(acell, tid, -1)
    return {"nucleus": g_nuc, "corona": g_ring, "acellular": g_acell}

def profile_ratio(d, umi):
    """S1.3c: UMI medio a 4-8 um dal nucleo piu' vicino / a 0-2 um (0 dentro il nucleo)."""
    a = umi[(d >= 4) & (d < 8)]; b = umi[(d >= 0) & (d < 2)]
    return float(a.mean() / b.mean()) if len(a) and len(b) and b.mean() > 0 else float("nan")
