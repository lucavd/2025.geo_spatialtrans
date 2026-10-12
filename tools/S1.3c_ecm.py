"""tools/S1.3c_ecm.py — S1.3c (BL-077): spazio extracellulare intra-regione dai conteggi dei bin di 2 um.

Pre-registrazione: results/S1.3c/S1.3c_preregistration.md (f341f71); definizioni operative: results/S1.3c/S1.3c_dev_notes.md.
Uso (root del repo):  .venv/bin/python tools/S1.3c_ecm.py <stage> [outdir]
  stage = umi | roi | ctrl | all ;  outdir default /mnt/micron/geo_spatialtrans/S1.3c  (C-6: rerun in un'altra cartella)
Stadi rilanciabili: saltano i file gia' scritti.
"""
import sys, json, hashlib
from pathlib import Path
import numpy as np, pandas as pd, h5py
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage.io import imread

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/mnt/micron/geo_spatialtrans/data_real")
R2M = Path("/mnt/micron/geo_spatialtrans/R2/masks")
R3R = Path("/mnt/micron/geo_spatialtrans/R3/real")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("/mnt/micron/geo_spatialtrans/S1.3c")
KS = (3, 5)                       # finestre 6 e 10 um
BG_MIN_UM = 50.0; BG_Q = 95.0; NPERM = 20
T_MED = {"A1": 61.3, "A2": 94.0, "A3": 167.7, "A4": 32.1, "A5": 819.1, "A6": 742.8}   # B-055..B-085 (um2)
PRIMARY = {"A1": ("cellpose_rgb", "cellpose_rgb"), "A2": ("cellpose_rgb", "cellpose_rgb"),
           "A3": ("cellpose_rgb", "cellpose_rgb"), "A4": ("spaceranger", "stardist_he"),
           "A5": ("cellpose_rgb", "cellpose_rgb"), "A6": ("spaceranger", "stardist_he")}   # (territori, maschera)
SENS = {"A1": "stardist_he", "A2": "stardist_he", "A3": "stardist_he", "A4": "cellpose_rgb",
        "A5": "stardist_he", "A6": "cellpose_rgb"}
E_THR, BA_A, BA_B, LARGE_FRAC = 0.05, 0.75, 0.65, 0.50

# ------------------------------------------------------------------------------------- funzioni pure (testate)
def umi_totals(indptr, data, indices=None, keep_rows=None):
    """Totale per colonna di una matrice CSC (colonne = barcode). keep_rows: maschera booleana sulle righe (feature)."""
    d = data.astype(np.int64)
    if keep_rows is not None and not keep_rows.all():
        d = np.where(keep_rows[indices], d, 0)
    cs = np.concatenate([[0], np.cumsum(d)])
    return cs[indptr[1:]] - cs[indptr[:-1]]

def window_ids(arow, acol, k):
    """Indici di finestra (tassellatura globale k x k sulla griglia dei bin)."""
    return arow // k, acol // k

def window_table(wr, wc, values, ok, k):
    """Somme per finestra e completezza: completa se contiene k^2 bin con ok=True. Ritorna DataFrame (wr, wc, n_ok, total)."""
    df = pd.DataFrame({"wr": wr, "wc": wc, "v": values, "ok": ok.astype(np.int64)})
    g = df.groupby(["wr", "wc"], sort=True).agg(n_ok=("ok", "sum"), n=("ok", "size"), total=("v", "sum")).reset_index()
    g["complete"] = (g.n_ok == k * k) & (g.n == k * k)
    return g

def permute_within(groups, values, rng):
    """Permuta values fra gli elementi con lo stesso gruppo (conserva le somme per gruppo)."""
    groups = np.asarray(groups); values = np.asarray(values)
    stable = np.argsort(groups, kind="stable")
    rnd = np.lexsort((rng.random(len(groups)), groups))
    out = np.empty_like(values)
    out[rnd] = values[stable]
    return out

def fit_model_a(dmin, empty_obs, step=0.25):
    """r* tale che mean(dmin > r*) ~ phi_obs (griglia step; a parita' il piu' piccolo); BA rispetto a empty_obs."""
    phi = empty_obs.mean()
    grid = np.arange(0, np.nanmax(dmin) + step, step)
    fr = np.array([(dmin > r).mean() for r in grid])
    r = grid[int(np.argmin(np.abs(fr - phi)))]
    pred = dmin > r
    return float(r), balanced_accuracy(pred, empty_obs)

def balanced_accuracy(pred, obs):
    pred = np.asarray(pred, bool); obs = np.asarray(obs, bool)
    tpr = (pred & obs).sum() / max(obs.sum(), 1); tnr = (~pred & ~obs).sum() / max((~obs).sum(), 1)
    if obs.all() or (~obs).all(): return float("nan")
    return float((tpr + tnr) / 2)

def components_large_frac(wr, wc, empty, s_um, t_med):
    """Componenti 4-connesse delle finestre vuote. Ritorna (area vuota, area in componenti > t_med, n componenti)."""
    if empty.sum() == 0: return 0.0, 0.0, 0
    r0, c0 = wr.min(), wc.min()
    G = np.zeros((wr.max() - r0 + 1, wc.max() - c0 + 1), bool)
    G[wr[empty] - r0, wc[empty] - c0] = True
    lab, n = ndimage.label(G, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    sizes = np.bincount(lab.ravel())[1:] * s_um ** 2
    return float(sizes.sum()), float(sizes[sizes > t_med].sum()), int(n)

def verdict(E, ba, large_frac):
    if not np.isfinite(E) or E < E_THR: return "nessuno"
    if np.isfinite(ba) and ba >= BA_A: return "A"
    if np.isfinite(ba) and ba < BA_B and large_frac >= LARGE_FRAC: return "B"
    return "indeterminato"

def profile_ratio(d, umi):
    a = umi[(d >= 4) & (d < 8)]; b = umi[(d >= 0) & (d < 2)]
    return float(a.mean() / b.mean()) if len(a) and len(b) and b.mean() > 0 else float("nan")

# ------------------------------------------------------------------------------------- I/O
def rois_table():
    r = pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")
    r["side_um"] = r.side_px * r.um_per_px
    return r

def sq2(ds): return DATA / ds / "binned_outputs_x/binned_outputs/square_002um"

def load_dataset_bins(ds, log):
    """Tabella dei bin del dataset: barcode, in_tissue, array_row/col, pixel, UMI grezzi (join per barcode)."""
    pos = pd.read_parquet(sq2(ds) / "spatial/tissue_positions.parquet")
    with h5py.File(sq2(ds) / "raw_feature_bc_matrix.h5") as f:
        m = f["matrix"]; ft = m["features"]["feature_type"][:].astype(str)
        keep = ft == "Gene Expression"
        log[f"{ds}_feature_types"] = {t: int((ft == t).sum()) for t in np.unique(ft)}
        tot = umi_totals(m["indptr"][:], m["data"][:], m["indices"][:] if not keep.all() else None, keep)
        bc = m["barcodes"][:].astype(str)
    u = pd.DataFrame({"barcode": bc, "umi": tot})
    b = pos.merge(u, on="barcode", how="left", validate="one_to_one")
    log[f"{ds}_n_bins"] = int(len(b)); log[f"{ds}_n_umi_missing"] = int(b.umi.isna().sum())
    b["umi"] = b.umi.fillna(0).astype(np.int64)
    return b

def check_filtered(ds, b):
    """C-3: totali della matrice filtrata = grezza sui bin in tessuto."""
    with h5py.File(sq2(ds) / "filtered_feature_bc_matrix.h5") as f:
        m = f["matrix"]; ft = m["features"]["feature_type"][:].astype(str); keep = ft == "Gene Expression"
        tot = umi_totals(m["indptr"][:], m["data"][:], m["indices"][:] if not keep.all() else None, keep)
        bc = m["barcodes"][:].astype(str)
    fdf = pd.DataFrame({"barcode": bc, "umi_f": tot})
    j = b[b.in_tissue == 1][["barcode", "umi"]].merge(fdf, on="barcode", how="outer", indicator=True)
    return {"n_in_tissue": int((b.in_tissue == 1).sum()), "n_filtered": int(len(fdf)),
            "n_only_one_side": int((j._merge != "both").sum()),
            "n_diff": int((j.umi != j.umi_f).sum()), "identical": bool((j._merge == "both").all() and (j.umi == j.umi_f).all())}

def background_thresholds(b):
    """Soglia di vuoto per k: q95 dei totali delle finestre di sfondo (tutti bin in_tissue=0 a >= 50 um dal tessuto)."""
    R, C = b.array_row.max() + 1, b.array_col.max() + 1
    T = np.zeros((R, C), bool); T[b.array_row, b.array_col] = b.in_tissue.values == 1
    pitch = 2.0
    dist = ndimage.distance_transform_edt(~T) * pitch
    okb = (b.in_tissue.values == 0) & (dist[b.array_row, b.array_col] >= BG_MIN_UM)
    out = {}
    for k in KS:
        wr, wc = window_ids(b.array_row.values, b.array_col.values, k)
        g = window_table(wr, wc, b.umi.values, okb, k)
        g = g[g.complete]
        out[k] = {"thr": float(np.percentile(g.total, BG_Q)) if len(g) else float("nan"), "n_windows": int(len(g)),
                  "mean": float(g.total.mean()) if len(g) else float("nan"), "frac_zero": float((g.total == 0).mean()) if len(g) else float("nan")}
    return out

def nuc_mask(A, roi, method):
    lab = np.load(R2M / f"{A}_{roi}_{method}_native.npz")["labels"]
    nuc = pd.read_parquet(ROOT / "results/R2/R2_nuclei_all.parquet", columns=["archetype", "roi_id", "method", "scale", "keep", "label"])
    keep = nuc[(nuc.archetype == A) & (nuc.roi_id == roi) & (nuc.method == method) & (nuc.scale == 1) & nuc.keep].label.values
    lut = np.zeros(int(lab.max()) + 1, bool); lut[keep[keep <= lab.max()].astype(np.int64)] = True
    return lut[lab]

def generators(A, roi, method, valid, px4):
    f = R3R / f"{A}_{roi}_{method}_gen.parquet"
    if f.exists():
        g = pd.read_parquet(f, columns=["x", "y"]); return g.x.values, g.y.values, "R3"
    nuc = pd.read_parquet(ROOT / "results/R2/R2_nuclei_all.parquet", columns=["archetype", "roi_id", "method", "scale", "keep", "x_um", "y_um"])
    d = nuc[(nuc.archetype == A) & (nuc.roi_id == roi) & (nuc.method == method) & (nuc.scale == 1) & nuc.keep]
    ci = np.floor(d.x_um.values / px4).astype(int); ri = np.floor(d.y_um.values / px4).astype(int)
    ok = (ci >= 0) & (ri >= 0) & (ci < valid.shape[1]) & (ri < valid.shape[0])
    ok[ok] = valid[ri[ok], ci[ok]]
    return d.x_um.values[ok], d.y_um.values[ok], "R2_keep_in_valid"

def roi_bins(b, r):
    s = b[(b.pxl_col_in_fullres >= r.c0) & (b.pxl_col_in_fullres < r.c0 + r.side_px) &
          (b.pxl_row_in_fullres >= r.r0) & (b.pxl_row_in_fullres < r.r0 + r.side_px)].copy()
    s["ci"] = np.floor(s.pxl_col_in_fullres - r.c0).astype(int).clip(0, r.side_px - 1)
    s["ri"] = np.floor(s.pxl_row_in_fullres - r.r0).astype(int).clip(0, r.side_px - 1)
    s["x"] = (s.pxl_col_in_fullres - r.c0) * r.um_per_px; s["y"] = (s.pxl_row_in_fullres - r.r0) * r.um_per_px
    return s

def nuc_fraction_shift(nm, ri, ci, dx_px, dy_px, h):
    """Frazione nucleare nel quadrato di lato 2h centrato su (ri+dy, ci+dx); NaN se il quadrato esce dal ROI."""
    S = np.pad(np.cumsum(np.cumsum(nm.astype(np.int32), 0), 1), ((1, 0), (1, 0)))
    r = np.round(ri + dy_px).astype(int); c = np.round(ci + dx_px).astype(int)
    r0, r1, c0, c1 = r - h, r + h, c - h, c + h
    ok = (r0 >= 0) & (c0 >= 0) & (r1 <= nm.shape[0]) & (c1 <= nm.shape[1])
    out = np.full(len(r), np.nan)
    out[ok] = (S[r1[ok], c1[ok]] - S[r0[ok], c1[ok]] - S[r1[ok], c0[ok]] + S[r0[ok], c0[ok]]) / (2 * h) ** 2
    return out

# ------------------------------------------------------------------------------------- analisi di un ROI
def analyse(s, valid_b, dnuc, terr, thr, k, A, rng_seed=1):
    """Tutte le grandezze della regola per una configurazione (s = bin del ROI, valid_b = bin in tessuto)."""
    umi = s.umi.values.astype(np.int64)
    wr, wc = window_ids(s.array_row.values, s.array_col.values, k)
    g = window_table(wr, wc, umi, valid_b, k); gt = g[g.complete]
    empty = (gt.total.values <= thr)
    phi_obs = float(empty.mean())
    # nullo
    idx = np.where(valid_b)[0]; phis = []
    for p in range(1, NPERM + 1):
        u2 = umi.copy(); u2[idx] = permute_within(terr[idx], umi[idx], np.random.default_rng(p))
        g2 = window_table(wr, wc, u2, valid_b, k); g2 = g2[g2.complete]
        phis.append(float((g2.total.values <= thr).mean()))
    phi_null = float(np.mean(phis))
    # modello A
    dmin = pd.DataFrame({"wr": wr, "wc": wc, "d": dnuc}).groupby(["wr", "wc"], sort=True).d.min().reset_index()
    dmin = gt[["wr", "wc"]].merge(dmin, on=["wr", "wc"], how="left").d.values
    r_star, ba = fit_model_a(dmin, empty)
    a_empty, a_large, ncomp = components_large_frac(gt.wr.values, gt.wc.values, empty, 2.0 * k, T_MED[A])
    E = phi_obs - phi_null
    return {"k": k, "s_um": 2 * k, "thr": thr, "n_tissue_windows": int(len(gt)), "n_empty": int(empty.sum()),
            "phi_obs": phi_obs, "phi_null": phi_null, "phi_null_min": float(min(phis)), "phi_null_max": float(max(phis)),
            "E": E, "r_star": r_star, "BA": ba, "empty_area": a_empty, "large_area": a_large, "n_comp": ncomp,
            "large_frac": a_large / a_empty if a_empty > 0 else float("nan"),
            "R_profile": profile_ratio(dnuc[valid_b], umi[valid_b].astype(float)),
            "verdict_roi": verdict(E, ba, a_large / a_empty if a_empty > 0 else float("nan"))}, gt.assign(empty=empty, dmin=dmin)

def load_roi_context(r, b, mask_method):
    A, roi = r.archetype, r.roi_id
    s = roi_bins(b, r)
    valid4 = imread(ROOT / f"results/R2/tissue_masks/{A}_{roi}_valid_ds4.png") > 127
    if valid4.ndim == 3: valid4 = valid4[..., 0]
    px4 = r.side_um / valid4.shape[1]
    vb = valid4[np.clip((s.y / px4).astype(int), 0, valid4.shape[0] - 1), np.clip((s.x / px4).astype(int), 0, valid4.shape[1] - 1)]
    nm = nuc_mask(A, roi, mask_method)
    dmap = ndimage.distance_transform_edt(~nm) * r.um_per_px
    dnuc = dmap[s.ri.values, s.ci.values]
    return s, valid4, px4, np.asarray(vb), nm, dmap, dnuc

def roi_stage_one(args):
    r, ds_npz, thr_all, variant = args
    A, roi = r.archetype, r.roi_id
    terr_m, mask_m = PRIMARY[A] if variant == "primary" else (SENS[A], SENS[A])
    f_out = OUT / "roi" / f"{A}_{roi}_{variant}.json"
    if f_out.exists(): return str(f_out)
    b = pd.read_parquet(ds_npz)
    s, valid4, px4, vb, nm, dmap, dnuc = load_roi_context(r, b, mask_m)
    gx, gy, gsrc = generators(A, roi, terr_m, valid4, px4)
    terr = cKDTree(np.c_[gx, gy]).query(np.c_[s.x.values, s.y.values])[1]
    res = {"archetype": A, "roi_id": roi, "variant": variant, "territories": terr_m, "mask": mask_m, "gen_source": gsrc,
           "n_gen": int(len(gx)), "n_bins_roi": int(len(s)), "n_bins_valid": int(vb.sum()),
           "n_bins_expected_C1": float(r.cov_bins * (r.side_um / 2.0) ** 2)}
    # C-4 (una volta, sui bin validi)
    idx = np.where(vb)[0]; u = s.umi.values
    up = permute_within(terr[idx], u[idx], np.random.default_rng(1))
    res["C4_sums_preserved"] = bool((np.bincount(terr[idx], up, minlength=len(gx)) == np.bincount(terr[idx], u[idx], minlength=len(gx))).all())
    um = np.random.default_rng(1).permutation(u[idx])
    res["C4_mutant_violates"] = bool(not (np.bincount(terr[idx], um, minlength=len(gx)) == np.bincount(terr[idx], u[idx], minlength=len(gx))).all())
    res["by_k"] = {}
    for k in KS:
        out, gt = analyse(s, vb, dnuc, terr, thr_all[str(k)]["thr"], k, A)
        res["by_k"][str(k)] = out
        if variant == "primary":
            gt.to_parquet(OUT / "roi" / f"{A}_{roi}_{variant}_k{k}_windows.parquet", index=False)
    if variant == "primary":
        res.update(roi_extras(r, s, vb, nm, dmap, dnuc, terr, thr_all, A))
        prof = pd.DataFrame({"d": dnuc[vb], "umi": s.umi.values[vb]})
        prof["d_bin"] = np.floor(prof.d).clip(upper=20)
        prof.groupby("d_bin").umi.agg(["mean", "size"]).reset_index().to_parquet(OUT / "roi" / f"{A}_{roi}_profile.parquet", index=False)
    f_out.write_text(json.dumps(res, indent=1, default=float))
    return str(f_out)

def roi_extras(r, s, vb, nm, dmap, dnuc, terr, thr_all, A):
    """C-2 registrazione, CP-1 placebo, CP-2 vuoto vero, descrittivo H&E."""
    out = {}
    upp = r.um_per_px; h = int(round(1.0 / upp)); u = s.umi.values.astype(float)
    best = None; surf = []
    for dx in range(-10, 11, 2):
        for dy in range(-10, 11, 2):
            nf = nuc_fraction_shift(nm, s.ri.values, s.ci.values, dx / upp, dy / upp, h)
            ok = vb & np.isfinite(nf)
            c = float(np.corrcoef(nf[ok], u[ok])[0, 1]) if ok.sum() > 10 and nf[ok].std() > 0 else float("nan")
            surf.append((dx, dy, c))
            if np.isfinite(c) and (best is None or abs(c) > abs(best[2])): best = (dx, dy, c)
    out["C2_best_dx_um"], out["C2_best_dy_um"], out["C2_best_corr"] = best
    out["C2_corr_at_0"] = [c for dx, dy, c in surf if dx == 0 and dy == 0][0]
    out["C2_max_abs_shift_bins"] = int(max(abs(best[0]), abs(best[1])) // 2)
    # CP-1 placebo
    rng = np.random.default_rng(20261013 + int(r.name))
    ang = rng.uniform(0, 2 * np.pi); mod = rng.uniform(30, 60)
    dy_px, dx_px = int(round(mod * np.sin(ang) / upp)), int(round(mod * np.cos(ang) / upp))
    nm_p = np.roll(nm, (dy_px, dx_px), axis=(0, 1))
    d_p = (ndimage.distance_transform_edt(~nm_p) * upp)[s.ri.values, s.ci.values]
    changed = float((d_p[vb] != dnuc[vb]).mean())
    out["CP1_shift_um"] = float(np.hypot(dy_px, dx_px) * upp); out["CP1_frac_changed"] = changed
    assert changed > 0.5, "CP-1: il placebo non ha cambiato d_nuc"
    out["CP1"] = {}
    for k in KS:
        wr, wc = window_ids(s.array_row.values, s.array_col.values, k)
        g = window_table(wr, wc, s.umi.values, vb, k); gt = g[g.complete]
        empty = gt.total.values <= thr_all[str(k)]["thr"]
        dm = pd.DataFrame({"wr": wr, "wc": wc, "d": d_p}).groupby(["wr", "wc"], sort=True).d.min().reset_index()
        dm = gt[["wr", "wc"]].merge(dm, on=["wr", "wc"], how="left").d.values
        rs, ba = fit_model_a(dm, empty)
        out["CP1"][str(k)] = {"r_star": rs, "BA": ba}
    out["CP1_R_profile"] = profile_ratio(d_p[vb], u[vb])
    # CP-2: finestre fuori tessuto R2 (non bolle) ma in_tissue = 1
    z = np.load(R2M / f"{r.archetype}_{r.roi_id}_masks.npz"); t4, b4 = z["tissue_ds4"], z["bubble_ds4"]
    px4 = r.side_um / t4.shape[1]
    yy = np.clip((s.y / px4).astype(int), 0, t4.shape[0] - 1); xx = np.clip((s.x / px4).astype(int), 0, t4.shape[1] - 1)
    okv = (~t4[yy, xx]) & (~b4[yy, xx]) & (s.in_tissue.values == 1)
    out["CP2"] = {}
    for k in KS:
        wr, wc = window_ids(s.array_row.values, s.array_col.values, k)
        g = window_table(wr, wc, s.umi.values, np.asarray(okv), k); g = g[g.complete]
        out["CP2"][str(k)] = {"n_windows": int(len(g)), "n_empty": int((g.total <= thr_all[str(k)]["thr"]).sum())}
    # descrittivo H&E sulla griglia ds4 (valida)
    valid4 = imread(ROOT / f"results/R2/tissue_masks/{r.archetype}_{r.roi_id}_valid_ds4.png") > 127
    if valid4.ndim == 3: valid4 = valid4[..., 0]
    f = r.side_px / valid4.shape[0]
    ii, jj = np.nonzero(valid4)
    dv = dmap[np.clip(((ii + 0.5) * f).astype(int), 0, r.side_px - 1), np.clip(((jj + 0.5) * f).astype(int), 0, r.side_px - 1)]
    out["HE_frac_beyond"] = {str(d): float((dv > d).mean()) for d in (2, 4, 6, 8)}
    return out

# ------------------------------------------------------------------------------------- controlli sintetici C-5
def ctrl_one(args):
    r, ds_npz, thr_all = args
    A, roi = r.archetype, r.roi_id
    f_out = OUT / "ctrl" / f"{A}_{roi}_C5.json"
    if f_out.exists(): return str(f_out)
    terr_m, mask_m = PRIMARY[A]; k = 3; thr = thr_all["3"]["thr"]
    b = pd.read_parquet(ds_npz)
    s, valid4, px4, vb, nm, dmap, dnuc = load_roi_context(r, b, mask_m)
    gx, gy, _ = generators(A, roi, terr_m, valid4, px4)
    terr = cKDTree(np.c_[gx, gy]).query(np.c_[s.x.values, s.y.values])[1]
    base, gt0 = analyse(s, vb, dnuc, terr, thr, k, A)
    res = {"archetype": A, "roi_id": roi, "base": base}
    # (b) tipo A
    r90 = float(np.percentile(dnuc[vb], 90))
    zA = vb & (dnuc > r90)
    sA = s.copy(); sA.loc[zA, "umi"] = 0
    res["typeA"], _ = analyse(sA, vb, dnuc, terr, thr, k, A); res["typeA"]["r90"] = r90; res["typeA"]["frac_zeroed"] = float(zA.sum() / vb.sum())
    # (c) tipo B: dischi
    rng = np.random.default_rng(20261012); R = np.sqrt(3 * T_MED[A] / np.pi)
    tree = cKDTree(np.c_[s.x.values, s.y.values]); inside = np.zeros(len(s), bool); n_disc = 0
    while (inside & vb).sum() < 0.10 * vb.sum():
        cx, cy = rng.uniform(0, r.side_um, 2)
        if not valid4[min(int(cy / px4), valid4.shape[0] - 1), min(int(cx / px4), valid4.shape[1] - 1)]: continue
        inside[tree.query_ball_point([cx, cy], R)] = True; n_disc += 1
    zB = vb & inside
    sB = s.copy(); sB.loc[zB, "umi"] = 0
    res["typeB"], gtB = analyse(sB, vb, dnuc, terr, thr, k, A)
    res["typeB"].update({"disc_radius_um": float(R), "n_discs": n_disc, "frac_zeroed": float(zB.sum() / vb.sum())})
    # (a) bilancio esatto su (c)
    wr, wc = window_ids(s.array_row.values, s.array_col.values, k)
    zz = pd.DataFrame({"wr": wr, "wc": wc, "z": zB.astype(int)}).groupby(["wr", "wc"], sort=True).z.sum().reset_index()
    m = gt0[["wr", "wc", "empty"]].merge(zz, on=["wr", "wc"], how="left").merge(gtB[["wr", "wc", "empty"]].rename(columns={"empty": "empty_B"}), on=["wr", "wc"])
    full = int(((m.z == k * k) & ~m.empty).sum())
    partial = int(((m.z > 0) & (m.z < k * k) & ~m.empty & m.empty_B).sum())
    delta = int(m.empty_B.sum() - m.empty.sum())
    res["C5a"] = {"delta_empty": delta, "full_geom": full, "partial_fallen": partial, "exact": delta == full + partial,
                  "geom_frac": full / len(m), "other_changes": int((m.empty & ~m.empty_B).sum())}
    f_out.write_text(json.dumps(res, indent=1, default=float))
    return str(f_out)

# ------------------------------------------------------------------------------------- stadi
def stage_umi():
    (OUT / "bins").mkdir(parents=True, exist_ok=True)
    rois = rois_table(); log = {}
    for ds, g in rois.groupby("dataset", sort=False):
        f_thr = OUT / "bins" / f"{ds}_bg.json"
        if f_thr.exists() and all((OUT / "bins" / f"{ds}_roi_bins.parquet").exists() for _ in [0]): continue
        b = load_dataset_bins(ds, log)
        c3 = check_filtered(ds, b); bg = background_thresholds(b)
        # bin di tutti i ROI del dataset (con margine nullo) -> un parquet per dataset
        keep = np.zeros(len(b), bool)
        for _, r in g.iterrows():
            keep |= ((b.pxl_col_in_fullres >= r.c0) & (b.pxl_col_in_fullres < r.c0 + r.side_px) &
                     (b.pxl_row_in_fullres >= r.r0) & (b.pxl_row_in_fullres < r.r0 + r.side_px)).values
        b[keep].to_parquet(OUT / "bins" / f"{ds}_roi_bins.parquet", index=False)
        json.dump({"C3": c3, "bg": {str(k): v for k, v in bg.items()}, "log": {k: v for k, v in log.items() if k.startswith(ds)},
                   "umi_total_dataset": int(b.umi.sum()), "umi_in_tissue": int(b.umi[b.in_tissue == 1].sum())},
                  open(f_thr, "w"), indent=1, default=float)
        print(ds, "C3", c3["identical"], "bg", {k: (v["thr"], v["n_windows"]) for k, v in bg.items()}, flush=True)

def stage_roi(ncores=15):
    from multiprocessing import Pool
    (OUT / "roi").mkdir(parents=True, exist_ok=True)
    rois = rois_table(); tasks = []
    for i, r in rois.iterrows():
        bg = json.load(open(OUT / "bins" / f"{r.dataset}_bg.json"))["bg"]
        for v in ("primary", "sens"):
            tasks.append((r, OUT / "bins" / f"{r.dataset}_roi_bins.parquet", bg, v))
    with Pool(ncores) as p:
        for f in p.imap_unordered(roi_stage_one, tasks): print("ok", f, flush=True)

def stage_ctrl(ncores=6):
    from multiprocessing import Pool
    (OUT / "ctrl").mkdir(parents=True, exist_ok=True)
    rois = rois_table(); first = rois.groupby("archetype", sort=False).head(1)
    tasks = [(r, OUT / "bins" / f"{r.dataset}_roi_bins.parquet", json.load(open(OUT / "bins" / f"{r.dataset}_bg.json"))["bg"]) for _, r in first.iterrows()]
    with Pool(ncores) as p:
        for f in p.imap_unordered(ctrl_one, tasks): print("ok", f, flush=True)

if __name__ == "__main__":
    st = sys.argv[1]
    if st in ("umi", "all"): stage_umi()
    if st in ("roi", "all"): stage_roi()
    if st in ("ctrl", "all"): stage_ctrl()
