"""tools/R4_roi.py — R4, stadio per ROI. Uso (root del repo): .venv/bin/python tools/R4_roi.py <KEY> <outdir>
KEY = A1 | A23 | A4 | A5 | A6. Scrive <outdir>/roi/<A>_<roi>_{bins.npz,cells.parquet,p2.npz,summary.json}; salta le ROI gia' scritte.
Pre-registrazione: docs/steps/R4.md (ff04e8f, e0c2341)."""
import sys, os, json, time
import numpy as np, pandas as pd, scipy.sparse as sp, pyarrow.parquet as pq
from scipy import ndimage
sys.path.insert(0, "tools")
import R4_common as C

M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"
KEY, OUT = sys.argv[1], sys.argv[2]
os.makedirs(f"{OUT}/roi", exist_ok=True)
DS = {"A1": ("Visium_HD_Mouse_Small_Intestine", ["A1"]), "A23": ("Visium_HD_Human_Colon_Cancer", ["A2", "A3"]),
      "A4": ("Visium_HD_Human_Lymph_Node_FFPE", ["A4"]), "A5": ("Visium_HD_6p5mm_Mouse_Brain", ["A5"]), "A6": ("Visium_HD_6p5mm_Human_Heart", ["A6"])}
SRC = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
ds, archs = DS[KEY]
d_exp = float(pd.read_csv("results/R6/R6_dexp.csv").query("dataset == 'pooled'")["d_exp_um"].iloc[0])
rois_all = pd.read_csv("results/R2/R2_rois_checked.csv")
rois = rois_all[rois_all.archetype.isin(archs)]
sq = lambda s: f"{D}/{ds}/binned_outputs_x/binned_outputs/{s}"

def p2_block(X, rng_b, rng_d):
    """Zeri in eccesso per i geni ammessi, con B repliche bootstrap delle cellule (pesi multinomiali), sui dati veri e con dropout."""
    out = {}
    for tag, Y in (("obs", X), ("drop", C.apply_dropout(X, C.DROPOUT, rng_d))):
        genes = C.eligible_genes(Y)
        ex = C.excess_zeros(Y, genes) if len(genes) else np.zeros(0)
        reps = np.zeros((C.B, len(genes)), np.float32)
        if len(genes):
            Yg = Y.tocsc()[:, genes].tocsr(); n = Yg.shape[0]
            W = rng_b.multinomial(n, np.ones(n) / n, size=C.B).astype(np.float64)   # B x cellule
            m1 = (Yg.T @ W.T).T / n                                                # B x geni
            m2 = (Yg.multiply(Yg).tocsr().T @ W.T).T / n
            nz = Yg.copy(); nz.data = (nz.data > 0).astype(float)
            z = 1.0 - (nz.T @ W.T).T / n
            var = np.maximum(m2 - m1 ** 2, 0.0)
            reps = (z - C.nb_zero_expected(m1.ravel(), var.ravel()).reshape(m1.shape)).astype(np.float32)
        out[tag] = (genes, ex, reps)
        # PH-2 (post hoc, dopo la revisione F06/F07): stessi geni e stessa matrice (stesso dropout), binomiale negativa con fattore di taglia
        lib = np.asarray(Y.sum(axis=1)).ravel().astype(float); sf = lib / lib.mean(); ok = sf > 0
        out[tag + "_sf"] = C.nb_sf_excess(Y[ok][:, genes], sf[ok]) if len(genes) else np.zeros(0)
        g05 = C.eligible_genes(Y, 0.5)                       # PH-3: soglia 0.5 UMI per cellula (linfonodo e tutti, sensibilita')
        out[tag + "_05"] = (g05, C.excess_zeros(Y, g05) if len(g05) else np.zeros(0), C.nb_sf_excess(Y[ok][:, g05], sf[ok]) if len(g05) else np.zeros(0))
    return out

t0 = time.time()
X, bc, gnames, gids = C.read_h5_csc(f"{sq('square_002um')}/filtered_feature_bc_matrix.h5")
tp = pq.read_table(f"{sq('square_002um')}/spatial/tissue_positions.parquet").to_pandas()
col_of = pd.Series(np.arange(len(bc)), index=bc)
X8, bc8, _, _ = C.read_h5_csc(f"{sq('square_008um')}/filtered_feature_bc_matrix.h5")
tp8 = pq.read_table(f"{sq('square_008um')}/spatial/tissue_positions.parquet").to_pandas()
col8 = pd.Series(np.arange(len(bc8)), index=bc8)
n_genes = X.shape[0]
print(f"[{KEY}] loaded 2um {X.shape} 8um {X8.shape} in {time.time()-t0:.0f}s", flush=True)

for _, R in rois.iterrows():
    A, roi, upp = R.archetype, R.roi_id, R.um_per_px
    tag = f"{A}_{roi}"; fsum = f"{OUT}/roi/{tag}_summary.json"
    if os.path.exists(fsum):
        print(f"[{KEY}] {tag} skip", flush=True); continue
    ridx = int(rois_all.index[(rois_all.archetype == A) & (rois_all.roi_id == roi)][0])
    rng_b = np.random.default_rng([C.SEED, ridx, 1]); rng_d = np.random.default_rng([C.SEED, ridx, 2])
    mask = np.load(f"{M}/{'R3b' if SRC[A]=='consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
    nuc = pd.read_parquet(f"{M}/R3b/nuclei/{tag}_consensus_native.parquet" if SRC[A] == "consensus" else f"results/R2/nuclei/{tag}_{SRC[A]}_native.parquet")
    keep = nuc.loc[nuc.keep, "label"].to_numpy()
    box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)].copy()
    box["rr"] = np.floor(box.pxl_row_in_fullres - R.r0).astype(int); box["cc"] = np.floor(box.pxl_col_in_fullres - R.c0).astype(int)
    filt = box[box.barcode.isin(col_of.index)]
    cols = col_of.reindex(filt.barcode).to_numpy(); Xb = X[:, cols]
    rr, cc = filt.rr.to_numpy(), filt.cc.to_numpy()
    umi, ngen = C.column_totals(Xb)
    # strati: maschera vera e spostata; corona 8 um (primaria) e 4 um (S-1)
    ms = C.shifted_mask(mask, upp)
    L = {}
    for mk, mm in (("", mask), ("_shift", ms)):
        L["nuc" + mk] = C.footprint_labels(mm, rr, cc, upp, d_exp, mode="nuc")
        L["exp" + mk] = C.footprint_labels(mm, rr, cc, upp, d_exp, mode="exp")
        L["exp4" + mk] = C.footprint_labels(mm, rr, cc, upp, 4.0, mode="exp")
    s = C.strata(L["nuc"], L["exp"], keep); s_sh = C.strata(L["nuc_shift"], L["exp_shift"], keep)
    s4 = C.strata(L["nuc"], L["exp4"], keep); s4_sh = C.strata(L["nuc_shift"], L["exp4_shift"], keep)
    # S-2: distanza del bin dal bordo del nucleo (dentro e fuori), per escludere i bin entro 2 um dal bordo
    def edge_far(mm, st):
        ins = (rr >= 0) & (rr < mm.shape[0]) & (cc >= 0) & (cc < mm.shape[1])
        d_in = ndimage.distance_transform_edt(mm > 0) * upp; d_out = ndimage.distance_transform_edt(mm == 0) * upp
        far = np.zeros(len(rr), bool)
        far[ins] = np.where(st[ins] == 1, d_in[rr[ins], cc[ins]] >= 2.0, d_out[rr[ins], cc[ins]] >= 2.0)
        return far
    far, far_sh = edge_far(mask, s), edge_far(ms, s_sh)
    # S-3 (linfonodo, cuore): solo nuclei con forma StarDist
    if SRC[A] == "consensus":
        cmap = pd.read_parquet(f"{M}/R3b/masks/{tag}_consensus_map.parquet")
        keep_sd = np.intersect1d(keep, cmap.loc[cmap.src == "sd", "label"].to_numpy())
        s_sd = C.strata(L["nuc"], L["exp"], keep_sd); s_sd_sh = C.strata(L["nuc_shift"], L["exp_shift"], keep_sd)
    else:
        keep_sd = None
    # G-3(c): impronte ricalcolate contro R6
    r6 = pd.read_parquet(f"{M}/R6/counts/{tag}_nuclei.parquet")
    g3c = {}
    agg = {}
    for fp, lab in (("exp", L["exp"]), ("nuc", L["nuc"]), ("shift", L["exp_shift"])):
        Cm, nb, um = C.aggregate_counts(Xb, lab, keep)
        agg[fp] = (Cm, nb, um)
        R6m = sp.load_npz(f"{M}/R6/counts/{tag}_{fp}.npz")
        same_mat = (Cm.shape == R6m.shape) and ((Cm.astype(np.int64) != R6m.astype(np.int64)).nnz == 0)
        same_tab = bool(np.array_equal(r6["label"].to_numpy(), keep) and np.array_equal(r6[f"n_bins_{fp}"].to_numpy(), nb)
                        and np.array_equal(r6[f"umi_{fp}"].to_numpy(), um))
        g3c[fp] = {"matrix_identical": bool(same_mat), "table_identical": same_tab}
    # cellule
    side = R.side_px * upp
    lab6 = pd.read_parquet(f"{M}/R6/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet", columns=["roi_id", "label", "cell_type"])
    lab6 = lab6[lab6.roi_id == roi].set_index("label").cell_type
    cells = nuc.loc[nuc.keep, ["label", "x_um", "y_um", "area_um2"]].reset_index(drop=True)
    cells["archetype"] = A; cells["roi_id"] = roi
    cells["interior"] = (cells.x_um >= C.INTERIOR_UM) & (cells.y_um >= C.INTERIOR_UM) & (cells.x_um <= side - C.INTERIOR_UM) & (cells.y_um <= side - C.INTERIOR_UM)
    cells["cell_type"] = lab6.reindex(cells.label).fillna("missing").to_numpy()
    cells["umi_exp"], cells["n_bins_exp"] = agg["exp"][2], agg["exp"][1]
    cells["umi_nuc"], cells["n_bins_nuc"] = agg["nuc"][2], agg["nuc"][1]
    if keep_sd is not None:
        cells["sd_shape"] = np.isin(cells.label, keep_sd)
    cells.to_parquet(f"{OUT}/roi/{tag}_cells.parquet", index=False)
    # P-2 per tipo, cellule interne
    Cexp = agg["exp"][0]
    p2 = {}; p2_meta = []
    for t, g in cells[cells.interior & ~cells.cell_type.isin(list(C.NON_TYPES) + ["missing"])].groupby("cell_type"):
        if len(g) < C.MIN_CELLS:
            p2_meta.append({"type": t, "n_cells": int(len(g)), "used": False}); continue
        res = p2_block(Cexp[g.index.to_numpy()], rng_b, rng_d)
        for k in ("obs", "drop"):
            genes, ex, reps = res[k]
            p2[f"{t}|{k}|genes"] = genes; p2[f"{t}|{k}|ex"] = ex; p2[f"{t}|{k}|reps"] = reps; p2[f"{t}|{k}|exsf"] = res[k + "_sf"]
            g05, e05, e05sf = res[k + "_05"]
            p2[f"{t}|{k}|genes05"] = g05; p2[f"{t}|{k}|ex05"] = e05; p2[f"{t}|{k}|ex05sf"] = e05sf
        p2_meta.append({"type": t, "n_cells": int(len(g)), "used": True, "n_genes_obs": int(len(res["obs"][0])), "n_genes_drop": int(len(res["drop"][0])),
                        "n_genes05_obs": int(len(res["obs_05"][0]))})
    np.savez_compressed(f"{OUT}/roi/{tag}_p2.npz", **p2)
    # 8 um nella ROI
    b8 = tp8[(tp8.pxl_row_in_fullres >= R.r0) & (tp8.pxl_row_in_fullres < R.r1) & (tp8.pxl_col_in_fullres >= R.c0) & (tp8.pxl_col_in_fullres < R.c1)]
    b8 = b8[b8.barcode.isin(col8.index)]
    umi8, ngen8 = C.column_totals(X8[:, col8.reindex(b8.barcode).to_numpy()])
    # D-5: spettro dei bin acellulari della ROI (tessere 32 x 32, almeno meta' acellulari)
    ar, ac = filt.array_row.to_numpy(), filt.array_col.to_numpy()
    img = np.zeros((ar.max() - ar.min() + 1, ac.max() - ac.min() + 1)); val = np.zeros(img.shape, bool)
    img[ar - ar.min(), ac - ac.min()] = umi; val[ar - ar.min(), ac - ac.min()] = s == 0
    S, ktiles = C.welch_spectrum(img, val, tile=32, min_valid=0.5)
    nanpk = {"ratio": float("nan"), "present": None, "fy": float("nan"), "fx": float("nan")}
    spec = C.spectral_peak(S) if ktiles >= 4 else dict(nanpk)
    spec_off = C.spectral_peak(S, off_axis=True) if ktiles >= 4 else dict(nanpk)
    spec.update({f"offaxis_{k}": v for k, v in spec_off.items()}); spec["n_tiles"] = int(ktiles)
    np.savez_compressed(f"{OUT}/roi/{tag}_bins.npz", umi=umi, ngen=ngen, s=s, s_shift=s_sh, s4=s4, s4_shift=s4_sh, far=far, far_shift=far_sh,
                        umi8=umi8, ngen8=ngen8, x_um=cc * upp, y_um=rr * upp, side_um=np.float64(R.side_px * upp), **({"s_sd": s_sd, "s_sd_shift": s_sd_sh} if keep_sd is not None else {}))
    summ = {"archetype": A, "roi_id": roi, "roi_index": ridx, "n_genes": int(n_genes), "n_bins_box": int(len(box)), "n_bins_tissue": int(len(filt)),
            "n_bins8_tissue": int(len(b8)), "n_nuclei_keep": int(len(keep)), "d_exp_um": d_exp, "G3c": g3c, "p2_meta": p2_meta, "spectrum": spec,
            "elapsed_s": round(time.time() - t0, 1)}
    json.dump(summ, open(fsum, "w"), indent=1)
    print(f"[{KEY}] {tag}: bins {len(filt)} 8um {len(b8)} nuclei {len(keep)} G3c {all(v['matrix_identical'] and v['table_identical'] for v in g3c.values())} "
          f"p2 {[(m['type'], m.get('n_genes_obs')) for m in p2_meta if m['used']]} {time.time()-t0:.0f}s", flush=True)
print(f"[{KEY}] done {time.time()-t0:.0f}s", flush=True)
