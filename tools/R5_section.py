"""tools/R5_section.py — R5, stadio per sezione: UMI totali per bin a 2 µm dalla matrice grezza (anche fuori dal tessuto), immagine 3350 x 3350,
trasformazione affine griglia -> immagine, controlli G-3 (a) grezza = filtrata sui bin in tessuto, (b) 8 µm = somma dei 16 bin a 2 µm.
Uso: .venv/bin/python tools/R5_section.py <KEY> <outdir>   (KEY = A1 | A23 | A4 | A5 | A6). Scrive <outdir>/section/<KEY>.npz e <KEY>_G3.json.
Pre-registrazione: docs/steps/R5.md (c5c9fce, ceca77b)."""
import sys, os, json, time
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq

D = "/mnt/micron/geo_spatialtrans/data_real"
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}

def feature_ids(path):
    with h5py.File(path, "r") as h:
        return h["matrix"]["features"]["id"][:].astype(str)

def totals(path, keep_ids=None):
    """UMI totali per barcode (tutte le feature sono 'Gene Expression': controllato). keep_ids: se dato, solo quelle feature
    (D-3: la grezza contiene feature fuori dal pannello della filtrata; si usa il pannello della filtrata, come S1.3c 'com')."""
    with h5py.File(path, "r") as h:
        m = h["matrix"]
        ft = set(m["features"]["feature_type"][:].astype(str))
        ip = m["indptr"][:]
        data = m["data"][:].astype(np.int64)
        if keep_ids is not None:
            ids = m["features"]["id"][:].astype(str)
            keep = np.isin(ids, keep_ids)
            data = data * keep[m["indices"][:]]
        bc = m["barcodes"][:].astype(str)
    cs = np.concatenate([[0], np.cumsum(data)])
    return pd.Series(cs[ip[1:]] - cs[ip[:-1]], index=bc), ft

def affine(tp):
    X = np.c_[tp.array_row.to_numpy(float), tp.array_col.to_numpy(float), np.ones(len(tp))]
    out = {}
    for k in ("pxl_row_in_fullres", "pxl_col_in_fullres"):
        coef, *_ = np.linalg.lstsq(X, tp[k].to_numpy(float), rcond=None)
        out[k] = coef; out[k + "_maxres"] = float(np.max(np.abs(X @ coef - tp[k].to_numpy(float))))
    M = np.array([[out["pxl_row_in_fullres"][0], out["pxl_row_in_fullres"][1]], [out["pxl_col_in_fullres"][0], out["pxl_col_in_fullres"][1]]])
    t = np.array([out["pxl_row_in_fullres"][2], out["pxl_col_in_fullres"][2]])
    return M, t, max(out["pxl_row_in_fullres_maxres"], out["pxl_col_in_fullres_maxres"])

if __name__ == "__main__":
    KEY, OUT = sys.argv[1], sys.argv[2]; os.makedirs(f"{OUT}/section", exist_ok=True)
    ds = DS[KEY]; sq = lambda s: f"{D}/{ds}/binned_outputs_x/binned_outputs/{s}"
    t0 = time.time()
    tp = pq.read_table(f"{sq('square_002um')}/spatial/tissue_positions.parquet").to_pandas()
    panel = feature_ids(f"{sq('square_002um')}/filtered_feature_bc_matrix.h5")
    raw_ids = feature_ids(f"{sq('square_002um')}/raw_feature_bc_matrix.h5")
    raw_all, ft_raw = totals(f"{sq('square_002um')}/raw_feature_bc_matrix.h5")
    raw, _ = totals(f"{sq('square_002um')}/raw_feature_bc_matrix.h5", panel)
    fil, ft_fil = totals(f"{sq('square_002um')}/filtered_feature_bc_matrix.h5")
    nr, nc = int(tp.array_row.max()) + 1, int(tp.array_col.max()) + 1
    U = np.zeros((nr, nc), np.int32); IT = np.zeros((nr, nc), bool)
    u = raw.reindex(tp.barcode).fillna(0).to_numpy(np.int64)
    U[tp.array_row.to_numpy(), tp.array_col.to_numpy()] = u
    IT[tp.array_row.to_numpy(), tp.array_col.to_numpy()] = tp.in_tissue.to_numpy().astype(bool)
    M, t, maxres = affine(tp)
    # G-3 (a): sui bin in tessuto la grezza e' uguale alla filtrata; i barcode della filtrata sono tutti in tessuto
    it_bc = tp.barcode[tp.in_tissue == 1]
    a_same = bool(np.array_equal(raw.reindex(it_bc).fillna(0).to_numpy(np.int64), fil.reindex(it_bc).fillna(0).to_numpy(np.int64)))
    a_cover = bool(fil.index.isin(it_bc).all())
    # G-3 (b): 8 µm grezzo = somma dei 16 bin a 2 µm (bin assenti dalla grezza = 0)
    tp8 = pq.read_table(f"{sq('square_008um')}/spatial/tissue_positions.parquet").to_pandas()
    raw8, ft8 = totals(f"{sq('square_008um')}/raw_feature_bc_matrix.h5", panel)
    n4r, n4c = -(-nr // 4), -(-nc // 4)
    Up = np.zeros((n4r * 4, n4c * 4), np.int64); Up[:nr, :nc] = U
    S = Up.reshape(n4r, 4, n4c, 4).sum(axis=(1, 3))
    u8 = raw8.reindex(tp8.barcode).fillna(0).to_numpy(np.int64)
    s8 = S[tp8.array_row.to_numpy(), tp8.array_col.to_numpy()]
    b_diff = int(np.sum(u8 != s8))
    it_all = raw_all.reindex(it_bc).fillna(0).to_numpy(np.int64); it_fil = fil.reindex(it_bc).fillna(0).to_numpy(np.int64)
    res = {"dataset": ds, "n_features_raw": int(len(raw_ids)), "n_features_filtered": int(len(panel)),
           "panel_subset_of_raw": bool(np.isin(panel, raw_ids).all()),
           "G3a_as_written_all_raw_features_equal": bool(np.array_equal(it_all, it_fil)), "G3a_as_written_n_bins_differing": int(np.sum(it_all != it_fil)),
           "G3a_as_written_umi_in_tissue_raw_all": int(it_all.sum()), "G3a_as_written_umi_in_tissue_filtered": int(it_fil.sum()), "grid": [nr, nc], "n_bins_positions": int(len(tp)), "n_barcodes_raw": int(len(raw)), "n_barcodes_filtered": int(len(fil)),
           "feature_types": sorted(ft_raw | ft_fil | ft8), "umi_raw_total": int(raw.sum()), "umi_grid_total": int(U.sum()),
           "raw_barcodes_not_in_positions": int((~raw.index.isin(tp.barcode)).sum()),
           "G3a_raw_equals_filtered_in_tissue": a_same, "G3a_filtered_all_in_tissue": a_cover, "G3b_n_8um_bins_differing": b_diff,
           "n_8um_bins": int(len(tp8)), "affine_M": M.tolist(), "affine_t": t.tolist(), "affine_max_residual_px": maxres,
           "elapsed_s": round(time.time() - t0, 1)}
    np.savez_compressed(f"{OUT}/section/{KEY}.npz", U=U, in_tissue=IT, M=M, t=t)
    json.dump(res, open(f"{OUT}/section/{KEY}_G3.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if not k.startswith("affine_M")}), flush=True)
