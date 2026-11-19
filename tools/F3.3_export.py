"""tools/F3.3_export.py - F3.3: aree nel formato comune dei due lati, come tools/F3.1_export_real.py: quadratini di Space Ranger con il centro
nel riquadro, conteggi dalla matrice grezza a 2 um ristretta al pannello (feature della filtrata) e, nel tumore del colon con 'common', ai
geni comuni a tutte le sezioni del tumore e al pannello del simulato (DEC-097, punto 4); nucleo sotto il centro del quadratino dalla
maschera dell'area (Cellpose RGB o consenso), solo nuclei keep; in_reg = centro corretto con lo scarto di registrazione della tabella
nelle regioni dell'area e a >= 25 um dal bordo. Tipi: nelle aree di taratura le etichette adottate di R6 (per il controllo G-1 d); nelle
aree di verifica nessuno (etichette trasferite dopo). Uso: .venv/bin/python tools/F3.3_export.py <panel|common> <aree|IN|holdout>"""
import sys, os, json, time, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow as pa, pyarrow.parquet as pq, shapely, h5py
sys.path.insert(0, "tools")
from R6_common import read_h5_csc
import S22_common as S
M = "/mnt/micron/geo_spatialtrans"; ROOT = os.environ.get("F33_ROOT", f"{M}/F3.3"); SEG = f"{ROOT}/seg"; RUNS = f"{ROOT}/runs"; FRAME_UM = 25.0
CG = "results/F3.3/F3.3_crc_common_genes.txt"
genes_mode, sel = sys.argv[1], sys.argv[2]
A = pd.read_csv("results/F3.3/F3.3_areas.csv")
A = A[A.level == "IN"] if sel == "IN" else A[A.level != "IN"] if sel == "holdout" else A[A.area.isin(sel.split(","))]
if genes_mode == "common": A = A[A.crc]
def panel(sq):
    with h5py.File(f"{sq}/square_002um/filtered_feature_bc_matrix.h5", "r") as h: return h["matrix"]["features"]["id"][:].astype(str)
if genes_mode == "common" and not os.path.exists(CG):
    T = pd.read_csv("results/F3.3/F3.3_areas.csv"); sets = [set(panel(sq)) for sq in T[T.crc].sq.unique()]
    sets.append(set(pd.read_csv(f"{M}/F3.1/main/sim/A2_r1_k0/genes.tsv", header=None)[0].astype(str)))
    pd.Series(sorted(set.intersection(*sets))).to_csv(CG, index=False, header=False)
keep_ids = set(pd.read_csv(CG, header=None)[0].astype(str)) if genes_mode == "common" else None
_cache = {}
def section(sq):
    if sq not in _cache:
        _cache.clear(); d2 = f"{sq}/square_002um"
        X, bc, gnames, gids = read_h5_csc(f"{d2}/raw_feature_bc_matrix.h5")
        g = np.isin(np.asarray(gids).astype(str), panel(sq))
        if keep_ids is not None: g &= np.isin(np.asarray(gids).astype(str), list(keep_ids))
        X = X[np.flatnonzero(g)]; gids = np.asarray(gids)[g]
        tp = pq.read_table(f"{d2}/spatial/tissue_positions.parquet").to_pandas()
        _cache[sq] = (X, pd.Series(np.arange(len(bc)), index=bc), gids, tp)
    return _cache[sq]
for _, r in A.iterrows():
    t0 = time.time(); name = r.area + ("_cg" if genes_mode == "common" else ""); out = f"{RUNS}/{name}"; os.makedirs(out, exist_ok=True)
    X, col_of, gids, tp = section(r.sq); upp = float(r.upp); side = float(r.side_px * upp)
    box = tp[(tp.pxl_row_in_fullres >= r.r0) & (tp.pxl_row_in_fullres < r.r1) & (tp.pxl_col_in_fullres >= r.c0) & (tp.pxl_col_in_fullres < r.c1)].copy()
    box = box.sort_values(["array_row", "array_col"]).reset_index(drop=True)
    rr_um = (box.pxl_row_in_fullres.to_numpy() - r.r0) * upp; cc_um = (box.pxl_col_in_fullres.to_numpy() - r.c0) * upp
    di = np.array([r.d_img_r, r.d_img_c], float); y_c, x_c = rr_um - di[0], cc_um - di[1]
    t = pd.read_csv(f"{SEG}/{r.area}/regions.tsv", sep="\t"); R = shapely.union_all(shapely.from_wkb([bytes.fromhex(h) for h in t.wkb]))
    fr = (y_c >= FRAME_UM) & (y_c <= side - FRAME_UM) & (x_c >= FRAME_UM) & (x_c <= side - FRAME_UM)
    in_reg = np.zeros(len(box), bool); in_reg[fr] = shapely.contains_xy(R, x_c[fr], y_c[fr])
    seg = r.seg; mask = np.load(f"{SEG}/{r.area}/{seg}_native.npz")["labels"]; nuc = pd.read_parquet(f"{SEG}/{r.area}/{seg}_native.parquet")
    keep = set(nuc.loc[nuc.keep, "label"].astype(int)); H, W = mask.shape
    def lab_at(r_px, c_px):
        r_px = np.floor(r_px).astype(np.int64); c_px = np.floor(c_px).astype(np.int64)
        ok = (r_px >= 0) & (r_px < H) & (c_px >= 0) & (c_px < W); lab = np.zeros(len(r_px), np.int64); lab[ok] = mask[r_px[ok], c_px[ok]]
        lab[~np.isin(lab, list(keep))] = 0; return lab
    label = lab_at(box.pxl_row_in_fullres.to_numpy() - r.r0, box.pxl_col_in_fullres.to_numpy() - r.c0); label_dimg = lab_at(y_c / upp, x_c / upp)
    cols = col_of.reindex(box.barcode).to_numpy(); has = ~np.isnan(cols)
    Xsub = X[:, cols[has].astype(np.int64)].tocsc(); T = Xsub.T.tocoo(); idx = np.flatnonzero(has)
    pq.write_table(pa.Table.from_pandas(pd.DataFrame({"bin": idx[T.row].astype(np.int64), "gene": T.col.astype(np.int64), "n": T.data.astype(np.int32)})), f"{out}/X.parquet")
    bins = pd.DataFrame({"bin_id": np.arange(len(box)), "row": box.array_row.to_numpy(), "col": box.array_col.to_numpy(), "x_um": cc_um, "y_um": rr_um,
                         "in_tissue": box.in_tissue.to_numpy() == 1, "in_filtered": has, "in_reg": in_reg, "label": label, "label_dimg": label_dimg})
    bins.to_parquet(f"{out}/bins.parquet", index=False)
    cells = pd.DataFrame({"label": sorted(keep)})
    if r.level == "IN":
        lab = pd.read_parquet(f"{M}/R6/labels/{r.A}_exp_umi20_seed20261022_nuclei_ph1.parquet"); lab = lab[lab.roi_id == r.roi_id]
        cells = cells.merge(lab[["label", "cell_type", "classifiable"]], on="label", how="left")
        cells["cell_type_r6"] = cells.cell_type; cells["cell_type"] = cells.cell_type.replace(S.MERGE)
    else:
        cells["cell_type"] = None; cells["classifiable"] = None; cells["cell_type_r6"] = None
    cells.to_parquet(f"{out}/cells.parquet", index=False)
    pd.Series(gids).to_csv(f"{out}/genes.tsv", index=False, header=False)
    info = {"tag": name, "area": r.area, "archetype": r.A, "level": r.level, "section": r.section, "genes": genes_mode, "side_um": side, "d_img_um": di.tolist(),
            "n_bins": int(len(box)), "n_in_filtered": int(has.sum()), "n_reg": int(in_reg.sum()), "n_reg_not_filtered": int((in_reg & ~has).sum()),
            "umi": int(T.data.sum()), "umi_reg": int(np.asarray(Xsub.sum(axis=0)).ravel()[in_reg[has]].sum()), "n_genes": int(len(gids)),
            "n_lab_bins": int((label > 0).sum()), "n_lab_bins_dimg": int((label_dimg > 0).sum()), "n_nuclei_keep": len(keep), "sec": round(time.time() - t0, 1)}
    json.dump(info, open(f"{out}/export.json", "w")); print(json.dumps(info), flush=True)
