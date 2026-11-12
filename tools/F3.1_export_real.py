"""tools/F3.1_export_real.py - F3.1: aree vere nel formato comune dei due lati (vedi tools/F3.1_export_sim.R).
Quadratini: tutti quelli di Space Ranger con il centro nel riquadro della ROI (results/R2/R2_rois_checked.csv); conteggi dalla matrice
grezza a 2 um ristretta ai geni del pannello (feature della filtrata), come S2.5 (in_filtered: quadratino nella filtrata).
label: nucleo della maschera adottata in R6 (Cellpose RGB; consenso in linfonodo e cuore) sotto il centro del quadratino, registrazione
di Space Ranger come arriva (R6, footprint 'nuc'); label_dimg: la stessa con lo scarto sistematico d_image di R5/S2.5 tolto (sensibilita').
Solo i nuclei tenuti (keep) di R2/R3b; gli altri -> 0. Tipi: etichette adottate di R6 (exp, umi20, seme 20261022, ph1), come S2.2.
in_reg: definizione reg di S2.5 (centro corretto con d_image nelle regioni dello Step 1 e a >= 25 um dal bordo).
Uso: .venv/bin/python tools/F3.1_export_real.py <out_root> <tag> [<tag> ...]"""
import sys, os, json, time, importlib.util, numpy as np, pandas as pd, scipy.sparse as sp, pyarrow.parquet as pq, shapely
sys.path.insert(0, "tools")
from R6_common import read_h5_csc
import S22_common as S   # vocabolario dei preset (MERGE), come S2.2
def _load(name, path):
    s = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A3": "Visium_HD_Human_Colon_Cancer",
      "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
SRC = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
FRAME_UM = 25.0
out_root = sys.argv[1]; tags = sys.argv[2:]
RL = _load("s23real", "tools/S2.3_real.py")
PJ = json.load(open(f"{RL.S20}/params_flat.json"))
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois["tag"] = rois.archetype + "_" + rois.roi_id
_cache = {}
def section(ds):
    if ds not in _cache:
        _cache.clear()
        sq = f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um"
        X, bc, gnames, gids = read_h5_csc(f"{sq}/raw_feature_bc_matrix.h5")          # grezza: tutti i quadratini, come S2.5
        with __import__("h5py").File(f"{sq}/filtered_feature_bc_matrix.h5", "r") as h: fid = h["matrix"]["features"]["id"][:].astype(str)
        g = np.isin(np.asarray(gids).astype(str), fid)                                  # solo i geni del pannello (feature della filtrata)
        X = X[np.flatnonzero(g)]; gnames = np.asarray(gnames)[g]; gids = np.asarray(gids)[g]
        tp = pq.read_table(f"{sq}/spatial/tissue_positions.parquet").to_pandas()
        _cache[ds] = (X, pd.Series(np.arange(len(bc)), index=bc), gnames, gids, tp)
    return _cache[ds]
for tag in tags:
    t0 = time.time(); Rw = rois.set_index("tag").loc[tag]; A = Rw.archetype; out = f"{out_root}/{tag}"; os.makedirs(out, exist_ok=True)
    X, col_of, gnames, gids, tp = section(DS[A]); upp = float(Rw.um_per_px); side = float(Rw.side_px * upp)
    box = tp[(tp.pxl_row_in_fullres >= Rw.r0) & (tp.pxl_row_in_fullres < Rw.r1) & (tp.pxl_col_in_fullres >= Rw.c0) & (tp.pxl_col_in_fullres < Rw.c1)].copy()
    box = box.sort_values(["array_row", "array_col"]).reset_index(drop=True)
    # posizione (um) nel riferimento della ROI, come S2.5: registrazione di Space Ranger, poi meno d_image
    rr_um = (box.pxl_row_in_fullres.to_numpy() - Rw.r0) * upp; cc_um = (box.pxl_col_in_fullres.to_numpy() - Rw.c0) * upp
    prm = PJ[RL.NC.KEY[A]]["C2"]; Dd = RL.NC.roi_data(Rw); di = np.asarray(RL.NC.d_image(Dd, prm), float)
    y_c, x_c = rr_um - di[0], cc_um - di[1]
    R = RL.region_geom(tag)
    fr = (y_c >= FRAME_UM) & (y_c <= side - FRAME_UM) & (x_c >= FRAME_UM) & (x_c <= side - FRAME_UM)
    in_reg = np.zeros(len(box), bool); in_reg[fr] = shapely.contains_xy(R, x_c[fr], y_c[fr])
    # nuclei
    mask = np.load(f"{M}/{'R3b' if SRC[A]=='consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
    nuc = pd.read_parquet(f"{M}/R3b/nuclei/{tag}_consensus_native.parquet" if SRC[A] == "consensus" else f"results/R2/nuclei/{tag}_{SRC[A]}_native.parquet")
    keep = set(nuc.loc[nuc.keep, "label"].astype(int))
    H, W = mask.shape
    def lab_at(r_px, c_px):
        r_px = np.floor(r_px).astype(np.int64); c_px = np.floor(c_px).astype(np.int64)
        ok = (r_px >= 0) & (r_px < H) & (c_px >= 0) & (c_px < W); lab = np.zeros(len(r_px), np.int64)
        lab[ok] = mask[r_px[ok], c_px[ok]]
        lab[~np.isin(lab, list(keep))] = 0
        return lab
    label = lab_at(box.pxl_row_in_fullres.to_numpy() - Rw.r0, box.pxl_col_in_fullres.to_numpy() - Rw.c0)   # = R6 footprint 'nuc'
    label_dimg = lab_at(y_c / upp, x_c / upp)
    # conteggi
    cols = col_of.reindex(box.barcode).to_numpy()
    has = ~np.isnan(cols)
    Xsub = X[:, cols[has].astype(np.int64)].tocsc()
    T = Xsub.T.tocoo(); idx = np.flatnonzero(has)
    pq.write_table(__import__("pyarrow").Table.from_pandas(pd.DataFrame({"bin": idx[T.row].astype(np.int64), "gene": T.col.astype(np.int64), "n": T.data.astype(np.int32)})), f"{out}/X.parquet")
    bins = pd.DataFrame({"bin_id": np.arange(len(box)), "row": box.array_row.to_numpy(), "col": box.array_col.to_numpy(),
                         "x_um": cc_um, "y_um": rr_um, "in_tissue": box.in_tissue.to_numpy() == 1, "in_filtered": has,
                         "in_reg": in_reg, "label": label, "label_dimg": label_dimg})
    bins.to_parquet(f"{out}/bins.parquet", index=False)
    lab = pd.read_parquet(f"{M}/R6/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet"); lab = lab[lab.roi_id == Rw.roi_id]
    cells = pd.DataFrame({"label": sorted(keep)}); cells = cells.merge(lab[["label", "cell_type", "classifiable"]], on="label", how="left")
    cells["cell_type_r6"] = cells.cell_type; cells["cell_type"] = cells.cell_type.replace(S.MERGE)
    cells.to_parquet(f"{out}/cells.parquet", index=False)
    pd.Series(gids).to_csv(f"{out}/genes.tsv", index=False, header=False)   # identificativi Ensembl (unici), come il simulato
    info = {"tag": tag, "archetype": A, "side_um": side, "d_img_um": di.tolist(), "n_bins": int(len(box)), "n_in_filtered": int(has.sum()),
            "n_reg": int(in_reg.sum()), "n_reg_not_filtered": int((in_reg & ~has).sum()), "umi": int(T.data.sum()),
            "umi_reg": int(np.asarray(Xsub.sum(axis=0)).ravel()[in_reg[has]].sum()), "n_lab_bins": int((label > 0).sum()),
            "n_lab_bins_dimg": int((label_dimg > 0).sum()), "n_nuclei_keep": len(keep), "sec": round(time.time() - t0, 1)}
    json.dump(info, open(f"{out}/export.json", "w"))
    print(json.dumps(info), flush=True)
