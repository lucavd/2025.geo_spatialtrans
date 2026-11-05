"""tools/S2.1_desc_real.py - S2.1, descrittivi D-1...D-3 sui vetrini veri (senza previsione). Definizioni: docs/steps/S2.1.md.
Uso (root del repo): .venv/bin/python tools/S2.1_desc_real.py [outcsv]
Per ogni area (30): quadratini di 2 um di Space Ranger con il centro nel riquadro, coordinate dell'area x = (colonna - c0) * um_per_px,
y = (riga - r0) * um_per_px (stesse dei livelli simulati); in tessuto = maschera valida ds4 di R2 al centro (la stessa che ha fatto le regioni
simulate). Blocchi di 8 um = (array_row // 4, array_col // 4) con tutti i 16 figli nel riquadro e almeno un figlio in tessuto.
Fonti: Space Ranger 4 (in_nucleus, cell_id di barcode_mappings; 4 sezioni) e il segmentatore dello Step 1 (maschere native di R2, consenso
R3b per linfonodo e cuore, nuclei 'keep' come R4): centro del quadratino nel nucleo; nuclei che toccano un blocco = etichette su 7 x 7 punti
per figlio."""
import sys, numpy as np, pandas as pd, pyarrow.parquet as pq
from PIL import Image
M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"
OUT = sys.argv[1] if len(sys.argv) > 1 else "results/S2.1/S2.1_desc_real_roi.csv"
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A3": "Visium_HD_Human_Colon_Cancer",
      "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
SRC = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
def cls(v, pre):
    v = np.asarray(v); return {f"{pre}_mean": v.mean(), f"{pre}_p0": np.mean(v == 0), f"{pre}_p1": np.mean(v == 1), f"{pre}_p2": np.mean(v == 2),
                               f"{pre}_p3": np.mean(v == 3), f"{pre}_p4": np.mean(v >= 4)}
rows = []
for ds in dict.fromkeys(DS.values()):
    tp = pq.read_table(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet").to_pandas()
    bm = None
    if ds != DS["A1"]:
        bm = pq.read_table(f"{D}/{ds}/{ds}_barcode_mappings.parquet", columns=["square_002um", "cell_id", "in_nucleus", "in_cell"]).to_pandas()
        bm = bm.set_index("square_002um")
    for _, R in rois[rois.dataset == ds].iterrows():
        A, roi, upp = R.archetype, R.roi_id, R.um_per_px; tag = f"{A}_{roi}"
        box = tp[(tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1) & (tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1)].copy()
        box["x"] = (box.pxl_col_in_fullres - R.c0) * upp; box["y"] = (box.pxl_row_in_fullres - R.r0) * upp
        m = np.array(Image.open(f"results/R2/tissue_masks/{tag}_valid_ds4.png")); m = (m[..., 0] if m.ndim == 3 else m) > (127 if m.max() > 1 else 0.5)
        px = R.side_px * upp / m.shape[1]
        ci = np.clip(np.floor(box.x / px).astype(int), 0, m.shape[1] - 1); ri = np.clip(np.floor(box.y / px).astype(int), 0, m.shape[0] - 1)
        box["tis"] = m[ri, ci]
        box["b8"] = (box.array_row // 4) * 100000 + box.array_col // 4
        n16 = box.groupby("b8").size(); anyt = box.groupby("b8").tis.any()
        good8 = n16.index[(n16 == 16) & anyt.reindex(n16.index).values]
        out = dict(tag=tag, archetype=A, n_bins_tissue=int(box.tis.sum()), n_8um_tissue=len(good8))
        # segmentatore dello Step 1
        lab = np.load(f"{M}/{'R3b' if SRC[A] == 'consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
        nuc = pd.read_parquet(f"{M}/R3b/nuclei/{tag}_consensus_native.parquet" if SRC[A] == "consensus" else f"results/R2/nuclei/{tag}_{SRC[A]}_native.parquet")
        keep = np.zeros(int(lab.max()) + 1, bool); kl = nuc.loc[nuc.keep, "label"].to_numpy().astype(int); keep[kl[kl <= lab.max()]] = True
        def lab_at(x, y):
            c = np.floor(x / upp).astype(int); r = np.floor(y / upp).astype(int)
            ok = (c >= 0) & (r >= 0) & (c < lab.shape[1]) & (r < lab.shape[0]); v = np.zeros(len(x), int); v[ok] = lab[r[ok], c[ok]]
            v[~keep[v]] = 0; return v
        cen = lab_at(box.x.to_numpy(), box.y.to_numpy())
        out["seg_d1_centre_in_nucleus"] = float(np.mean(cen[box.tis.to_numpy()] > 0))
        off = (np.arange(7) - 3) * (2.0 / 7.0)
        dx, dy = np.meshgrid(off, off); dx = dx.ravel(); dy = dy.ravel()
        sub = box[box.b8.isin(good8)]
        X = (sub.x.to_numpy()[:, None] + dx[None, :]).ravel(); Y = (sub.y.to_numpy()[:, None] + dy[None, :]).ravel()
        L = lab_at(X, Y); B = np.repeat(sub.b8.to_numpy(), len(dx))
        pairs = pd.DataFrame({"b8": B, "l": L}); pairs = pairs[pairs.l > 0].drop_duplicates()
        nN = pairs.groupby("b8").size().reindex(good8, fill_value=0).to_numpy()
        out.update(cls(nN, "seg_d2_nuclei"))
        if bm is not None:
            j = bm.reindex(box.barcode)
            inn = j.in_nucleus.fillna(False).to_numpy().astype(bool)
            out["sr_d1_centre_in_nucleus"] = float(np.mean(inn[box.tis.to_numpy()]))
            sb = pd.DataFrame({"b8": box.b8.to_numpy(), "cell": j.cell_id.to_numpy(), "inn": inn})
            sb = sb[sb.b8.isin(good8)]
            nn = sb[sb.inn & sb.cell.notna()].drop_duplicates(["b8", "cell"]).groupby("b8").size().reindex(good8, fill_value=0).to_numpy()
            nc = sb[sb.cell.notna()].drop_duplicates(["b8", "cell"]).groupby("b8").size().reindex(good8, fill_value=0).to_numpy()
            out.update(cls(nn, "sr_d2_nuclei")); out.update(cls(nc, "sr_d3_cells"))
            out["sr_frac_bins_in_cell"] = float(np.mean(j.cell_id.notna().to_numpy()[box.tis.to_numpy()]))
        rows.append(out); print(tag, {k: round(v, 3) for k, v in out.items() if isinstance(v, float) and ("mean" in k or "d1" in k)}, flush=True)
pd.DataFrame(rows).to_csv(OUT, index=False)
