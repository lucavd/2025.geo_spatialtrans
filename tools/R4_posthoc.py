"""tools/R4_posthoc.py — R4, analisi post hoc (dichiarate dopo i risultati, non previsioni).
  roi <KEY> <outdir> : PH-1 rapporto nucleo/corona pesato per cellula (rapporto fra le medie per cellula; maschera vera e spostata) e mediana
                       dei rapporti per cellula (revisione F08); PH-4 parita' di riga e colonna
                       della griglia a 2 um (ROI e sezione intera, bin in tessuto). -> <outdir>/posthoc/<A>_<roi>_cells.parquet, <KEY>_parity.csv
  (PH-2 e PH-3, zeri in eccesso con fattore di taglia e soglia 0.5, sono calcolati in tools/R4_roi.py sulla stessa matrice e sullo stesso
   dropout di P-2/CP-2, revisione F07; qui solo il riassunto, bootstrap sulle sole ROI)
  summary <outdir>   : tabelle <outdir>/tables/R4_PH*.csv"""
import sys, os, json, glob
import numpy as np, pandas as pd, scipy.sparse as sp, pyarrow.parquet as pq
sys.path.insert(0, "tools")
import R4_common as C

M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"
MODE, ARGS = sys.argv[1], sys.argv[2:]
DS = {"A1": ("Visium_HD_Mouse_Small_Intestine", ["A1"]), "A23": ("Visium_HD_Human_Colon_Cancer", ["A2", "A3"]),
      "A4": ("Visium_HD_Human_Lymph_Node_FFPE", ["A4"]), "A5": ("Visium_HD_6p5mm_Mouse_Brain", ["A5"]), "A6": ("Visium_HD_6p5mm_Human_Heart", ["A6"])}
SRC = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]

def per_label(lab_nuc, lab_exp, keep, u):
    """Somme di UMI e numero di bin nel nucleo e nella corona di ogni nucleo tenuto."""
    pos = pd.Series(np.arange(len(keep)), index=keep)
    out = {}
    for nm, sel, lab in (("nuc", np.isin(lab_nuc, keep), lab_nuc), ("ring", (lab_nuc == 0) & np.isin(lab_exp, keep), lab_exp)):
        idx = pos.reindex(lab[sel]).to_numpy()
        out[f"umi_{nm}"] = np.bincount(idx, weights=u[sel], minlength=len(keep)); out[f"n_{nm}"] = np.bincount(idx, minlength=len(keep))
    return out

if MODE == "roi":
    KEY, OUT = ARGS; os.makedirs(f"{OUT}/posthoc", exist_ok=True)
    ds, archs = DS[KEY]; sq = lambda s: f"{D}/{ds}/binned_outputs_x/binned_outputs/{s}"
    d_exp = float(pd.read_csv("results/R6/R6_dexp.csv").query("dataset == 'pooled'")["d_exp_um"].iloc[0])
    rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois = rois[rois.archetype.isin(archs)]
    X, bc, _, _ = C.read_h5_csc(f"{sq('square_002um')}/filtered_feature_bc_matrix.h5")
    tp = pq.read_table(f"{sq('square_002um')}/spatial/tissue_positions.parquet").to_pandas()
    umi_all, _ = C.column_totals(X)
    sec = pd.DataFrame({"barcode": bc, "umi": umi_all}).merge(tp[["barcode", "array_row", "array_col"]], on="barcode")
    par = [{"unit": "section", "dataset": ds, "row_par": rp, "col_par": cp, "n_bins": len(g), "mean_umi": float(g.umi.mean())}
           for (rp, cp), g in sec.groupby([sec.array_row % 2, sec.array_col % 2])]
    col_of = pd.Series(np.arange(len(bc)), index=bc)
    for _, R in rois.iterrows():
        A, roi, upp = R.archetype, R.roi_id, R.um_per_px; tag = f"{A}_{roi}"
        mask = np.load(f"{M}/{'R3b' if SRC[A]=='consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
        nuc = pd.read_parquet(f"{M}/R3b/nuclei/{tag}_consensus_native.parquet" if SRC[A] == "consensus" else f"results/R2/nuclei/{tag}_{SRC[A]}_native.parquet")
        keep = nuc.loc[nuc.keep, "label"].to_numpy()
        box = tp[(tp.pxl_row_in_fullres >= R.r0) & (tp.pxl_row_in_fullres < R.r1) & (tp.pxl_col_in_fullres >= R.c0) & (tp.pxl_col_in_fullres < R.c1)].copy()
        box["rr"] = np.floor(box.pxl_row_in_fullres - R.r0).astype(int); box["cc"] = np.floor(box.pxl_col_in_fullres - R.c0).astype(int)
        filt = box[box.barcode.isin(col_of.index)]; u = umi_all[col_of.reindex(filt.barcode).to_numpy()].astype(float)
        rr, cc = filt.rr.to_numpy(), filt.cc.to_numpy()
        ms = C.shifted_mask(mask, upp)
        t = nuc.loc[nuc.keep, ["label", "x_um", "y_um"]].reset_index(drop=True)
        for mk, mm in (("", mask), ("_shift", ms)):
            ln = C.footprint_labels(mm, rr, cc, upp, d_exp, mode="nuc"); le = C.footprint_labels(mm, rr, cc, upp, d_exp, mode="exp")
            for k, v in per_label(ln, le, keep, u).items(): t[k + mk] = v
        t["archetype"] = A; t["roi_id"] = roi
        t.to_parquet(f"{OUT}/posthoc/{tag}_cells.parquet", index=False)
        f2 = filt.assign(umi=u)
        par += [{"unit": tag, "dataset": ds, "row_par": rp, "col_par": cp, "n_bins": len(g), "mean_umi": float(g.umi.mean())}
                for (rp, cp), g in f2.groupby([f2.array_row % 2, f2.array_col % 2])]
        print(f"[{KEY}] {tag} done", flush=True)
    pd.DataFrame(par).to_csv(f"{OUT}/posthoc/{KEY}_parity.csv", index=False)

elif MODE == "summary":
    OUT, = ARGS; TB = f"{OUT}/tables"; os.makedirs(TB, exist_ok=True)
    rng_for = lambda a, q: np.random.default_rng([C.SEED, ARCH.index(a), 1000 + q])
    # PH-1: rapporto pesato per cellula (media per cellula della densita' nel nucleo / media per cellula della densita' nella corona) e
    # mediana dei rapporti per cellula (cellule con UMI > 0 nella corona); entrambi con la maschera spostata come nullo
    rows, roi_rows = [], []
    for f in sorted(glob.glob(f"{OUT}/posthoc/A?_*_cells.parquet")):
        t = pd.read_parquet(f); A, roi = t.archetype.iloc[0], t.roi_id.iloc[0]
        r = {"archetype": A, "roi_id": roi}
        for mk in ("", "_shift"):
            ok = (t[f"n_nuc{mk}"] > 0) & (t[f"n_ring{mk}"] > 0)
            dn = t.loc[ok, f"umi_nuc{mk}"] / t.loc[ok, f"n_nuc{mk}"]; dr = t.loc[ok, f"umi_ring{mk}"] / t.loc[ok, f"n_ring{mk}"]
            r[f"cellweighted{mk}"] = float(dn.mean() / dr.mean()); r[f"n_cells{mk}"] = int(ok.sum())
            pos = dr > 0
            r[f"median_cell_ratio{mk}"] = float(np.median(dn[pos] / dr[pos])); r[f"n_cells_ratio{mk}"] = int(pos.sum())
            r[f"pooled{mk}"] = float(t[f"umi_nuc{mk}"].sum() / t[f"n_nuc{mk}"].sum() / (t[f"umi_ring{mk}"].sum() / t[f"n_ring{mk}"].sum()))
        r["cellweighted_star"] = r["cellweighted"] / r["cellweighted_shift"]; r["median_cell_ratio_star"] = r["median_cell_ratio"] / r["median_cell_ratio_shift"]
        roi_rows.append(r)
    PR = pd.DataFrame(roi_rows); PR.to_csv(f"{TB}/R4_PH1_roi.csv", index=False)
    for a in ARCH:
        Ra = PR[PR.archetype == a]
        for qi, col in enumerate(("cellweighted", "cellweighted_shift", "cellweighted_star", "median_cell_ratio", "median_cell_ratio_shift",
                                  "median_cell_ratio_star", "pooled", "pooled_shift")):
            e, lo, hi, _ = C.boot_one_stage(Ra[col].values, rng=rng_for(a, qi))
            rows.append({"archetype": a, "quantity": col, "estimate": e, "lo": lo, "hi": hi, "class_vs_0.909_1.1": C.verdict(lo, hi, 0.909, 1.1),
                         "class_vs_0.95_1.05": C.verdict(lo, hi, 0.95, 1.05)})
    pd.DataFrame(rows).to_csv(f"{TB}/R4_PH1.csv", index=False)
    # PH-2 / PH-3 dai file p2 dello stadio per ROI (stessi geni e stesso dropout di P-2/CP-2). Bootstrap sulle sole ROI: le classi NON sono
    # confrontabili con i verdetti pre-registrati di P-2/CP-2 (bootstrap a due stadi); 'non valutabile' con meno di 10 geni, come in P-2.
    out = []
    for a in ARCH:
        tags = sorted(os.path.basename(f)[:-7] for f in glob.glob(f"{OUT}/roi/{a}_*_p2.npz"))
        for kind in ("obs", "drop"):
            for model, gkey, ekey in (("NB with cell size factor, mean >= 2", "genes", "exsf"), ("NB one mean, mean >= 0.5", "genes05", "ex05"),
                                      ("NB with cell size factor, mean >= 0.5", "genes05", "ex05sf")):
                per, genes = [], set()
                for t_ in tags:
                    z = np.load(f"{OUT}/roi/{t_}_p2.npz"); xs = []
                    for k in z.files:
                        ty, kd, wh = k.split("|")
                        if kd == kind and wh == ekey: xs.append(z[k])
                        if kd == kind and wh == gkey: genes.update(z[k].tolist())
                    per.append(np.concatenate(xs) if xs else np.zeros(0))
                n_pairs = int(sum(len(x) for x in per))
                if len(genes) < 10 or n_pairs == 0:
                    out.append({"archetype": a, "kind": kind, "model": model, "n_genes": len(genes), "n_pairs": n_pairs, "estimate": np.nan, "lo": np.nan, "hi": np.nan,
                                "class_vs_pm0.02_roi_bootstrap": "non valutabile"}); continue
                q = 10 + ["obs", "drop"].index(kind) * 5 + ["genes", "genes05"].index(gkey) * 2 + (ekey.endswith("sf"))
                e, lo, hi, _ = C.boot_one_stage(np.arange(len(per)), stat=lambda ix: np.median(np.concatenate([per[int(i)] for i in ix])), rng=rng_for(a, q))
                out.append({"archetype": a, "kind": kind, "model": model, "n_genes": len(genes), "n_pairs": n_pairs, "estimate": e, "lo": lo, "hi": hi,
                            "class_vs_pm0.02_roi_bootstrap": C.verdict(lo, hi, -0.02, 0.02)})
    pd.DataFrame(out).to_csv(f"{TB}/R4_PH2.csv", index=False)
    # PH-4
    par = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{OUT}/posthoc/*_parity.csv"))], ignore_index=True)
    w = []
    for (unit, ds), g in par.groupby(["unit", "dataset"]):
        m = {(int(r.row_par), int(r.col_par)): (r.mean_umi, r.n_bins) for r in g.itertuples()}
        def pooled(keys): return sum(m[k][0] * m[k][1] for k in keys) / sum(m[k][1] for k in keys)
        w.append({"unit": unit, "dataset": ds, "col_odd_over_even": pooled([(0, 1), (1, 1)]) / pooled([(0, 0), (1, 0)]),          # media raggruppata (RV-2)
                  "row_odd_over_even": pooled([(1, 0), (1, 1)]) / pooled([(0, 0), (0, 1)]), "mean_umi_per_bin": pooled(list(m)),
                  **{f"mean_r{k[0]}c{k[1]}": v[0] for k, v in m.items()}, **{f"n_r{k[0]}c{k[1]}": int(v[1]) for k, v in m.items()}})
    pd.DataFrame(w).to_csv(f"{TB}/R4_PH4_parity.csv", index=False)
    print("posthoc summary written")
