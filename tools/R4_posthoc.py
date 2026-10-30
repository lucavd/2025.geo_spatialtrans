"""tools/R4_posthoc.py — R4, analisi post hoc (dichiarate dopo i risultati, non previsioni).
  roi <KEY> <outdir> : PH-1 rapporto nucleo/corona appaiato dentro la cellula (maschera vera e spostata); PH-4 parita' di riga e colonna
                       della griglia a 2 um (ROI e sezione intera, bin in tessuto). -> <outdir>/posthoc/<A>_<roi>_cells.parquet, <KEY>_parity.csv
  p2 <outdir>        : PH-2 zeri in eccesso con binomiale negativa a fattore di taglia per cellula; PH-3 linfonodo con soglia 0.5 UMI per cellula
                       (e dropout di controllo). Dalle matrici per cellula di R6 (identiche a quelle ricalcolate, cancello G-3c).
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

def nb_sf_excess(X, s):
    """PH-2: per gene, P(0) osservata meno P(0) di una binomiale negativa con media s_i * mu e theta dai momenti con fattore di taglia."""
    X = X.tocsc(); n = X.shape[0]
    mu = np.asarray(X.sum(axis=0)).ravel() / s.sum()
    x2 = np.asarray(X.multiply(X).T @ np.ones(n)).ravel(); xs = np.asarray(X.T @ s).ravel()
    resid = x2 - 2 * mu * xs + mu ** 2 * (s ** 2).sum() - mu * s.sum()          # sum (x - s mu)^2 - s mu
    den = mu ** 2 * (s ** 2).sum()
    theta = np.where(resid > 0, den / np.maximum(resid, 1e-300), np.inf)
    nz = X.copy(); nz.data = (nz.data > 0).astype(float); z_obs = 1 - np.asarray(nz.sum(axis=0)).ravel() / n
    lam = np.outer(s, mu)                                                         # cellule x geni
    with np.errstate(over="ignore", invalid="ignore"):
        p0 = np.where(np.isinf(theta)[None, :], np.exp(-lam), np.power(1 + lam / np.where(np.isinf(theta), 1, theta)[None, :], -np.where(np.isinf(theta), 1, theta)[None, :]))
    return z_obs - p0.mean(axis=0)

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

elif MODE == "p2":
    OUT, = ARGS; os.makedirs(f"{OUT}/posthoc", exist_ok=True); rows = []
    for f in sorted(glob.glob(f"{OUT}/roi/*_cells.parquet")):
        ce = pd.read_parquet(f); A, roi = ce.archetype.iloc[0], ce.roi_id.iloc[0]; tag = f"{A}_{roi}"
        Mx = sp.load_npz(f"{M}/R6/counts/{tag}_exp.npz").tocsr()
        ridx = int(pd.read_csv("results/R2/R2_rois_checked.csv").reset_index().query("archetype == @A and roi_id == @roi")["index"].iloc[0])
        rng = np.random.default_rng([C.SEED, ridx, 3])
        for t, g in ce[ce.interior & ~ce.cell_type.isin(list(C.NON_TYPES) + ["missing"])].groupby("cell_type"):
            if len(g) < C.MIN_CELLS: continue
            Xt = Mx[g.index.to_numpy()]
            for thr in ((C.MIN_MEAN, 0.5) if A == "A4" else (C.MIN_MEAN,)):
                for kind, Y in (("obs", Xt), ("drop", C.apply_dropout(Xt, C.DROPOUT, rng))):
                    genes = C.eligible_genes(Y, thr)
                    if not len(genes): continue
                    lib = np.asarray(Y.sum(axis=1)).ravel().astype(float); s = lib / lib.mean()
                    ok = s > 0
                    ex_sf = nb_sf_excess(Y[ok][:, genes], s[ok]); ex_plain = C.excess_zeros(Y, genes)
                    for gi, a_, b_ in zip(genes, ex_plain, ex_sf):
                        rows.append({"archetype": A, "roi_id": roi, "cell_type": t, "min_mean": thr, "kind": kind, "gene_idx": int(gi),
                                     "excess_nb": float(a_), "excess_nb_sizefactor": float(b_)})
        print(f"p2 {tag}", flush=True)
    pd.DataFrame(rows).to_parquet(f"{OUT}/posthoc/p2_sizefactor.parquet", index=False)

elif MODE == "summary":
    OUT, = ARGS; TB = f"{OUT}/tables"; os.makedirs(TB, exist_ok=True)
    rng_for = lambda a, q: np.random.default_rng([C.SEED, ARCH.index(a), 1000 + q])
    # PH-1
    rows, roi_rows = [], []
    for f in sorted(glob.glob(f"{OUT}/posthoc/A?_*_cells.parquet")):
        t = pd.read_parquet(f); A, roi = t.archetype.iloc[0], t.roi_id.iloc[0]
        r = {"archetype": A, "roi_id": roi}
        for mk in ("", "_shift"):
            ok = (t[f"n_nuc{mk}"] > 0) & (t[f"n_ring{mk}"] > 0)
            dn = t.loc[ok, f"umi_nuc{mk}"] / t.loc[ok, f"n_nuc{mk}"]; dr = t.loc[ok, f"umi_ring{mk}"] / t.loc[ok, f"n_ring{mk}"]
            r[f"paired{mk}"] = float(dn.mean() / dr.mean()); r[f"n_cells{mk}"] = int(ok.sum())
            r[f"pooled{mk}"] = float(t[f"umi_nuc{mk}"].sum() / t[f"n_nuc{mk}"].sum() / (t[f"umi_ring{mk}"].sum() / t[f"n_ring{mk}"].sum()))
        r["paired_star"] = r["paired"] / r["paired_shift"]; roi_rows.append(r)
    PR = pd.DataFrame(roi_rows); PR.to_csv(f"{TB}/R4_PH1_roi.csv", index=False)
    for a in ARCH:
        Ra = PR[PR.archetype == a]
        for qi, col in enumerate(("paired", "paired_shift", "paired_star", "pooled", "pooled_shift")):
            e, lo, hi, _ = C.boot_one_stage(Ra[col].values, rng=rng_for(a, qi))
            rows.append({"archetype": a, "quantity": col, "estimate": e, "lo": lo, "hi": hi, "class_vs_0.909_1.1": C.verdict(lo, hi, 0.909, 1.1),
                         "class_vs_0.95_1.05": C.verdict(lo, hi, 0.95, 1.05)})
    pd.DataFrame(rows).to_csv(f"{TB}/R4_PH1.csv", index=False)
    # PH-2 / PH-3
    p2 = pd.read_parquet(f"{OUT}/posthoc/p2_sizefactor.parquet"); out = []
    for (a, thr, kind), g in p2.groupby(["archetype", "min_mean", "kind"]):
        rois_ = sorted(g.roi_id.unique())
        for qi, col in enumerate(("excess_nb", "excess_nb_sizefactor")):
            per = [g[g.roi_id == r_][col].to_numpy() for r_ in rois_]
            e, lo, hi, _ = C.boot_one_stage(np.arange(len(per)), stat=lambda ix: np.median(np.concatenate([per[int(i)] for i in ix])), rng=rng_for(a, 10 + qi + (5 if kind == "drop" else 0) + (20 if thr < 1 else 0)))
            out.append({"archetype": a, "min_mean": thr, "kind": kind, "model": col, "n_genes": int(g.gene_idx.nunique()), "n_pairs": int(len(g)),
                        "estimate": e, "lo": lo, "hi": hi, "class_vs_pm0.02": C.verdict(lo, hi, -0.02, 0.02)})
    pd.DataFrame(out).to_csv(f"{TB}/R4_PH2.csv", index=False)
    # PH-4
    par = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{OUT}/posthoc/*_parity.csv"))], ignore_index=True)
    w = []
    for (unit, ds), g in par.groupby(["unit", "dataset"]):
        m = {(int(r.row_par), int(r.col_par)): r.mean_umi for r in g.itertuples()}
        tot = np.mean(list(m.values()))
        w.append({"unit": unit, "dataset": ds, "col_odd_over_even": (m[(0, 1)] + m[(1, 1)]) / (m[(0, 0)] + m[(1, 0)]),
                  "row_odd_over_even": (m[(1, 0)] + m[(1, 1)]) / (m[(0, 0)] + m[(0, 1)]), "mean_umi": tot, **{f"mean_r{k[0]}c{k[1]}": v for k, v in m.items()}})
    pd.DataFrame(w).to_csv(f"{TB}/R4_PH4_parity.csv", index=False)
    print("posthoc summary written")
