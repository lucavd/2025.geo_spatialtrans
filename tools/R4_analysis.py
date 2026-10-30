"""tools/R4_analysis.py — R4: aggregazione, bootstrap per blocchi sulle ROI, verdetti pre-registrati, tabelle.
Uso: .venv/bin/python tools/R4_analysis.py <outdir>   -> <outdir>/tables/R4_*.csv"""
import sys, os, json, glob
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, "tools")
import R4_common as C

OUT = sys.argv[1]; TB = f"{OUT}/tables"; os.makedirs(TB, exist_ok=True)
M = "/mnt/micron/geo_spatialtrans"
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; AIDX = {a: i for i, a in enumerate(ARCH)}
KEYOF = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
rng_for = lambda a, q: np.random.default_rng([C.SEED, AIDX[a], q])
INF = np.inf
R1LAB = {"confermata": "uniforme", "smentita": "due densita", "indeterminata": "indeterminata (uniforme per difetto)", "non valutabile": "non valutabile"}
def rule1(lo, hi): return R1LAB[C.verdict(lo, hi, 0.909, 1.1)]
P1_RANGE = {a: ((1.1, INF) if a in ("A5", "A6") else (0.909, 1.1)) for a in ARCH}
P4_EXPECT = {"A5": (1.25, INF), "A6": (1.25, INF)}

def q(x, p): return float(np.percentile(x, p)) if len(x) else float("nan")
def cv(x): x = np.asarray(x, float); return float(x.std() / x.mean()) if len(x) and x.mean() > 0 else float("nan")

# ------------------------------------------------------------------ per ROI
rows, cells_all, p2_all = [], [], {}
for f in sorted(glob.glob(f"{OUT}/roi/*_summary.json")):
    S = json.load(open(f)); A, roi = S["archetype"], S["roi_id"]; tag = f"{A}_{roi}"
    b = np.load(f"{OUT}/roi/{tag}_bins.npz"); u, s, ssh = b["umi"], b["s"], b["s_shift"]
    r = {"archetype": A, "roi_id": roi, "n_bins2": len(u), "n_bins8": len(b["umi8"]), "n_nuclei": S["n_nuclei_keep"]}
    for lab, x in (("2", u), ("8", b["umi8"])):
        r.update({f"mean_umi{lab}": float(x.mean()), f"median_umi{lab}": float(np.median(x)), f"q90_umi{lab}": q(x, 90), f"q99_umi{lab}": q(x, 99),
                  f"cv_umi{lab}": cv(x), f"frac_empty{lab}": float((x == 0).mean())})
    r["mean_genes2"] = float(b["ngen"].mean()); r["mean_genes8"] = float(b["ngen8"].mean())
    r["sparsity2"] = 1 - r["mean_genes2"] / S["n_genes"]; r["sparsity8"] = 1 - r["mean_genes8"] / S["n_genes"]
    for code, nm in C.STRATA.items():
        m = s == code
        r[f"share_{nm}"] = float(m.mean()); r[f"mean_umi2_{nm}"] = float(u[m].mean()) if m.any() else float("nan")
        r[f"frac_empty2_{nm}"] = float((u[m] == 0).mean()) if m.any() else float("nan")
    r["rho_obs"] = C.density_ratio(u, s); r["rho_shift"] = C.density_ratio(u, ssh); r["rho_star"] = r["rho_obs"] / r["rho_shift"]
    r["rho4_star"] = C.density_ratio(u, b["s4"]) / C.density_ratio(u, b["s4_shift"])
    r["rho_far_star"] = C.density_ratio(u[b["far"]], s[b["far"]]) / C.density_ratio(u[b["far_shift"]], ssh[b["far_shift"]])
    eok = (b["x_um"] >= 8) & (b["y_um"] >= 8) & (b["x_um"] <= b["side_um"] - 8) & (b["y_um"] <= b["side_um"] - 8)   # S-4 (revisione F15)
    r["rho_edge_star"] = C.density_ratio(u[eok], s[eok]) / C.density_ratio(u[eok], ssh[eok])
    r["rho_sd_star"] = (C.density_ratio(u, b["s_sd"]) / C.density_ratio(u, b["s_sd_shift"])) if "s_sd" in b.files else float("nan")
    r.update({f"spec_{k}": v for k, v in S["spectrum"].items()})
    r["G3c_all_identical"] = all(v["matrix_identical"] and v["table_identical"] for v in S["G3c"].values())
    rows.append(r)
    ce = pd.read_parquet(f"{OUT}/roi/{tag}_cells.parquet"); cells_all.append(ce)
    z = np.load(f"{OUT}/roi/{tag}_p2.npz"); d = {}
    for k in z.files:
        t, kind, what = k.split("|"); d.setdefault(kind, {}).setdefault(t, {})[what] = z[k]
    p2_all[tag] = d
R = pd.DataFrame(rows).sort_values(["archetype", "roi_id"]).reset_index(drop=True)
CE = pd.concat(cells_all, ignore_index=True)
R.to_csv(f"{TB}/R4_roi.csv", index=False)

# ------------------------------------------------------------------ G-3
g3 = []
for f in sorted(glob.glob(f"{OUT}/section/*.json")):
    S = json.load(open(f))
    for sc, v in S["G3a"].items():
        g3.append({"check": "G-3a", "dataset": S["dataset"], "unit": sc, "ours_mean_umi": v["ours"]["mean_umi"], "ref_mean_umi": v["ref"]["mean_umi"],
                   "rel_diff_umi": v["rel_diff_umi"], "ours_mean_genes": v["ours"]["mean_genes"], "ref_mean_genes": v["ref"]["mean_genes"],
                   "rel_diff_genes": v["rel_diff_genes"], "n_bins_ours": v["ours"]["n_bins"], "n_bins_ref": v["ref"]["n_bins"], "pass": v["pass"]})
    gb = S["G3b"]
    g3.append({"check": "G-3b", "dataset": S["dataset"], "unit": "raw 8um = sum of 16 raw 2um", "n_bins_ours": gb["n_bins8"], "n_diff": gb["n_diff"],
               "max_abs_diff": gb["max_abs_diff"], "n_only_one_side_nonzero": gb["n_only_one_side_nonzero"], "pass": gb["pass"]})
for r in rows:
    g3.append({"check": "G-3c", "dataset": f"{r['archetype']}_{r['roi_id']}", "unit": "exp, nuc, shift vs R6", "pass": r["G3c_all_identical"]})
pd.DataFrame(g3).to_csv(f"{TB}/R4_G3.csv", index=False)

# ------------------------------------------------------------------ previsioni
pred = []
def add(pid, a, est, lo, hi, rng_, verdict, **kw):
    pred.append({"id": pid, "archetype": a, "estimate": est, "lo": lo, "hi": hi,
                 "range_lo": rng_[0] if rng_ else np.nan, "range_hi": rng_[1] if rng_ else np.nan, "verdict": verdict, **kw})
refs = []
def ref(a, name, est, lo, hi, method, **kw): refs.append({"archetype": a, "quantity": name, "estimate": est, "lo": lo, "hi": hi, "method": method, **kw})

types_rows, sens_rows, p2_gene_rows = [], [], []
for a in ARCH:
    Ra = R[R.archetype == a]
    # P-1 e CP-1 (uno stadio sulle ROI)
    e, lo, hi, _ = C.boot_one_stage(Ra.rho_star.values, rng=rng_for(a, 1))
    add("P-1", a, e, lo, hi, P1_RANGE[a], C.verdict(lo, hi, *P1_RANGE[a]), rule_R1=rule1(lo, hi), n_roi=len(Ra))
    from scipy import stats as _st
    v = Ra.rho_star.values; se = v.std(ddof=1) / np.sqrt(len(v))
    for nm, mult in (("P-1 interval: estimate +- 2 SE (post hoc, review F09)", 2.0), ("P-1 interval: estimate +- t(4) SE (post hoc, review F09)", float(_st.t.ppf(0.975, len(v) - 1)))):
        sens_rows.append({"archetype": a, "sensitivity": nm, "estimate": float(v.mean()), "lo": float(v.mean() - mult * se), "hi": float(v.mean() + mult * se),
                          "verdict_as_P1": C.verdict(v.mean() - mult * se, v.mean() + mult * se, *P1_RANGE[a]), "rule_R1": rule1(v.mean() - mult * se, v.mean() + mult * se)})
    e, lo, hi, _ = C.boot_one_stage(Ra.rho_shift.values, rng=rng_for(a, 2))
    add("CP-1", a, e, lo, hi, (0.95, 1.05), "superata" if C.verdict(lo, hi, 0.95, 1.05) == "confermata" else "non superata", n_roi=len(Ra))
    e, lo, hi, _ = C.boot_one_stage(Ra.rho_obs.values, rng=rng_for(a, 3))
    ref(a, "rho_obs nucleus/ring UMI per bin (uncorrected)", e, lo, hi, "mean of ROIs, bootstrap ROIs")
    for qi, (col, nm) in enumerate((("rho4_star", "S-1 ring 4 um"), ("rho_far_star", "S-2 bins >= 2 um from nuclear edge"), ("rho_sd_star", "S-3 StarDist-shaped nuclei only"),
                                   ("rho_edge_star", "S-4 bins >= 8 um from the ROI edge (post hoc, review F15)"))):
        v = Ra[col].values
        if np.all(np.isnan(v)): continue
        e, lo, hi, _ = C.boot_one_stage(v, rng=rng_for(a, 10 + qi))
        sens_rows.append({"archetype": a, "sensitivity": nm, "estimate": e, "lo": lo, "hi": hi, "verdict_as_P1": C.verdict(lo, hi, *P1_RANGE[a]),
                          "rule_R1": rule1(lo, hi)})
    # descrittivi per bin (uno stadio)
    for qi, col in enumerate(["mean_umi2", "median_umi2", "q90_umi2", "q99_umi2", "cv_umi2", "frac_empty2", "mean_umi8", "median_umi8", "q90_umi8",
                              "q99_umi8", "cv_umi8", "frac_empty8", "mean_genes2", "mean_genes8", "sparsity2", "sparsity8",
                              "share_nucleus", "share_ring", "share_acellular", "share_other_nucleus", "mean_umi2_nucleus", "mean_umi2_ring",
                              "mean_umi2_acellular", "frac_empty2_nucleus", "frac_empty2_ring", "frac_empty2_acellular"]):
        v = Ra[col].values
        e, lo, hi, _ = C.boot_one_stage(v, rng=rng_for(a, 100 + qi))
        ref(a, col, e, lo, hi, "mean of ROIs, bootstrap ROIs", roi_min=float(np.nanmin(v)), roi_max=float(np.nanmax(v)))
    # cellule interne
    Ca = CE[(CE.archetype == a) & CE.interior]
    groups = [g.umi_exp.to_numpy(float) for _, g in Ca.groupby("roi_id")]
    for qi, (nm, fn) in enumerate((("cell_umi_median", np.median), ("cell_umi_q01", lambda x: np.percentile(x, 1)), ("cell_umi_q99", lambda x: np.percentile(x, 99)),
                                   ("cell_umi_cv", cv))):
        e, lo, hi, _ = C.boot_two_stage(groups, fn, rng=rng_for(a, 200 + qi))
        ref(a, nm, e, lo, hi, "two-stage bootstrap (ROI, cells)", n_cells=int(len(Ca)))
    dens = [(g.umi_exp / g.n_bins_exp.clip(lower=1)).to_numpy(float) for _, g in Ca.groupby("roi_id")]
    e, lo, hi, _ = C.boot_two_stage(dens, np.median, rng=rng_for(a, 210))
    ref(a, "cell_umi_per_bin_median", e, lo, hi, "two-stage bootstrap (ROI, cells)", n_cells=int(len(Ca)))
    ref(a, "cell_umi_max", float(Ca.umi_exp.max()), np.nan, np.nan, "max over interior cells")
    # P-3 (tutte le cellule tenute, anche al bordo)
    Call = CE[CE.archetype == a]
    n_over = int((Call.umi_exp > 50000).sum())
    add("P-3", a, n_over, n_over, n_over, (0, 0), "confermata" if n_over == 0 else "smentita", max_cell_umi=float(Call.umi_exp.max()), n_cells=int(len(Call)))
    # P-4 (due stadi: ROI, poi cellule). Tipi con >= 5 % dei nuclei interni (denominatore: tutti i nuclei interni, come pre-registrato;
    # correzione RV-1/F03 della revisione); sensibilita': quota fra le sole cellule tipizzate (versione del codice prima della revisione).
    Ct = Ca[~Ca.cell_type.isin(list(C.NON_TYPES) + ["missing"])]
    share_all = Ca.cell_type.value_counts(normalize=True); share_typed = Ct.cell_type.value_counts(normalize=True)
    for t in sorted(share_typed.index):
        gt = [g[g.cell_type == t].umi_exp.to_numpy(float) for _, g in Ct.groupby("roi_id")]
        gt = [x for x in gt if len(x)]
        e, lo, hi, _ = C.boot_two_stage(gt, np.median, B=C.B, rng=rng_for(a, 300 + sorted(share_typed.index).index(t)))
        dt = [(g[g.cell_type == t].umi_exp / g[g.cell_type == t].n_bins_exp.clip(lower=1)).to_numpy(float) for _, g in Ct.groupby("roi_id")]
        dt = [x for x in dt if len(x)]
        ed, lod, hid, _ = C.boot_two_stage(dt, np.median, B=C.B, rng=rng_for(a, 400 + sorted(share_typed.index).index(t)))
        types_rows.append({"archetype": a, "cell_type": t, "share": float(share_all[t]), "share_typed": float(share_typed[t]), "n_cells": int((Ct.cell_type == t).sum()),
                           "in_P4": bool(share_all[t] >= C.MIN_TYPE_FRAC), "median_umi_cell": e, "lo": lo, "hi": hi,
                           "median_umi_per_bin": ed, "lo_per_bin": lod, "hi_per_bin": hid})
    def p4(tsel, q):
        per_roi = [(g.umi_exp.to_numpy(float), g.cell_type.to_numpy()) for _, g in Ct[Ct.cell_type.isin(tsel)].groupby("roi_id")]
        def ratio_stat(u, t):
            med = [np.median(u[t == k]) for k in tsel if (t == k).any()]
            return max(med) / min(med) if len(med) == len(tsel) and min(med) > 0 else np.nan
        est = ratio_stat(np.concatenate([x[0] for x in per_roi]), np.concatenate([x[1] for x in per_roi]))
        rg = rng_for(a, q); reps = np.empty(C.B)
        for bb in range(C.B):
            pick = rg.integers(0, len(per_roi), len(per_roi)); us, ts = [], []
            for j in pick:
                uu, tt = per_roi[j]; ii = rg.integers(0, len(uu), len(uu)); us.append(uu[ii]); ts.append(tt[ii])
            reps[bb] = ratio_stat(np.concatenate(us), np.concatenate(ts))
        return est, float(np.nanpercentile(reps, 2.5)), float(np.nanpercentile(reps, 97.5)), int(np.isnan(reps).sum())
    R3 = {"confermata": "sopra 1.25", "smentita": "sotto 1.25", "indeterminata": "a cavallo di 1.25"}
    tsel = sorted(t for t in share_typed.index if share_all[t] >= C.MIN_TYPE_FRAC)
    est4, lo4, hi4, nn = p4(tsel, 4)
    add("P-4", a, est4, lo4, hi4, P4_EXPECT.get(a), C.verdict(lo4, hi4, *P4_EXPECT[a]) if a in P4_EXPECT else "nessuna attesa",
        rule_R3=R3.get(C.verdict(lo4, hi4, 1.25, INF), "non valutabile"), types=";".join(tsel), n_nan_reps=nn)
    tsel_t = sorted(share_typed[share_typed >= C.MIN_TYPE_FRAC].index)
    if tsel_t != tsel:
        e_, l_, h_, _ = p4(tsel_t, 7)
        sens_rows.append({"archetype": a, "sensitivity": "P-4 with the 5 % share among typed cells only (code before review, RV-1/F03)", "estimate": e_, "lo": l_, "hi": h_,
                          "verdict_as_P1": "", "rule_R1": "", "rule_R3": R3.get(C.verdict(l_, h_, 1.25, INF), "non valutabile"), "types": ";".join(tsel_t)})
    # P-2 e CP-2 (due stadi: ROI, poi cellule tramite le repliche salvate)
    tags = sorted(t for t in p2_all if t.startswith(a + "_"))
    for kind, pid in (("obs", "P-2"), ("drop", "CP-2")):
        obs = [np.concatenate([v["ex"] for v in p2_all[t].get(kind, {}).values()] or [np.zeros(0)]) for t in tags]
        reps_t = [[v["reps"] for v in p2_all[t].get(kind, {}).values()] for t in tags]
        genes = set()
        for t in tags:
            for tt, v in p2_all[t].get(kind, {}).items(): genes.update(v["genes"].tolist())
        allobs = np.concatenate(obs) if obs else np.zeros(0)
        if len(genes) < 10 or len(allobs) == 0:
            add(pid, a, np.nan, np.nan, np.nan, (-0.02, 0.02), "non valutabile", n_genes=len(genes), n_pairs=int(len(allobs))); continue
        rg = rng_for(a, 5 if kind == "obs" else 6); rep = np.empty(C.B)
        for bb in range(C.B):
            pick = rg.integers(0, len(tags), len(tags)); parts = []
            for p_, j in enumerate(pick):
                for M_ in reps_t[j]:
                    if M_.shape[1]: parts.append(M_[(bb + 7 * p_) % C.B])
            rep[bb] = np.median(np.concatenate(parts)) if parts else np.nan
        lo2, hi2 = float(np.nanpercentile(rep, 2.5)), float(np.nanpercentile(rep, 97.5))
        if pid == "P-2":
            add(pid, a, float(np.median(allobs)), lo2, hi2, (-0.02, 0.02), C.verdict(lo2, hi2, -0.02, 0.02), n_genes=len(genes), n_pairs=int(len(allobs)),
                frac_pairs_excess_gt_0p05=float((allobs > 0.05).mean()))
        else:
            add(pid, a, float(np.median(allobs)), lo2, hi2, (0.02, INF), "superata" if lo2 > 0.02 else "non superata", n_genes=len(genes), n_pairs=int(len(allobs)))
    gn = pd.read_csv(f"{M}/R6/counts/{KEYOF[a]}_genes.tsv", sep="\t")
    for t in tags:
        for tt, v in p2_all[t].get("obs", {}).items():
            for gi, ex in zip(v["genes"], v["ex"]):
                p2_gene_rows.append({"archetype": a, "roi": t, "cell_type": tt, "gene": gn.gene_name.iloc[int(gi)], "excess_zero": float(ex)})
P = pd.DataFrame(pred); P.to_csv(f"{TB}/R4_predictions.csv", index=False)
pd.DataFrame(refs).to_csv(f"{TB}/R4_refs.csv", index=False)
pd.DataFrame(types_rows).to_csv(f"{TB}/R4_types.csv", index=False)
pd.DataFrame(sens_rows).to_csv(f"{TB}/R4_sensitivity.csv", index=False)
pd.DataFrame(p2_gene_rows).to_csv(f"{TB}/R4_p2_genes.csv", index=False)

# ------------------------------------------------------------------ D-4, D-5, D-6
d4, d5 = [], []
for f in sorted(glob.glob(f"{OUT}/section/*.json")):
    S = json.load(open(f)); key = S["key"]
    for a in [x for x in ARCH if KEYOF[x] == key]:
        Ra = R[R.archetype == a]
        for sc, col in (("square_002um", "mean_umi2"), ("square_008um", "mean_umi8")):
            sec = S["G3a"][sc]["ours"]["mean_umi"]
            d4.append({"archetype": a, "unit": sc, "section_mean_umi": sec, "roi_mean_umi": float(Ra[col].mean()), "ratio_roi_to_section": float(Ra[col].mean() / sec),
                       "roi_min_ratio": float(Ra[col].min() / sec), "roi_max_ratio": float(Ra[col].max() / sec)})
    for where in ("D5_background", "D5_tissue_tiles"):
        d5.append({"dataset": S["dataset"], "where": where, **S[where]})
for r in rows:
    d5.append({"dataset": f"{r['archetype']}_{r['roi_id']}", "where": "ROI acellular bins",
               **{k[5:]: r[k] for k in r if k.startswith("spec_")}})
pd.DataFrame(d4).to_csv(f"{TB}/R4_D4.csv", index=False)
pd.DataFrame(d5).to_csv(f"{TB}/R4_D5.csv", index=False)
# D-6: vecchia pipeline (full_test.R, 06e, 06f) accanto al vero; sparsita' per cellula dalle matrici di R6 (identiche, G-3c)
d6 = []
for a in ARCH:
    Ca = CE[(CE.archetype == a) & CE.interior]; nz = ncell = 0; ng = None; nz5 = 0
    for roi, g in Ca.groupby("roi_id"):
        Mx = sp.load_npz(f"{M}/R6/counts/{a}_{roi}_exp.npz").tocsr()
        allk = CE[(CE.archetype == a) & (CE.roi_id == roi)].reset_index(drop=True)
        rowsel = allk.index[allk.interior].to_numpy(); Mx = Mx[rowsel]
        nz += Mx.nnz; ncell += Mx.shape[0]; ng = Mx.shape[1]
    rr = pd.DataFrame(refs); g = lambda nm: float(rr[(rr.archetype == a) & (rr.quantity == nm)].estimate.iloc[0])
    pp = P[(P.archetype == a)]
    d6.append({"archetype": a, "old_mean_library_per_cell": 8000, "real_median_umi_per_cell": g("cell_umi_median"),
               "old_library_cv": 0.22, "real_cell_umi_cv": g("cell_umi_cv"), "old_dropout_range": "0.45-0.65",
               "real_excess_zero_median_P2": float(pp[pp.id == "P-2"].estimate.iloc[0]),
               "old_sparsity_band_pct": "30-92 (5000 genes)", "real_cell_sparsity_pct_all_genes": 100 * (1 - nz / (ncell * ng)), "n_genes": ng,
               "old_cap_umi": 50000, "real_max_cell_umi": float(pp[pp.id == "P-3"].max_cell_umi.iloc[0])})
pd.DataFrame(d6).to_csv(f"{TB}/R4_D6_old_vs_real.csv", index=False)
print(P[["id", "archetype", "estimate", "lo", "hi", "verdict"]].to_string())
print("tables:", sorted(os.listdir(TB)))
