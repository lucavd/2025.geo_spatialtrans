"""tools/S1.6_explore_refit.py - S1.6, esplorazione E5 (BL-169, revisione RC-8 e RC-9 di S1.5e), prima della pre-registrazione, nessuna semina nuova.
Rifa' la scelta della famiglia RV di S1.5d (vuoti che spostano le cellule; e' la tabella dei preset e i parametri per piega dei bracci RV e Q) dalle
stesse metriche gia' calcolate (/mnt/micron/geo_spatialtrans/S1.5d/calpos_met/all_units.csv), con la varianza dei bersagli calcolata sulle sole ROI
della piega (le 4 ROI di addestramento; «pooled»: le 5, invariata), e segnala anche i bordi bassi della griglia di void_frac.
Confronta gli ottimi con results/S1.5d/S1.5d_het_RV.csv. Uscita: results/S1.6/S1.6_explore_refit_RV.csv."""
import glob, os, numpy as np, pandas as pd
RES = "results/S1.5d"; MET = "/mnt/micron/geo_spatialtrans/S1.5d/calpos_met"
import importlib.util
s = importlib.util.spec_from_file_location("met", "tools/S1.5d_metrics.py"); met = importlib.util.module_from_spec(s); s.loader.exec_module(met)
real = pd.read_csv(f"{RES}/S1.5d_het_real.csv"); U = pd.read_csv(f"{MET}/all_units.csv")
G = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{RES}/S1.5d_het_grid*.csv"))]).drop_duplicates("combo")
T = met.CAL_TARGETS
S = U.groupby(["combo", "archetype", "roi_id"]).agg({**{t: "mean" for t in T}, "n_failed": "max", "rep": "count"}).reset_index()
old = pd.read_csv(f"{RES}/S1.5d_het_RV.csv")
rows = []
for A, ra in real.groupby("archetype"):
    rlist = sorted(ra.roi_id)
    for fold in [f"loo_{r}" for r in rlist] + ["pooled"]:
        use = rlist if fold == "pooled" else [r for r in rlist if r != fold[4:]]
        for vmode in ["all5", "train"]:
            var = ra[T].var(ddof=1) if vmode == "all5" else ra[ra.roi_id.isin(use)][T].var(ddof=1)
            s_ = S[(S.archetype == A) & S.roi_id.isin(use)].merge(ra[["roi_id"] + T], on="roi_id", suffixes=("", "_real"))
            s_["L"] = sum((s_[t] - s_[t + "_real"]) ** 2 / var[t] for t in T)
            L = s_.groupby("combo").agg(L=("L", "sum"), n_roi=("roi_id", "nunique"), n_rep=("rep", "min"), failed=("n_failed", "max")).reset_index().merge(G, on="combo")
            L = L[(L.n_roi == len(use)) & (L.n_rep == 2) & (L.failed == 0) & (L.dens_sigma == 0)].sort_values(["L", "void_frac", "dens_sigma", "combo"])
            b = L.iloc[0]; pv, lv = b.void_frac, b.void_ell_um
            pos = sorted(set(L.void_frac[L.void_frac > 0]))
            edge = []
            if pv > 0 and pv >= max(L.void_frac): edge.append("void_frac_hi")
            if pv > 0 and pv <= min(pos): edge.append("void_frac_lo_pos")
            if pv > 0 and lv >= L.void_ell_um.max(): edge.append("void_ell_hi")
            if pv > 0 and lv <= L.void_ell_um.min(): edge.append("void_ell_lo")
            rows.append(dict(archetype=A, fold=fold, var=vmode, combo=b.combo, void_frac=pv, void_ell_um=lv if pv > 0 else np.nan, loss=b.L, edge="|".join(edge)))
R = pd.DataFrame(rows); os.makedirs("results/S1.6", exist_ok=True); R.to_csv("results/S1.6/S1.6_explore_refit_RV.csv", index=False)
a5 = R[R["var"] == "all5"].merge(old[["archetype", "fold", "combo"]], on=["archetype", "fold"], suffixes=("", "_s15d"))
print("E5 controllo: con la varianza di S1.5d gli ottimi sono quelli di S1.5d_het_RV.csv:", int((a5.combo == a5.combo_s15d).sum()), "su", len(a5))
w = R.pivot_table(index=["archetype", "fold"], columns="var", values="combo", aggfunc="first").reset_index()
w["changed"] = w["all5"] != w["train"]
print("E5: pieghe lascia-fuori con ottimo diverso usando la varianza delle sole ROI della piega:", int(w[w.fold != "pooled"].changed.sum()), "su", int((w.fold != "pooled").sum()))
print(w[w.changed].to_string(index=False))
print("E5 bordi (varianza delle ROI della piega):"); print(R[(R["var"] == "train") & (R.edge != "")][["archetype", "fold", "combo", "edge"]].to_string(index=False))
