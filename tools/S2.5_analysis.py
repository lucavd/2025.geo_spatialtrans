"""tools/S2.5_analysis.py — S2.5: esiti delle previsioni, regola R-2, controprove e misure (docs/steps/S2.5.md).
Lato vero: <RD>/S2.5_real_areas.csv, S2.5_real_refs.csv (tools/S2.5_real.py); lato simulato: <PD>/S2.5_sim_<variante>.csv (tools/S2.5_stats.py)
e i JSON del pannello (tempi, memoria, blocchi, M). Uso:
  .venv/bin/python tools/S2.5_analysis.py verdicts <RD> <PD> <variante> <outdir>   (main, main_r2, depth_off)
  .venv/bin/python tools/S2.5_analysis.py r2 <RD> <PD> <outdir>                    (regola R-2 dal pannello main)"""
import sys, os, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, "tools")
import S2_5_common as CM
SEED = 20261121
TIS = {"A1": "intestino", "A2": "tumore", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
PRED = [("P-1", "mean_umi2"), ("P-2", "frac_empty2"), ("P-3", "mean_genes2"), ("P-4", "mean_umi8"), ("P-5", "frac_empty8"), ("P-6", "mean_genes8"),
        ("P-7", "star"), ("P-9", "aniso8")]
ADOPT = {"mean_umi2": "mean UMI per 2 um bin", "frac_empty2": "fraction of empty 2 um bins", "mean_genes2": "genes detected per 2 um bin",
         "mean_umi8": "mean UMI per 8 um bin", "frac_empty8": "fraction of empty 8 um bins", "mean_genes8": "genes detected per 8 um bin", "star": "nucleus/ring"}

def sim_table(PD, VAR):
    S = pd.read_csv(f"{PD}/S2.5_sim_{VAR}.csv")
    J = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{PD}/{VAR}/*_k?.json"))])
    return S, J

def verdicts(RD, PD, VAR, OUT):
    os.makedirs(OUT, exist_ok=True)
    RA = pd.read_csv(f"{RD}/S2.5_real_areas.csv"); RF = pd.read_csv(f"{RD}/S2.5_real_refs.csv")
    S, J = sim_table(PD, VAR)
    cols = [f"reg_{m}" for _, m in PRED]
    Sa = S.groupby(["tag", "archetype"])[cols].mean().reset_index()        # media dei semi nell'area
    AD = pd.read_csv("results/R7/Mc/R7_refs_adopted.csv"); AD = AD[AD.is_S25_target.astype(str) == "True"]
    rows, sens = [], []
    for ai, A in enumerate(sorted(TIS)):
        sA = Sa[Sa.archetype == A]; rA = RA[RA.archetype == A]
        n8 = float(rA.reg_n8.mean())
        for qi, (pid, m) in enumerate(PRED):
            r = RF[(RF.archetype == A) & (RF.metric == m) & (RF.definition == "reg")].iloc[0]
            lo_r, hi_r, rule = float(r.lo), float(r.hi), "intervallo corretto (DEC-061)"
            if hi_r == 0: lo_r, hi_r, rule = 0.0, 1.0 / n8, "vero sempre 0: meno di un blocco vuoto per area"
            v = sA[f"reg_{m}"].to_numpy(float)
            ci = CM.corrected_interval(v, A, np.random.default_rng([SEED, 7, ai, qi, ["main", "main_r2", "depth_off"].index(VAR)]), side="simulated")
            rows.append({"id": pid, "archetype": A, "tessuto": TIS[A], "metric": m, "real_est": float(r.est), "rel_lo": lo_r, "rel_hi": hi_r, "rel_rule": rule,
                         "sim_est": ci["est"], "sim_lo": ci["lo"], "sim_hi": ci["hi"], "ratio": ci["est"] / float(r.est) if float(r.est) > 0 else np.nan,
                         "n_areas": ci["n"], "esito": CM.verdict(ci["lo"], ci["hi"], lo_r, hi_r)})
            if m in ADOPT:   # sensibilita': tabella adottata di R7 (definizione sr), nessun esito pre-registrato
                ad = AD[(AD.archetype == A) & AD.metric.str.contains(ADOPT[m], regex=False)]
                if len(ad):
                    a = ad.iloc[0]; lo_a, hi_a = float(a.adopted_lo), float(a.adopted_hi)
                    if hi_a == 0: lo_a, hi_a = 0.0, 1.0 / n8
                    sens.append({"archetype": A, "metric": m, "adopted_id": a.id, "adopted_est": float(a.est), "adopted_lo": lo_a, "adopted_hi": hi_a,
                                 "sim_lo": ci["lo"], "sim_hi": ci["hi"], "esito_sensibilita": CM.verdict(ci["lo"], ci["hi"], lo_a, hi_a)})
            if m == "aniso8":   # sensibilita' (revisione G-5): DEC-060 punto 6, fattore delle 'altre quantita'' (kappa_boot_box) invece di kappa_boot_bin
                cb = CM.corrected_interval(rA["reg_aniso8"].to_numpy(float), A, np.random.default_rng([SEED, ai, qi, 0]), q="box")
                sens.append({"archetype": A, "metric": "aniso8_kappa_box", "adopted_id": f"kappa_box {cb['kappa']:.3f}", "adopted_est": cb["est"], "adopted_lo": cb["lo"], "adopted_hi": cb["hi"],
                             "sim_lo": ci["lo"], "sim_hi": ci["hi"], "esito_sensibilita": CM.verdict(ci["lo"], ci["hi"], cb["lo"], cb["hi"])})
            if m == "star":   # sensibilita': banda di +-10 % di S2.3 sul rapporto simulato / vero
                sens.append({"archetype": A, "metric": "star_ratio_pm10", "adopted_id": "S2.3 P-5", "adopted_est": 1.0, "adopted_lo": 0.909, "adopted_hi": 1.1,
                             "sim_lo": ci["lo"] / float(r.est), "sim_hi": ci["hi"] / float(r.est), "esito_sensibilita": CM.verdict(ci["lo"] / float(r.est), ci["hi"] / float(r.est), 0.909, 1.1)})
        # P-8 tempo e memoria per area (media dei semi), media delle 5 aree +- 2 errori standard (come S2.4 P-7)
        jA = J[J.archetype == A].groupby("tag")[["elapsed_sim_s", "peak_mb"]].mean()
        for q, hi in (("elapsed_sim_s", 300.0), ("peak_mb", 8192.0)):
            v = jA[q].to_numpy(float); e = float(v.mean()); se = float(v.std(ddof=1) / np.sqrt(len(v)))
            rows.append({"id": "P-8", "archetype": A, "tessuto": TIS[A], "metric": q, "real_est": np.nan, "rel_lo": 0.0, "rel_hi": hi, "rel_rule": "bilancio del design, sezione 12",
                         "sim_est": e, "sim_lo": e - 2 * se, "sim_hi": e + 2 * se, "ratio": np.nan, "n_areas": len(v), "esito": CM.verdict(e - 2 * se, e + 2 * se, 0.0, hi),
                         "max_value": float(J[J.archetype == A][q].max())})
    P = pd.DataFrame(rows); P.to_csv(f"{OUT}/S2.5_verdicts_{VAR}.csv", index=False)
    pd.DataFrame(sens).to_csv(f"{OUT}/S2.5_sensitivity_{VAR}.csv", index=False)
    cnt = P.groupby("id").esito.value_counts().unstack(fill_value=0); cnt.to_csv(f"{OUT}/S2.5_verdict_counts_{VAR}.csv")
    # misure
    Jm = J.groupby(["tag", "archetype"]).agg(M=("M", "mean"), chunk_size=("chunk_size", "mean"), n_in_tissue=("n_in_tissue", "first"),
                                             n_tissue_no_territory=("n_tissue_no_territory", "first"), depth_factor=("depth_factor", "first"),
                                             captured_over_target=("n_captured", "sum"), target=("n_target_R4", "sum")).reset_index()
    Jm["captured_over_target"] = Jm.captured_over_target / Jm.target; Jm["share_tissue_no_territory"] = Jm.n_tissue_no_territory / Jm.n_in_tissue
    Jm.drop(columns="target").to_csv(f"{OUT}/S2.5_measures_{VAR}.csv", index=False)
    pd.set_option("display.width", 260); pd.set_option("display.max_rows", 300)
    print(P.round(5).to_string(index=False)); print(cnt); print(pd.DataFrame(sens).round(5).to_string(index=False)); print(Jm.round(5).to_string(index=False))

def rule_r2(RD, PD, OUT):
    RF = pd.read_csv(f"{RD}/S2.5_real_refs.csv"); S, J = sim_table(PD, "main")
    Sa = S.groupby(["tag", "archetype"])["reg_mean_umi2"].mean().reset_index()
    rows = []
    for A in sorted(TIS):
        Ur = float(RF[(RF.archetype == A) & (RF.metric == "mean_umi2") & (RF.definition == "reg")].est.iloc[0])
        Us = float(Sa[Sa.archetype == A].reg_mean_umi2.mean())
        jA = J[J.archetype == A]; f = float(jA.depth_factor.iloc[0]); b = float(jA.b.iloc[0])
        assert jA.depth_factor.nunique() == 1 and jA.b.nunique() == 1
        rows.append({"archetype": A, "U_real": Ur, "U_sim": Us, "b": b, "f_old": f, "f_new": f * (Ur - b) / (Us - b), "ratio": (Ur - b) / (Us - b)})
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/S2.5_R2_depth.csv", index=False); print(R.to_string(index=False))

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "verdicts": verdicts(*sys.argv[2:6])
    elif mode == "r2": rule_r2(*sys.argv[2:5])
