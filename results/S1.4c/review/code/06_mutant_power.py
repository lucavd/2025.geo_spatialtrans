#!/usr/bin/env python
"""review/code/06_mutant_power.py — potere dei mutanti e dei test di decisione; controlli sul file di bootstrap riportato.
 (a) per ogni mutante, quali asserzioni falliscono (motivo del rilevamento), dai log ufficiali;
 (b) MU5 nel banco: B-7 con il Voronoi dei generatori contro il Voronoi dei centroidi (fold 1, s = 1, K-c), per archetipo;
 (c) mutanti dello strato di decisione costruiti qui (copie in scratch di tools/S1.4c_decision.R) contro la suite
     R/testing/test_S1.4c_decision.R: rilevati o no; effetto sui dati reali delle mutazioni non rilevate (codice proprio);
 (d) S1.4c_bootstrap.csv riportato: P(K > K-c), P(K-2 > K-1), pareggi.
Uso (root del repo): .venv/bin/python results/S1.4c/review/code/06_mutant_power.py
"""
import os, re, subprocess, shutil, json, importlib.util
import numpy as np, pandas as pd
S = "/mnt/micron/geo_spatialtrans/S1.4c"; SC = "/mnt/micron/geo_spatialtrans/S1.4c_review/code"; OUTD = "results/S1.4c/review/code"
res = {}
# (a)
fails = {}
for m, f in [("MU1", "test_R_MU1.log"), ("MU2", "test_R_MU2.log"), ("MU5", "test_R_MU5.log"), ("MU3", "test_py_MU3.log"), ("MU4", "test_py_MU4.log")]:
    L = open(f"{S}/{f}").read().splitlines(); fails[m] = [l for l in L if l.startswith("FAIL")] + [l for l in L if l.startswith("RESULT")]
res["mutant_failures"] = fails
# (b)
mo = pd.read_csv(f"{SC}/exp/mut_MU5_obs.csv"); ob = pd.read_csv("results/S1.4c/S1.4c_obs_raw.csv")
b = ob[(ob.fold == 1) & (ob.seed_idx == 1) & (ob.variant == "K-c")].set_index("archetype")
res["MU5_bench"] = [dict(archetype=r.archetype, b7_none=float(b.loc[r.archetype, "b7"]), b7_MU5=float(r.b7), B7_none=float(b.loc[r.archetype, "B7"]), B7_MU5=float(r.B7))
                    for r in mo.itertuples()]
mg = pd.read_csv(f"{SC}/exp/mut_MU12_checks.csv")
res["MU12_bench"] = mg[mg.variant != "K-c"].groupby(["mutant", "archetype", "variant"])[["c3_viol", "c3_sf_viol", "cno_pairs", "c_area_rel_max", "n_capped"]].max().reset_index().to_dict("records")
# (c) mutanti dello strato di decisione
src = open("tools/S1.4c_decision.R").read()
MUTS = {
 "DM1 bootstrap senza reinserimento": ("sample.int(nseed, nseed, replace = TRUE)", "sample.int(nseed, nseed, replace = FALSE)"),
 "DM2 bootstrap non appaiato (estrazione per variante)": ("sel <- if (is.null(seeds_of)) z[[m]] else z[[m]][match(seeds_of(A, k), z$seed_idx)]",
      "sel <- if (is.null(seeds_of)) z[[m]] else z[[m]][sample.int(10, 10, replace = TRUE)]"),
 "DM3 verdetto di fold sulla media del metro": ("fold_verdict <- function(d, yard) ifelse(d <= stats::median(yard)", "fold_verdict <- function(d, yard) ifelse(d <= mean(yard)"),
 "DM4 verdetto di fold: PASS se d < mediana (stretto)": ("ifelse(d <= stats::median(yard)", "ifelse(d < stats::median(yard)"),
 "DM5 distanza di fold = media dei semi": ("fold_verdict(stats::median(sel)", "fold_verdict(mean(sel)"),
 "DM6 archetipo PASS con >= 3 fold PASS": ('if (sum(v == "PASS") >= 4)', 'if (sum(v == "PASS") >= 3)'),
 "DM7 margine stretto (>)": ("as.numeric(score[cand]) >= as.numeric(score[[\"K-c\"]]) + margin", "as.numeric(score[cand]) > as.numeric(score[[\"K-c\"]]) + margin"),
 "DM8 CP-2 contro la mediana del metro": ("cp2_exceeds <- function(dseg, yards) sum(dseg > vapply(yards, max, 0)) >= 3", "cp2_exceeds <- function(dseg, yards) sum(dseg > vapply(yards, median, 0)) >= 3"),
 "DM9 CP-2 basta 2/5 ROI": ("vapply(yards, max, 0)) >= 3", "vapply(yards, max, 0)) >= 2"),
 "DM10 parita' -> la piu' complessa": ("w <- w[order(-w$score, match(w$variant, order_simple)), ]", "w <- w[order(-w$score, -match(w$variant, order_simple)), ]"),
}
det = []
for name, (a, bb) in MUTS.items():
    assert src.count(a) == 1, name
    d = f"{SC}/dm/{name.split()[0]}"; os.makedirs(d + "/tools", exist_ok=True); os.makedirs(d + "/R/testing", exist_ok=True)
    open(d + "/tools/S1.4c_decision.R", "w").write(src.replace(a, bb)); shutil.copy("R/testing/test_S1.4c_decision.R", d + "/R/testing/")
    p = subprocess.run(["Rscript", "--vanilla", "R/testing/test_S1.4c_decision.R"], cwd=d, capture_output=True, text=True,
                       env=dict(os.environ, R_LIBS_SITE="/nonexistent", R_LIBS_USER="/nonexistent"))
    out = p.stdout.splitlines(); det.append(dict(mutant=name, detected=any(l.startswith("FAIL") for l in out),
                                                failing=[l for l in out if l.startswith("FAIL")], result=[l for l in out if l.startswith("RESULT")]))
res["decision_mutants"] = det
# effetto sui dati reali delle mutazioni non rilevate (DM3, DM4, DM5, DM6) con codice proprio
O3 = importlib.util.spec_from_file_location("x", f"{OUTD}/03_decision_recompute.py")
# ricostruzione locale dei metri e delle distanze (senza rieseguire tutto 03): metro dai file del banco ricalcolati in 03 non salvati;
# si usano le colonne yard_median / yard_max del fold_verdicts riportato, verificate uguali in 03 attraverso i verdetti (126/126)
fv = pd.read_csv("results/S1.4c/S1.4c_fold_verdicts.csv"); fv = fv[fv.metric.isin(["B1", "B2", "B7", "B8"])]
obs = pd.read_csv("results/S1.4c/S1.4c_obs_raw.csv")
def arch(v, need=4): return "PASS" if v.count("PASS") >= need else ("FAIL" if v.count("FAIL") >= 2 else "WARN")
def score(rule):
    s = {}
    for V in ["K-c", "K-1", "K-2"]:
        n = 0
        for A in sorted(fv.archetype.unique()):
            for m in ["B1", "B2", "B7", "B8"]:
                z = fv[(fv.archetype == A) & (fv.variant == V) & (fv.metric == m)].sort_values("fold"); v = []
                for r in z.itertuples():
                    d = obs[(obs.archetype == A) & (obs.variant == V) & (obs.fold == r.fold)][m]
                    dd = d.mean() if rule == "mean_seeds" else d.median()
                    thr = r.yard_median
                    if rule == "strict": v.append("PASS" if dd < thr else ("WARN" if dd <= r.yard_max else "FAIL"))
                    else: v.append("PASS" if dd <= thr else ("WARN" if dd <= r.yard_max else "FAIL"))
                n += arch(v, 3 if rule == "arch3" else 4) == "PASS"
        s[V] = n
    return s
res["undetected_effect_scores"] = {r: score(r) for r in ["base", "mean_seeds", "strict", "arch3"]}
# DM3 (media del metro) richiede i metri completi: rinviato a 03 (scenario 'stretto' gia' calcolato li'); qui solo base/mean/strict/arch3
# (d)
bs = pd.read_csv("results/S1.4c/S1.4c_bootstrap.csv")
res["reported_bootstrap"] = dict(B=len(bs), p_K1_gt_Kc=float(np.mean(bs["K-1"] > bs["K-c"])), p_K2_gt_Kc=float(np.mean(bs["K-2"] > bs["K-c"])),
    p_K2_gt_K1=float(np.mean(bs["K-2"] > bs["K-1"])), p_tie=float(np.mean(bs["K-2"] == bs["K-1"])), p_K1_gt_K2=float(np.mean(bs["K-1"] > bs["K-2"])),
    p_K2_ge_Kc3=float(np.mean(bs["K-2"] >= bs["K-c"] + 3)), p_K1_ge_Kc3=float(np.mean(bs["K-1"] >= bs["K-c"] + 3)))
json.dump(res, open(f"{OUTD}/06_mutant_power.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
