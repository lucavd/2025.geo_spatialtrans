"""tools/R7_refs_adopted.py — R7, decisione di Luca DEC-061: intervalli dei riferimenti adottati da qui in avanti.
Regole: intervallo corretto (lo_new, hi_new); se il CV fra aree supera 0.45, l'intervallo con il fattore per l'asimmetria (lo_cv, hi_cv;
DEC-060 punto 7); per le finestre contate a mano con due intervalli (R2c, B-127...B-132) il più largo fra quello a due stadi registrato
e quello a uno stadio corretto; nel cervello le componenti dei nuclei con il fattore della struttura misurata, segnate come incerte. Limiti inferiori negativi portati a 0 (quantità non negative).
-> results/R7/Mc/R7_refs_adopted.csv"""
import numpy as np, pandas as pd
M = pd.read_csv("results/R7/Mc/R7_Mc_refs.csv")
M["adopted_lo"] = M.lo_new; M["adopted_hi"] = M.hi_new; M["adopted_rule"] = "corretto"
hc = (M.cv_between_areas > 0.45) & M.lo_cv.notna()
M.loc[hc, "adopted_lo"] = M.loc[hc, "lo_cv"]; M.loc[hc, "adopted_hi"] = M.loc[hc, "hi_cv"]; M.loc[hc, "adopted_rule"] = "corretto con il fattore per l'asimmetria (CV > 0.45)"
w = M.interval_construction.isin(["window_one_stage", "window_two_stage"]) & M.id.isin([f"B-{k}" for k in range(127, 133)])
for i, g in M[w].groupby("id"):
    k = (g.adopted_hi - g.adopted_lo).idxmax()
    M.loc[g.index, "adopted_lo"] = M.loc[k, "adopted_lo"]; M.loc[g.index, "adopted_hi"] = M.loc[k, "adopted_hi"]
    M.loc[g.index, "adopted_rule"] = f"finestre: il più largo ({M.loc[k, 'interval_construction']})"
neg = M.adopted_lo < 0; M["adopted_lo_clipped_at_0"] = neg; M.loc[neg, "adopted_lo"] = 0.0   # tutte le quantità sono non negative
M["uncertain"] = (M.archetype == "A5") & (M.qtype == "nuc") & ~M.interval_construction.isin(["loao", "fixed", "window_one_stage", "window_two_stage"])
M.to_csv("results/R7/Mc/R7_refs_adopted.csv", index=False)
print("limiti inferiori portati a 0:", int(neg.sum())); print(M.adopted_rule.value_counts().to_dict(), "incerte:", int(M.uncertain.sum()), "componenti:", len(M))
print(M[w][["id", "archetype", "interval_construction", "lo", "hi", "lo_new", "hi_new", "adopted_lo", "adopted_hi"]].round(3).to_string(index=False))
