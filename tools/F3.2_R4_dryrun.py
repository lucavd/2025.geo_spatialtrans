#!/usr/bin/env python3
"""tools/F3.2_R4_dryrun.py - F3.2: prova a secco della regola R-4 (obiettivo come intervallo) sui soli dati di taratura di F3.1.
Per tessuto e misura (densita', UMI per quadratino, M-1...M-6 di F3.1): (a) quante delle 5 aree vere cadono nell'intervallo [min, max] delle 5 aree
simulate (versione piu' tollerante del 90 % centrale, con 5 valori); (b) la media simulata cade in [min, max] delle aree vere?
Lato simulato: semina uniforme (predefinito, 'sim') e semina a compartimenti ('simq6'). Nessun dato di verifica.
Uso: python3 tools/F3.2_R4_dryrun.py -> results/F3.2/F3.2_R4_dryrun.csv
"""
import pandas as pd, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/F3.1/F3.1_PH6_areas.csv")
Q = ["density_per_mm2", "umi_per_bin", "M-1", "M-2", "M-3", "M-4", "M-5", "M-6"]
rows = []
for arm in ["sim", "simq6"]:
    for A, g in d.groupby("A"):
        re_ = g[g.side == "real"]; si = g[g.side == arm]
        for q in Q:
            rv, sv = re_[q].dropna().values, si[q].dropna().values
            if len(rv) == 0 or len(sv) == 0: continue
            a_in = int(((rv >= sv.min()) & (rv <= sv.max())).sum())
            b_in = bool(rv.min() <= sv.mean() <= rv.max())
            cv_r = rv.std(ddof=1) / abs(rv.mean()) if len(rv) > 1 else np.nan; cv_s = sv.std(ddof=1) / abs(sv.mean()) if len(sv) > 1 else np.nan
            rows.append(dict(arm=arm, archetype=A, quantity=q, n_real=len(rv), n_sim=len(sv), real_in_sim_range=a_in, sim_mean_in_real_range=b_in,
                             a_ok=a_in == len(rv), cv_real=cv_r, cv_sim=cv_s, sd_ratio_sim_over_real=(sv.std(ddof=1) / rv.std(ddof=1)) if len(rv) > 1 and len(sv) > 1 else np.nan))
o = pd.DataFrame(rows); o.to_csv(ROOT / "results/F3.2/F3.2_R4_dryrun.csv", index=False)
for arm, g in o.groupby("arm"):
    print(arm, "misure x tessuti:", len(g), "| (a) tutte le aree vere dentro il simulato:", int(g.a_ok.sum()), "| aree vere dentro:", int(g.real_in_sim_range.sum()), "/", int(g.n_real.sum()),
          "| (b) media simulata dentro il vero:", int(g.sim_mean_in_real_range.sum()), "| (a) e (b):", int((g.a_ok & g.sim_mean_in_real_range).sum()),
          "| mediana del rapporto delle deviazioni standard simulato/vero:", round(float(g.sd_ratio_sim_over_real.median()), 3))
print(o[o.arm == "sim"].pivot(index="quantity", columns="archetype", values="real_in_sim_range").to_string())
