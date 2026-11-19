"""tools/F3.3_sens_insd.py - F3.3, sensibilita' segnalata dalla revisione: giudizi dentro campione con la deviazione standard delle sole aree di
taratura (invece di quella combinata con gli insiemi di verifica), e i loro effetti su P-3 (celle non smentite), P-4 e P-5.
Uso: .venv/bin/python tools/F3.3_sens_insd.py"""
import sys, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "tools"); import F3_3_stats as F
C = pd.read_csv("results/F3.3/F3.3_cells.csv"); U = pd.read_csv("results/F3.3/F3.3_area_values.csv")
MID = {f"M-{i}": f"M{i}" for i in range(1, 7)}; rows = []
for (A, m), g in C.groupby(["A", "measure"]):
    rin = U[(U.A == A) & (U.src == "IN")][MID[m]].to_numpy(float); sd = rin.std(ddof=1); n = len(rin); t = stats.t.ppf(0.975, n - 1)
    sim_units = U[(U.A == A) & (U.src == "SIM")].groupby("u2")[MID[m]].mean().to_numpy(float); se_s = sim_units.std(ddof=1) / np.sqrt(len(sim_units))
    g0 = g.iloc[0]; w = t * np.sqrt(se_s ** 2 + sd ** 2 / n); v = F.vrule(g0.gap_in, w, g0.h)
    for _, r in g.iterrows(): rows.append({**r.to_dict(), "verdict_in_insd": v})
D = pd.DataFrame(rows)
out = []
for lev in ("L1", "L2"):
    d = D[D.level == lev]
    for col in ("verdict_in", "verdict_in_insd"):
        s = d[d[col] == "smentita"]; k = int(((s.verdict_out == "smentita") & s.same_sign).sum())
        c = d[d[col] == "confermata"]; bad = c[c.verdict_out == "smentita"]; pm = bad.groupby("measure").set.nunique()
        good = d[d[col] != "smentita"]
        qg = (((good.gap_out ** 2 - good.sd_pooled ** 2 / good.n_out) - (good.gap_in ** 2 - good.sd_pooled ** 2 / good.n_in)) / good.h ** 2).mean()
        out.append({"level": lev, "in_verdict_from": col, "P4_k": k, "P4_n": len(s), "P4_share": k / len(s), "P5_max_sections": int(pm.max()) if len(pm) else 0,
                    "P5_cells": "; ".join(bad.A + " " + bad.measure + " " + bad.set), "n_good_cells": good[["A", "measure"]].drop_duplicates().shape[0], "QG_point": F.ssqrt(qg)})
O = pd.DataFrame(out); O.to_csv("results/F3.3/F3.3_sens_insd.csv", index=False); print(O.round(3).to_string(index=False))
print(D[(D.verdict_in != D.verdict_in_insd)][["A", "measure", "verdict_in", "verdict_in_insd"]].drop_duplicates().to_string(index=False))
