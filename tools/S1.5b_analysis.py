"""tools/S1.5b_analysis.py - S1.5b: verdetti di P-1...P-8, regola del braccio, controprova CP-1 (definizioni: docs/steps/S1.5b.md).
Unita' indipendente: ROI (valore del ROI = media dei 10 semi); stima = media dei 5 ROI +- 2 errori standard (sd fra ROI / sqrt(5)).
Esiti: confermata (stima +- 2 ES dentro l'intervallo), smentita (tutta fuori), indeterminata (a cavallo).
Uso: .venv/bin/python tools/S1.5b_analysis.py   Uscite: results/S1.5b/S1.5b_verdicts.csv, S1.5b_arm_rule.csv"""
import numpy as np, pandas as pd
RES = "results/S1.5b"
U = pd.read_csv(f"{RES}/S1.5b_panel_units.csv"); T = pd.read_csv(f"{RES}/S1.5b_panel_types.csv")
RR = pd.read_csv("results/S1.5/S1.5_real_refs_roi.csv")            # intervalli di P-1 e P-2 di S1.5
RT = pd.read_csv(f"{RES}/S1.5b_real_refs_roi_type.csv"); RT = RT[RT.level == "type"]
R5 = pd.read_csv(f"{RES}/S1.5b_real_refs_roi.csv")
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
def est(v):
    v = np.asarray(v, float); m = v.mean(); se = v.std(ddof=1) / np.sqrt(len(v)); return m, se, len(v)
def verdict(m, se, lo, hi, clip=False):
    a, b = m - 2 * se, m + 2 * se
    if clip: a, b = max(a, 0.0), min(b, 1.0)
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"
def roi_mean(df, col, by=("archetype", "roi_id")): return df.groupby(list(by))[col].mean().reset_index()
# colonne dei riferimenti reali di S1.5
ra, rn = "median_area_nuc", "median_nc"
rows = []
for arm in ["group", "recentred", "single"]:
    u = U[U.arm == arm]
    for A in sorted(u.archetype.unique()):
        ua = u[u.archetype == A]; ra_ = RR[RR.archetype == A]
        for pid, col, rcol, clip in [("P-1", "median_area_nuc", ra, False), ("P-8", "median_nc", rn, False)]:
            m, se, n = est(roi_mean(ua, col)[col]); lo, hi = ra_[rcol].min(), ra_[rcol].max()
            rows.append(dict(id=pid, archetype=A, unit="tessuto", arm=arm, est=m, se=se, n_roi=n, lo=lo, hi=hi, verdict=verdict(m, se, lo, hi)))
        if A == "A2":
            m, se, n = est(roi_mean(ua, "frac_nc095")["frac_nc095"]); rows.append(dict(id="P-5", archetype=A, unit="tessuto", arm=arm, est=m, se=se, n_roi=n, lo=0, hi=0.0287, verdict=verdict(m, se, 0, 0.0287, True)))
        m, se, n = est(roi_mean(ua, "frac_trunc")["frac_trunc"]); rows.append(dict(id="P-6", archetype=A, unit="tessuto", arm=arm, est=m, se=se, n_roi=n, lo=0, hi=0.05, verdict=verdict(m, se, 0, 0.05, True)))
        if A == "A5":
            r5 = R5[R5.archetype == "A5"]
            lo, hi = max(2.1, r5.ratio_neuron_glia.min()), r5.ratio_neuron_glia.max()
            m, se, n = est(roi_mean(ua, "ratio_neuron_glia")["ratio_neuron_glia"]); rows.append(dict(id="P-3", archetype=A, unit="neuroni/glia", arm=arm, est=m, se=se, n_roi=n, lo=lo, hi=hi, verdict=verdict(m, se, lo, hi)))
            lo, hi = r5.frac_band.min(), r5.frac_band.max()
            m, se, n = est(roi_mean(ua, "frac_band")["frac_band"]); rows.append(dict(id="P-4", archetype=A, unit="banda", arm=arm, est=m, se=se, n_roi=n, lo=lo, hi=hi, verdict=verdict(m, se, lo, hi, True)))
        tt = T[(T.arm == arm) & (T.archetype == A)]
        for t in MAIN[A]:
            x = tt[tt.cell_type == t]; rt = RT[(RT.archetype == A) & (RT.unit == t)]
            for pid, col, rcol in [("P-2", "median_area_nuc", "median_area_nuc"), ("P-7", "median_nc", "median_nt")]:
                v = roi_mean(x, col)[col]; m, se, n = est(v); lo, hi = rt[rcol].min(), rt[rcol].max()
                rows.append(dict(id=pid, archetype=A, unit=t, arm=arm, est=m, se=se, n_roi=n, lo=lo, hi=hi, verdict=verdict(m, se, lo, hi)))
V = pd.DataFrame(rows); V.to_csv(f"{RES}/S1.5b_verdicts.csv", index=False)
cnt = lambda pid, arm: int(((V.id == pid) & (V.arm == arm) & (V.verdict == "confermata")).sum())
rule = pd.DataFrame([dict(arm=a, p2_confirmed=cnt("P-2", a), p2_n=int(((V.id == "P-2") & (V.arm == a)).sum()), p1_confirmed=cnt("P-1", a)) for a in ["group", "recentred", "single"]])
g, r = rule.set_index("arm").loc["group"], rule.set_index("arm").loc["recentred"]
win = "recentred" if r.p2_confirmed > g.p2_confirmed else "group" if g.p2_confirmed > r.p2_confirmed else ("recentred" if r.p1_confirmed > g.p1_confirmed else "group")
step = "P-2" if r.p2_confirmed != g.p2_confirmed else ("P-1" if r.p1_confirmed != g.p1_confirmed else "pareggio -> (a)")
rule["winner"] = win; rule["decided_by"] = step; rule.to_csv(f"{RES}/S1.5b_arm_rule.csv", index=False)
pd.set_option("display.width", 250); print(rule.to_string(index=False))
print(V[V.arm.isin([win, "single"])].round(4).to_string(index=False))
