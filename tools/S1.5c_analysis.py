"""tools/S1.5c_analysis.py - S1.5c: verdetti di P-1...P-10 per braccio, regole delle tre scelte, controprova CP-1, descrittive (definizioni: docs/steps/S1.5c.md).
Unita' indipendente: ROI (valore del ROI = media dei 10 semi; stima = media dei 5 ROI +- 2 errori standard, sd fra ROI / sqrt(5)); P-1: seme (media delle 4
finestre per seme, poi i 10 semi). Esiti: confermata (stima +- 2 ES dentro l'intervallo), smentita (tutta fuori), indeterminata (a cavallo); frazioni tagliate a [0, 1].
Ingressi: results/S1.5c/S1.5c_panel_units.csv, _panel_types.csv, _unit_metrics.csv, _windows.csv, _real_intervals.csv; results/S1.5/S1.5_real_refs_roi.csv;
results/S1.5b/S1.5b_real_refs_roi_type.csv. Uso: .venv/bin/python tools/S1.5c_analysis.py   Uscite: S1.5c_verdicts.csv, S1.5c_rules.csv, S1.5c_cp1.csv, S1.5c_desc.csv"""
import numpy as np, pandas as pd
RES = "results/S1.5c"
ARMS = ["U", "C-a", "C-b", "N-a", "N-b", "G"]
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
CONTACT = {"A1": (0.366, 0.650), "A2": (0.217, 0.567), "A3": (0.083, 0.283), "A4": (0.350, 0.750), "A5": (0.000, 0.217), "A6": (0.000, 0.070)}  # B-127...B-132
U = pd.read_csv(f"{RES}/S1.5c_panel_units.csv"); T = pd.read_csv(f"{RES}/S1.5c_panel_types.csv")
UM = pd.read_csv(f"{RES}/S1.5c_unit_metrics.csv", keep_default_na=False, na_values=[""]); UM["bin"] = UM["bin"].astype(str).replace("nan", "")
W = pd.read_csv(f"{RES}/S1.5c_windows.csv")
RI = pd.read_csv(f"{RES}/S1.5c_real_intervals.csv", keep_default_na=False, na_values=[""]); RI["bin"] = RI["bin"].astype(str).replace("nan", ""); RI["cell_type"] = RI["cell_type"].fillna("")
RR = pd.read_csv("results/S1.5/S1.5_real_refs_roi.csv")
RT = pd.read_csv("results/S1.5b/S1.5b_real_refs_roi_type.csv"); RT = RT[RT.level == "type"]
assert (U.groupby(["arm", "archetype", "roi_id"]).seed.count() == 10).all() and U.arm.nunique() == 6, "servono 10 semi per ROI e braccio"
def est(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; return v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v)
def verdict(m, se, lo, hi, clip=False):
    a, b = m - 2 * se, m + 2 * se
    if clip: a, b = max(a, 0.0), min(b, 1.0)
    if not np.isfinite(m): return "non valutabile"
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"
rows = []
def add(pid, A, unit, arm, v, lo, hi, clip=False, src=""):
    m, se, n = est(v); rows.append(dict(id=pid, archetype=A, unit=unit, arm=arm, est=m, se=se, n_units=n, lo=lo, hi=hi, verdict=verdict(m, se, lo, hi, clip), source=src))
ps = W.groupby(["arm", "archetype", "seed"]).s1.mean().reset_index()
for arm in ARMS:
    u = U[U.arm == arm].groupby(["archetype", "roi_id"])[["median_area_nuc", "median_nc", "frac_nc095", "frac_trunc", "cv_terr", "cv_terr_nofrag", "median_ecc_T"]].mean().reset_index()
    u = u.merge(RR, on=["archetype", "roi_id"], suffixes=("", "_real"))
    um = UM[UM.arm == arm].groupby(["archetype", "roi_id", "metric", "cell_type", "bin"]).value.mean().reset_index()
    tt = T[T.arm == arm].groupby(["archetype", "roi_id", "cell_type"])[["median_area_nuc", "median_nc"]].mean().reset_index()
    for A in sorted(MAIN):
        s = u[u.archetype == A]; r = RR[RR.archetype == A]
        lo, hi = CONTACT[A]; add("P-1", A, "tessuto", arm, ps[(ps.arm == arm) & (ps.archetype == A)].s1, lo, hi, True, "B-127...B-132")
        add("P-2", A, "tessuto", arm, s.cv_terr / s.cv_terr_real - 1, -0.10, 0.10, False, "B-2 di S1.3")
        add("P-3", A, "tessuto", arm, s.median_nc, r.median_nc.min(), r.median_nc.max(), True, "P-2 di S1.5")
        x = um[(um.archetype == A) & (um.metric == "cv_local4")]; ri = RI[(RI.archetype == A) & (RI.metric == "cv_local4")].iloc[0]
        add("P-7", A, "tessuto", arm, x.value, ri["min"], ri["max"], False, "fra ROI reali")
        add("P-8a", A, "tessuto", arm, s.median_area_nuc, r.median_area_nuc.min(), r.median_area_nuc.max(), False, "P-1 di S1.5b")
        if A == "A2": add("P-8c", A, "tessuto", arm, s.frac_nc095, 0, 0.0287, True, "P-5 di S1.5b")
        add("P-9", A, "tessuto", arm, s.median_ecc_T - s.median_ecc_T_real, -0.05, 0.05, False, "B-4 di S1.3")
        add("P-10", A, "tessuto", arm, s.frac_trunc, 0, 0.05, True, "S1.4")
        for t in MAIN[A]:
            rt = RT[(RT.archetype == A) & (RT.unit == t)]; x = tt[(tt.archetype == A) & (tt.cell_type == t)]
            add("P-4", A, t, arm, x.median_nc, rt.median_nt.min(), rt.median_nt.max(), True, "P-7 di S1.5b")
            add("P-8b", A, t, arm, x.median_area_nuc, rt.median_area_nuc.min(), rt.median_area_nuc.max(), False, "P-2 di S1.5b")
            ri = RI[(RI.archetype == A) & (RI.metric == "nn_rel") & (RI.cell_type == t)]
            if len(ri): add("P-5", A, t, arm, um[(um.archetype == A) & (um.metric == "nn_rel") & (um.cell_type == t)].value, ri["min"].iloc[0], ri["max"].iloc[0], False, "fra ROI reali")
            for b in ["10-20", "40-80"]:
                ri = RI[(RI.archetype == A) & (RI.metric == "markconn") & (RI.cell_type == t) & (RI.bin == b)]
                if len(ri): add("P-6", A, f"{t} {b}", arm, um[(um.archetype == A) & (um.metric == "markconn") & (um.cell_type == t) & (um.bin == b)].value,
                                ri["min"].iloc[0], ri["max"].iloc[0], False, "fra ROI reali")
V = pd.DataFrame(rows); V.to_csv(f"{RES}/S1.5c_verdicts.csv", index=False)
conf = lambda arm, ids: int(((V.arm == arm) & V.id.isin(ids) & (V.verdict == "confermata")).sum())
SPATIAL = ["P-1", "P-2", "P-3", "P-4", "P-5", "P-6", "P-7"]; SIZE = ["P-8a", "P-8b", "P-8c"]
# regola 1: distanza minima (C-a contro N-a): P-1 + P-5; pareggio -> tessuto
r1c, r1n = conf("C-a", ["P-1", "P-5"]), conf("N-a", ["P-1", "P-5"]); dist = "N" if r1n > r1c else "C"
# regola 2: modello di taglia sulla distanza vincente: P-8b; poi P-8a; poi gruppo
ga, gb = f"{dist}-a", f"{dist}-b"; b8a, b8b = conf(ga, ["P-8b"]), conf(gb, ["P-8b"]); a8a, a8b = conf(ga, ["P-8a"]), conf(gb, ["P-8a"])
win = gb if (b8b > b8a or (b8b == b8a and a8b > a8a)) else ga
# regola 3: predefinito se piu' conferme spaziali di U e almeno altrettante nelle taglie
sp_w, sp_u, sz_w, sz_u = conf(win, SPATIAL), conf("U", SPATIAL), conf(win, SIZE), conf("U", SIZE)
adopt = sp_w > sp_u and sz_w >= sz_u
R = pd.DataFrame([dict(rule="1 distanza minima", a="C-a", a_value=r1c, b="N-a", b_value=r1n, counted="P-1 + P-5 confermate", result=dist + " (" + ("per nucleo" if dist == "N" else "tessuto") + ")"),
                  dict(rule="2 modello di taglia", a=ga, a_value=b8a, b=gb, b_value=b8b, counted=f"P-8b confermate (P-8a: {a8a} contro {a8b})", result=win),
                  dict(rule="3 predefinito", a=win, a_value=sp_w, b="U", b_value=sp_u, counted=f"spaziali P-1...P-7 confermate; taglie P-8 {sz_w} contro {sz_u}", result="adottata" if adopt else "non adottata")])
R.to_csv(f"{RES}/S1.5c_rules.csv", index=False)
cp = pd.DataFrame([dict(arm=a, **{p: conf(a, [p]) for p in SPATIAL + SIZE}, spatial=conf(a, SPATIAL), size=conf(a, SIZE)) for a in ARMS])
cp["n_items"] = [int(((V.arm == a) & V.id.isin(SPATIAL)).sum()) for a in ARMS]; cp.to_csv(f"{RES}/S1.5c_cp1.csv", index=False)
# descrittive: F-1 (cuore senza territori con frammenti ricevuti), composizione, distanza / somma dei raggi, tempi
d = []
for arm in ARMS:
    u = U[U.arm == arm]
    for A in sorted(MAIN):
        ua = u[u.archetype == A]; s = ua.groupby("roi_id")[["cv_terr_nofrag"]].mean().reset_index().merge(RR[RR.archetype == A], on="roi_id")
        m, se, n = est(s.cv_terr_nofrag / s.cv_terr - 1); d.append(dict(desc="D-desc-1 CV senza frammenti ricevuti, sim/reale - 1", archetype=A, arm=arm, est=m, se=se))
        d.append(dict(desc="D-desc-4 tempo mediano per unita' (s)", archetype=A, arm=arm, est=ua.elapsed_s.median(), se=np.nan))
        d.append(dict(desc="WARN medi per unita'", archetype=A, arm=arm, est=ua.n_warn.mean(), se=np.nan))
        for mtr in ["gapratio_q05", "gapratio_lt1"]:
            x = UM[(UM.arm == arm) & (UM.archetype == A) & (UM.metric == mtr)].groupby("roi_id").value.mean(); m, se, n = est(x)
            d.append(dict(desc=f"D-desc-3 {mtr}", archetype=A, arm=arm, est=m, se=se))
pd.DataFrame(d).to_csv(f"{RES}/S1.5c_desc.csv", index=False)
pd.set_option("display.width", 250); print(R.to_string(index=False)); print(cp.to_string(index=False))
