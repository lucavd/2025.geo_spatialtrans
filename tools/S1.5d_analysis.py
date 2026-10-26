"""tools/S1.5d_analysis.py - S1.5d: verdetti P-1...P-14 per braccio, regole A, C e B, controprove CP-1...CP-3, descrittive (definizioni: docs/steps/S1.5d.md).
Unita' indipendente: ROI (valore della ROI = media dei 10 semi; stima = media delle 5 ROI +- 2 errori standard, sd fra ROI / sqrt(5)); P-1: seme (media delle 4
finestre per seme, poi i 10 semi). Esiti: confermata (stima +- 2 ES dentro l'intervallo), smentita (tutta fuori), indeterminata (a cavallo); frazioni tagliate a [0, 1].
«Entro 2 errori standard» nelle controprove: differenze appaiate per ROI, media +- 2 ES contiene 0 (nota D-8).
Uso: .venv/bin/python tools/S1.5d_analysis.py phase1|phase2|final   Uscite: results/S1.5d/S1.5d_verdicts.csv, S1.5d_ruleA.csv, S1.5d_ruleC.csv, S1.5d_ruleB.csv,
S1.5d_cp.csv, S1.5d_desc.csv"""
import sys, numpy as np, pandas as pd
RES = "results/S1.5d"; STAGE = sys.argv[1]
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
CONTACT = {"A1": (0.366, 0.650), "A2": (0.217, 0.567), "A3": (0.083, 0.283), "A4": (0.350, 0.750), "A5": (0.000, 0.217), "A6": (0.000, 0.070)}  # B-127...B-132
rd = lambda f: pd.read_csv(f, keep_default_na=False, na_values=[""])
U = pd.read_csv(f"{RES}/S1.5d_panel_units.csv"); T = pd.read_csv(f"{RES}/S1.5d_panel_types.csv")
UM = rd(f"{RES}/S1.5d_unit_metrics.csv"); UM[["cell_type", "bin"]] = UM[["cell_type", "bin"]].fillna("")
W = pd.read_csv(f"{RES}/S1.5d_windows.csv")
RI = rd("results/S1.5c/S1.5c_real_intervals.csv"); RI[["cell_type", "bin"]] = RI[["cell_type", "bin"]].fillna("")
RD_ = pd.read_csv(f"{RES}/S1.5d_real_intervals.csv"); Q90 = pd.read_csv(f"{RES}/S1.5d_real_q90.csv")
RR = pd.read_csv("results/S1.5/S1.5_real_refs_roi.csv"); RT = pd.read_csv("results/S1.5b/S1.5b_real_refs_roi_type.csv"); RT = RT[RT.level == "type"]
ARMS = [a for a in ["U", "N", "R", "RV", "RD", "RVD", "VF", "W3", "U3", "W3p", "Q", "U3c"] if a in set(U.arm)]
assert (U.groupby(["arm", "archetype", "roi_id"]).seed.count() == 10).all(), "servono 10 semi per ROI e braccio"
for a in ARMS: assert {"cv_local4", "nn_rel", "H_rel_3.0", "xcv2_50", "frac_gt3", "q90_an"} <= set(UM[UM.arm == a].metric), f"metriche mancanti nel braccio {a}"
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
ROI_U = {}
for arm in ARMS:
    u = U[U.arm == arm].groupby(["archetype", "roi_id"])[["median_area_nuc", "median_nc", "frac_nc095", "frac_trunc", "cv_terr", "cv_terr_nofrag", "median_ecc_T"]].mean().reset_index()
    u = u.merge(RR, on=["archetype", "roi_id"], suffixes=("", "_real")); ROI_U[arm] = u
    um = UM[UM.arm == arm].groupby(["archetype", "roi_id", "metric", "cell_type", "bin"]).value.mean().reset_index()
    tt = T[T.arm == arm].groupby(["archetype", "roi_id", "cell_type"])[["median_area_nuc", "median_nc"]].mean().reset_index()
    for A in sorted(MAIN):
        s = u[u.archetype == A]; r = RR[RR.archetype == A]
        lo, hi = CONTACT[A]; add("P-1", A, "tessuto", arm, ps[(ps.arm == arm) & (ps.archetype == A)].s1, lo, hi, True, "B-127...B-132")
        add("P-2", A, "tessuto", arm, s.cv_terr / s.cv_terr_real - 1, -0.10, 0.10, False, "B-2 di S1.3")
        add("P-3", A, "tessuto", arm, s.median_nc, r.median_nc.min(), r.median_nc.max(), True, "P-2 di S1.5")
        x = um[(um.archetype == A) & (um.metric == "cv_local4")]; ri = RI[(RI.archetype == A) & (RI.metric == "cv_local4")].iloc[0]
        add("P-7", A, "tessuto", arm, x.value, ri["min"], ri["max"], False, "fra ROI reali (S1.5c)")
        add("P-8a", A, "tessuto", arm, s.median_area_nuc, r.median_area_nuc.min(), r.median_area_nuc.max(), False, "P-1 di S1.5b")
        if A == "A2": add("P-8c", A, "tessuto", arm, s.frac_nc095, 0, 0.0287, True, "P-5 di S1.5b")
        add("P-9", A, "tessuto", arm, s.median_ecc_T - s.median_ecc_T_real, -0.05, 0.05, False, "B-4 di S1.3")
        add("P-10", A, "tessuto", arm, s.frac_trunc, 0, 0.05, True, "S1.4")
        for pid, mt, clip in [("P-11", "H_rel_3.0", True), ("P-12", "xcv2_50", False), ("P-13", "frac_gt3", True)]:
            ri = RD_[(RD_.archetype == A) & (RD_.metric == mt)].iloc[0]
            add(pid, A, "tessuto", arm, um[(um.archetype == A) & (um.metric == mt)].value, ri["min"], ri["max"], clip, "fra ROI vere (S1.5d)")
        for t in MAIN[A]:
            rt = RT[(RT.archetype == A) & (RT.unit == t)]; x = tt[(tt.archetype == A) & (tt.cell_type == t)]
            add("P-4", A, t, arm, x.median_nc, rt.median_nt.min(), rt.median_nt.max(), True, "P-7 di S1.5b")
            add("P-8b", A, t, arm, x.median_area_nuc, rt.median_area_nuc.min(), rt.median_area_nuc.max(), False, "P-2 di S1.5b")
            q = Q90[(Q90.archetype == A) & (Q90.cell_type == t)]
            add("P-14", A, t, arm, um[(um.archetype == A) & (um.metric == "q90_an") & (um.cell_type == t)].value, q.q90_min.iloc[0], q.q90_max.iloc[0], False, "S1.5d_real_q90.csv")
            ri = RI[(RI.archetype == A) & (RI.metric == "nn_rel") & (RI.cell_type == t)]
            if len(ri): add("P-5", A, t, arm, um[(um.archetype == A) & (um.metric == "nn_rel") & (um.cell_type == t)].value, ri["min"].iloc[0], ri["max"].iloc[0], False, "fra ROI reali (S1.5c)")
            for b in ["10-20", "40-80"]:
                ri = RI[(RI.archetype == A) & (RI.metric == "markconn") & (RI.cell_type == t) & (RI.bin == b)]
                if len(ri): add("P-6", A, f"{t} {b}", arm, um[(um.archetype == A) & (um.metric == "markconn") & (um.cell_type == t) & (um.bin == b)].value,
                                ri["min"].iloc[0], ri["max"].iloc[0], False, "fra ROI reali (S1.5c)")
V = pd.DataFrame(rows); V.to_csv(f"{RES}/S1.5d_verdicts.csv", index=False)
assert (V.n_units > 0).all(), V[V.n_units == 0][["id", "archetype", "unit", "arm"]].drop_duplicates().to_string()
conf = lambda arm, ids: int(((V.arm == arm) & V.id.isin(ids) & (V.verdict == "confermata")).sum())
SPATIAL = ["P-1", "P-2", "P-3", "P-4", "P-5", "P-6", "P-7"]; SIZE = ["P-8a", "P-8b", "P-8c"]; SIZE_C = ["P-8a", "P-8b", "P-8c", "P-10", "P-14"]; RULE_A = ["P-1", "P-2", "P-7", "P-13"]
# regola A: famiglia dell'eterogeneita' (RV, RD, RVD): conferme in P-1, P-2, P-7, P-13; pari -> meno parametri (RV, RD: 2; RVD: 4); RV contro RD pari -> RV
cA = {f: conf(f, RULE_A) for f in ["RV", "RD", "RVD"]}; npar = {"RV": 2, "RD": 2, "RVD": 4}; order = {"RV": 0, "RD": 1, "RVD": 2}
winA = sorted(cA, key=lambda f: (-cA[f], npar[f], order[f]))[0]
pd.DataFrame([dict(winner=winA, RV=cA["RV"], RD=cA["RD"], RVD=cA["RVD"], counted="P-1, P-2, P-7, P-13 confermate (24 verdetti)")]).to_csv(f"{RES}/S1.5d_ruleA.csv", index=False)
print("regola A:", cA, "->", winA)
if STAGE in ("phase2", "final"):
    def ruleC(base, emp):
        zb, ze = conf(base, SIZE_C), conf(emp, SIZE_C); sb, se_ = conf(base, SPATIAL), conf(emp, SPATIAL)
        return dict(base=base, emp=emp, size_base=zb, size_emp=ze, spatial_base=sb, spatial_emp=se_, winner="empirical" if (ze > zb and se_ >= sb) else "gaussian")
    C = pd.DataFrame([ruleC(winA, "W3"), ruleC("U", "U3")]); C.to_csv(f"{RES}/S1.5d_ruleC.csv", index=False); print(C.to_string(index=False))
    if "U3c" in ARMS:                                                  # sensibilita' (revisione RC-03), non sostituisce la regola
        Cs = pd.DataFrame([ruleC("U", "U3c")]); Cs.to_csv(f"{RES}/S1.5d_ruleC_sensitivity_U3c.csv", index=False); print("sensibilita' U3c:"); print(Cs.to_string(index=False))
if STAGE == "final":
    C = pd.read_csv(f"{RES}/S1.5d_ruleC.csv"); C = C[C.emp.isin(["W3", "U3"])].reset_index(drop=True); fin = "W3" if C.winner.iloc[0] == "empirical" else winA; ref = "U3" if C.winner.iloc[1] == "empirical" else "U"
    sp_f, sp_r, sz_f, sz_r = conf(fin, SPATIAL), conf(ref, SPATIAL), conf(fin, SIZE), conf(ref, SIZE)
    B = pd.DataFrame([dict(final=fin, reference=ref, spatial_final=sp_f, spatial_ref=sp_r, size_final=sz_f, size_ref=sz_r, result="adottata" if (sp_f > sp_r and sz_f >= sz_r) else "non adottata")])
    B.to_csv(f"{RES}/S1.5d_ruleB.csv", index=False); print(B.to_string(index=False))
    def paired(a, b, A, col=None, metric=None):
        if metric is None: x = ROI_U[a][ROI_U[a].archetype == A].set_index("roi_id"); y = ROI_U[b][ROI_U[b].archetype == A].set_index("roi_id"); d = (x[col] - y[col]) if col != "p2" else (x.cv_terr / x.cv_terr_real - y.cv_terr / y.cv_terr_real)
        else:
            g = lambda arm: UM[(UM.arm == arm) & (UM.archetype == A) & (UM.metric == metric)].groupby("roi_id").value.mean()
            d = g(a) - g(b)
        m, se, n = est(d); return m, se, abs(m) <= 2 * se
    cp = []
    H = {f: pd.read_csv(f"{RES}/S1.5d_het_{f}.csv") for f in ["RV", "RD", "RVD"]}
    for f in H:
        h = H[f][(H[f].archetype == "A6") & (H[f].fold == "pooled")].iloc[0]
        cp.append(dict(cp="CP-1 cuore, parametri di tutte le ROI", item=f, value=f"void_frac {h.void_frac:g}, dens_sigma {h.dens_sigma:g}", pass_=bool(h.void_frac <= 0.025 and h.dens_sigma <= 0.3), judged=f == winA))
    for nm, col, mt in [("P-2", "p2", None), ("P-7", None, "cv_local4"), ("P-11", None, "H_rel_3.0"), ("P-13", None, "frac_gt3")]:
        m, se, ok = paired(winA, "R", "A6", col, mt); cp.append(dict(cp="CP-1 cuore, vincente - R appaiato", item=nm, value=f"{m:.4g} +- {2*se:.2g}", pass_=bool(ok), judged=True))
    rv11 = conf("RV", ["P-11"]); vf11 = conf("VF", ["P-11"]); rv13 = conf("RV", ["P-13"]); vf13 = conf("VF", ["P-13"])
    a1 = V[(V.arm == "VF") & (V.id == "P-11") & (V.archetype == "A1")].verdict.iloc[0]
    cp += [dict(cp="CP-2 vuoti fini", item="P-11 confermate VF < RV", value=f"{vf11} contro {rv11}", pass_=vf11 < rv11, judged=True),
           dict(cp="CP-2 vuoti fini", item="P-13 confermate VF < RV", value=f"{vf13} contro {rv13}", pass_=vf13 < rv13, judged=True),
           dict(cp="CP-2 vuoti fini", item="P-11 intestino smentita in VF", value=a1, pass_=a1 == "smentita", judged=True)]
    if "W3p" in ARMS:
        p8c = V[(V.arm == "W3p") & (V.id == "P-8c")].verdict.iloc[0]
        cp.append(dict(cp="CP-3 residui senza decimi", item="P-8c' confermata in W3p", value=p8c, pass_=p8c == "confermata", judged=True))
        for A in sorted(MAIN):
            m, se, ok = paired("W3p", "W3", A, "median_area_nuc"); cp.append(dict(cp="CP-3 residui senza decimi", item=f"P-8a {A} W3p - W3", value=f"{m:.4g} +- {2*se:.2g}", pass_=bool(ok), judged=True))
        tw = lambda arm: T[T.arm == arm].groupby(["archetype", "roi_id", "cell_type"]).median_area_nuc.mean()
        dd = (tw("W3p") - tw("W3")).reset_index()
        for A in sorted(MAIN):
            for t in MAIN[A]:
                m, se, n = est(dd[(dd.archetype == A) & (dd.cell_type == t)].median_area_nuc)
                cp.append(dict(cp="CP-3 residui senza decimi", item=f"P-8b {A} {t} W3p - W3", value=f"{m:.4g} +- {2*se:.2g}", pass_=bool(abs(m) <= 2 * se), judged=True))
    pd.DataFrame(cp).rename(columns={"pass_": "pass"}).to_csv(f"{RES}/S1.5d_cp.csv", index=False)
    d = []
    for arm in ARMS:
        u = U[U.arm == arm]
        for A in sorted(MAIN):
            ua = u[u.archetype == A]
            d.append(dict(desc="D-desc-1 rifiuti per quota / proposte accettate dal diradamento", archetype=A, arm=arm, est=(ua.n_quota_reject / ua.n_thin_accept).median(), se=np.nan))
            x = UM[(UM.arm == arm) & (UM.archetype == A) & (UM.metric == "nn_cv")].groupby("roi_id").value.mean(); m, se, n = est(x)
            d.append(dict(desc="D-desc-2 CV della distanza dal vicino", archetype=A, arm=arm, est=m, se=se))
            d.append(dict(desc="D-desc-4 tempo mediano per unita' (s)", archetype=A, arm=arm, est=ua.elapsed_s.median(), se=np.nan))
            for mtr in ["gapratio_q05", "gapratio_lt1"]:
                x = UM[(UM.arm == arm) & (UM.archetype == A) & (UM.metric == mtr)].groupby("roi_id").value.mean(); m, se, n = est(x)
                d.append(dict(desc=f"D-desc-7 (non pre-registrata, come D-desc-3 di S1.5c) {mtr}", archetype=A, arm=arm, est=m, se=se))
    pd.DataFrame(d).to_csv(f"{RES}/S1.5d_desc.csv", index=False)
    cnt = pd.DataFrame([dict(arm=a, **{p: conf(a, [p]) for p in SPATIAL + SIZE_C + ["P-9", "P-11", "P-12", "P-13"]}, spatial=conf(a, SPATIAL), size=conf(a, SIZE), size_C=conf(a, SIZE_C)) for a in ARMS])
    cnt.to_csv(f"{RES}/S1.5d_counts.csv", index=False); pd.set_option("display.width", 250); print(cnt.to_string(index=False))
