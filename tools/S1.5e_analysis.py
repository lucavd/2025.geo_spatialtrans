"""tools/S1.5e_analysis.py - S1.5e: verdetti P-1...P-17 per braccio, regole A, Q e B, controprove CP-1...CP-3, descrittive (definizioni: docs/steps/S1.5e.md).
Unita' indipendente: ROI (valore della ROI = media dei 10 semi; stima = media delle ROI +- 2 errori standard, sd fra ROI / sqrt(n)); P-1: seme (media delle 4
finestre per seme, poi i 10 semi); P-15: i 4 ROI con una finestra manuale. Esiti: confermata (stima +- 2 ES dentro l'intervallo), smentita (tutta fuori),
indeterminata (a cavallo); frazioni tagliate a [0, 1]. «Entro 2 errori standard»: differenze appaiate (per ROI o per seme), media +- 2 ES contiene 0.
Uso: .venv/bin/python tools/S1.5e_analysis.py phase1|final   (phase1: regola A; final: regole Q e B, controprove, descrittive, conteggi)
Uscite: results/S1.5e/S1.5e_verdicts.csv, S1.5e_ruleA.csv, S1.5e_ruleQ.csv, S1.5e_ruleB.csv, S1.5e_cp.csv, S1.5e_desc.csv, S1.5e_counts.csv"""
import sys, numpy as np, pandas as pd
RES = "results/S1.5e"; STAGE = sys.argv[1]
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
CONTACT = {"A1": (0.366, 0.650), "A2": (0.217, 0.567), "A3": (0.083, 0.283), "A4": (0.350, 0.750), "A5": (0.000, 0.217), "A6": (0.000, 0.070)}  # B-127...B-132
rd = lambda f: pd.read_csv(f, keep_default_na=False, na_values=[""])
U = pd.read_csv(f"{RES}/S1.5e_panel_units.csv"); T = pd.read_csv(f"{RES}/S1.5e_panel_types.csv")
UM = rd(f"{RES}/S1.5e_unit_metrics.csv"); UM[["cell_type", "bin"]] = UM[["cell_type", "bin"]].fillna("")
W = pd.read_csv(f"{RES}/S1.5e_windows.csv"); WN = pd.read_csv(f"{RES}/S1.5e_windows_nn.csv")
RI = rd("results/S1.5c/S1.5c_real_intervals.csv"); RI[["cell_type", "bin"]] = RI[["cell_type", "bin"]].fillna("")
RD_ = pd.read_csv("results/S1.5d/S1.5d_real_intervals.csv"); Q90 = pd.read_csv("results/S1.5d/S1.5d_real_q90.csv")
RR = pd.read_csv("results/S1.5/S1.5_real_refs_roi.csv"); RT = pd.read_csv("results/S1.5b/S1.5b_real_refs_roi_type.csv"); RT = RT[RT.level == "type"]
DR = pd.read_csv(f"{RES}/S1.5e_explore_density_ref.csv"); DR["est"] = DR.dens_seg_roi * DR.dens_manual_window / DR.dens_seg_window
VS = pd.read_csv(f"{RES}/S1.5e_real_void_shape.csv")
ARMS = [a for a in ["U", "RV", "X", "XA", "XAD", "Q"] if a in set(U.arm)]
assert (U.groupby(["arm", "archetype", "roi_id"]).seed.count() == 10).all(), "servono 10 semi per ROI e braccio"
for a in ARMS: assert {"cv_local4", "nn_rel", "H_rel_3.0", "xcv2_50", "frac_gt3", "q90_an", "roi_density", "frac_elong_gt3", "n_comp"} <= set(UM[UM.arm == a].metric), f"metriche mancanti nel braccio {a}"
def est(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; return (v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v)) if len(v) > 1 else (np.nan, np.nan, len(v))
def verdict(m, se, lo, hi, clip=False):
    if not np.isfinite(m): return "non valutabile"
    a, b = m - 2 * se, m + 2 * se
    if clip: a, b = max(a, 0.0), min(b, 1.0)
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
        s = u[u.archetype == A]; r = RR[RR.archetype == A]; ua = um[um.archetype == A]
        lo, hi = CONTACT[A]; add("P-1", A, "tessuto", arm, ps[(ps.arm == arm) & (ps.archetype == A)].s1, lo, hi, True, "B-127...B-132")
        add("P-2", A, "tessuto", arm, s.cv_terr / s.cv_terr_real - 1, -0.10, 0.10, False, "B-2 di S1.3")
        add("P-3", A, "tessuto", arm, s.median_nc, r.median_nc.min(), r.median_nc.max(), True, "P-2 di S1.5")
        ri = RI[(RI.archetype == A) & (RI.metric == "cv_local4")].iloc[0]
        add("P-7", A, "tessuto", arm, ua[ua.metric == "cv_local4"].value, ri["min"], ri["max"], False, "fra ROI reali (S1.5c)")
        add("P-8a", A, "tessuto", arm, s.median_area_nuc, r.median_area_nuc.min(), r.median_area_nuc.max(), False, "P-1 di S1.5b")
        if A == "A2": add("P-8c", A, "tessuto", arm, s.frac_nc095, 0, 0.0287, True, "P-5 di S1.5b")
        add("P-9", A, "tessuto", arm, s.median_ecc_T - s.median_ecc_T_real, -0.05, 0.05, False, "B-4 di S1.3")
        add("P-10", A, "tessuto", arm, s.frac_trunc, 0, 0.05, True, "S1.4")
        for pid, mt, clip in [("P-11", "H_rel_3.0", True), ("P-12", "xcv2_50", False), ("P-13", "frac_gt3", True)]:
            ri = RD_[(RD_.archetype == A) & (RD_.metric == mt)].iloc[0]
            add(pid, A, "tessuto", arm, ua[ua.metric == mt].value, ri["min"], ri["max"], clip, "fra ROI vere (S1.5d)")
        d = DR[DR.archetype == A]; x = ua[(ua.metric == "roi_density") & ua.roi_id.isin(d.roi_id)]
        add("P-15", A, "tessuto", arm, x.value, d.est.min(), d.est.max(), False, "fra ROI vere (E1)")
        if A != "A4":
            v = VS[VS.archetype == A]
            add("P-16", A, "tessuto", arm, ua[ua.metric == "frac_elong_gt3"].value, v.frac_elong_gt3.min(), v.frac_elong_gt3.max(), True, "fra ROI vere (E4)")
            add("P-17", A, "tessuto", arm, ua[ua.metric == "n_comp"].value, v.n_comp.min(), v.n_comp.max(), False, "fra ROI vere (E4)")
        for t in MAIN[A]:
            rt = RT[(RT.archetype == A) & (RT.unit == t)]; x = tt[(tt.archetype == A) & (tt.cell_type == t)]
            add("P-4", A, t, arm, x.median_nc, rt.median_nt.min(), rt.median_nt.max(), True, "P-7 di S1.5b")
            add("P-8b", A, t, arm, x.median_area_nuc, rt.median_area_nuc.min(), rt.median_area_nuc.max(), False, "P-2 di S1.5b")
            q = Q90[(Q90.archetype == A) & (Q90.cell_type == t)]
            add("P-14", A, t, arm, ua[(ua.metric == "q90_an") & (ua.cell_type == t)].value, q.q90_min.iloc[0], q.q90_max.iloc[0], False, "S1.5d_real_q90.csv")
            ri = RI[(RI.archetype == A) & (RI.metric == "nn_rel") & (RI.cell_type == t)]
            if len(ri): add("P-5", A, t, arm, ua[(ua.metric == "nn_rel") & (ua.cell_type == t)].value, ri["min"].iloc[0], ri["max"].iloc[0], False, "fra ROI reali (S1.5c)")
            for b in ["10-20", "40-80"]:
                ri = RI[(RI.archetype == A) & (RI.metric == "markconn") & (RI.cell_type == t) & (RI.bin == b)]
                if len(ri): add("P-6", A, f"{t} {b}", arm, ua[(ua.metric == "markconn") & (ua.cell_type == t) & (ua.bin == b)].value,
                                ri["min"].iloc[0], ri["max"].iloc[0], False, "fra ROI reali (S1.5c)")
V = pd.DataFrame(rows); V.to_csv(f"{RES}/S1.5e_verdicts.csv", index=False)
assert (V.n_units > 1).all(), V[V.n_units <= 1][["id", "archetype", "unit", "arm"]].drop_duplicates().to_string()
conf = lambda arm, ids: int(((V.arm == arm) & V.id.isin(ids) & (V.verdict == "confermata")).sum())
SPATIAL = ["P-1", "P-2", "P-3", "P-4", "P-5", "P-6", "P-7"]; SIZE = ["P-8a", "P-8b", "P-8c"]; TYPES = ["P-4", "P-5", "P-6"]
RULE_A = ["P-1", "P-2", "P-7", "P-13", "P-15", "P-16", "P-17"]
CAND = [a for a in ["RV", "X", "XA", "XAD"] if a in ARMS]; npar = {"RV": 2, "X": 2, "XA": 3, "XAD": 5}; order = {"RV": 0, "X": 1, "XA": 2, "XAD": 3}
cA = {f: conf(f, RULE_A) for f in CAND}; winA = sorted(cA, key=lambda f: (-cA[f], npar[f], order[f]))[0]
nA = int((V.arm == CAND[0]).values[V.id.isin(RULE_A).values].sum())
pd.DataFrame([dict(winner=winA, **cA, n_verdicts=nA, counted="P-1, P-2, P-7, P-13, P-15, P-16, P-17 confermate")]).to_csv(f"{RES}/S1.5e_ruleA.csv", index=False)
print("regola A:", cA, "su", nA, "verdetti ->", winA)
if STAGE == "final":
    if "Q" in ARMS:
        sw, sq, tw, tq = conf(winA, SPATIAL), conf("Q", SPATIAL), conf(winA, TYPES), conf("Q", TYPES)
        adopt = sq >= sw and tq >= tw and (sq > sw or tq > tw)
        pd.DataFrame([dict(winner=winA, spatial_quotas=sw, spatial_noquotas=sq, types_quotas=tw, types_noquotas=tq,
                           result="senza quote" if adopt else "quote")]).to_csv(f"{RES}/S1.5e_ruleQ.csv", index=False)
        print(pd.read_csv(f"{RES}/S1.5e_ruleQ.csv").to_string(index=False))
    sp_f, sp_r, sz_f, sz_r = conf(winA, SPATIAL), conf("U", SPATIAL), conf(winA, SIZE), conf("U", SIZE)
    B = pd.DataFrame([dict(final=winA, reference="U", spatial_final=sp_f, spatial_ref=sp_r, size_final=sz_f, size_ref=sz_r,
                           result="adottata" if (sp_f > sp_r and sz_f >= sz_r) else "non adottata")])
    B.to_csv(f"{RES}/S1.5e_ruleB.csv", index=False); print(B.to_string(index=False))
    def paired(a, b, A, col=None, metric=None):
        if metric is None:
            x = ROI_U[a][ROI_U[a].archetype == A].set_index("roi_id"); y = ROI_U[b][ROI_U[b].archetype == A].set_index("roi_id")
            d = (x[col] - y[col]) if col != "p2" else (x.cv_terr / x.cv_terr_real - y.cv_terr / y.cv_terr_real)
        else:
            g = lambda arm: UM[(UM.arm == arm) & (UM.archetype == A) & (UM.metric == metric)].groupby("roi_id").value.mean(); d = g(a) - g(b)
        m, se, n = est(d); return m, se, bool(abs(m) <= 2 * se)
    cp = []
    XA = pd.read_csv(f"{RES}/S1.5e_het_XA.csv"); pool = XA[XA.fold == "pooled"].set_index("archetype")
    h = pool.loc["A6"]; cp.append(dict(cp="CP-1 cuore, parametri di XA su tutte le ROI", item="void_frac <= 0.05 e kappa <= 2", value=f"void_frac {h.void_frac:g}, kappa {h.void_aspect:g}",
                                      pass_=bool(h.void_frac <= 0.05 and h.void_aspect <= 2), judged=True))
    for nm, col, mt in [("P-2", "p2", None), ("P-7", None, "cv_local4"), ("P-11", None, "H_rel_3.0"), ("P-13", None, "frac_gt3")]:
        m, se, ok = paired(winA, "RV", "A6", col, mt); cp.append(dict(cp="CP-1 cuore, vincente - RV appaiato per ROI", item=nm, value=f"{m:.4g} +- {2*se:.2g}", pass_=ok, judged=True))
    for A, rule in [("A2", "<= 2"), ("A3", "<= 2"), ("A6", "<= 2"), ("A5", ">= 4")]:
        k = pool.loc[A].void_aspect; ok = (k <= 2) if rule == "<= 2" else (k >= 4)
        cp.append(dict(cp="CP-2 allungamento scelto dalla calibrazione di XA (tutte le ROI)", item=f"{A} kappa {rule}", value=f"kappa {k:g}", pass_=bool(ok), judged=True))
    s1rv = WN[WN.arm == "RV"].groupby(["archetype", "bin"])[["n", "n_contact"]].sum(); s1rv = (s1rv.n_contact / s1rv.n).rename("s1_rv").reset_index()
    for arm in [a for a in ["X", "XA", "XAD", "RV"] if a in ARMS]:
        x = WN[WN.arm == arm].groupby(["archetype", "seed", "bin"])[["n", "n_contact"]].sum().reset_index().merge(s1rv, on=["archetype", "bin"], how="left")
        x["s1_rv"] = x.s1_rv.fillna(0)
        for A in sorted(MAIN):
            g = x[x.archetype == A].groupby("seed").apply(lambda z: pd.Series(dict(obs=z.n_contact.sum() / z.n.sum(), pred=(z.n * z.s1_rv).sum() / z.n.sum())))
            m, se, n = est(g.obs - g.pred)
            cp.append(dict(cp="CP-3 i contatti seguono le distanze (osservato - previsto da RV per classe, semi)", item=f"{arm} {A}", value=f"{m:.4f} +- {2*se:.2g} (oss. {g.obs.mean():.3f})",
                           pass_=bool(abs(m) <= 2 * se), judged=arm == winA))
    pd.DataFrame(cp).rename(columns={"pass_": "pass"}).to_csv(f"{RES}/S1.5e_cp.csv", index=False)
    d = []
    for arm in ARMS:
        u = U[U.arm == arm]
        for A in sorted(MAIN):
            ua = u[u.archetype == A]
            d.append(dict(desc="D-desc-1 rifiuti per quota / proposte accettate dal diradamento", archetype=A, arm=arm, est=(ua.n_quota_reject / ua.n_thin_accept).median(), se=np.nan))
            d.append(dict(desc="D-desc-1b cellule tolte dai vuoti / cellule del piano", archetype=A, arm=arm, est=(ua.n_removed_void / (ua.n + ua.n_removed_void)).mean(), se=np.nan))
            d.append(dict(desc="D-desc-4 tempo mediano per unita' (s)", archetype=A, arm=arm, est=ua.elapsed_s.median(), se=np.nan))
            z = WN[(WN.arm == arm) & (WN.archetype == A)].groupby("bin")[["n", "n_contact"]].sum()
            for b, r_ in z.iterrows():
                d.append(dict(desc=f"D-desc-3 contatto nella classe di distanza {b}", archetype=A, arm=arm, est=r_.n_contact / r_.n if r_.n else np.nan, se=np.nan))
                d.append(dict(desc=f"D-desc-3 frazione dei focali nella classe {b}", archetype=A, arm=arm, est=r_.n / z.n.sum(), se=np.nan))
    pd.DataFrame(d).to_csv(f"{RES}/S1.5e_desc.csv", index=False)
    ALL = SPATIAL + SIZE + ["P-9", "P-10", "P-11", "P-12", "P-13", "P-14", "P-15", "P-16", "P-17"]
    cnt = pd.DataFrame([dict(arm=a, **{p: conf(a, [p]) for p in ALL}, spatial=conf(a, SPATIAL), size=conf(a, SIZE), types=conf(a, TYPES), ruleA=conf(a, RULE_A)) for a in ARMS])
    cnt.to_csv(f"{RES}/S1.5e_counts.csv", index=False); pd.set_option("display.width", 250); print(cnt.to_string(index=False))
