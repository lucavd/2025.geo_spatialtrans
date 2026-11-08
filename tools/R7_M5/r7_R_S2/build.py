import numpy as np, pandas as pd, re, json, sys
R = "/home/user/2025.geo_spatialtrans/results/"
OUT = "/tmp/r7_R_S2/"
INF = np.inf
# ---- kappa table R7_kappa_final_v4.csv (copied verbatim). keys: boot <q>, se_<q>, r<q>
def _k(bb, bn, bx, sb, sn, sx, rb, rn, rx):
    return dict(bin=bb, nuc=bn, box=bx, se_bin=sb, se_nuc=sn, se_box=sx, rbin=rb, rnuc=rn, rbox=rx)
K = {
 "A1": _k(1.856554, 1.854074, 1.837434, 1.586026, 1.583085, 1.569508, 1.032707, 1.032584, 1.029498),
 "A2": _k(1.655263, 1.666445, 1.662134, 1.413347, 1.422352, 1.420388, 1.004158, 1.004919, 1.004163),
 "A3": _k(3.192814, 2.817212, 3.192814, 2.739571, 2.416476, 2.739571, 1.253434, 1.188013, 1.253434),
 "A4": _k(2.570965, 2.965472, 2.965472, 2.206106, 2.549470, 2.549470, 1.147244, 1.216181, 1.216181),
 "A5": _k(1.818286, 3.865405, 3.865405, 1.555139, 3.324402, 3.324402, 1.028976, 1.372253, 1.372253),
 "A6": _k(1.624284, 1.626878, 1.625283, 1.387311, 1.386960, 1.386430, 1.000000, 1.000000, 1.000000),
}
GEOA5 = _k(1.760, 1.751, 1.743, 1.503, 1.496, 1.491, 1.019, 1.020, 1.018)   # geometric A5 factors (sensitivity verdict_geom_A5)
SIM = dict(boot=1.632431, se=1.390531)
WIN = dict(boot=1.979389, se=1.594398)
KWIN = WIN["boot"]
STRESS = 1.343
def verdict(lo, hi, a, b):
    if not (np.isfinite(lo) and np.isfinite(hi)): return "non valutabile"
    if pd.isna(a) and pd.isna(b): return "non valutabile"
    if pd.isna(a): a = -np.inf
    if pd.isna(b): b = np.inf
    if lo >= a and hi <= b: return "confermata"
    if hi < a or lo > b: return "smentita"
    return "indeterminata"
def tissues(arch):
    return [t for t in re.split(r"[;,/ ]+", str(arch)) if t in K]
def kfac(unit, con, side, m, arch, qtype, stress=False, geo=False):
    """con: boot | se (area, window) | two (window two-stage). side: real | simulated."""
    ts = tissues(arch)
    def kk(t, key): return (GEOA5 if (geo and t == "A5") else K[t])[key]
    if unit == "area":
        if side == "simulated":
            return SIM["boot"] if con == "boot" else SIM["se"] * 2 / m
        if not ts: return np.nan
        k = max(kk(t, qtype) if con == "boot" else kk(t, "se_" + qtype) * 2 / m for t in ts)
        return k * STRESS if stress else k
    if unit == "range":
        if not ts: return np.nan
        return max(kk(t, "r" + qtype) for t in ts)
    if unit == "window":
        if con == "two": return 1.0
        k = WIN["boot"] if con == "boot" else WIN["se"] * 2 / m
        return k * STRESS if stress else k
    return 1.0
def label(unit, con, m):
    if unit == "area": return "boot" if con == "boot" else f"se_m (m={m:g})"
    if unit == "window": return "window_two_stage" if con == "two" else ("window_one_stage (boot)" if con == "boot" else f"window_one_stage (se_m, m={m:g})")
    return unit
def scale(lo, hi, c, k):
    if not np.isfinite(k) or k == 1.0: return lo, hi
    nlo = c - (c - lo) * k if np.isfinite(lo) else lo
    nhi = c + (hi - c) * k if np.isfinite(hi) else hi
    return nlo, nhi
ROWS = []
def f(x):
    try: return float(x)
    except Exception: return np.nan
CORR = ("area", "range", "window")
def add(step, pid, arch, qty, est, lo, hi, a, b, eu, mu, qt, rec, rec_map, touches="none", notes="",
        mcenter=None, a_new=None, b_new=None, a_new_s=None, b_new_s=None, clip=None, a_new_geo=None,
        ec="boot", es="real", em=2, mc="boot", ms="real", mm=2, force_nv=False, alt_sim=False):
    est, lo, hi, a, b = map(f, (est, lo, hi, a, b))
    old = "non valutabile" if force_nv else verdict(lo, hi, a, b)
    reproduces = "NA" if rec_map is None else str(old == rec_map).upper()
    out = {}
    for tag, stress, geo in (("n", False, False), ("s", True, False), ("g", False, True)):
        ke = kfac(eu, ec, es, em, arch, qt, stress, geo)
        ce = est if np.isfinite(est) else (lo + hi) / 2
        if eu == "range": ce = (lo + hi) / 2
        nlo, nhi = scale(lo, hi, ce, ke) if eu in CORR else (lo, hi)
        km = kfac(mu, mc, ms, mm, arch, qt, stress, geo)
        if a_new is not None:
            na, nb = (a_new, b_new) if not stress else (a_new_s, b_new_s)
            if geo and a_new_geo is not None: na = a_new_geo
        elif mu in CORR:
            cm = mcenter if mcenter is not None else (a + b) / 2
            na, nb = scale(a, b, cm, km)
        else:
            na, nb = a, b
        if clip is not None:
            nlo, nhi, na, nb = [min(max(v, clip[0]), clip[1]) if np.isfinite(v) else v for v in (nlo, nhi, na, nb)]
        out[tag] = (nlo, nhi, na, nb, ke if eu in CORR else 1.0, km if mu in CORR else 1.0)
    nlo, nhi, na, nb, ke, km = out["n"]
    vnew = "non valutabile" if force_nv else verdict(nlo, nhi, na, nb)
    s_ = out["s"]; vst = "non valutabile" if force_nv else verdict(*s_[:4])
    g_ = out["g"]; vg = (("non valutabile" if force_nv else verdict(*g_[:4])) if tissues(arch) == ["A5"] and len(str(arch)) == 2 else "")
    corr_e = eu in CORR and ke != 1.0; corr_m = mu in CORR and km != 1.0
    affected = corr_e or corr_m
    trans = "unaffected" if not affected else ("unchanged" if vnew == old else f"{old}->{vnew}")
    note = notes
    if len(tissues(arch)) > 1 and affected:
        note = (note + "; " if note else "") + f"panel-wide row ({arch}): largest kappa among tissues used"
    if eu == "loao" or mu == "loao":
        note = (note + "; " if note else "") + "LOAO component not corrected: held-out area shares up to 60% of nuclei with training areas in A3/A4"
    if (eu == "window" and ec == "two") or (mu == "window" and mc == "two"):
        note = (note + "; " if note else "") + "two-stage window interval kept unchanged (rule d: two-stage is 1.1-4x wider than one-stage; not widened) - flagged"
    if alt_sim and eu == "area":
        k_alt = kfac(eu, ec, "simulated", em, arch, qt)
        l_alt, h_alt = scale(lo, hi, est if np.isfinite(est) else (lo + hi) / 2, k_alt)
        if clip is not None: l_alt, h_alt = max(l_alt, clip[0]), min(h_alt, clip[1])
        note = (note + "; " if note else "") + f"paired sim-vs-real per area on overlapping real areas: side taken as real (tissue factor); with sim_5 factor ({k_alt:.3f}) the verdict would be {verdict(l_alt, h_alt, na, nb)}"
    if (eu == "area" and es == "simulated"):
        note = (note + "; " if note else "") + "estimate interval across simulated areas: sim_5 factor (not stressed)"
    ROWS.append(dict(step=step, prediction_id=pid, archetype=arch, quantity=qty, est=est, lo=lo, hi=hi, a=a, b=b,
        est_unit=eu, margin_unit=mu, qtype=qt, verdict_recorded=rec, reproduces=reproduces,
        lo_new=nlo, hi_new=nhi, a_new=na, b_new=nb, kappa_est=ke, kappa_margin=km,
        verdict_new=vnew, verdict_new_stress=vst, verdict_geom_A5=vg, transition=trans, touches_conclusion=touches, notes=note,
        interval_construction=label(eu, ec, em), side=(es if eu == "area" else ("real" if eu in ("range", "window") else "")),
        margin_construction=label(mu, mc, mm), margin_side=(ms if mu == "area" else ("real" if mu in ("range", "window") else "")),
        _old=old))
MAPV = {"confermata": "confermata", "smentita": "smentita", "indeterminata": "indeterminata",
        "non valutabile": "non valutabile", "PASS": "confermata", "FAIL": "smentita",
        "superata": "confermata", "smentita per eccesso": "smentita", "smentita per difetto": "smentita",
        "non calcolabile come scritto": "non valutabile", "non calcolabile": "non valutabile"}
# ======================= R2b (protocol v1, check B/CP on point estimates) =======================
def r2b():
    S = pd.read_csv(R + "R2b/R2b_archetype_summary.csv"); S = S[S.rater == "Luca"].set_index("archetype")
    T = pd.read_csv(R + "R2b/R2b_test_results.csv").set_index("id")
    st = "R2b"; NOTE = "protocol v1 check on a point estimate (no interval used in the verdict); window-level t interval (BL-041) not used"
    tc_dens = "density reference = interval [evident, total] of manual count (B-043...B-048); segmenter choice for R3 (ROADMAP s6, 2026-10-07)"
    det = T.loc["B-R2b.1", "detail"]; fail1 = re.findall(r"A\d", det.split(":")[-1])
    for A in S.index:
        cp, sd = S.loc[A, "density_cellpose_rgb_windows"], S.loc[A, "density_stardist_he_windows"]; m = S.loc[A, "density_manual"]
        rec = "FAIL" if A in fail1 else "PASS"
        add(st, "B-R2b.1", A, "manual nuclear density (/mm2) within [0.9*min(CP,SD), 1.1*max(CP,SD)] on same windows", m, m, m,
            0.9 * min(cp, sd), 1.1 * max(cp, sd), "fixed", "fixed", "nuc", f"{rec} (panel {T.loc['B-R2b.1','result']} {det})", MAPV[rec], tc_dens, NOTE)
    det = T.loc["B-R2b.2", "detail"]; pct = dict(re.findall(r"(A\d) ([+-]?\d+)%", det))
    for A in S.index:
        thr = 0.15 if A in ("A1", "A2", "A3", "A4") else 0.20
        e = S.loc[A, "consensus_R2_5roi"] / S.loc[A, "density_manual"] - 1
        rec = "PASS" if abs(int(pct[A])) / 100 <= thr else "FAIL"
        add(st, "B-R2b.2", A, "consensus R2 (5 ROI) / manual density - 1", e, e, e, -thr, thr, "fixed", "fixed", "nuc",
            f"{rec} (panel {T.loc['B-R2b.2','result']}; recorded {pct[A]}%)", MAPV[rec], tc_dens, NOTE + f"; recorded rounded value {pct[A]}% vs computed {100*e:.1f}%")
    det = T.loc["B-R2b.3a", "detail"]; ok = dict(re.findall(r"(A\d):(ok|NO)", det))
    for A in S.index:
        rc, rs = S.loc[A, "recall_cellpose_rgb"], S.loc[A, "recall_stardist_he"]
        if A in ("A1", "A2", "A3"): e, a, b, q = rc - rs, 0, INF, "recall CP - recall SD (>0)"
        elif A in ("A5", "A6"): e, a, b, q = rs - rc, 0, INF, "recall SD - recall CP (>0)"
        else: e, a, b, q = min(rc, rs), 0.85, 1, "min(recall CP, recall SD) >= 0.85"
        rec = "PASS" if ok[A] == "ok" else "FAIL"
        add(st, "B-R2b.3a", A, q, e, e, e, a, b, "fixed", "fixed", "nuc", f"{rec} (panel FAIL; {det})", MAPV[rec], "segmenter choice per tissue (R2b -> R3, D-R3.1)", NOTE)
    for A, thr in (("A5", 0.85), ("A6", 0.70)):
        p = S.loc[A, "precision_cellpose_rgb"]
        add(st, "B-R2b.3b", A, "Cellpose precision (< thr expected)", p, p, p, 0, thr, "fixed", "fixed", "nuc",
            f"{T.loc['B-R2b.3b','result']} ({T.loc['B-R2b.3b','detail']})", "smentita", "BL-032/033 reformulated", NOTE)
    order = ["A4", "A1", "A2", "A3", "A5", "A6"]; dens = S.density_manual
    n_ok = sum(dens[order[i]] > dens[order[i + 1]] for i in range(5))
    add(st, "B-R2b.4", "A1;A2;A3;A4;A5;A6", "ordering of manual densities A4>A1>A2>A3>A5>A6 (adjacent pairs in order, of 5)", n_ok, n_ok, n_ok, 5, 5,
        "fixed", "fixed", "nuc", T.loc["B-R2b.4", "result"], MAPV[T.loc["B-R2b.4", "result"]], "none", NOTE)
    add(st, "B-R2b.5", "A2;A3;A4;A5", "nuclear area from manual contours", np.nan, np.nan, np.nan, np.nan, np.nan, "other", "fixed", "nuc",
        "N/A (not executed)", "non valutabile", "none", "phase 2 not run (BL-025)")
    IR = pd.read_csv(R + "R2b/R2b_interrater.csv"); IR = IR[(IR.rater_a == "Federica") & (IR.rater_b == "Luca")]
    for A, g in IR.groupby("archetype"):
        e = abs(g.n_a.sum() - g.n_b.sum()) / g.n_b.sum()
        add(st, "CP-R2b.1a", A, "|count Federica - count Luca| / count Luca (2 windows)", e, e, e, 0, 0.10, "fixed", "fixed", "nuc",
            f"{T.loc['CP-R2b.1a','result']} (panel; {T.loc['CP-R2b.1a','detail']})", None, "none", NOTE + "; recorded per-archetype outcome only as '5/6 (max 13% A4)': reproduces NA")
    e = abs(IR.n_a.sum() - IR.n_b.sum()) / IR.n_b.sum()
    add(st, "CP-R2b.1a", "A1;A2;A3;A4;A5;A6", "|count F - count L| / count L, total over 12 windows", e, e, e, 0, 0.10, "fixed", "fixed", "nuc",
        f"{T.loc['CP-R2b.1a','result']} (PASS on total)", "confermata", "none", NOTE)
    f1 = 2 * IR.tp_3um.sum() / (IR.n_a.sum() + IR.n_b.sum())
    add(st, "CP-R2b.1b", "A1;A2;A3;A4;A5;A6", "point-to-point F1 (3 um) Luca vs Federica, pooled", f1, f1, f1, 0.90, 1, "fixed", "fixed", "nuc",
        T.loc["CP-R2b.1b", "result"], MAPV[T.loc["CP-R2b.1b", "result"]], "none", NOTE)
    for A in S.index:
        e = S.loc[A, "common_fn_frac_cp_sd"]
        add(st, "CP-R2b.2", A, "manual points missed by both CP and SD", e, e, e, 0, 0.05, "fixed", "fixed", "nuc",
            T.loc["CP-R2b.2", "result"], "smentita", "density reference as interval (BL-035)", NOTE)
    for A, thr in (("A5", 0.15), ("A6", 0.30)):
        e = S.loc[A, "ghost_cellpose_rgb"]
        add(st, "CP-R2b.3", A, "Cellpose objects without manual point (ghosts)", e, e, e, thr, 1, "fixed", "fixed", "nuc",
            T.loc["CP-R2b.3", "result"], "smentita", "BL-032/033", NOTE)
    NU = pd.read_csv(R + "R2b/R2b_null_summary.csv"); NU = NU[NU.rater == "Luca"]
    for A, g in NU.groupby("archetype"):
        e = g.f1_gain.min()
        add(st, "CP-R2b.4", A, "min over methods of F1(Luca) - F1(null)", e, e, e, 0.20, INF, "fixed", "fixed", "nuc",
            T.loc["CP-R2b.4", "result"], "confermata", "none", NOTE)
    OD = pd.read_csv(R + "R2b/R2b_od_strata.csv").set_index("archetype")
    for A in OD.index:
        e = OD.loc[A, "ratio_missed_over_matched"]
        add(st, "CP-R2b.5", A, "OD missed / OD matched (Luca), <= 0.75", e, e, e, 0, 0.75, "fixed", "fixed", "nuc",
            T.loc["CP-R2b.5", "result"], "confermata", "reference redefined as interval [evident, total] (Luca, 2026-09-21)", NOTE)
    for vid in ("VLM-1", "VLM-2"):
        add(st, vid, "A1;A2;A3;A4;A5;A6", T.loc[vid, "check"], np.nan, np.nan, np.nan, np.nan, np.nan, "fixed", "fixed", "nuc",
            f"{T.loc[vid,'result']} ({T.loc[vid,'detail']})", None, "none", "VLM raters, descriptive extra; no interval")
# ======================= R2c =======================
def bdec(s1, lo, hi, cp, sd):
    inc = lo <= cp <= hi; ins = lo <= sd <= hi
    if inc and ins:
        dc, ds = abs(cp - s1), abs(sd - s1)
        if abs(dc - ds) < 0.02: return "indeciso"
        return "cellpose" if dc < ds else "StarDist"
    if inc: return "cellpose"
    if ins: return "StarDist"
    return "nessuno"
def r2c():
    st = "R2c"
    S = pd.read_csv(R + "R2c/R2c_s1.csv").set_index("archetype")
    P = pd.read_csv(R + "R2c/R2c_predictions.csv").set_index("prediction")
    TC = "S-1 reference for the simulator = manual number with honest margin (Luca, 2026-10-20, ROADMAP s6; B-127...B-132)"
    yes = dict(re.findall(r"(A\d) (sì|no)", P.loc["P-1", "detail"]))
    for A in S.index:
        s1 = S.loc[A, "s1_primary"]; a, b = (0.30, 1) if A in ("A1", "A2", "A3", "A4") else (0, 0.15)
        rec = "confermata" if yes[A] == "sì" else "smentita"
        add(st, "P-1", A, "manual S-1 (fraction of focal nuclei touching)", s1, s1, s1, a, b, "fixed", "fixed", "nuc",
            f"{rec} per tissue ({yes[A]}); panel {P.loc['P-1','verdict']} {P.loc['P-1','n_ok']}/{P.loc['P-1','n']}", rec, "none",
            "point estimate vs fixed S1.4c threshold; manual window CI not used in this prediction")
    yes2 = dict(re.findall(r"(A\d) (sì|no)", P.loc["P-2", "detail"]))
    for A in ("A1", "A2", "A3", "A4"):
        s1 = S.loc[A, "s1_primary"]; cp, sd = S.loc[A, "seg_cellpose_rgb"], S.loc[A, "seg_stardist_he"]
        rec = "confermata" if yes2[A] == "sì" else "smentita"
        add(st, "P-2", A, "manual S-1 strictly between StarDist and Cellpose S-1", s1, s1, s1, min(cp, sd), max(cp, sd), "fixed", "fixed", "nuc",
            f"{rec} per tissue; panel {P.loc['P-2','verdict']} {P.loc['P-2','n_ok']}/{P.loc['P-2','n']}", rec, "none", "point estimates")
    for A in S.index:
        add(st, "P-3", A, "sensitivity/specificity CP vs SD on contact pairs", np.nan, np.nan, np.nan, np.nan, np.nan, "fixed", "fixed", "nuc",
            P.loc["P-3", "verdict"], "non valutabile", "none", "power: 1-12 touching pairs per tissue (BL-132)")
    newdec = {}
    for A in S.index:
        s1, lo, hi = S.loc[A, ["s1_primary", "lo_primary", "hi_primary"]]
        cp, sd = S.loc[A, "seg_cellpose_rgb"], S.loc[A, "seg_stardist_he"]
        rb = S.loc[A, "b_dec"]
        nlo, nhi = scale(lo, hi, s1, KWIN); nlo, nhi = max(nlo, 0), min(nhi, 1)
        slo, shi = scale(lo, hi, s1, KWIN * STRESS); slo, shi = max(slo, 0), min(shi, 1)
        nd = bdec(s1, nlo, nhi, cp, sd); sdd = bdec(s1, slo, shi, cp, sd); newdec[A] = nd
        rep = bdec(s1, lo, hi, cp, sd)
        for seg, v in (("cellpose", cp), ("StarDist", sd)):
            recm = {"nessuno": "smentita", "indeciso": "confermata"}.get(rb, "confermata" if rb == seg else None)
            add(st, f"B-dec ({seg})", A, f"{seg} S-1 inside manual 95% CI (cluster bootstrap over windows)", v, v, v, lo, hi, "fixed", "window", "nuc",
                f"b_dec={rb}", recm, TC,
                f"segmenter S-1 is a point; manual one-stage cluster-percentile CI over 4 windows (pre-registered) widened by window_4 kappa_boot around s1={s1:.4f}, clipped to [0,1]; "
                f"b_dec recomputed old={rep}, new={nd}, stress={sdd}", mcenter=s1, clip=(0, 1))
    for A in S.index:
        s1, lo, hi, k2, w, rk = S.loc[A, ["s1_primary", "lo_primary", "hi_primary", "k2_s1", "widen", "k2_verdict"]]
        def cls(L, H):
            if L <= k2 <= H: return "PASS"
            if L - w <= k2 <= H + w: return "WARN"
            return "FAIL"
        nlo, nhi = scale(lo, hi, s1, KWIN); nlo, nhi = max(nlo, 0), min(nhi, 1)
        slo, shi = scale(lo, hi, s1, KWIN * STRESS); slo, shi = max(slo, 0), min(shi, 1)
        add(st, "B-K2", A, "K-2 simulated S-1 (median of 10 seeds) inside manual 95% CI", k2, k2, k2, lo, hi, "fixed", "window", "nuc",
            f"{rk}", {"PASS": "confermata", "WARN": "smentita", "FAIL": "smentita"}[rk],
            "K-2 nucleus model adopted in derive_cell_geometry() (S1.4d)",
            f"3-level rule: PASS in CI, WARN in CI +/- |S1 Luca - S1 Federica| ({w:.4f}, not corrected), FAIL otherwise; old recomputed {cls(lo,hi)}, new {cls(nlo,nhi)}, stress {cls(slo,shi)}; verdict columns use CI only (WARN maps to smentita)",
            mcenter=s1, clip=(0, 1))
    F = pd.read_csv(R + "R2c/R2c_s1_federica.csv"); F = F[F["set"] == "Federica"].set_index("archetype")
    C5 = pd.read_csv(R + "R2c/R2c_cp5.csv").set_index("archetype")
    nsame = 0; nsame_new = 0
    for A in ("A1", "A2", "A3", "A4"):
        s1, lo, hi = F.loc[A, ["s1_primary", "lo_primary", "hi_primary"]]; cp, sd = F.loc[A, "seg_cellpose_rgb"], F.loc[A, "seg_stardist_he"]
        e_new = None
        old = bdec(s1, lo, hi, cp, sd); nlo, nhi = scale(lo, hi, s1, KWIN); nlo, nhi = max(nlo, 0), min(nhi, 1)
        new = bdec(s1, nlo, nhi, cp, sd)
        same_old = old == S.loc[A, "b_dec"]; same_new = new == newdec[A]; nsame += same_old; nsame_new += same_new
        e_old = float(same_old)
        recm = "confermata" if bool(C5.loc[A, "same"]) else "smentita"
        ROWS.append(dict(step=st, prediction_id="CP-5", archetype=A, quantity="B-dec with Federica's answers equals B-dec with Luca's (1=same)",
            est=e_old, lo=e_old, hi=e_old, a=1, b=1, est_unit="window", margin_unit="fixed", qtype="nuc",
            verdict_recorded=f"same={C5.loc[A,'same']} (Federica {C5.loc[A,'b_dec_Federica']})", reproduces=str((verdict(e_old,e_old,1,1))==recm).upper(),
            lo_new=e_old, hi_new=e_old, a_new=1, b_new=1, kappa_est=1.0, kappa_margin=1.0,
            verdict_new=verdict(e_old, e_old, 1, 1), verdict_new_stress=verdict(e_old, e_old, 1, 1), transition="unaffected", touches_conclusion="none",
            interval_construction="other (2 windows)", side="", margin_construction="fixed", margin_side="", verdict_geom_A5="",
            notes=f"component over 2 windows only (Federica): not correctable, kept unchanged and flagged (rule d). For information only: with Luca's one-stage CI x window kappa, Luca b_dec={newdec[A]}; Federica CI x kappa would give {new}", _old=verdict(e_old,e_old,1,1)))
    CP = pd.read_csv(R + "R2c/R2c_counterproofs.csv")
    thr = {"CP-1": (0.60, 1), "CP-2": (0.70, 1), "CP-3": (0.95, 1), "D-3": (0, 0.05), "CP-4": (0.80, 1)}
    for _, r in CP.iterrows():
        if r.check not in thr: 
            add(st, r.check, r.rater if r.check != "CP-1 per archetipo" else r.rater, f"{r.check} (descriptive)", r.value, r.value, r.value, np.nan, np.nan,
                "fixed", "fixed", "nuc", r.verdict, None, "none", "descriptive, no threshold"); continue
        arch = r.rater if str(r.rater).startswith("A") else "A1;A2;A3;A4;A5;A6"
        e = r.upper if r.check == "D-3" else r.value
        a, b = thr[r.check]
        add(st, r.check, arch, {"CP-1": "Cohen kappa Luca-Federica", "CP-2": "intra-rater kappa", "CP-3": "trap pairs judged separate",
            "D-3": "Clopper-Pearson upper limit of contacts in band", "CP-4": "AUC"}[r.check] + f" ({r.rater})", e, e, e, a, b, "fixed", "fixed", "nuc",
            r.verdict, MAPV.get(r.verdict) if r.verdict in ("PASS",) else None, "none", "point statistic vs convention threshold (Landis & Koch)")
    add(st, "CP-6", "A1;A2;A3;A4;A5;A6", "re-review of fast 'separate' answers (Fisher one-sided)", np.nan, np.nan, np.nan, np.nan, np.nan, "fixed", "fixed", "nuc",
        "non smentita (non sostenuta)", None, "none", "post hoc; test, no interval")
# ======================= R3 =======================
def parse_thr(t):
    t = str(t)
    m = re.match(r"\[([\d.]+), ([\d.]+)\]", t)
    if m: return float(m.group(1)), float(m.group(2))
    m = re.match(r"([\d.]+) ± (\d+)%", t)
    if m: v = float(m.group(1)); p = float(m.group(2)) / 100; return v * (1 - p), v * (1 + p)
    m = re.match(r">=? ?([\d.]+)/(\d+)", t)
    if m: return float(m.group(1)), float(m.group(2))
    m = re.match(r"<= ([\d.]+)/(\d+)", t)
    if m: return 0, float(m.group(1))
    m = re.match(r">=? ?([\d.]+)", t)
    if m: return float(m.group(1)), INF
    m = re.match(r"<=? ?([\d.]+)", t)
    if m: return 0 if float(m.group(1)) < 1 else -INF, float(m.group(1))
    return np.nan, np.nan
def r3():
    st = "R3"; V = pd.read_csv(R + "R3/R3_verdicts.csv")
    for _, r in V.iterrows():
        a, b = parse_thr(r.threshold)
        if r.id == "B-R3.3b": a, b = -INF, 0.222
        if r.id == "CP-4" and r.archetype != "tutti": a, b = 0, INF
        rec = f"{r.assertion} (predicted {r.predicted}; prediction_confirmed={r.prediction_confirmed})"
        recm = MAPV.get(r.assertion)
        tc = "none"
        if r.id in ("B-R3.1", "B-R3.2", "B-R3.3a", "B-R3.3b"): tc = "territory references B-055...B-090 (R3 -> S1.3)"
        if r.id == "CP-2a" or r.id == "CP-2b": tc = "A6 without cell-shape reference (BL-064 -> R3b)"
        add(st, r.id, "A1;A2;A3;A4;A5;A6" if r.archetype == "tutti" else r.archetype, r.statistic, r.value, r.value, r.value, a, b,
            "fixed", "fixed", "nuc", rec, recm, tc,
            "protocol v1 check on a point estimate pooled over 5 ROI (no interval); threshold " + str(r.threshold) + ("; count of ROI out of 5, not an interval" if "/" in str(r.threshold) else ""))
# ======================= R3b =======================
def r3b():
    st = "R3b"; P = pd.read_csv(R + "R3b/R3b_predictions.csv")
    S1 = pd.read_csv(R + "R2c/R2c_s1.csv").set_index("archetype")
    TC = "mask CO (Space Ranger nuclei with StarDist shape) adopted for A4/A6 via decision rule (Luca, R3b closure 2026-10-21)"
    res = {}
    for _, r in P.iterrows():
        A = r.archetype; MK = r["mask"]
        if r.id == "P-1": eu, mu, q, tc, nt = "window", "other", "nuclei of mask / Luca points in same windows", TC, "est = mean of 4 window ratios +- 2 SE (SD/sqrt 4): one-stage window interval, window_4 kappa_se; margin B-046/B-048 Poisson CI of manual total (not across areas or windows): not corrected"
        elif r.id == "P-2": eu, mu, q, tc, nt = "window", "window", "fraction of Luca focal nuclei touching in mask", TC, "est: two-stage bootstrap windows+focal; margin B-130/B-132 two-stage CI over windows. Both two-stage: unchanged (rule d)"
        elif r.id in ("P-3", "CP-1"): eu, mu, q, tc, nt = "seed", "window", "simulated contact fraction on R2c windows" + (" (Luca points as centres)" if r.id == "CP-1" else ""), TC if r.id == "P-3" else "explanation of A4 contact excess (S1.4d) rejected", "est across 10 seeds (size model M2 from leave-one-ROI-out fold, not an interval component); margin B-130/B-132 two-stage CI over windows: unchanged (rule d)"
        else: eu, mu, q, tc, nt = "loao", "fixed", "fraction of simulated nuclei truncated by phi cap (mean of 10 seeds per ROI)", "phi cap and size model coherent (S1.4d criterion)", "est: SE across the 5 R2 ROI, each simulated with the size model M2 of the fold in which it is held out (leave-one-ROI-out)"
        mcen = S1.loc[A, "s1_primary"] if mu == "window" else None
        ecx = "se" if r.id == "P-1" else ("two" if r.id == "P-2" else "boot")
        add(st, f"{r.id} ({MK})", A, q, r.est, r.lo, r.hi, r.int_lo, r.int_hi, eu, mu, "nuc", r.verdict, r.verdict, tc, nt + "; fractions clipped to [0,1] as pre-registered",
            mcenter=mcen, clip=(0, 1) if r.id != "P-1" else None, ec=ecx, em=2, mc="two")
        res[(r.id, A, MK)] = (ROWS[-1]["_old"], ROWS[-1]["verdict_new"], ROWS[-1]["verdict_new_stress"], r.est)
    for A in ("A4", "A6"):
        for k, lab in ((0, "old"), (1, "new"), (2, "stress")):
            elig = [m for m in ("CO", "SD") if res[("P-1", A, m)][k] != "smentita" and res[("P-2", A, m)][k] != "smentita"]
            if len(elig) == 1: ch = elig[0]
            elif len(elig) == 0: ch = "none (to Luca)"
            else:
                cpok = res[("CP-1", A, "manuale+CO")][k] != "smentita"
                p3 = [m for m in elig if res[("P-3", A, m)][k] == "confermata"]
                if cpok and len(p3) == 1: ch = p3[0]
                else: ch = min(elig, key=lambda m: abs(np.log(res[("P-1", A, m)][3])))
            res[("rule", A, lab)] = ch
    for A in ("A4", "A6"):
        o, n, s = res[("rule", A, "old")], res[("rule", A, "new")], res[("rule", A, "stress")]
        ROWS.append(dict(step=st, prediction_id="decision rule (mask)", archetype=A, quantity="mask chosen by pre-registered rule from P-1,P-2,P-3,CP-1",
            est=np.nan, lo=np.nan, hi=np.nan, a=np.nan, b=np.nan, est_unit="window", margin_unit="window", qtype="nuc", verdict_recorded="CO",
            reproduces=str(o == "CO").upper(), lo_new=np.nan, hi_new=np.nan, a_new=np.nan, b_new=np.nan, kappa_est=np.nan, kappa_margin=np.nan,
            verdict_new=n, verdict_new_stress=s, transition="unchanged" if n == o else f"{o}->{n}", touches_conclusion=TC,
            interval_construction="decision rule", side="", margin_construction="decision rule", margin_side="", verdict_geom_A5="",
            notes="decision outcome (mask), not a three-way verdict; recomputed from the corrected verdicts of this table", _old=o))
# ======================= R4 =======================
def r4():
    st = "R4"; P = pd.read_csv(R + "R4/R4_predictions.csv")
    for _, r in P.iterrows():
        A = r.archetype
        if r.id == "P-1":
            add(st, "P-1", A, "UMI/bin nucleus/ring ratio over shifted-mask null", r.estimate, r.lo, r.hi, r.range_lo, r.range_hi, "area", "fixed", "bin",
                r.verdict, r.verdict, f"rule R-1: {r.rule_R1} (Step 2 cell density model; reference later = PH-1 weighted ratio, DEC-034)", "bootstrap over 5 ROI; margin Luca 2026-10-30")
        elif r.id == "CP-1":
            add(st, "CP-1", A, "nucleus/ring ratio with mask shifted 20 um (geometry null)", r.estimate, r.lo, r.hi, r.range_lo, r.range_hi, "area", "fixed", "bin",
                r.verdict, None, "reserve on P-1 reading (explained by PH-1)", "recorded superata/non superata: superata=confermata; non superata covers smentita or indeterminata, reproduces checked as old!=confermata in notes")
            ROWS[-1]["reproduces"] = str((ROWS[-1]["_old"] == "confermata") == (r.verdict == "superata")).upper()
        elif r.id == "P-2":
            add(st, "P-2", A, "median excess zeros per gene within type vs NB", r.estimate, r.lo, r.hi, r.range_lo, r.range_hi, "area", "fixed", "nuc",
                r.verdict, r.verdict, "rule R-2: counts with NB only, no dropout step (Step 2 design)", "cell-level quantity; bootstrap over 5 ROI")
        elif r.id == "CP-2":
            add(st, "CP-2", A, "P-2 statistic after 30% random dropout", r.estimate, r.lo, r.hi, r.range_lo, r.range_hi, "area", "fixed", "nuc",
                r.verdict, "confermata" if r.verdict == "superata" else ("non valutabile" if r.verdict == "non valutabile" else None), "power of P-2", "superata = whole interval above 0.02")
        elif r.id == "P-3":
            add(st, "P-3", A, "number of cells with > 50,000 UMI (exact count)", r.estimate, r.lo, r.hi, 0, 0, "fixed", "fixed", "nuc",
                r.verdict, r.verdict, "none (BL-011 cap)", f"exact count, max {r.max_cell_umi}")
        elif r.id == "P-4":
            rec = r.verdict if r.verdict != "nessuna attesa" else f"nessuna attesa; rule R-3 class '{r.rule_R3}'"
            recm = r.verdict if r.verdict != "nessuna attesa" else {"sopra 1.25": "confermata", "a cavallo di 1.25": "indeterminata"}.get(r.rule_R3)
            add(st, "P-4", A, "max/min median UMI per cell among types >= 5% nuclei", r.estimate, r.lo, r.hi, 1.25, INF, "area", "fixed", "nuc",
                rec, recm, "rule R-3: library size per cell type in preset (Luca chose per type in all six tissues, R4 closure)", "A1-A4 judged against rule R-3 (> 1.25) as recorded")
# ======================= R5 =======================
def r5():
    st = "R5"; E = pd.read_csv(R + "R5/R5_predictions_edges.csv")
    tcs = {"P-1": "rule R-2: sigma per tissue", "P-2": "rule R-1: nucleus = gaussian + exponential tail", "P-3": "rule R-3: isotropic", "P-4": "rule R-4: shift sensitivity in report"}
    for _, r in E.iterrows():
        A = "A2;A3" if r.key == "A23" else ("A1;A2;A3;A4;A5;A6" if r.key == "all" else r.key)
        m = re.findall(r"[\d.]+", str(r.interval)); a, b = (float(m[0]), float(m[1])) if len(m) >= 2 else (float(m[0]), INF)
        if r.prediction == "P-2": a, b = 0.01, INF
        rec = str(r.verdict)
        recm = {"senza attesa (indeterminata)": "indeterminata", "sigma per tessuto": "smentita" if False else None}.get(rec, rec if rec in MAPV else None)
        if r.prediction == "R-2": a, b = 1.25, INF; recm = "confermata"; rec = "sigma per tessuto (whole interval above 1.25)"
        add(st, r.prediction, A, r.quantity, r.estimate, r.lo, r.hi, a, b, "other", "fixed", "bin", rec, recm, tcs.get(r.prediction, "rule R-2: sigma per tissue"),
            "unit = edge tile (bootstrap over tiles), not R2 areas" + ("; P-2 coded as 'tail present' = interval above 0.01" if r.prediction == "P-2" else ""))
    N = pd.read_csv(R + "R5/R5_predictions_nc.csv")
    for _, r in N.iterrows():
        A = r.archetype
        if r.prediction == "P-5":
            ex = r.expected; a, b = {"uniforme": (0.909, 1.1), "due densita (sotto)": (-INF, 0.909), "due densita (sopra)": (1.1, INF)}[ex]
            q = f"r_true nucleus/ring before diffusion, class vs expected '{ex}'"; tc = "rule R-5: preset uses r_true (superseded: CP-1 failed; S2.0 internal check)"
        else:
            a, b = 0.909, 1.1; q = "predicted/observed weighted ratio without bins within 2 um of nucleus"; tc = "reliability of r_true (CP-1 failed outside brain)"
        rec = str(r.verdict); recm = MAPV.get(rec, None)
        add(st, f"{r.prediction} [{r.variant}]", A, q, r.estimate, r.lo, r.hi, a, b, "area", "fixed", "bin", rec, recm, tc,
            "two-source bootstrap (tiles + ROI); ROI = R2 areas")
    C2 = pd.read_csv(R + "R5/R5_CP2_lumens.csv")
    for _, r in C2.iterrows():
        A = "A2;A3" if r.key == "A23" else r.key
        add(st, f"CP-2 [{r.band_um} um]", A, "lumen UMI observed/predicted by band", r.ratio, r.lo, r.hi, 0.8, 1.25, "other", "fixed", "bin", r.verdict, r.verdict,
            "none (BL-077/BL-201)", "bootstrap over tiles")
    C3 = pd.read_csv(R + "R5/R5_CP3.csv")
    for _, r in C3.iterrows():
        A = "A2;A3" if r.key == "A23" else r.key
        for comp, (m, lo, hi, tgt) in {"rows": (r.d_dr_median, r.d_dr_lo, r.d_dr_hi, -2), "cols": (r.d_dc_median, r.d_dc_lo, r.d_dc_hi, 2)}.items():
            add(st, f"CP-3 ({comp})", A, "change of estimated shift when mask moved by 2 um (target +-2 within 0.25)", m, m, m, tgt - 0.25, tgt + 0.25, "other", "fixed", "bin",
                f"FAIL as written (pass={r['pass']})", None, "none", "judged on median change (tile-based); recorded FAIL in 5/5")
# ======================= R6 =======================
def ratio_range(r, k):
    return ((1 + r) + k * (r - 1)) / ((1 + r) - k * (r - 1))
def r6():
    st = "R6"
    for fn, tag in (("R6_predictions.csv", ""), ("R6_predictions_ph1.csv", " [PH-1 labels]")):
        P = pd.read_csv(R + "R6/" + fn)
        for _, r in P.iterrows():
            A = r.archetype; lo_m, up_m = r.lower, r.upper
            kw = {}
            if r.id == "P-1" or r.id == "P-2":
                mu = "range"; rr = lo_m
                kr = K[A]["rnuc"]; kw = dict(a_new=ratio_range(rr, kr), b_new=INF, a_new_s=ratio_range(rr, kr), b_new_s=INF, a_new_geo=ratio_range(rr, GEOA5["rnuc"]) if A == "A5" else None)
                b = INF; nt = "margin = max/min of per-ROI median nuclear area (S1.5 refs): range; corrected on linear scale around midpoint (ratio' = ((1+r)+k(r-1))/((1+r)-k(r-1)))"
                tc = "type-specific nuclear sizes in presets (Luca: one size for lymphocytes in A4; PH-1 labels adopted)"
            elif r.id == "P-3":
                mu = "range"; kr = K[A]["rnuc"]; kw = dict(a_new=lo_m * kr, b_new=INF, a_new_s=lo_m * kr, b_new_s=INF, a_new_geo=lo_m * GEOA5["rnuc"] if A == "A5" else None); b = INF
                nt = "margin = range of per-ROI median eccentricity (0.874-0.770): width scaled by kappa_range_nuc"; tc = "none (cardiomyocyte nuclei not more elongated)"
            elif r.id == "P-4":
                mu = "fixed"; b = INF; nt = "margin Luca's choice"; tc = "single composition for lymph node preset (B/T references provisional)"
            elif r.id == "P-5":
                mu = "fixed"; b = up_m; nt = "margin +-0.15 around design value (Luca)"; tc = "measured type fractions replace design hypotheses in presets"
            else:
                mu = "fixed"; b = INF; nt = "margin Luca's choice"; tc = "labelling criterion validated (decision 1: PH-1 labels)"
            add(st, r.id + tag, A, r.quantity, r.estimate, r.lo, r.hi, lo_m, b, "area", mu, "nuc", r.verdict, r.verdict, tc,
                nt + "; est: bootstrap over 5 ROI (two-stage where stated)" + ("; post hoc labels (PH-1), verdicts of record are pre-registered criterion" if tag else ""), **kw)
# ======================= S2.0 =======================
def s20():
    st = "S2.0"; P = pd.read_csv(R + "S2.0/S2.0_predictions.csv")
    for _, r in P.iterrows():
        A = r.unit.replace("A23", "A2;A3").replace(", ", ";")
        if r.id == "CP-1": continue
        m = re.findall(r"[\d.]+", r.interval); a, b = float(m[0]), float(m[1])
        if r.id in ("P-1", "P-2"):
            eu = "other"; nt = "unit: blocks of 25 rows/cols (P-1) or 300 um tiles (P-2), not R2 areas"
            tc = "rule R-B: random row/column factors added to stripe layer" if r.id == "P-1" else "rule R-B: one parity factor per section (except smentita sections: smooth field)"
        else:
            eu = "area"; nt = "bootstrap over 5 ROI"
            tc = "rule R-A: internal check admits two-population rule C-1 (DEC-038)" if r.id == "P-3" else "rule R-C: multiplicative gamma noise layer per bin with v per tissue"
        add(st, r.id, A, r.quantity, r.estimate, r.lo, r.hi, a, b, eu, "fixed", "bin", r.verdict, MAPV.get(r.verdict), tc, nt + f"; extra: {r.extra}")
        if r.id == "P-3" and r.unit == "A2":
            mm = re.findall(r"([\d.]+) \[([\d.]+), ([\d.]+)\]", r.extra)[0]
            add(st, "P-3 [sensitivity: successful ROI]", "A2", r.quantity, mm[0], mm[1], mm[2], a, b, "area", "fixed", "bin",
                "not judged (as written: non calcolabile; read as not smentita in R-A, deviation 8)", None, "rule R-A (DEC-038)", "sensitivity in 'extra' column")
    L = pd.read_csv(R + "S2.0/S2.0_CP1_lumens.csv"); L = L[L.judged_CP1 == True]
    for _, r in L.iterrows():
        A = "A2;A3" if r.key == "A23" else r.key
        add(st, f"CP-1 [{r.band} um]", A, "inner lumens observed/predicted with rule C2", r.obs_over_pred, r.lo, r.hi, 0.8, 1.25, "other", "fixed", "bin",
            r.verdict, r.verdict, "generality of rule C2 inside tissue", "bootstrap over tiles; CP-1 passes if no cell smentita (recorded superata, 0/8)")
# ======================= S2.1 =======================
def s21():
    st = "S2.1"; C = pd.read_csv(R + "S2.1/S2.1_CP1_roi.csv")
    for A, g in C.groupby("archetype"):
        n = g[g.grid == "none"]
        add(st, "CP-1a", A, "max over 5 areas of fraction of bin-nucleus pairs with fraction in [0.1,0.9]", n.frac_mid.max(), n.frac_mid.max(), n.frac_mid.max(), 0, 0.001,
            "fixed", "fixed", "box", "confermata (CP-1)", "confermata", "build_hd_grid() geometry accepted (S2.1 closure)", "deterministic per-area count; no interval")
        add(st, "CP-1b", A, "min over 5 areas of Jaccard (our bins >= 0.5 vs Space Ranger)", n.jaccard.min(), n.jaccard.min(), n.jaccard.min(), 0.999, 1,
            "fixed", "fixed", "box", "confermata (CP-1)", "confermata", "build_hd_grid() geometry accepted (S2.1 closure)", "deterministic per-area count; no interval")
        for gname in ("half", "transpose", "reverse"):
            w = g[g.grid == gname]
            add(st, f"CP-1c ({gname})", A, f"max Jaccard with wrong grid '{gname}' (must be < 0.9)", w.jaccard.max(), w.jaccard.max(), w.jaccard.max(), 0, 0.9,
                "fixed", "fixed", "box", "confermata (CP-1)", "confermata", "power of CP-1", "deterministic; no interval")
# ======================= S2.2 =======================
def s22():
    st = "S2.2"; P = pd.read_csv(R + "S2.2/S2.2_predictions.csv")
    TCW = "winner arm used for P-1...P-4 (DEC-045); verdicts not re-run on another arm"
    for _, r in P.iterrows():
        A = r.archetype
        add(st, r.id, A, f"{r.quantity} relative error sim vs real, arm {r.arm}", r.estimate, r.lo, r.hi, -r.relevance, r.relevance, "loao", "range", "bin",
            r.verdict, r.verdict, "none" if r.id != "P-3" else "none (P-3 A1 smentita explained by PH-1: spatial position of subtypes)",
            "est: mean over the 5 held-out areas (profiles from the other 4: leave-one-area-out folds) +-2 SE between folds; margin = max difference between two real areas (range), centred on 0", mcenter=0.0)
    p5 = pd.read_csv(R + "S2.2/S2.2_P5_CP2.csv").set_index("arm")
    r = p5.loc["exp_G"]
    add(st, "P-5", "A1", "excess zeros index, simulated - real (intestine)", r["diff"], r.diff_lo, r.diff_hi, -0.02, 0.02, "loao", "fixed", "nuc", "smentita", "smentita",
        "motivates default 'one total distribution per subtype' (DEC-046, with PH-3)", "bootstrap over held-out areas of the race bench (leave-one-area-out folds)")
    r = p5.loc["exp_G_nosub"]
    add(st, "CP-2", "A1", "excess zeros index with subtypes off (simulated); must be < 0.02", r.sim, r.sim_lo, r.sim_hi, -INF, 0.02, "loao", "fixed", "nuc", "confermata", "confermata",
        "none", "CP-2 judged on the simulated index (0.004) as written; interval sim_lo/sim_hi")
    D = pd.read_csv(R + "S2.2/S2.2_depth_factor.csv")
    for _, r in D.iterrows():
        a, b = (1, 2) if r.archetype != "A6" else (np.nan, np.nan)
        rec = {"A1": "smentita", "A2": "indeterminata", "A3": "confermata", "A4": "smentita", "A5": "confermata", "A6": "nessun intervallo (come atteso, circa 3)"}[r.archetype]
        add(st, "P-6", r.archetype, "depth factor f", r.f, r.lo, r.hi, a, b, "area", "fixed", "bin", rec, MAPV.get(rec, "non valutabile"),
            "depth factor per tissue (R-3, DEC-045) - value used regardless of verdict", "percentile bootstrap over 5 areas (flagged optimistic in step, deviation 5); real UMI / simulated cells per area", alt_sim=True)
    C = pd.read_csv(R + "S2.2/S2.2_CP1.csv")
    for _, r in C.iterrows():
        add(st, f"CP-1 ({r.null})", r.archetype, "paired score difference null - winner (must be > 2 SE)", r["diff"], r["diff"] - 2 * r.se, r["diff"] + 2 * r.se, 0, INF,
            "loao", "fixed", "box", f"null_loses={r.null_loses}", "confermata" if r.null_loses else None, "power of the race; A6 winner read as order choice",
            "SE between 5 areas")
    S = pd.read_csv(R + "S2.2/S2.2_race_scores.csv"); W = pd.read_csv(R + "S2.2/S2.2_race_winners.csv").set_index("archetype")
    srank = {"atlas": 0, "nuc": 1, "exp": 2}; mrank = {"M": 0, "G": 1}
    def order(a): s_, m_ = a.split("_"); return (srank[s_], mrank[m_])
    for A, g in S.groupby("archetype"):
        best = W.loc[A, "best_score"]; ties = {"old": [best], "new": [best], "stress": [best]}
        for _, r in g.iterrows():
            if r.arm == best or not r.is_second: continue
            add(st, f"R-1 race ({r.arm} vs {best})", A, "paired score difference arm - best (tie if < 2 SE)", r.diff_vs_best, r.diff_vs_best - 2 * r.se_diff,
                r.diff_vs_best + 2 * r.se_diff, 0, INF, "loao", "fixed", "box", f"tie_with_best={r.tie_with_best}",
                "indeterminata" if r.tie_with_best else "confermata", "profile source and count model per tissue (R-1, DEC-045)",
                "decision rule mapped to verdict (tie rule applies to second vs first only): confermata = clear loss, indeterminata = tie (then order rule)")
            if ROWS[-1]["_old"] != "confermata": ties["old"].append(r.arm)
            if ROWS[-1]["verdict_new"] != "confermata": ties["new"].append(r.arm)
            if ROWS[-1]["verdict_new_stress"] != "confermata": ties["stress"].append(r.arm)
        win = {k: (v[0] if len(v) == 1 else min(v, key=order)) for k, v in ties.items()}
        ROWS.append(dict(step=st, prediction_id="R-1 winner", archetype=A, quantity="winning arm (source_model) by pre-registered race rule",
            est=np.nan, lo=np.nan, hi=np.nan, a=np.nan, b=np.nan, est_unit="loao", margin_unit="fixed", qtype="box", verdict_recorded=W.loc[A, "winner"],
            reproduces=str(win["old"] == W.loc[A, "winner"]).upper(), lo_new=np.nan, hi_new=np.nan, a_new=np.nan, b_new=np.nan, kappa_est=1.0, kappa_margin=1.0,
            verdict_new=win["new"], verdict_new_stress=win["stress"], transition="unchanged" if win["new"] == win["old"] else f"{win['old']}->{win['new']}",
            touches_conclusion="profile source and count model per tissue (DEC-045)", interval_construction="loao", side="", margin_construction="fixed", margin_side="", verdict_geom_A5="" if A != "A5" else win["old"], notes=f"LOAO component not corrected (held-out area shares up to 60% of nuclei with training areas in A3/A4; S2.2 deviation 5, PH-5); tie sets old={ties['old']}, new={ties['new']}, stress={ties['stress']}; order atlas>nuc>exp, M>G; atlas not run in A4/A5 (deviation 1)", _old=win["old"]))
# ======================= S2.3 =======================
def s23():
    st = "S2.3"; P = pd.read_csv(R + "S2.3/S2.3_P123_bench.csv")
    for _, r in P.iterrows():
        pid = {"q": "P-1", "e1": "P-2", "e2": "P-3"}[r.measure]
        mu = "fixed" if pid == "P-1" else "range"
        if pid == "P-1":
            tc = "R-1 calibration of the field on nucleus (DEC-049)" if r.group == "nucleus" else "R-4: one field outside nuclei, limit declared (DEC-050)"
        else: tc = "R-4 context (DEC-050)" if r.group != "nucleus" else "none"
        add(st, f"{pid} ({r.group})", r.archetype, f"bench {r.measure} sim vs obs, group {r.group}", r.estimate, r.lo2se, r.hi2se, r.int_lo, r.int_hi, "area", mu, "bin",
            r.verdict, MAPV.get(r.verdict), tc, "est +-2 SE between 5 areas" + ("; margin = max difference between two real areas (range) centred on 0" if mu == "range" else "; margin Luca 2026-11-03"),
            mcenter=0.0 if mu == "range" else None, ec="se", em=2, alt_sim=True)
    G = pd.read_csv(R + "S2.3/S2.3_P4_verdicts.csv"); G = G[G.variant == "come scritto"]
    for _, r in G.iterrows():
        A = "A2;A3" if r.key == "A23" else r.key
        add(st, f"P-4 [{r.band} um]", A, "glass observed/predicted by distance band", r.ratio, r.lo, r.hi, 0.8, 1.25, "other", "fixed", "bin", r.verdict, MAPV.get(r.verdict),
            "R-2: second tail term in background layer (DEC-051)", "bootstrap over tiles")
    names = {"panel_star": ("P-5", "contrast nucleus/corona sim/real", "none"), "panel_star_far": ("P-6", "internal control sim/real", "none"),
             "panel_radial_R": ("P-7", "radial ratio R sim/real", "R-3: BL-087 open if smentita in A5 or A6 (radial gradient option)"),
             "panel_nucleus_psi": ("P-8 (nucleus)", "psi nucleus sim/real", "none (explained by depth, PH-2)"),
             "panel_corona_psi": ("P-8 (corona)", "psi corona sim/real", "none"),
             "lumen_band_2_4": ("P-9 [2-4 um]", "lumen band 2-4 um sim/real", "none"), "lumen_band_4_8": ("P-9 [4-8 um]", "lumen band 4-8 um sim/real", "none")}
    Q = pd.read_csv(R + "S2.3/S2.3_P5_P9.csv")
    for _, r in Q.iterrows():
        pid, q, tc = names[r.stat]
        add(st, pid, r.archetype, q, r.estimate, r.lo2se, r.hi2se, r.int_lo, r.int_hi, "area", "fixed", "bin", r.verdict, MAPV.get(r.verdict), tc, "est +-2 SE between 5 areas (panel sim/real per area)", ec="se", em=2, alt_sim=True)
    C = pd.read_csv(R + "S2.3/S2.3_CP.csv")
    for _, r in C.iterrows():
        add(st, f"{r.cp} ({r.stat})", r.archetype, f"{r.stat} sim/real in control panel", r.estimate, r.lo2se, r.hi2se, r.int_lo, r.int_hi, "area", "fixed", "bin",
            r.verdict, MAPV.get(r.verdict), "power of P-8/P-9/P-7", "counter-proof passes when smentita (expected out of interval); est +-2 SE between 5 areas", ec="se", em=2, alt_sim=True)
# ======================= S2.4 =======================
def s24():
    st = "S2.4"; V = pd.read_csv(R + "S2.4/S2.4_verdicts.csv")
    for _, r in V.iterrows():
        A = "A2;A3" if r.A == "A23" else r.A
        if r.id in ("P-1", "P-2", "P-3"):
            add(st, f"{r.id} ({r.quantita})", A, f"{r.quantita}, sim - real", r.stima, r.stima - 2 * r.se, r.stima + 2 * r.se, r.lo, r.hi, "area", "range", "bin", r.esito, MAPV.get(r.esito),
                "stripe generator parameters (parity per section, fine SD, lag-1 from R-1, DEC-055)", "est +-2 SE between areas (sim - real per area); margin = max difference between two real areas (range) centred on 0", mcenter=0.0, ec="se", em=2, alt_sim=True)
        elif r.id in ("P-4", "P-6", "P-7"):
            qt = "box" if r.id == "P-7" else "bin"
            tc = {"P-4": "second tail values (DEC-053); heart = negative control", "P-6": "none (stripes do not create empty bins: declared limit)", "P-7": "none (resource budget, design s12)"}[r.id]
            add(st, f"{r.id} ({r.quantita})", A, r.quantita, r.stima, r.stima - 2 * r.se, r.stima + 2 * r.se, r.lo, r.hi, "area" if r.n > 1 else "other", "fixed", qt, r.esito,
                MAPV.get(r.esito) if r.esito in MAPV else None, tc, f"est +-2 SE between {int(r.n) if r.n==r.n else 0} areas" + ("; descriptive per pre-registration" if r.esito == "descrittivo" else "") +
                ("; computational quantity of the simulation only" if r.id == "P-7" else ""), ec="se", em=2,
                es=("simulated" if r.id == "P-7" else "real"), alt_sim=(r.id != "P-7"))
        elif r.id == "P-5":
            add(st, f"P-5 ({r.quantita})", A, "glass observed/predicted after R-3", r.stima, r.ci_lo, r.ci_hi, r.lo, r.hi, "other", "fixed", "bin", r.esito, MAPV.get(r.esito),
                "R-3 background + second tail (DEC-053); tumour far-glass limit (DEC-059)", "bootstrap over 300 um tiles" + ("; single tile: step rule declares non valutabile (kept)" if r.n == 1 else ""), force_nv=(r.n == 1))
    C = pd.read_csv(R + "S2.4/S2.4_CP12.csv")
    for _, r in C.iterrows():
        if r.id == "CP-1":
            add(st, "CP-1 (sd fine)", r.A, "fine SD with stripes off (< 0.05)", r.sd_fine_max, r.sd_fine_max, r.sd_fine_max, 0, 0.05, "fixed", "fixed", "bin", f"ok={r.ok}", "confermata" if r.ok else None, "power of P-1...P-3", "mean over areas, point")
            add(st, "CP-1 (parity)", r.A, "|parity - 1| with stripes off (<= 0.01)", r.parity_dev_max, r.parity_dev_max, r.parity_dev_max, 0, 0.01, "fixed", "fixed", "bin", f"ok={r.ok}", "confermata" if r.ok else None, "power of P-1...P-3", "mean over areas, point")
        else:
            add(st, f"CP-2 [{r.band}]", "A2", "lumen band with tail - without tail (must be > 2 SE)", r.main_minus_off, r.main_minus_off - 2 * r.se, r.main_minus_off + 2 * r.se, 0, INF,
                "area", "fixed", "bin", f"ok={r.ok}", "confermata" if r.ok else None, "power of P-4 on second tail", "SE = SD of bootstrap over simulated areas (tags) of A2, judged as > 2 SE (tools/S2.4_analysis.py l.90); both terms simulated", ec="se", em=2, es="simulated")
    C3 = pd.read_csv(R + "S2.4/S2.4_CP3.csv")
    for _, r in C3.iterrows():
        A = "A2;A3" if r.key == "A23" else r.key
        exp = r.control == "positive"
        add(st, f"CP-3 ({r.control}, {r.axis})", A, "slow part identified on synthetic glass (positive: yes; negative: no)", float(r.identified), float(r.identified), float(r.identified),
            1 if exp else 0, 1 if exp else 0, "other", "fixed", "bin", f"identified={r.identified}", None, "R-2 slow part not added (DEC-056)", "synthetic replicate of each section; recorded CP-3 FAIL as written (A4 rows, A6)")
def main():
    for fn in (r2b, r2c, r3, r3b, r4, r5, r6, s20, s21, s22, s23, s24):
        try: fn()
        except Exception as e:
            import traceback; traceback.print_exc(); print("FAILED", fn.__name__)
    D = pd.DataFrame(ROWS)
    cols = ["step","prediction_id","archetype","quantity","est","lo","hi","a","b","est_unit","margin_unit","qtype","verdict_recorded","reproduces",
            "lo_new","hi_new","a_new","b_new","kappa_est","kappa_margin","verdict_new","verdict_new_stress","verdict_geom_A5","transition","touches_conclusion","interval_construction","side","margin_construction","margin_side","notes"]
    for cc in ("verdict_geom_A5", "interval_construction", "side", "margin_construction", "margin_side"):
        D[cc] = D[cc].fillna("")
    D["reproduces"] = D["reproduces"].astype(str).str.upper().replace({"NAN": "NA"})
    D[cols].to_csv(OUT + "R7_M5_R_S2.csv", index=False)
    print(D.groupby("step").size().to_dict())
    print("rows", len(D), "affected", (D.transition != "unaffected").sum())
    print(D.transition.value_counts().to_dict())
    print("not reproduced:", D[D.reproduces == "FALSE"][["step","prediction_id","archetype","verdict_recorded"]].to_string())
main()
