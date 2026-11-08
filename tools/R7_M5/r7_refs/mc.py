# R7 M-c: corrected reference intervals (read-only on repo; writes /tmp/r7_refs)
import re, numpy as np, pandas as pd
R = "/home/user/2025.geo_spatialtrans/results/"
W = "/tmp/r7_refs/"
bio = pd.read_csv(W + "bio_all.csv", dtype=str).set_index("id")
sel = pd.read_csv(W + "sel.csv", dtype=str)
K = pd.read_csv(W + "kappa4.csv"); CVF = pd.read_csv(W + "cvf.csv")
kap = K.set_index(["archetype", "unit"])
T975_3 = 3.182446305284263  # t quantile 0.975, df 3 (as in R2b)
rows = []; checks = []
selu = sel.set_index("id")

def num(x):
    return float(str(x).replace("\u2212", "-"))

def rx(pat, s, n=2):
    m = re.search(pat, s)
    if not m: return None
    return [num(g) for g in m.groups()[:n]]

GEO5 = {"nuc": (1.751, 1.496, 1.020), "bin": (1.760, 1.503, 1.019), "box": (1.743, 1.491, 1.018)}  # boot, se, range (geometric A5)
STRESS = 1.343

def kappa_for(a, constr, q, m=None, side="real", geo=False):
    """constr: boot / se_m / range / window_one_stage_boot / window_one_stage_se / others -> NaN (unchanged)."""
    if geo:
        kb, ks, kr = GEO5[q]
    else:
        row = ("all", "sim_5") if side == "simulated" else (a, "area")
        kb, ks, kr = (float(kap.loc[row, f"kappa_boot_{q}"]), float(kap.loc[row, f"kappa_se_{q}"]), float(kap.loc[row, f"kappa_range_{q}"]))
    w = ("all", "window_4")
    if constr == "boot": return kb
    if constr == "se_m": return ks * 2.0 / m
    if constr == "range": return kr
    if constr == "window_one_stage_boot": return float(kap.loc[w, f"kappa_boot_{q}"])
    if constr == "window_one_stage_se": return float(kap.loc[w, f"kappa_se_{q}"]) * 2.0 / m
    return np.nan

def cvfac(cv):
    if cv > 1.0: return float(CVF.cv_factor.iloc[-1]), "CV>1.0: cv_factor capped at the CV=1.0 value"
    return float(np.interp(cv, CVF.cv.values, CVF.cv_factor.values)), ""

def dec(x):
    s = str(x)
    return len(s.split(".")[1]) if "." in s else 0

def add(bid, metric, kind, q, est, lo, hi, source, per_area=None, per_src="", s25=False, notes="", text_lohi=None, tol=None,
        constr=None, m=None, side="real"):
    if constr is None:
        constr = {"ci": "boot", "range": "range", "loao": "loao", "other": "fixed"}[kind]
    b = bio.loc[bid]; a = b.archetype; n = [notes] if notes else []
    lo = float(lo); hi = float(hi); est = float(est) if est is not None else np.nan
    if text_lohi is not None:
        tl, th = text_lohi
        t = tol if tol is not None else 0.51 * 10 ** (-max(dec(tl), dec(th)))
        ok = abs(lo - tl) <= t + 1e-9 and abs(hi - th) <= t + 1e-9
        checks.append((bid, metric, lo, hi, tl, th, ok))
        if not ok: n.append(f"file lo/hi {lo:.4g}/{hi:.4g} differ from row text {tl}/{th}")
    elif source != "row text":
        checks.append((bid, metric, lo, hi, np.nan, np.nan, None))
    corrected = constr in ("boot", "se_m", "range", "window_one_stage_boot", "window_one_stage_se")
    k = kappa_for(a, constr, q, m, side)
    centre = np.nan
    if constr == "range":
        centre = (lo + hi) / 2
    elif np.isfinite(est):
        centre = est
    else:
        centre = (lo + hi) / 2
        if corrected: n.append("no point estimate: scaled around midpoint")
    sc = lambda kk: (centre - kk * (centre - lo), centre + kk * (hi - centre))
    lo_new, hi_new = sc(k) if corrected else (lo, hi)
    if constr == "loao": n.append("LOAO: not corrected (rule 4); held-out area shares up to 60% of its nuclei with the training areas in A3/A4")
    if constr == "window_two_stage": n.append("two-stage window interval: not widened (R-3 does not hold for windows), unchanged and flagged")
    if constr in ("se_m", "window_one_stage_se"): n.append(f"+-m SE/t interval, m = {m:.4f}: factor = kappa_se x 2/m")
    wr = (hi_new - lo_new) / (hi - lo) if (hi - lo) > 0 else np.nan
    if (hi - lo) == 0: n.append("zero-width interval: correction has no effect")
    cv = np.nan; lo_cv = hi_cv = np.nan
    if per_area is not None:
        v = np.asarray(per_area, float); v = v[np.isfinite(v)]
        if len(v) >= 2 and np.mean(v) != 0:
            cv = float(np.std(v, ddof=1) / abs(np.mean(v)))
            n.append(f"per-area n={len(v)} from {per_src}")
    if np.isfinite(cv) and cv > 0.45 and constr in ("boot", "se_m", "range"):
        f, fl = cvfac(cv)
        lo_cv, hi_cv = sc(k * f)
        n.append(f"cv_factor {f:.4f}")
        if fl: n.append(fl)
        if lo_cv < 0 <= lo: n.append("lo_cv < 0 (not clipped)")
    if lo_new < 0 <= lo: n.append("lo_new < 0 for a non-negative quantity (literal rule, not clipped)")
    # stress: area (boot, se_m) and window factors x 1.343; range and unchanged kinds stay as lo_new
    if constr in ("boot", "se_m", "window_one_stage_boot", "window_one_stage_se"):
        lo_st, hi_st = sc(k * STRESS)
        if lo_st < 0 <= lo: n.append("lo_stress < 0 (not clipped)")
    else:
        lo_st, hi_st = lo_new, hi_new
    # geometric A5 factors, A5 rows only
    lo_g = hi_g = np.nan
    if a == "A5":
        if constr in ("boot", "se_m", "range"):
            lo_g, hi_g = sc(kappa_for(a, constr, q, m, side, geo=True))
            if lo_g < 0 <= lo: n.append("lo_geomA5 < 0 (not clipped)")
        else:
            lo_g, hi_g = lo_new, hi_new
    ic = {"window_one_stage_boot": "window_one_stage", "window_one_stage_se": "window_one_stage"}.get(constr, constr)
    rows.append(dict(id=bid, archetype=a, metric=metric, session=b.session, unit=selu.loc[bid, "unit"],
                     interval_kind=kind, interval_construction=ic, side=side, m_se=m if m is not None else np.nan, qtype=q,
                     est=est, lo=lo, hi=hi, kappa=k, lo_new=lo_new, hi_new=hi_new, width_ratio=wr,
                     cv_between_areas=cv, lo_cv=lo_cv, hi_cv=hi_cv, lo_stress=lo_st, hi_stress=hi_st, lo_geomA5=lo_g, hi_geomA5=hi_g,
                     source=source, is_S25_target=s25, notes="; ".join(n)))

def ids_where(session, pat):
    s = sel[sel.session == session]
    return [i for i in s.id if re.search(pat, bio.loc[i, "metric"])]

def bold(bid):
    return re.search(r"\*\*(.+?)\*\*", bio.loc[bid, "metric"]).group(1)

# ---------- R2: nuclear area, Cellpose-RGB scale 1
r2 = pd.read_csv(R + "R2/R2_roi_summary.csv"); r2 = r2[(r2.method == "cellpose_rgb") & (r2.scale == 1)]
for bid in ids_where("R2", "area nucleare mediana"):
    a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]; v = r2[r2.archetype == a].area_median.values
    add(bid, "median nuclear area (Cellpose-RGB), min-max of ROI medians", "range", "nuc", num(bio.loc[bid, "value"]), v.min(), v.max(),
        "results/R2/R2_roi_summary.csv:area_median (cellpose_rgb, scale 1)", v, "R2_roi_summary.csv:area_median",
        text_lohi=rx(r"mediane ROI: ([\d.]+)\u2013([\d.]+)", t), notes=f"est = row value; median of ROI medians in file {np.median(v):.2f}; n ROI {len(v)}")

# ---------- R2b: manual density, t CI across 4 windows
s2b = pd.read_csv(R + "R2b/R2b_archetype_summary.csv"); s2b = s2b[s2b.rater == "Luca"].set_index("archetype")
mw = pd.read_csv(R + "R2b/R2b_manual_windows.csv"); mw = mw[mw.rater == "Luca"]
for bid in sel[sel.session == "R2b"].id:
    a = bio.loc[bid, "archetype"]; r = s2b.loc[a]; m, sd, nw = r.density_mean_windows, r.density_sd_windows, int(r.n_windows)
    hw = T975_3 * sd / np.sqrt(nw); dw = mw[mw.archetype == a].density_manual_net.values
    cvw = np.std(dw, ddof=1) / np.mean(dw)
    add(bid, "manual nuclear density, t 95% CI across the 4 windows", "window", "nuc", m, m - hw, m + hw,
        "results/R2b/R2b_archetype_summary.csv:density_mean_windows,density_sd_windows (rater Luca; t df=3)",
        text_lohi=rx(r"IC t 95 % fra le 4 finestre ([\d.]+)\u2013([\d.]+)", bio.loc[bid, "interval"]), tol=1.0,
        constr="window_one_stage_se", m=T975_3,
        notes=f"CV between the {len(dw)} windows (Luca) {cvw:.3f}: windows, not areas, no cv_factor applied; the row's reference interval [evidenti, totale] and the Poisson CI are not window intervals and stay unchanged")

# ---------- R3: territories (cellpose_rgb)
r3 = pd.read_csv(R + "R3/R3_roi_summary.csv"); r3 = r3[r3.method == "cellpose_rgb"]
for bid in sel[sel.session == "R3"].id:
    a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]; met = bio.loc[bid, "metric"]; d = r3[r3.archetype == a]
    if met.startswith("area mediana del territorio"):
        v = d.median_area.values
        add(bid, "median interior territory area, min-max of ROI medians", "range", "nuc", num(bio.loc[bid, "value"]), v.min(), v.max(),
            "results/R3/R3_roi_summary.csv:median_area (cellpose_rgb)", v, "R3_roi_summary.csv:median_area",
            text_lohi=rx(r"mediane ROI ([\d.]+)\u2013([\d.]+)", t), notes="the bracket [CP-4 corrected, value] is a fixed segmenter correction, unchanged")
    elif met.startswith("eq_radius"):
        lo, hi = rx(r"\[([\d.]+), ([\d.]+)\]", t); v = d.median_eq_r.values
        add(bid, "median territory eq_radius, bracket [value / CP-4 factor, value]", "other", "nuc", num(bio.loc[bid, "value"]), lo, hi,
            "row text", v, "R3_roi_summary.csv:median_eq_r",
            notes="RECLASSIFIED window->other: the bracket is a fixed correction by the CP-4 factor estimated on the 4 windows, not an interval across windows; no window CI of the factor is reported, nothing to scale")
    else:
        v = d.cv_loc.values
        add(bid, "local CV of territory area, min-max across ROIs", "range", "nuc", num(bio.loc[bid, "value"]), v.min(), v.max(),
            "results/R3/R3_roi_summary.csv:cv_loc (cellpose_rgb)", v, "R3_roi_summary.csv:cv_loc",
            text_lohi=rx(r"min\u2013max ROI ([\d.]+)\u2013([\d.]+)", t))

# ---------- R3b (A4, A6)
for bid in sel[sel.session == "R3b"].id:
    t = bio.loc[bid, "interval"]; met = bio.loc[bid, "metric"]
    if met.startswith("area mediana"):
        lo, hi = rx(r"mediane ROI ([\d.]+)-([\d.]+)", t); lab = "median interior territory area (R3b CO centres), min-max of ROI medians"
    else:
        lo, hi = rx(r"min-max ROI ([\d.]+)-([\d.]+)", t); lab = "local CV of territory area (R3b), min-max across ROIs"
    add(bid, lab, "range", "nuc", num(bio.loc[bid, "value"]), lo, hi, "row text", notes="per-ROI values not in a machine-readable file (R3b_refs.csv has only min/max)")

# ---------- S1.4: C8 interval (LOAO power check)
for bid in sel[sel.session == "S1.4"].id:
    lo, hi = rx(r"\[([\d.]+), ([\d.]+)\]", bio.loc[bid, "value"]); est = rx(r"mediana ([\d.]+)", bio.loc[bid, "interval"], 1)[0]
    add(bid, "C8 interval [q05, q95] of interior eq_radius (pooled cells); power checked leaving one ROI out", "loao", "nuc", est, lo, hi, "row text",
        notes="RECLASSIFIED area_R2->loao: the interval is a pooled-cell quantile range; its only area-based element is the leave-one-ROI-out power check")

# ---------- S1.4b: contact fractions
fb = pd.read_csv(R + "S1.4b/S1.4b_b8_b9.csv")
for bid in sel[sel.session == "S1.4b"].id:
    a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]; e1, e2 = [num(x) for x in bio.loc[bid, "value"].split(" / ")]
    lo, hi = rx(r"cerchio ([\d.]+)\u2013([\d.]+)", t); v = fb[fb.archetype == a].b8_circle.values
    use = len(v) == 5 and abs(v.min() - lo) < 0.0015 and abs(v.max() - hi) < 0.0015
    add(bid, "fraction of nuclei touching a neighbour, equal-area circles, min-max per ROI", "range", "nuc", e1, lo, hi, "row text",
        v if use else None, "S1.4b_b8_b9.csv:b8_circle (fold = held-out ROI; min/max match row text)" if use else "",
        notes="" if use else "per-ROI values not matched to a file")
    lo, hi = rx(r"ellisse ([\d.]+)\u2013([\d.]+)", t)
    add(bid, "fraction of nuclei touching a neighbour, equal-area moment ellipses, min-max per ROI", "range", "nuc", e2, lo, hi, "row text",
        notes="per-ROI values not in a machine-readable file")

# ---------- S1.4c: contact in masks
rc = pd.read_csv(R + "S1.4c/S1.4c_cp2.csv"); rc = rc[rc.metric == "S1"].set_index("archetype")
for bid in sel[sel.session == "S1.4c"].id:
    a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]; e1, e2 = [num(x) for x in bio.loc[bid, "value"].split(" / ")]
    r = rc.loc[a]
    col = {r.primary: "primary_value", r.other: "other_value"}
    for seg, e, pat in (("cellpose_rgb", e1, r"cellpose ([\d.]+)\u2013([\d.]+)"), ("stardist_he", e2, r"StarDist ([\d.]+)\u2013([\d.]+)")):
        v = np.array([float(x) for x in r[col[seg]].split(";")])
        add(bid, f"fraction of nuclei in contact in masks ({seg}), min-max per ROI", "range", "nuc", e, v.min(), v.max(),
            f"results/S1.4c/S1.4c_cp2.csv:{col[seg]} (metric S1, {seg}; 3 significant digits)", v, f"S1.4c_cp2.csv:{col[seg]}",
            text_lohi=rx(pat, t))

# ---------- R2c: S-1 manual contact fraction, two-stage bootstrap (windows, focal nuclei)
pc = pd.read_csv(R + "R2c/R2c_posthoc_ci.csv"); pc = pc[pc.ci == "due stadi"].set_index("archetype")
s1 = pd.read_csv(R + "R2c/R2c_s1.csv"); s1 = s1[s1.set == "Luca"].set_index("archetype")
for bid in sel[sel.session == "R2c"].id:
    a = bio.loc[bid, "archetype"]; r = s1.loc[a]
    add(bid, "S-1 manual fraction of focal nuclei with a touching neighbour, one-stage window-cluster percentile bootstrap (pre-registered)", "window", "nuc",
        r.s1_primary, r.lo_primary, r.hi_primary, "results/R2c/R2c_s1.csv:s1_primary,lo_primary,hi_primary (set Luca)",
        constr="window_one_stage_boot", notes=f"one-stage interval over the 4 windows (rule d); n focal {int(r.n_focal)}")
    r = pc.loc[a]
    add(bid, "S-1 manual fraction of focal nuclei with a touching neighbour, two-stage bootstrap (windows, focal nuclei), as cited in the row", "window", "nuc",
        r.est, r.lo, r.hi, "results/R2c/R2c_posthoc_ci.csv:est,lo,hi (ci = due stadi)",
        text_lohi=rx(r"due stadi ([\d.]+)\u2013([\d.]+)", bio.loc[bid, "interval"]), tol=0.006, constr="window_two_stage")

# ---------- R6: composition and size per type (PH-1 labels)
t6 = pd.read_csv(R + "R6/R6_refs_type_ph1.csv"); p6 = pd.read_csv(R + "R6/R6_refs_roi_type_ph1.csv")
for bid in sel[sel.session == "R6"].id:
    a = bio.loc[bid, "archetype"]; ty = bold(bid); r = t6[(t6.archetype == a) & (t6.cell_type == ty)].iloc[0]
    pr = p6[(p6.archetype == a) & (p6.cell_type == ty)]; t = bio.loc[bid, "interval"]
    add(bid, f"{ty}: fraction of nuclei, bootstrap over ROIs", "ci", "nuc", r.fraction, r.fraction_lo, r.fraction_hi,
        "results/R6/R6_refs_type_ph1.csv:fraction,fraction_lo,fraction_hi", pr.fraction.values, "R6_refs_roi_type_ph1.csv:fraction",
        text_lohi=rx(r"bootstrap ROI \[([\d.]+), ([\d.]+)\]", t))
    add(bid, f"{ty}: nuclear density per ROI, min-max", "range", "nuc", r.density_mm2_mean_roi, r.density_mm2_min_roi, r.density_mm2_max_roi,
        "results/R6/R6_refs_type_ph1.csv:density_mm2_min_roi,density_mm2_max_roi", pr.density_mm2.values, "R6_refs_roi_type_ph1.csv:density_mm2",
        text_lohi=rx(r"min-max ([\d.]+)-([\d.]+)\)", t), tol=1.0)
    add(bid, f"{ty}: median nuclear area, bootstrap over ROIs", "ci", "nuc", r.median_area_nuc, r.area_lo, r.area_hi,
        "results/R6/R6_refs_type_ph1.csv:median_area_nuc,area_lo,area_hi", pr.median_area_nuc.values, "R6_refs_roi_type_ph1.csv:median_area_nuc",
        text_lohi=rx(r"\u00b5m\u00b2 \[([\d.]+), ([\d.]+)\]", t))

# ---------- S1.5d: leave-one-out folds
for bid in sel[sel.session == "S1.5d"].id:
    t = bio.loc[bid, "interval"]; v = bio.loc[bid, "value"]
    if "void_frac" in bio.loc[bid, "metric"]:
        lo, hi = rx(r"void_frac ([\d.]+)-([\d.]+)", t); est = num(v.split(" / ")[0]); lab = "seeding heterogeneity void_frac, range over leave-one-out folds"; extra = ""
    else:
        lo, hi = rx(r"s ([\d.]+)-([\d.]+)", t); est = num(v); lab = f"{bold(bid)}: latent field amplitude s, range over leave-one-out folds"
        extra = "; ell folds " + re.search(r"ell ([\d.]+-[\d.]+)", t).group(1) + " um, also LOAO, unchanged"
    add(bid, lab, "loao", "nuc", est, lo, hi, "row text",
        notes="RECLASSIFIED area_R2->loao: 'pieghe lascia-fuori' are leave-one-ROI-out folds" + extra)

# ---------- S1.6: nuclear size per type
s6 = pd.read_csv(R + "S1.6/S1.6_real_refs_roi_type.csv"); s6 = s6[s6.level == "type"]
q6 = pd.read_csv(R + "S1.6/S1.6_real_q90_roi.csv")
for bid in sel[sel.session == "S1.6"].id:
    a = bio.loc[bid, "archetype"]; ty = bold(bid); t = bio.loc[bid, "interval"]
    d = s6[(s6.archetype == a) & (s6.unit == ty)]; dq = q6[(q6.archetype == a) & (q6.cell_type == ty)]
    v = d.median_area_nuc.values
    add(bid, f"{ty}: median nuclear area per ROI, min-max", "range", "nuc", num(re.match(r"[\d.]+", bio.loc[bid, "value"]).group()), v.min(), v.max(),
        "results/S1.6/S1.6_real_refs_roi_type.csv:median_area_nuc", v, "S1.6_real_refs_roi_type.csv:median_area_nuc",
        text_lohi=rx(r"area mediana fra ROI ([\d.]+)-([\d.]+)", t), tol=0.006, notes=f"file mean of ROIs {v.mean():.2f}")
    m = re.search(r"nucleo/territorio ([\d.]+) \(([\d.]+)-([\d.]+)\)", t); v = d.median_nt.values
    add(bid, f"{ty}: median nucleus/territory ratio per ROI, min-max", "range", "nuc", num(m.group(1)), v.min(), v.max(),
        "results/S1.6/S1.6_real_refs_roi_type.csv:median_nt", v, "S1.6_real_refs_roi_type.csv:median_nt",
        text_lohi=(num(m.group(2)), num(m.group(3))), tol=0.00006)
    m = re.search(r"percentile ([\d.]+) \(([\d.]+)-([\d.]+)\)", t); v = dq.q90_an.values
    add(bid, f"{ty}: 90th percentile of nuclear area per ROI, min-max", "range", "nuc", num(m.group(1)), v.min(), v.max(),
        "results/S1.6/S1.6_real_q90_roi.csv:q90_an", v, "S1.6_real_q90_roi.csv:q90_an",
        text_lohi=(num(m.group(2)), num(m.group(3))), tol=0.006)

# ---------- R4
r4 = pd.read_csv(R + "R4/R4_refs.csv"); r4r = pd.read_csv(R + "R4/R4_roi.csv"); r4t = pd.read_csv(R + "R4/R4_types.csv")
ph1 = pd.read_csv(R + "R4/R4_PH1.csv"); ph1r = pd.read_csv(R + "R4/R4_PH1_roi.csv"); r4p = pd.read_csv(R + "R4/R4_predictions.csv")
par = pd.read_csv(R + "R4/R4_PH4_parity.csv")
def r4add(bid, q, lab, qt, s25, pattern):
    a = bio.loc[bid, "archetype"]; r = r4[(r4.archetype == a) & (r4.quantity == q)].iloc[0]
    pa = r4r[r4r.archetype == a][q].values if q in r4r.columns else None
    add(bid, lab, "ci", qt, r.estimate, r.lo, r.hi, f"results/R4/R4_refs.csv:{q} ({r.method})", pa, f"R4_roi.csv:{q}" if pa is not None else "",
        s25=s25, text_lohi=rx(pattern, bio.loc[bid, "interval"]),
        notes=("" if pa is not None else "per-ROI values not in R4_roi.csv"))
for bid in sel[sel.session == "R4"].id:
    met = bio.loc[bid, "metric"]; a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]
    if met.startswith("UMI per bin a 2"):
        r4add(bid, "mean_umi2", "mean UMI per 2 um bin in tissue", "bin", True, r"media [\d.]+ \[([\d.]+), ([\d.]+)\]")
        r4add(bid, "median_umi2", "median UMI per 2 um bin in tissue", "bin", False, r"mediana [\d.]+ \[([\d.]+), ([\d.]+)\]")
        r4add(bid, "frac_empty2", "fraction of empty 2 um bins in tissue", "bin", True, r"vuoti [\d.]+ \[([\d.]+), ([\d.]+)\]")
    elif met.startswith("UMI per bin a 8"):
        r4add(bid, "mean_umi8", "mean UMI per 8 um bin in tissue", "bin", True, r"media [\d.]+ \[([\d.]+), ([\d.]+)\]")
        r4add(bid, "median_umi8", "median UMI per 8 um bin in tissue", "bin", False, r"mediana [\d.]+ \[([\d.]+), ([\d.]+)\]")
        r4add(bid, "frac_empty8", "fraction of empty 8 um bins in tissue", "bin", True, r"vuoti [\d.]+ \[([\d.]+), ([\d.]+)\]")
    elif met.startswith("geni rilevati"):
        r4add(bid, "mean_genes2", "genes detected per 2 um bin", "bin", True, r"2 \u00b5m [\d.]+ \[([\d.]+), ([\d.]+)\]")
        r4add(bid, "mean_genes8", "genes detected per 8 um bin", "bin", True, r"8 \u00b5m [\d.]+ \[([\d.]+), ([\d.]+)\]")
    elif met.startswith("rapporto fra UMI"):
        r = ph1[(ph1.archetype == a) & (ph1.quantity == "cellweighted_star")].iloc[0]; v = ph1r[ph1r.archetype == a].cellweighted_star.values
        add(bid, "nucleus/ring UMI-per-bin ratio, cell-weighted, divided by the 20 um shifted-mask ratio (Step 2 reference)", "ci", "bin",
            r.estimate, r.lo, r.hi, "results/R4/R4_PH1.csv:cellweighted_star", v, "R4_PH1_roi.csv:cellweighted_star", s25=True,
            text_lohi=rx(r"^\[([\d.]+), ([\d.]+)\]", t))
        r = r4p[(r4p.archetype == a) & (r4p.id == "P-1")].iloc[0]; v = r4r[r4r.archetype == a].rho_star.values
        add(bid, "nucleus/ring ratio, pre-registered pooled-bin version (P-1)", "ci", "bin", r.estimate, r.lo, r.hi,
            "results/R4/R4_predictions.csv:P-1 estimate,lo,hi", v, "R4_roi.csv:rho_star",
            text_lohi=rx(r"bin raggruppati\) [\d.]+ \[([\d.]+), ([\d.]+)\]", t), notes=f"mean of R4_roi rho_star {np.mean(v):.4f}")
    elif met.startswith("UMI per cellula (nucleo"):
        r4add(bid, "cell_umi_median", "median UMI per cell (nucleus + 8 um ring), two-stage bootstrap (ROI, cells)", "nuc", False, r"mediana [\d.]+ \[([\d.]+), ([\d.]+)\]")
        r4add(bid, "cell_umi_cv", "CV of UMI per cell, two-stage bootstrap (ROI, cells)", "nuc", False, r"CV [\d.]+ \[([\d.]+), ([\d.]+)\]")
        r4add(bid, "cell_umi_per_bin_median", "median UMI per bin of the cell, two-stage bootstrap (ROI, cells)", "nuc", False, r"della cellula [\d.]+ \[([\d.]+), ([\d.]+)\]")
    elif met.startswith("UMI per cellula per tipo"):
        ty = bold(bid); r = r4t[(r4t.archetype == a) & (r4t.cell_type == ty)].iloc[0]
        add(bid, f"{ty}: median UMI per cell", "ci", "nuc", r.median_umi_cell, r.lo, r.hi, "results/R4/R4_types.csv:median_umi_cell,lo,hi",
            text_lohi=rx(r"^\[([\d.]+), ([\d.]+)\]", t), notes="per-ROI values not in file")
        add(bid, f"{ty}: median UMI per bin of the cell", "ci", "nuc", r.median_umi_per_bin, r.lo_per_bin, r.hi_per_bin,
            "results/R4/R4_types.csv:median_umi_per_bin,lo_per_bin,hi_per_bin", text_lohi=rx(r"per bin [\d.]+ \[([\d.]+), ([\d.]+)\]", t),
            notes="per-cell quantity, qtype nuc; per-ROI values not in file")
    elif met.startswith("alternanza"):
        d = par[par.unit.str.startswith(a + "_")]; e1, e2 = [num(x) for x in bio.loc[bid, "value"].split(" / ")]
        m = re.search(r"colonne ([\d.]+)-([\d.]+), righe ([\d.]+)-([\d.]+)", t)
        v = d.col_odd_over_even.values
        add(bid, "2 um capture alternation, odd/even columns: min-max across ROIs (est = section value)", "range", "bin", e1, v.min(), v.max(),
            "results/R4/R4_PH4_parity.csv:col_odd_over_even (ROI units)", v, "R4_PH4_parity.csv:col_odd_over_even", text_lohi=(num(m.group(1)), num(m.group(2))))
        v = d.row_odd_over_even.values
        add(bid, "2 um capture alternation, odd/even rows: min-max across ROIs (est = section value)", "range", "bin", e2, v.min(), v.max(),
            "results/R4/R4_PH4_parity.csv:row_odd_over_even (ROI units)", v, "R4_PH4_parity.csv:row_odd_over_even", text_lohi=(num(m.group(3)), num(m.group(4))))
    else:
        print("UNHANDLED R4", bid, met[:60])

# ---------- R5
ph2 = pd.read_csv(R + "R5/R5_PH2_inner_summary.csv"); ph2["arch"] = ph2.tag.str.split("_").str[0]
pn = pd.read_csv(R + "R5/R5_predictions_nc.csv")
for bid in sel[sel.session == "R5"].id:
    a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]
    if "nucleo dei bordi" in bio.loc[bid, "metric"]:
        r = pn[(pn.archetype == a) & (pn.prediction == "P-5") & pn.variant.str.startswith("as written")].iloc[0]
        add(bid, "nucleus/ring ratio before diffusion with the edge kernel (P-5), bootstrap over ROIs x 20 draws", "ci", "bin", r.estimate, r.lo, r.hi,
            "results/R5/R5_predictions_nc.csv:P-5 as written estimate,lo,hi", text_lohi=rx(r"^\[([\d.]+), ([\d.]+)\]", t),
            notes="two-stage: ROIs resampled (area) and simulation draws; treated as area")
        continue
    d = ph2[(ph2.arch == a) & ph2.sigma_in_star.notna()]
    v = d.sigma_in_star.values
    add(bid, f"in-tissue gaussian blur sigma_in, min-max across ROIs with a solution ({len(v)}/5)", "range", "bin", None,
        v.min(), v.max(), "results/R5/R5_PH2_inner_summary.csv:sigma_in_star", v if len(v) >= 2 else None, "R5_PH2_inner_summary.csv:sigma_in_star",
        text_lohi=rx(r"^([\d.]+)-([\d.]+)", bio.loc[bid, "value"]), tol=0.006,
        notes=("no point estimate (row value is the range)" + ("; single solved ROI" if len(v) == 1 else "")))
    v = d.r_at_star.values
    add(bid, f"nucleus/ring ratio before diffusion at sigma_in, min-max across solved ROIs ({len(v)}/5)", "range", "bin", None,
        v.min(), v.max(), "results/R5/R5_PH2_inner_summary.csv:r_at_star", v if len(v) >= 2 else None, "R5_PH2_inner_summary.csv:r_at_star",
        text_lohi=rx(r"diffusione ([\d.]+)-([\d.]+)", t), tol=0.006, notes="no point estimate")

# ---------- S2.0
ir = pd.read_csv(R + "S2.0/S2.0_inner_roi.csv"); ir = ir[ir.rule == "C2"]
isum = pd.read_csv(R + "S2.0/S2.0_inner_summary.csv"); isum = isum[isum.rule == "C2"].set_index("archetype")
dsp = pd.read_csv(R + "S2.0/S2.0_PH23_disp.csv").set_index("archetype")
for bid in sel[sel.session == "S2.0"].id:
    a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]
    if bio.loc[bid, "metric"].startswith("rapporto nucleo/citoplasma"):
        d = ir[ir.archetype == a]; v = d.r.values
        add(bid, "nucleus/cytoplasm ratio before displacement (rule C2), min-max across ROIs", "range", "bin", float(np.median(v)), v.min(), v.max(),
            "results/S2.0/S2.0_inner_roi.csv:r (rule C2)", v, "S2.0_inner_roi.csv:r", text_lohi=rx(r"ROI ([\d.]+)-([\d.]+)", t), tol=0.006)
        r = isum.loc[a]; v = d.ratio_far.values
        add(bid, "internal check predicted/observed without bins within 2 um of the nucleus edge, bootstrap over ROIs", "ci", "bin",
            r.ratio_mean, r.lo, r.hi, "results/S2.0/S2.0_inner_summary.csv:ratio_mean,lo,hi (rule C2)", v, "S2.0_inner_roi.csv:ratio_far",
            text_lohi=rx(r"interno [\d.]+ \[([\d.]+), ([\d.]+)\]", t), notes=f"recorded verdict vs [0.909, 1.1]: {r.verdict}")
    else:
        r = dsp.loc[a]
        add(bid, "extra multiplicative variance v inside the nucleus, bootstrap over ROIs", "ci", "bin", r.v, r.v_lo, r.v_hi,
            "results/S2.0/S2.0_PH23_disp.csv:v,v_lo,v_hi", text_lohi=rx(r"^\[([\d.]+), ([\d.]+)\]", t), notes="per-ROI values not in file")
        add(bid, "extra residual correlation between adjacent 2 um bins, bootstrap over ROIs", "ci", "bin", r.exc_lag1, r.exc_lag1_lo, r.exc_lag1_hi,
            "results/S2.0/S2.0_PH23_disp.csv:exc_lag1", text_lohi=rx(r"a 2 \u00b5m [\d.]+ \[([\d.]+), ([\d.]+)\]", t), notes="per-ROI values not in file")
        add(bid, "extra residual correlation at 4 um, bootstrap over ROIs", "ci", "bin", r.exc_lag2, r.exc_lag2_lo, r.exc_lag2_hi,
            "results/S2.0/S2.0_PH23_disp.csv:exc_lag2", text_lohi=rx(r"a 4 \u00b5m -?[\d.]+ \[(-?[\d.]+), (-?[\d.]+)\]", t), notes="per-ROI values not in file")

# ---------- S2.2
e1 = pd.read_csv(R + "S2.2/S2.2_E1_real.csv"); df_ = pd.read_csv(R + "S2.2/S2.2_depth_factor.csv").set_index("archetype")
dpr = pd.read_csv(R + "S2.2/S2.2_depth_per_roi.csv")
for bid in sel[sel.session == "S2.2"].id:
    a = bio.loc[bid, "archetype"]; t = bio.loc[bid, "interval"]; met = bio.loc[bid, "metric"]
    if met.startswith("fattore di profondit"):
        r = df_.loc[a]; v = dpr[dpr.archetype == a].f.values
        add(bid, "simulator depth factor, bootstrap over the 5 ROIs", "ci", "bin", r.f, r.lo, r.hi, "results/S2.2/S2.2_depth_factor.csv:f,lo,hi",
            v, "S2.2_depth_per_roi.csv:f", text_lohi=rx(r"5 ROI: ([\d.]+)\u2013([\d.]+)", t))
        add(bid, "simulator depth factor, min-max across ROIs", "range", "bin", r.f, r.f_min, r.f_max, "results/S2.2/S2.2_depth_factor.csv:f_min,f_max",
            v, "S2.2_depth_per_roi.csv:f", text_lohi=rx(r"min\u2013max ([\d.]+)\u2013([\d.]+)", t))
    else:
        if met.startswith("somiglianza"):
            col, lab = "m3", "mean cosine between adjacent 8 um blocks in typed-nucleus squares"
        else:
            col, lab = "m4_ratio", "same-cell / different-cell cosine ratio of squares (Daniele's measure)"
        v = e1[e1.archetype == a][col].values
        add(bid, lab + ", min-max across 5 ROI", "range", "bin", num(bio.loc[bid, "value"]), v.min(), v.max(), f"results/S2.2/S2.2_E1_real.csv:{col}",
            v, f"S2.2_E1_real.csv:{col}", text_lohi=rx(r"5 ROI: ([\d.]+)\u2013([\d.]+)", t), notes="row note: overlapping ROIs (BL-215), interval optimistic")

out = pd.DataFrame(rows)
done = set(out.id); miss = [i for i in sel.id if i not in done]
print("components", len(out), "ids", len(done), "missing", miss)
print("multi-archetype rows:", [i for i in done if not re.fullmatch(r"A\d", bio.loc[i, "archetype"])])
ck = pd.DataFrame(checks, columns=["id", "metric", "lo", "hi", "text_lo", "text_hi", "ok"])
print("text checks", int(ck.ok.notna().sum()), "mismatch", int((ck.ok == False).sum()), "no text pattern (file only)", int(ck.ok.isna().sum()))
pd.set_option("display.width", 250)
print(ck[ck.ok != True][["id", "metric", "lo", "hi", "text_lo", "text_hi", "ok"]].to_string(max_colwidth=45))
out.to_csv(W + "R7_Mc_refs.csv", index=False); ck.to_csv(W + "R7_Mc_textcheck.csv", index=False)
