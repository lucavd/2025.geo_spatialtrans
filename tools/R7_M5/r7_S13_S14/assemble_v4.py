
import sys, re, math, json
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
REPO = "/home/user/2025.geo_spatialtrans"; RES = REPO + "/results"; T = "/tmp/r7_S13_S14"
sys.path.insert(0, REPO + "/tools")
try:
    from R4_common import verdict
    VERDICT_SRC = "imported from tools/R4_common.py"
except Exception as e:
    VERDICT_SRC = "IMPORT FAILED (%s); identical copy used" % e
    def verdict(lo, hi, a, b):
        if not (np.isfinite(lo) and np.isfinite(hi)): return "non valutabile"
        if lo >= a and hi <= b: return "confermata"
        if hi < a or lo > b: return "smentita"
        return "indeterminata"
KT = {'A1': {'kappa_boot_bin': 1.856554021148876, 'kappa_boot_nuc': 1.8540743155566048, 'kappa_boot_box': 1.8374341234262737, 'kappa_se_bin': 1.5860264117842129, 'kappa_se_nuc': 1.583085172267185, 'kappa_se_box': 1.5695075595462076, 'kappa_range_bin': 1.0327069651572496, 'kappa_range_nuc': 1.0325837249509415, 'kappa_range_box': 1.0294979967840427}, 'A2': {'kappa_boot_bin': 1.6552632135756324, 'kappa_boot_nuc': 1.666445140393112, 'kappa_boot_box': 1.6621336397193407, 'kappa_se_bin': 1.4133471972021925, 'kappa_se_nuc': 1.4223516561672318, 'kappa_se_box': 1.42038754038451, 'kappa_range_bin': 1.0041579432733128, 'kappa_range_nuc': 1.0049188480678313, 'kappa_range_box': 1.0041629918993775}, 'A3': {'kappa_boot_bin': 3.192814308027333, 'kappa_boot_nuc': 2.817212095699824, 'kappa_boot_box': 3.192814308027333, 'kappa_se_bin': 2.739570848619368, 'kappa_se_nuc': 2.416476005411261, 'kappa_se_box': 2.739570848619368, 'kappa_range_bin': 1.2534342414882282, 'kappa_range_nuc': 1.1880128870565405, 'kappa_range_box': 1.2534342414882282}, 'A4': {'kappa_boot_bin': 2.570964962130436, 'kappa_boot_nuc': 2.965471527739436, 'kappa_boot_box': 2.965471527739436, 'kappa_se_bin': 2.2061061528365347, 'kappa_se_nuc': 2.549470325025491, 'kappa_se_box': 2.549470325025491, 'kappa_range_bin': 1.1472439207231864, 'kappa_range_nuc': 1.216180536882715, 'kappa_range_box': 1.216180536882715}, 'A5': {'kappa_boot_bin': 1.8182864626070927, 'kappa_boot_nuc': 3.865405212182495, 'kappa_boot_box': 3.865405212182495, 'kappa_se_bin': 1.5551385711316092, 'kappa_se_nuc': 3.324402032276403, 'kappa_se_box': 3.324402032276403, 'kappa_range_bin': 1.02897639105081, 'kappa_range_nuc': 1.372253355074593, 'kappa_range_box': 1.372253355074593}, 'A6': {'kappa_boot_bin': 1.624284044489665, 'kappa_boot_nuc': 1.6268777698741623, 'kappa_boot_box': 1.6252827285381342, 'kappa_se_bin': 1.3873110269132694, 'kappa_se_nuc': 1.386960324713831, 'kappa_se_box': 1.386430485194469, 'kappa_range_bin': 1.0, 'kappa_range_nuc': 1.0, 'kappa_range_box': 1.0}, 'sim_5': {'kappa_boot_bin': 1.6324309910704722, 'kappa_boot_nuc': 1.6324309910704722, 'kappa_boot_box': 1.6324309910704722, 'kappa_se_bin': 1.390530688688309, 'kappa_se_nuc': 1.390530688688309, 'kappa_se_box': 1.390530688688309, 'kappa_range_bin': 1.0, 'kappa_range_nuc': 1.0, 'kappa_range_box': 1.0}, 'window_4': {'kappa_boot_bin': 1.979389361922006, 'kappa_boot_nuc': 1.979389361922006, 'kappa_boot_box': 1.979389361922006, 'kappa_se_bin': 1.5943981227415345, 'kappa_se_nuc': 1.5943981227415345, 'kappa_se_box': 1.5943981227415345, 'kappa_range_bin': 1.0, 'kappa_range_nuc': 1.0, 'kappa_range_box': 1.0}}
GEO = {"nuc": {"boot": 1.751, "se": 1.496, "range": 1.020}, "bin": {"boot": 1.760, "se": 1.503, "range": 1.019}, "box": {"boot": 1.743, "se": 1.491, "range": 1.018}}
# compatibility dict for notes/aggregates: r<t> = kappa_range_<t>, <t> = kappa_se_<t>
K = {a_: {"rnuc": KT[a_]["kappa_range_nuc"], "rbin": KT[a_]["kappa_range_bin"], "rbox": KT[a_]["kappa_range_box"],
          "nuc": KT[a_]["kappa_se_nuc"], "bin": KT[a_]["kappa_se_bin"], "box": KT[a_]["kappa_se_box"]} for a_ in ["A1","A2","A3","A4","A5","A6"]}
KWB = KT["window_4"]["kappa_boot_nuc"]; KWS = KT["window_4"]["kappa_se_nuc"]; KSIM = KT["sim_5"]["kappa_se_nuc"]
KW = KWB
STRESS = 1.343
INF = float("inf")
ARCHS = ["A1","A2","A3","A4","A5","A6"]
CORR_IC = {"se_m", "boot", "range", "window_one_stage"}
COLS = ["step","prediction_id","archetype","quantity","est","lo","hi","a","b","est_unit","margin_unit","qtype","interval_construction","side","verdict_recorded","reproduces",
        "lo_new","hi_new","a_new","b_new","kappa_est","kappa_margin","verdict_new","verdict_new_stress","verdict_geom_A5","transition","touches_conclusion","notes"]
UNIT_IC = {"range": "range", "loao": "loao", "seed": "seed", "fixed": "fixed", "other": "other", "window": "window_one_stage", "area": "se_m"}
ROWS = []
def add(**kw):
    r = {c: kw.get(c, np.nan) for c in COLS}
    for c in ["lo_new","hi_new","a_new","b_new"]:
        if c not in kw: r[c] = r[c.replace("_new","")]
    for c in ["kappa_est","kappa_margin"]:
        if c not in kw: r[c] = 1.0
    if "verdict_new" not in kw: r["verdict_new"] = r["verdict_recorded"]
    if "verdict_new_stress" not in kw: r["verdict_new_stress"] = r["verdict_new"]
    if r["archetype"] != "A5": r["verdict_geom_A5"] = ""
    elif "verdict_geom_A5" not in kw or pd.isna(kw.get("verdict_geom_A5")): r["verdict_geom_A5"] = r["verdict_new"]
    ice = kw.get("_ice", UNIT_IC.get(r["est_unit"], str(r["est_unit"]))); icm = kw.get("_icm", UNIT_IC.get(r["margin_unit"], str(r["margin_unit"])))
    sde = kw.get("_sde", "real" if ice in CORR_IC else "-"); sdm = kw.get("_sdm", "real" if icm in CORR_IC else "-")
    r["interval_construction"] = "est=%s; margin=%s" % (ice, icm); r["side"] = "est=%s; margin=%s" % (sde, sdm)
    aff = kw.get("_aff", (ice in CORR_IC) or (icm in CORR_IC))
    if "transition" not in kw:
        r["transition"] = ("unchanged" if str(r["verdict_recorded"]) == str(r["verdict_new"]) else "%s->%s" % (r["verdict_recorded"], r["verdict_new"])) if aff else "unaffected"
    elif aff and "->" in kw["transition"] and kw["transition"].split("->")[0].split()[-1] == kw["transition"].split("->")[1].split()[0]: r["transition"] = "unchanged"
    if aff and str(r["verdict_recorded"]) == str(r["verdict_new"]): r["transition"] = "unchanged"
    rp = str(kw.get("reproduces", "NA")); r["reproduces"] = {"True": "TRUE", "False": "FALSE"}.get(rp, "NA" if rp in ("NA", "nan") else rp)
    if "touches_conclusion" not in kw: r["touches_conclusion"] = "none"
    ROWS.append(r)
def add2(base, **kw):
    d = dict(base); d.update(kw); add(**d)
def comp(kind, A, q, m=2.0):
    """factor for one component: f (new), fs (stress), fg (geometric A5), ic, side, unit"""
    if kind == "se_real":
        f = KT[A]["kappa_se_" + q] * 2.0 / m; g = GEO[q]["se"] * 2.0 / m if A == "A5" else f
        return dict(f=f, fs=f * STRESS, fg=g, ic="se_m", side="real", unit="area")
    if kind == "boot_real":
        f = KT[A]["kappa_boot_" + q]; g = GEO[q]["boot"] if A == "A5" else f
        return dict(f=f, fs=f * STRESS, fg=g, ic="boot", side="real", unit="area")
    if kind == "se_sim":
        f = KT["sim_5"]["kappa_se_" + q] * 2.0 / m
        return dict(f=f, fs=f * STRESS, fg=f, ic="se_m", side="simulated", unit="area")
    if kind == "range":
        f = KT[A]["kappa_range_" + q]; g = GEO[q]["range"] if A == "A5" else f
        return dict(f=f, fs=f, fg=g, ic="range", side="real", unit="range")
    if kind == "win_boot":
        f = KT["window_4"]["kappa_boot_" + q]
        return dict(f=f, fs=f * STRESS, fg=f, ic="window_one_stage", side="real", unit="window")
    if kind == "win_se":
        f = KT["window_4"]["kappa_se_" + q] * 2.0 / m
        return dict(f=f, fs=f * STRESS, fg=f, ic="window_one_stage", side="real", unit="window")
    if kind == "win_two":
        return dict(f=1.0, fs=1.0, fg=1.0, ic="window_two_stage", side="real (not corrected)", unit="window")
    return dict(f=1.0, fs=1.0, fg=1.0, ic=kind, side="-", unit=kind)
def resc(c, lo, hi, f):
    lo2 = lo if not np.isfinite(lo) else c - f * (c - lo)
    hi2 = hi if not np.isfinite(hi) else c + f * (hi - c)
    return lo2, hi2
def judge(est, lo, hi, a, b, ek, mk, q, A, mcen=None, clip01=False, m=2.0, mbase=None):
    """ek/mk = component kinds; mbase = (a0, b0) one-stage window interval to rescale for the margin (rule d)."""
    E = comp(ek, A, q, m); M = comp(mk, A, q)
    ce = est if np.isfinite(est) else (lo + hi) / 2
    a0, b0 = (mbase if mbase is not None else (a, b))
    cm = (a0 + b0) / 2 if mk == "range" else (mcen if mcen is not None else (a0 + b0) / 2)
    def go(f1, f2):
        l, h = resc(ce, lo, hi, f1)
        if clip01: l, h = max(0.0, l), min(1.0, h)
        if M["ic"] in CORR_IC:
            if mk == "range" and not (np.isfinite(a0) and np.isfinite(b0)): raise ValueError("range margin needs finite ends")
            aa, bb = resc(cm, a0, b0, f2)
        else: aa, bb = a, b
        return l, h, aa, bb, verdict(l, h, aa, bb)
    l, h, aa, bb, v = go(E["f"], M["f"]); _, _, _, _, vs = go(E["fs"], M["fs"])
    out = dict(lo_new=l, hi_new=h, a_new=aa, b_new=bb, verdict_new=v, verdict_new_stress=vs, kappa_est=E["f"], kappa_margin=M["f"],
               est_unit=E["unit"], margin_unit=M["unit"], _ice=E["ic"], _icm=M["ic"], _sde=E["side"], _sdm=M["side"])
    if A == "A5":
        lg, hg, ag, bg, vg = go(E["fg"], M["fg"]); out["verdict_geom_A5"] = vg; out["_g"] = (lg, hg, ag, bg)
    return out
def num(x):
    try: return float(x)
    except Exception: return np.nan
def rd(p): return pd.read_csv(p, keep_default_na=True)
LOG = {"not_reproduced": [], "unresolved": [], "info": [VERDICT_SRC]}

# =========================== S1.3 ===========================
v = rd(RES + "/S1.3/S1.3_verdicts.csv")
roi = rd(RES + "/S1.3/S1.3_checkB_roi.csv"); pk = rd(RES + "/S1.3/S1.3_checkB_roi_pkgref.csv")
GM = {"G1": "CSR", "G2": "RSA", "G3": "RSArule"}
TOL = {"B1": 0.05, "B2": 0.10, "B3": 0.10, "B4": 0.05}
for _, r in v.iterrows():
    sec, pid, A, mdl, met = r["section"], str(r["id"]), str(r["archetype"]), str(r["model"]), str(r["metric"])
    base = dict(step="S1.3", prediction_id="%s %s %s" % (sec, pid, mdl), archetype=A, quantity=met, verdict_recorded=r["outcome"])
    pc = "prediction=%s, confirmed=%s" % (r["prediction"], r["confirmed"])
    if sec in ("B", "Bsens") and pid in TOL:
        tab = roi if sec == "B" else pk
        z = tab[(tab.archetype == A) & (tab.model == GM[mdl])]
        dev = z[pid].to_numpy(float); tol = TOL[pid]
        ins = (np.abs(dev) < tol) if pid == "B4" else (np.abs(dev) <= tol)
        n_in = int(ins.sum()); rec = "PASS" if n_in >= 4 else "FAIL"
        kr = K[A]["rnuc"]; mid = (dev.min() + dev.max()) / 2; d2 = mid + kr * (dev - mid)
        ins2 = (np.abs(d2) < tol) if pid == "B4" else (np.abs(d2) <= tol); n2 = int(ins2.sum())
        rep = (rec == r["outcome"])
        if not rep: LOG["not_reproduced"].append("S1.3 %s %s %s %s: recount %d/5 -> %s vs recorded %s" % (sec, pid, A, mdl, n_in, rec, r["outcome"]))
        add2(base, est=float(np.median(dev)), lo=float(dev.min()), hi=float(dev.max()), a=-tol, b=tol, est_unit="other", margin_unit="fixed", qtype="nuc",
            reproduces=str(rep), touches_conclusion="none (plausibility check of Voronoi territories; failures attributed upstream, BL-071)",
            notes="count rule: PASS if >=4/5 R2 areas (ROI) within fixed tolerance; recount %d/5; no interval across areas, unit 'other' -> unchanged by R7 rule. Sensitivity (not the rule): per-ROI deviations stretched around the midpoint of their range by kappa_range_nuc=%.4f -> %d/5 -> %s. %s" % (n_in, kr, n2, "PASS" if n2 >= 4 else "FAIL", pc))
    elif sec == "B" and pid.startswith("ordine"):
        col = "cv_sim" if pid == "ordine cv" else "cv_loc_sim"
        z = roi[roi.archetype == A].pivot(index="roi_id", columns="model", values=col)
        ok = (z["RSArule"] > np.minimum(z["RSA"], z["CSR"])) & (z["RSArule"] < np.maximum(z["RSA"], z["CSR"])) if A == "A4" else (z["RSArule"] < z["RSA"])
        n_in = int(ok.sum()); rec = "PASS" if n_in >= 4 else "FAIL"; rep = rec == r["outcome"]
        if not rep: LOG["not_reproduced"].append("S1.3 %s %s %s" % (pid, A, mdl))
        add2(base, est_unit="other", margin_unit="fixed", qtype="nuc", reproduces=str(rep), notes="count rule over 5 R2 areas (simulated ordering G3 vs G2/CSR), recount %d/5; %s" % (n_in, pc))
    elif sec == "B":
        add2(base, est=num(str(r["value"]).split(" ")[0]), est_unit="seed", margin_unit="fixed", qtype="nuc", notes="descriptive (INFO); noise = D between replicates 1 and 2 (seed)")
    elif sec == "CP":
        add2(base, est_unit="seed", margin_unit="fixed", qtype="nuc", notes="counterproof on synthetic/seeded or real-centroid geometry; counts, no interval across areas; %s" % pc)
    elif sec == "D":
        add2(base, est_unit="seed", margin_unit="fixed", qtype="nuc", notes="descriptive / prediction on seeded synthetic squares; %s" % pc)
    else:
        tc = "engine GEOS for tessellate_voronoi (Luca's choice, rule had no eligible engine)" if pid in ("scelta",) else "none"
        add2(base, est_unit="fixed", margin_unit="fixed", qtype="box", touches_conclusion=tc, notes="software gate / engine comparison (no biological interval); %s" % pc)

# =========================== S1.3b ===========================
v = rd(RES + "/S1.3b/S1.3b_verdicts.csv")
for _, r in v.iterrows():
    tc = "overlay variant V2 (GEOS snap-rounding 1e-6 um) adopted" if r["section"] == "DEC" else "none"
    add(step="S1.3b", prediction_id="%s %s %s" % (r["section"], r["id"], r["variant"]), archetype="A1-A6 (nulls)", quantity=str(r["metric"]),
        verdict_recorded=r["outcome"], est_unit="fixed", margin_unit="fixed", qtype="box", touches_conclusion=tc,
        notes="numerical/software check on simulated nulls and real tessellations; tolerance numeric, no area interval; prediction=%s confirmed=%s" % (r["prediction"], r["confirmed"]))

# =========================== S1.3c ===========================
ar = rd(RES + "/S1.3c/S1.3c_archetype.csv")
for _, r in ar.iterrows():
    A = r["archetype"]; e = float(r["E_median"]); rec = r["verdict"]
    rep = ((e < 0.05) == (rec == "nessuno"))
    if not rep: LOG["not_reproduced"].append("S1.3c E %s %s %s" % (r["variant"], r["s_um"], A))
    mid = (r["E_min"] + r["E_max"]) / 2; mx2 = mid + K[A]["rbin"] * (r["E_max"] - mid)
    add(step="S1.3c", prediction_id="E<0.05 (%s, k=%d, s=%d um)" % (r["variant"], r["k"], r["s_um"]), archetype=A, quantity="E = phi_obs - phi_null, extracellular fraction (median of 5 R2 ROIs)",
        est=e, lo=e, hi=e, a=-INF, b=0.05, est_unit="other", margin_unit="fixed", qtype="bin", verdict_recorded=rec, reproduces=str(rep),
        touches_conclusion="design decision 1 (no extracellular space inside regions) kept by default; BL-077 left open (step not closed)",
        notes="archetype verdict uses the median of the 5 R2 ROIs as a point (min-max %.4f-%.4f reported only descriptively) vs conventional threshold 0.05 -> unit 'other', unchanged. Sensitivity: max stretched by kappa_range_bin=%.4f -> %.4f (<0.05: %s)" % (r["E_min"], r["E_max"], K[A]["rbin"], mx2, mx2 < 0.05))
pr = rd(RES + "/S1.3c/S1.3c_predictions.csv")
prim = ar[(ar.variant == "primary") & (ar.k == 3)].set_index("archetype")["E_median"]
chk = {"P1": "non valutabile" if (prim["A4"] == prim.min() and (prim == prim.min()).sum() > 1) else ("PASS" if prim["A4"] == prim.min() and prim["A4"] < 0.05 else "FAIL"),
       "P2": "PASS" if prim["A3"] >= 0.05 else "FAIL", "P3": "PASS" if (prim["A1"] <= prim["A2"] < prim["A3"]) else "FAIL",
       "P4": "non valutabile", "P5": "PASS" if prim["A6"] < 0.05 else "FAIL"}
for _, r in pr.iterrows():
    rep = chk[r["id"]] == r["esito"]
    if not rep: LOG["not_reproduced"].append("S1.3c %s" % r["id"])
    add(step="S1.3c", prediction_id=r["id"], archetype="A1-A6", quantity=r["previsione"], verdict_recorded=r["esito"], est_unit="other", margin_unit="fixed", qtype="bin",
        reproduces=str(rep), notes="comparison of archetype medians (points) / threshold 0.05 convention; observed: %s" % r["osservato"])
cp = rd(RES + "/S1.3c/S1.3c_controprove.csv")
for _, r in cp.iterrows():
    add(step="S1.3c", prediction_id="%s (k=%s)" % (r["cp"], r["k"]), archetype=r["archetype"], quantity=str(r["value"])[:120], verdict_recorded=r["status"],
        est_unit="other", margin_unit="fixed", qtype="bin", notes="counterproof: per-ROI counts / window counts vs fixed criteria; no interval across areas")
cc = rd(RES + "/S1.3c/S1.3c_checkC.csv")
for _, r in cc.iterrows():
    add(step="S1.3c", prediction_id=str(r["check"]), archetype=str(r["unit"]).split("_")[0], quantity="%s: %s" % (r["unit"], str(r["value"])[:100]), verdict_recorded=r["status"],
        est_unit="fixed", margin_unit="fixed", qtype="bin", notes="check C (software/data gate, per ROI)")

# =========================== S1.3d ===========================
v = rd(RES + "/S1.3d/S1.3d_verdicts.csv")
for _, r in v.iterrows():
    tc = "overlay F1 for O1/O2/O3 (rule)" if str(r["id"]).startswith("P-F1") else "none"
    add(step="S1.3d", prediction_id="%s %s %s" % (r["section"], r["id"], r["cand"]), archetype="synthetic/nulls", quantity=str(r["metric"]), verdict_recorded=r["outcome"],
        est_unit="fixed", margin_unit="fixed", qtype="box", touches_conclusion=tc, notes="software/numerical check vs exact arbiter; prediction=%s confirmed=%s" % (r["prediction"], r["confirmed"]))
b1 = rd(RES + "/S1.3d/S1.3d_b1.csv")
for _, r in b1.iterrows():
    add(step="S1.3d", prediction_id="B-1 %s" % r["cand"], archetype="A1-A6", quantity="62 B verdicts of S1.3 unchanged per candidate", verdict_recorded="PASS" if r["pred_same"] == r["n_pred"] else "FAIL",
        est_unit="fixed", margin_unit="fixed", qtype="box", notes="invariance of S1.3 check-B outcomes (%d/%d)" % (r["pred_same"], r["n_pred"]))

# =========================== S1.4 / S1.4b (yardstick rules) ===========================
F = rd(T + "/folds_S14_S14b.csv")
PRIM = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "spaceranger", "A5": "cellpose_rgb", "A6": "spaceranger"}
def arch_v(o):
    o = list(o)
    return "PASS" if o.count("PASS") >= 4 else ("FAIL" if o.count("FAIL") >= 2 else "WARN")
def pat(o): return "".join(x[0] for x in o)
def yrule(step, A, variant, metric, method=None):
    z = F[(F.step == step) & (F.archetype == A) & (F.variant == variant) & (F.metric == metric) & (F.method == (method or PRIM[A]))].sort_values("fold")
    assert len(z) == 5, (step, A, variant, metric, len(z))
    return dict(old=arch_v(z.v_old), new=arch_v(z.v_new), alt=arch_v(z.v_alt), a5=arch_v(z.v_a5), p5=pat(z.v_a5), po=pat(z.v_old), pn=pat(z.v_new), pa=pat(z.v_alt),
                dmin=z.d.min(), dmax=z.d.max(), dmed=z.d.median(), b=z.ymax.mean(), bn=z.ymax_new.mean(), k=z.kappa.iloc[0])
YNOTE = ("estimate = per-fold distance of the held-out R2 area (leave-one-area-out, NOT corrected, flagged: held-out area shares up to 60%% of nuclei with training areas in A3/A4); "
         "margin = inter-ROI yardstick D(real k, real j) over the 4 other R2 areas, PASS<=median_j, WARN<=max_j (3-tier, archetype PASS if >=4/5 folds PASS, FAIL if >=2 FAIL; "
         "not the verdict() rule, mapped as PASS/WARN/FAIL); correction: range endpoints min_j,max_j stretched around their midpoint by kappa_range_nuc, median_j kept (interior point); "
         "a=0, b=mean over folds of max_j. folds %s -> %s; sensitivity with the whole yardstick vector stretched (median moves too): %s -> %s")
v = rd(RES + "/S1.4/S1.4_verdicts.csv")
MAP14 = {"B-1": "ks_nc", "B-2": "ks_an", "B-3": "d_sp"}
TC14 = "nucleus size model M2 (beta, sigma; D-S1.4.1) in use; step closed with reservations (BL-098)"
def b_from_rule(step, A, variant, metric):
    y = yrule(step, A, variant, metric)
    return y
for _, r in v.iterrows():
    sec, pid, A, met, val, out = r["section"], r["id"], str(r["archetype"]), str(r["metric"]), str(r["value"]), r["outcome"]
    base = dict(step="S1.4", prediction_id="%s %s" % (sec, pid), archetype=A, quantity=met, verdict_recorded=out)
    if sec == "B" and pid in MAP14:
        y = yrule("S1.4", A, "M2", MAP14[pid]); rep = (y["old"] == out) and ("fold: %s" % y["po"]) in met
        if not rep: LOG["not_reproduced"].append("S1.4 %s %s: %s (%s) vs %s" % (pid, A, y["old"], y["po"], out))
        add2(base, est=y["dmed"], lo=y["dmin"], hi=y["dmax"], a=0.0, b=y["b"], a_new=0.0, b_new=y["bn"], est_unit="loao", margin_unit="range", qtype="nuc",
            reproduces=str(rep), verdict_new=y["new"], verdict_new_stress=y["new"], verdict_geom_A5=y["a5"], kappa_est=1.0, kappa_margin=y["k"], touches_conclusion=TC14,
            notes="model M2, primary segmenter; " + YNOTE % (y["po"], y["pn"], y["pa"], y["alt"]) + "; stress unchanged (range factor not stressed)")
    elif sec == "B" and pid == "B-4":
        e = num(val); rec = "PASS" if e <= 0.05 else "FAIL"; rep = rec == out
        add2(base, est=e, lo=e, hi=e, a=0.0, b=0.05, est_unit="loao", margin_unit="fixed", qtype="nuc", reproduces=str(rep),
            touches_conclusion="led to BL-098 (containment vs real nuclei); superseded by K-2 clip model (S1.4c/S1.4d)", notes="mean clamped fraction over 5 LOO folds (point) vs 5%% (ROADMAP); loao flagged, not corrected")
    elif sec == "CP" and pid == "CP-2":
        b_ = num(val.split(" ")[0]); rec = "PASS" if abs(b_) <= 0.02 else "FAIL"
        add2(base, est=b_, est_unit="other", margin_unit="fixed", qtype="nuc", reproduces=str(rec == out), notes="beta with permuted A_T; OLS SE over cells (not areas); point vs |b|<=0.02")
    elif sec == "CP" and pid == "CP-1" and A != "A1-A6":
        z = F[(F.step == "S1.4") & (F.archetype == A) & (F.method == PRIM[A]) & (F.metric.isin(["ks_nc", "ks_an"]))]
        S = z.groupby(["variant", "fold"]).d.sum().unstack(0)
        n = int((S["M2"] < S["M0"]).sum()); rec = "PASS" if n >= 4 else "FAIL"
        add2(base, est_unit="loao", margin_unit="fixed", qtype="nuc", reproduces=str(rec == out and val.startswith("%d/5" % n)), touches_conclusion=TC14,
            notes="M2 beats M0 in %d/5 LOO folds (S = KS N/C + KS A_N); loao flagged, not corrected" % n)
    elif sec == "CP" and pid == "CP-1":
        add2(base, est_unit="loao", margin_unit="fixed", qtype="nuc", touches_conclusion=TC14, notes="6/6 archetypes; loao flagged")
    elif sec == "CP" and pid == "CP-3":
        meth = "spaceranger" if "spaceranger" in met else "stardist_he"
        y1 = yrule("S1.4", A, "M2", "ks_nc", meth); y2 = yrule("S1.4", A, "M2", "ks_an", meth)
        ro = "B-1 %s, B-2 %s" % (y1["old"], y2["old"]); rn = "B-1 %s, B-2 %s" % (y1["new"], y2["new"])
        rep = val.endswith("%s %s" % (y1["old"], y2["old"]))
        add2(base, verdict_recorded=ro, verdict_new=rn, est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(rep), kappa_margin=y1["k"],
            notes="INFO row (segmenter sensitivity, %s); recorded outcome INFO; " % meth + YNOTE % (y1["po"] + "/" + y2["po"], y1["pn"] + "/" + y2["pn"], y1["pa"] + "/" + y2["pa"], y1["alt"] + "/" + y2["alt"]))
    elif sec == "P" and pid in ("P2", "P3"):
        res = {}
        for lab in ("old", "new", "alt", "a5"):
            if pid == "P2":
                n = sum(yrule("S1.4", a_, "M0", "ks_nc")[lab] == "FAIL" for a_ in ARCHS); res[lab] = ("PASS" if n == 6 else "FAIL", n)
            else:
                n = sum((yrule("S1.4", a_, "M2", "ks_nc")[lab] == "PASS") and (yrule("S1.4", a_, "M2", "ks_an")[lab] == "PASS") for a_ in ARCHS); res[lab] = ("PASS" if n >= 4 else "FAIL", n)
        add2(base, est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(res["old"][0] == out), verdict_new=res["new"][0], verdict_new_stress=res["new"][0],
            kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS), touches_conclusion=TC14,
            notes="panel-wide count recomputed from per-tissue corrected archetype verdicts (each tissue its own kappa_range_nuc; kappa_margin = largest); count old %d, new %d, full-stretch %d; geometric-A5 sensitivity %d -> %s" % (res["old"][1], res["new"][1], res["alt"][1], res["a5"][1], res["a5"][0]))
    elif sec == "P":
        add2(base, est_unit="other", margin_unit="fixed", qtype="nuc", touches_conclusion=TC14 if pid == "P1" else "none", notes="point estimates vs fixed thresholds (beta in (0,1), clamp 5%, C8 detection count)")
    elif sec == "B":
        add2(base, est_unit="seed" if pid == "B-seed" else "other", margin_unit="fixed", qtype="nuc", notes="descriptive (INFO)")
    else:
        add2(base, est_unit="fixed", margin_unit="fixed", qtype="nuc", notes="check C (software gate)")

v = rd(RES + "/S1.4b/S1.4b_verdicts.csv")
MAP4b = {"B-1": "ks_nc", "B-2": "ks_an", "B-3": "d_sp", "B-7": "b7"}
TC4b = "decision (c) containment with declared truncation (rule tie 4 vs 4); superseded in use by the K-2 'clip' nucleus (S1.4c/S1.4d, DEC-003)"
def arch_tab(lab):
    t = {}
    for a_ in ARCHS:
        for var in ("V-a/V-b", "V-c"):
            for bid, m in MAP4b.items(): t[(a_, var, bid)] = yrule("S1.4b", a_, var, m)[lab]
    return t
AT = {lab: arch_tab(lab) for lab in ("old", "new", "alt", "a5")}
for _, r in v.iterrows():
    sec, pid, A, met, val, out = r["section"], r["id"], str(r["archetype"]), str(r["metric"]), str(r["value"]), r["outcome"]
    base = dict(step="S1.4b", prediction_id="%s %s" % (sec, pid), archetype=A, quantity=met, verdict_recorded=out)
    var = "V-a/V-b" if met.startswith("V-a/V-b") else ("V-c" if met.startswith("V-c") else None)
    if sec == "B" and pid in MAP4b and var:
        y = yrule("S1.4b", A, var, MAP4b[pid]); rep = (y["old"] == out) and ("fold %s" % y["po"]) in met
        if not rep: LOG["not_reproduced"].append("S1.4b %s %s %s: %s (%s) vs %s" % (pid, A, var, y["old"], y["po"], out))
        add2(base, prediction_id="B %s %s" % (pid, var), est=y["dmed"], lo=y["dmin"], hi=y["dmax"], a=0.0, b=y["b"], a_new=0.0, b_new=y["bn"], est_unit="loao", margin_unit="range", qtype="nuc",
            reproduces=str(rep), verdict_new=y["new"], verdict_new_stress=y["new"], verdict_geom_A5=y["a5"], kappa_margin=y["k"], touches_conclusion=TC4b,
            notes="distance = median over 10 seeds per fold; " + YNOTE % (y["po"], y["pn"], y["pa"], y["alt"]))
    elif sec == "B" and pid == "B-4":
        e = num(val.split(" ")[0]); rec = "PASS" if e <= 0.05 else "FAIL"
        add2(base, prediction_id="B B-4 %s" % var, est=e, lo=e, hi=e, a=0.0, b=0.05, est_unit="loao", margin_unit="fixed", qtype="nuc", reproduces=str(rec == out), touches_conclusion=TC4b,
            notes="mean truncated fraction over LOO folds (point) vs 5%; loao flagged")
    elif sec == "P" and pid in ("P3", "P4", "P6"):
        res = {}
        for lab in ("old", "new", "alt", "a5"):
            t = AT[lab]
            if pid == "P3":
                n = sum(t[(a_, "V-a/V-b", "B-1")] == "PASS" and t[(a_, "V-a/V-b", "B-2")] == "PASS" for a_ in ARCHS); res[lab] = ("PASS" if n >= 4 else "FAIL", "%d/6" % n)
            elif pid == "P4":
                n1 = sum(t[(a_, "V-a/V-b", "B-7")] in ("PASS", "WARN") for a_ in ARCHS); n2 = sum(t[(a_, "V-c", "B-7")] == "FAIL" for a_ in ARCHS[:5])
                res[lab] = ("PASS" if (n1 >= 5 and n2 == 5) else "FAIL", "V-a/b %d/6; V-c FAIL %d/5" % (n1, n2))
            else:
                ok = t[("A5", "V-a/V-b", "B-1")] != "PASS" and t[("A5", "V-c", "B-1")] != "PASS"; res[lab] = ("PASS" if ok else "FAIL", "%s/%s" % (t[("A5", "V-c", "B-1")], t[("A5", "V-a/V-b", "B-1")]))
        add2(base, est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(res["old"][0] == out), verdict_new=res["new"][0], verdict_new_stress=res["new"][0],
            kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS), touches_conclusion="none (prediction)",
            notes="recomputed from per-tissue corrected archetype verdicts (kappa_margin = largest kappa_range_nuc); old %s, new %s, full-stretch %s; geometric-A5 sensitivity %s -> %s" % (res["old"][1], res["new"][1], res["alt"][1], res["a5"][1], res["a5"][0]))
    elif sec == "D":
        res = {}
        for lab in ("old", "new", "alt", "a5"):
            t = AT[lab]; s_ab = sum(t[(a_, "V-a/V-b", b_)] == "PASS" for a_ in ARCHS for b_ in MAP4b); s_c = sum(t[(a_, "V-c", b_)] == "PASS" for a_ in ARCHS for b_ in MAP4b)
            ch = "(b) nucleo non contenuto" if s_ab > s_c else "(c) contenimento con troncamento dichiarato"   # eligible a FALSE, b TRUE, c TRUE (recorded)
            res[lab] = (ch, s_ab, s_c)
        rep = res["old"][0] == out and ("%d vs %d" % (res["old"][1], res["old"][2])) in val
        add2(base, est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(rep), verdict_new=res["new"][0], verdict_new_stress=res["new"][0],
            kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS), touches_conclusion=TC4b,
            notes="decision rule: PASS count {B-1,B-2,B-3,B-7} {a,b} vs c (tie -> c); eligibility a FALSE, b TRUE, c TRUE taken as recorded; scores old %d vs %d, new %d vs %d, full-stretch %d vs %d; geometric-A5 sensitivity %d vs %d -> %s" % (res["old"][1], res["old"][2], res["new"][1], res["new"][2], res["alt"][1], res["alt"][2], res["a5"][1], res["a5"][2], res["a5"][0]))
    elif sec == "CP":
        add2(base, est_unit="other", margin_unit="fixed", qtype="nuc", notes="Laguerre vs 20 permutations per ROI (count of ROIs); no interval across areas")
    elif sec == "P":
        add2(base, est_unit="other" if pid != "P2" else "loao", margin_unit="fixed", qtype="nuc", notes="point estimates vs fixed criteria")
    elif sec == "B":
        add2(base, est_unit="other", margin_unit="fixed", qtype="nuc", notes="descriptive (INFO)")
    else:
        add2(base, est_unit="fixed", margin_unit="fixed", qtype="nuc", notes="check C (software gate)")

# =========================== S1.4c ===========================
vo = rd(T + "/S14c_verdicts_oldnew.csv"); vrec = rd(RES + "/S1.4c/S1.4c_verdicts.csv")
TC4c = "nucleus model K-2 (clip + re-inflate to target area) chosen by rule, implemented as default nucleus_model='clip' in S1.4d (DEC-003), in use"
for _, r in vo.iterrows():
    A = r["archetype"]; rr = vrec[(vrec.archetype == A) & (vrec.variant == r["variant"]) & (vrec.metric == r["metric"])].iloc[0]
    rep = (rr["verdict"] == r["verdict_old"]) and int(rr["n_pass"]) == int(r["n_pass_old"]) and int(rr["n_fail"]) == int(r["n_fail_old"])
    if not rep: LOG["not_reproduced"].append("S1.4c %s %s %s" % (A, r["variant"], r["metric"]))
    add(step="S1.4c", prediction_id="B %s %s" % (r["metric"], r["variant"]), archetype=A, quantity="%s, %s (in score panel: %s)" % (r["metric"], r["variant"], r["in_panel_old"]),
        verdict_recorded=rr["verdict"], est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(rep), verdict_new=r["verdict_new"], verdict_new_stress=r["verdict_new"], verdict_geom_A5=r["verdict_a5"],
        kappa_margin=K[A]["rnuc"], touches_conclusion=TC4c if r["in_panel_old"] else "none (outside score panel)",
        notes="fold verdict on median of 10 seeds; P/W/F folds old %d/%d/%d, new %d/%d/%d, full-stretch %d/%d/%d (%s). " % (r["n_pass_old"], r["n_warn_old"], r["n_fail_old"], r["n_pass_new"], r["n_warn_new"], r["n_fail_new"], r["n_pass_alt"], r["n_warn_alt"], r["n_fail_alt"], r["verdict_alt"]) + "estimate leave-one-area-out (flagged, not corrected); yardstick endpoints stretched by kappa_range_nuc, median kept")
dec = rd(T + "/S14c_decision_oldnew.csv"); drec = rd(RES + "/S1.4c/S1.4c_decision.csv")
for V in ("K-1", "K-2"):
    o = dec[(dec.variant == V) & (dec.which == "old")].iloc[0]; n = dec[(dec.variant == V) & (dec.which == "new")].iloc[0]; al = dec[(dec.variant == V) & (dec.which == "alt")].iloc[0]; a5r = dec[(dec.variant == V) & (dec.which == "a5")].iloc[0]
    rr = drec[drec.variant == V].iloc[0]
    rep = (o["score"] == rr["score"]) and (o["choice"] == rr["choice"]) and abs(o["p_boot"] - rr["p_boot"]) < 1e-12
    add(step="S1.4c", prediction_id="Regola (decisione) %s" % V, archetype="A1-A6", quantity="score %s >= score K-c + 3 and P(boot over seeds) >= 0.9; choice" % V,
        verdict_recorded="choice %s (score %d vs K-c %d, p_boot %.3f)" % (rr["choice"], rr["score"], rr["score_Kc"], rr["p_boot"]),
        verdict_new="choice %s (score %d vs K-c %d, p_boot %.3f)" % (n["choice"], n["score"], n["score_Kc"], n["p_boot"]),
        verdict_new_stress="choice %s (score %d vs K-c %d, p_boot %.3f)" % (n["choice"], n["score"], n["score_Kc"], n["p_boot"]),
        est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(rep), kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS), touches_conclusion=TC4c,
        transition="choice %s->%s" % (rr["choice"], n["choice"]),
        notes="rule rerun in R with tools/S1.4c_decision.R (read-only) on corrected yardsticks, same bootstrap seed 20261015 (bootstrap over seeds: unit seed, not corrected); panel old %s, new %s; full-stretch sensitivity: choice %s, score %d vs %d, p_boot %.3f; geometric-A5 sensitivity: choice %s, score %d vs %d, p_boot %.3f" % (o["panel"], n["panel"], al["choice"], al["score"], al["score_Kc"], al["p_boot"], a5r["choice"], a5r["score"], a5r["score_Kc"], a5r["p_boot"]))
c2 = rd(T + "/S14c_cp2.csv"); c2r = rd(RES + "/S1.4c/S1.4c_cp2.csv")
for _, r in c2.iterrows():
    rr = c2r[(c2r.metric == r["metric"]) & (c2r.archetype == r["archetype"])].iloc[0]
    rep = (bool(rr["exceeds"]) == bool(r["exceeds_old"])) and int(rr["n_roi_above"]) == int(r["n_roi_above_old"])
    f = lambda x: "exceeds (>=3/5 ROI above max of yardstick)" if bool(x) else "does not exceed"
    add(step="S1.4c", prediction_id="CP-2 %s" % r["metric"], archetype=r["archetype"], quantity="segmenter distance vs max of inter-ROI yardstick (n ROI above: old %d, new %d)" % (r["n_roi_above_old"], r["n_roi_above_new"]),
        verdict_recorded=f(rr["exceeds"]), verdict_new=f(r["exceeds_new"]), verdict_new_stress=f(r["exceeds_new"]), est_unit="other", margin_unit="range", qtype="nuc", reproduces=str(rep),
        verdict_geom_A5=f(r["exceeds_a5"]), kappa_margin=r["kappa"], touches_conclusion="demotion of S-1/S-2/S-3 from the score panel (3/6 rule) -> S1.4c decision", notes="full-stretch sensitivity: %s" % f(r["exceeds_alt"]))
cp1 = rd(T + "/S14c_cp1_oldnew.csv").set_index("which"); cp1r = rd(RES + "/S1.4c/S1.4c_cp1.csv").iloc[0]
g = lambda w: "refuted (gain A5-A6 %d >= half of gain A1-A4 %d)" % (cp1.loc[w, "g56"], cp1.loc[w, "g14"]) if bool(cp1.loc[w, "refuted"]) else "not refuted (gain A5-A6 %d < half of %d)" % (cp1.loc[w, "g56"], cp1.loc[w, "g14"])
add(step="S1.4c", prediction_id="CP-1", archetype="A1-A6", quantity="specificity of K-2 gain over K-c to dense tissues (PASS counts, 7 metrics)", verdict_recorded=g("old"), verdict_new=g("new"), verdict_new_stress=g("new"),
    est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(int(cp1.loc["old", "g14"]) == int(cp1r["gain_A1_A4"]) and bool(cp1.loc["old", "refuted"]) == bool(cp1r["refuted"])),
    kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS), touches_conclusion=TC4c, notes="full-stretch sensitivity: " + g("alt") + "; geometric-A5 sensitivity: " + g("a5"))
pr = rd(RES + "/S1.4c/S1.4c_predictions.csv")
for _, r in pr.iterrows():
    pid = r["id"]; base = dict(step="S1.4c", prediction_id=pid, archetype="A1-A6", quantity=r["text"], verdict_recorded="confermata" if r["confirmed"] else "smentita")
    if pid == "P3":
        res = {w: all(vo[(vo.archetype == a_) & (vo.variant == "K-1") & (vo.metric == "B2")]["verdict_" + w].iloc[0] == "FAIL" for a_ in ARCHS[:4]) for w in ("old", "new", "alt")}
        add2(base, est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str(res["old"] == bool(r["confirmed"])), verdict_new="confermata" if res["new"] else "smentita",
            kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS[:4]), notes="K-1 B-2 verdicts A1-A4 recomputed; full-stretch: %s" % res["alt"])
    elif pid == "P4":
        o = dec[(dec.variant == "K-2") & (dec.which == "old")].iloc[0]; n = dec[(dec.variant == "K-2") & (dec.which == "new")].iloc[0]
        add2(base, est_unit="loao", margin_unit="range", qtype="nuc", reproduces=str((o["score"] >= o["score_Kc"] + 3) == bool(r["confirmed"])),
            verdict_new="confermata" if n["score"] >= n["score_Kc"] + 3 else "smentita", kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS), touches_conclusion=TC4c,
            notes="score K-2 %d vs K-c %d (new %d vs %d)" % (o["score"], o["score_Kc"], n["score"], n["score_Kc"]))
    else:
        add2(base, est_unit="other", margin_unit="fixed", qtype="nuc", notes="real or simulated point values vs fixed conventional thresholds; observed: %s" % str(r["observed"])[:150])
for _, r in rd(RES + "/S1.4c/S1.4c_b4.csv").iterrows():
    add(step="S1.4c", prediction_id="B-4 %s" % r["variant"], archetype=r["archetype"], quantity="truncated fraction (mean over LOO folds of seed medians)", est=r["trunc_all"], lo=r["trunc_all"], hi=r["trunc_all"], a=0.0, b=0.05,
        verdict_recorded=r["verdict"], est_unit="loao", margin_unit="fixed", qtype="nuc", reproduces=str((r["trunc_all"] <= 0.05) == (r["verdict"] == "PASS")), notes="gate reported outside score; loao flagged")
for _, r in rd(RES + "/S1.4c/S1.4c_cp0.csv").iterrows():
    add(step="S1.4c", prediction_id="CP-0", archetype=r["archetype"], quantity="AUC of contact premise on real neighbour pairs", est=r["auc"], lo=r["auc"], hi=r["auc"], a=0.75, b=1.0,
        verdict_recorded="AUC ok" if r["auc_ok"] else "AUC not ok", est_unit="other", margin_unit="fixed", qtype="nuc", reproduces=str((r["auc"] >= 0.75) == bool(r["auc_ok"])), notes="pairs pooled over 5 ROIs, point vs 0.75 convention")
for _, r in rd(RES + "/S1.4c/S1.4c_checks.csv").iterrows():
    add(step="S1.4c", prediction_id=str(r["check"]), archetype="A1-A6", quantity=str(r["what"]), verdict_recorded=r["verdict"], est_unit="fixed", margin_unit="fixed", qtype="nuc", notes="check C (software gate)")

# =========================== S1.4d ===========================
s1 = rd(RES + "/R2c/R2c_s1.csv"); s1 = s1[s1.set == "Luca"].set_index("archetype")["s1_primary"]
TC4d = {"P-1": "nucleus_model='clip' (K-2) default in derive_cell_geometry (DEC-003); step closed",
        "P-4": "choice of centres in the lymph node deferred to R3b (BL-130/BL-121)",
        "P-2": "per-tissue cap phi (nucleus_nc_max preset values 0.757/1/0.947/0.816/0.474/0.180)",
        "P-3": "DEC-005: no cap on nucleus area in colon tumour; size-model tail to be fixed with cell types (BL-138)",
        "CP-1": "validity of P-1 as discriminating test of K-2 vs scaled nucleus K-c"}
s1all = rd(RES + "/R2c/R2c_s1.csv"); s1all = s1all[s1all.set == "Luca"].set_index("archetype")
ONE = {a_: (float(s1all.loc[a_, "lo_primary"]), float(s1all.loc[a_, "hi_primary"])) for a_ in ARCHS}
def p14d(r, idlab, ek, mk, tc, note):
    A = r["archetype"]; rec = r["verdict"]
    if A == "S":
        add(step="S1.4d", prediction_id=idlab, archetype="S (section)", quantity=idlab, est=r["est"], lo=r["lo2"], hi=r["hi2"], a=r["int_lo"], b=r["int_hi"], verdict_recorded=rec,
            est_unit="seed", margin_unit="fixed", qtype="box", reproduces=str(verdict(r["lo2"], r["hi2"], r["int_lo"], r["int_hi"]) == rec), touches_conclusion="none (performance)", notes="3 repetitions; threshold from C-perf (Luca)")
        return
    ov = verdict(r["lo2"], r["hi2"], r["int_lo"], r["int_hi"]); rep = ov == rec
    if not rep: LOG["not_reproduced"].append("S1.4d %s %s" % (idlab, A))
    mcen = float(s1[A]) if mk == "win_boot" else None; mb = ONE[A] if mk == "win_boot" else None
    j = judge(r["est"], r["est"] - 2 * r["se"], r["est"] + 2 * r["se"], r["int_lo"], r["int_hi"], ek, mk, "nuc", A, mcen=mcen, clip01=True, mbase=mb)
    if mk == "win_boot":
        nt = ("estimate: mean over 10 seeds +- 2 SE (seed, unchanged), clipped to [0,1] (D-5); margin: recorded = two-stage interval of B-127..B-132 (NOT widened, rule d); new = one-stage percentile cluster bootstrap over the 4 windows "
              "(R2c_s1.csv lo_primary/hi_primary = %.4f-%.4f) widened by window kappa_boot %.4f around the manual value %.4f (stress x1.343); a_new<0 or b_new>1 left as computed" % (mb[0], mb[1], KWB, mcen))
    elif ek == "se_real":
        nt = "estimate: mean of 5 real R2 ROIs +- 2 SE across ROIs (simulated nuclei on real tessellations, real side) -> kappa_se_nuc %.4f x 2/2 (stress x1.343), clipped to [0,1] (D-5); margin fixed 0-0.05 (ROADMAP, Luca)" % K[A]["nuc"]
    else:
        nt = "seed/fixed components, unchanged"
    add(step="S1.4d", prediction_id=idlab, archetype=A, quantity=note, est=r["est"], lo=r["lo2"], hi=r["hi2"], a=r["int_lo"], b=r["int_hi"], qtype="nuc",
        verdict_recorded=rec, reproduces=str(rep), touches_conclusion=tc, notes=nt, **j)
pp = rd(RES + "/S1.4d/S1.4d_predictions.csv"); p3r = rd(RES + "/S1.4d/S1.4d_p3_roi.csv")
for _, r in pp.iterrows():
    pid = r["id"]
    if pid in ("P-1", "P-4"):
        p14d(r, pid, "seed", "win_boot", TC4d[pid], "fraction of simulated nuclei touching a neighbour, R2c windows" + (" (StarDist centres)" if pid == "P-4" else ""))
    elif pid == "P-2":
        p14d(r, pid, "se_real", "fixed", TC4d[pid], "fraction of nuclei truncated by cap phi, all cells")
    elif pid == "P-3":
        A = r["archetype"]; rv = p3r[p3r.archetype == A]["real_nc95"].to_numpy(float)
        rec = r["verdict"]; ov = verdict(r["lo2"], r["hi2"], r["int_lo"], r["int_hi"]); rep = ov == rec
        E = comp("se_real", A, "nuc"); M = comp("range", A, "nuc"); mid = (rv.min() + rv.max()) / 2
        bn = mid + M["f"] * (rv.max() - mid); bg = mid + M["fg"] * (rv.max() - mid)
        iv = lambda fe: (max(0.0, r["est"] - fe * 2 * r["se"]), min(1.0, r["est"] + fe * 2 * r["se"]))
        ln, hn = iv(E["f"]); ls, hs = iv(E["fs"]); lg, hg = iv(E["fg"])
        add(step="S1.4d", prediction_id=pid, archetype=A, quantity="fraction of cells with nucleus >= 0.95 territory", est=r["est"], lo=r["lo2"], hi=r["hi2"], a=r["int_lo"], b=r["int_hi"],
            lo_new=ln, hi_new=hn, a_new=0.0, b_new=bn, kappa_est=E["f"], kappa_margin=M["f"], est_unit="area", margin_unit="range", _ice="se_m", _icm="range (upper bound; a=0 fixed)", _sde="real", _sdm="real", _aff=True,
            verdict_new=verdict(ln, hn, 0.0, bn), verdict_new_stress=verdict(ls, hs, 0.0, bn), verdict_geom_A5=verdict(lg, hg, 0.0, bg) if A == "A5" else np.nan,
            qtype="nuc", verdict_recorded=rec, reproduces=str(rep), touches_conclusion=TC4d["P-3"] if A == "A2" else "none",
            notes="estimate mean of 5 real ROIs +- 2 SE (real side, kappa_se_nuc %.4f x 2/2, clipped [0,1]); margin [0, max of 5 real ROIs]: a=0 fixed by convention, b = range endpoint stretched around (min+max)/2=%.5f by kappa_range_nuc %.4f" % (E["f"], mid, M["f"]))
    else:
        p14d(r, pid, "seed", "fixed", "none", pid)
cpd = rd(RES + "/S1.4d/S1.4d_counterproof.csv")
for _, r in cpd.iterrows():
    lab = "CP-1" if r["id"] == "CP-1" else "P-4 (descrittiva)"
    p14d(r, lab, "seed", "win_boot", TC4d["CP-1"] if lab == "CP-1" else "none", "contact fraction with the old scaled nucleus K-c (counterproof)" if lab == "CP-1" else "contact fraction, StarDist centres, heart (descriptive)")
cpn = [x for x in ROWS if x["step"] == "S1.4d" and x["prediction_id"] == "CP-1" and x["archetype"] in ("A1", "A2", "A3", "A4")]
add(step="S1.4d", prediction_id="CP-1 (esito)", archetype="A1-A4", quantity="K-c confirmed or indeterminate in >=2 of 4 dense tissues -> P-1 not discriminating",
    verdict_recorded="superata" if sum(x["verdict_recorded"] != "smentita" for x in cpn) < 2 else "fallita", verdict_new="superata" if sum(x["verdict_new"] != "smentita" for x in cpn) < 2 else "fallita",
    verdict_new_stress="superata" if sum(x["verdict_new_stress"] != "smentita" for x in cpn) < 2 else "fallita", est_unit="seed", margin_unit="window", qtype="nuc", kappa_margin=KW,
    reproduces="True", touches_conclusion=TC4d["CP-1"], notes="aggregate of the 4 per-tissue CP-1 rows above; K-c smentita new %d/4, stress %d/4" % (sum(x["verdict_new"] == "smentita" for x in cpn), sum(x["verdict_new_stress"] == "smentita" for x in cpn)))
for _, r in rd(RES + "/S1.4d/S1.4d_gates.csv").iterrows():
    add(step="S1.4d", prediction_id=str(r["gate"]), archetype="A1-A6", quantity=str(r["check"]), verdict_recorded="PASS" if bool(r["pass"]) else "FAIL", est_unit="fixed", margin_unit="fixed", qtype="nuc",
        touches_conclusion="DEC-004 (G-3 failure = threshold error)" if not bool(r["pass"]) else "none", notes="gate (software)")

# =========================== S1.4e / S1.4f (shared) ===========================
def refs(path):
    x = rd(path); return x
R1 = rd(RES + "/S1.4e/S1.4e_refs.csv"); R1c = {r["archetype"]: r["est"] for _, r in R1[(R1.ref == "R-1") & (R1.cls == "lt6")].iterrows()}
KINDS = {"P-1": ("seed", "win_boot"), "P-2": ("seed", "win_two"), "P-3": ("se_sim", "range"), "P-4": ("se_sim", "range"), "P-5": ("se_sim", "range")}
QD = {"P-1": "fraction of focal nuclei touching a neighbour (R2c windows)", "P-2": "fraction of pairs <6 um in contact", "P-3": "alignment <5 um (share of elongated nuclei transverse to neighbour)",
      "P-4": "median nuclear eccentricity", "P-5": "alignment share by distance class"}
def pred_row(step, r, tc, extra_id=""):
    pid = r["id"]; A = r["archetype"]; ek, mk = KINDS[pid]
    lo, hi = r["est"] - 2 * r["se"], r["est"] + 2 * r["se"]; rec = r["verdict"]
    ov = verdict(lo, hi, r["lo"], r["hi"])
    rep = "NA" if rec == "descrittiva" else str(ov == rec)
    if rep == "False": LOG["not_reproduced"].append("%s %s %s %s" % (step, pid, A, r["arm"]))
    mcen = float(s1[A]) if pid == "P-1" else (float(R1c[A]) if pid == "P-2" else None)
    j = judge(r["est"], lo, hi, r["lo"], r["hi"], ek, mk, "nuc", A, mcen=mcen, mbase=ONE[A] if pid == "P-1" else None)
    dnote = ""
    if rec == "descrittiva":
        dnote = "; verdict() on corrected intervals would give %s (stress %s), kept descriptive as preregistered" % (j["verdict_new"], j["verdict_new_stress"])
        j["verdict_new"] = j["verdict_new_stress"] = "descrittiva"
        if "verdict_geom_A5" in j: j["verdict_geom_A5"] = "descrittiva"
    cls = ("" if ("cls" not in r or pd.isna(r.get("cls"))) else " " + str(r["cls"]))
    if pid == "P-1":
        nt = ("margin: recorded = two-stage interval B-127..B-132 (rounded in the table, NOT widened, rule d); new = one-stage percentile cluster bootstrap over the 4 windows (R2c_s1 %.4f-%.4f) widened by window kappa_boot %.4f around %.4f (stress x1.343)" % (ONE[A][0], ONE[A][1], KWB, mcen))
    elif pid == "P-2":
        nt = "margin R-1 = two-stage bootstrap (windows, then pairs); no one-stage interval exists -> kept unchanged and flagged (rule d)"
    else:
        nt = "margin %s = [min, max] of the 5 real R2 ROIs -> stretched around midpoint by kappa_range_nuc %.4f (not stressed)" % (r["source"], K[A]["rnuc"])
    ne = "estimate mean over 10 seeds +- 2 SE (seed, unchanged)" if ek == "seed" else "estimate mean of 5 SIMULATED ROI units (bench B = S1.5d RV units: own simulation and seed per ROI, no shared cells) +- 2 SE -> sim_5 kappa_se_nuc %.4f x 2/2 (stress x1.343)" % KSIM
    add(step=step, prediction_id="%s %s%s%s" % (pid, r["arm"], cls, extra_id), archetype=A, quantity=QD[pid] + cls, est=r["est"], lo=lo, hi=hi, a=r["lo"], b=r["hi"], qtype="nuc",
        verdict_recorded=rec, reproduces=rep, touches_conclusion=tc(pid, r["arm"]), notes=ne + "; " + nt + ("; recorded 'descrittiva' by preregistration (D-5)" + dnote if rec == "descrittiva" else ""), **j)
    return ROWS[-1]
def tc_e(pid, arm):
    if pid in ("P-1", "P-2"): return "rule N -> preset nucleus kept as disc K2 (DEC-026/DEC-027)"
    if pid == "P-3": return "ET over-aligns -> next attempt partial orientation (BL-171, S1.4f)"
    return "none (weak test by construction)"
ve = rd(RES + "/S1.4e/S1.4e_verdicts.csv"); E = {}
for _, r in ve.iterrows(): E[(r["id"], r["archetype"], r["arm"])] = pred_row("S1.4e", r, tc_e)
def rule_score(D, arm, key):
    return sum(D[("P-1", a_, arm)][key] == "confermata" for a_ in ARCHS) + sum(D[("P-2", a_, arm)][key] == "confermata" for a_ in ARCHS[:4])
def rule_N(D, key, arms, cond_c=False):
    sc = {a: rule_score(D, a, key) for a in arms}; p3c = {a: sum(D[("P-3", a_, a)][key] == "confermata" for a_ in ["A1", "A2", "A3", "A4", "A6"]) for a in arms}
    p3r = {a: sum(D[("P-3", a_, a)][key] == "smentita" for a_ in ["A1", "A2", "A3", "A4", "A6"]) for a in arms}
    ok = []
    for a in arms:
        if a == "K2": continue
        ca = sc[a] >= sc["K2"] + 2
        cb = not any(D[("P-1", a_, "K2")][key] == "confermata" and D[("P-1", a_, a)][key] == "smentita" for a_ in ARCHS)
        cc = (p3r[a] <= 2) if cond_c else True
        if ca and cb and cc: ok.append(a)
    if not ok: w = "K2"
    else:
        ok.sort(key=lambda a: (-sc[a], -p3c[a])); w = ok[0]
        if len(ok) > 1 and sc[ok[0]] == sc[ok[1]] and p3c[ok[0]] == p3c[ok[1]]: w = "K2"
    return w, sc
for key in ():
    pass
wo, so = rule_N(E, "verdict_recorded", ["K2", "ER", "ET"]); wn, sn = rule_N(E, "verdict_new", ["K2", "ER", "ET"]); ws, ss = rule_N(E, "verdict_new_stress", ["K2", "ER", "ET"])
rn = rd(RES + "/S1.4e/S1.4e_ruleN.csv")
add(step="S1.4e", prediction_id="Regola N", archetype="A1-A6", quantity="arm proposed for presets (confirmed P-1 x6 + P-2 x4; >= K2+2; no K2-confirmed P-1 tissue refuted)",
    verdict_recorded="%s (scores %s)" % (wo, so), verdict_new="%s (scores %s)" % (wn, sn), verdict_new_stress="%s (scores %s)" % (ws, ss), transition="winner %s->%s" % (wo, wn),
    est_unit="seed", margin_unit="window", qtype="nuc", kappa_margin=KW, reproduces=str(wo == rn.winner.iloc[0] and so == dict(zip(rn.arm, rn.score_P1_P2))),
    touches_conclusion="preset nucleus shape = disc K2 (DEC-026, DEC-027)", notes="recomputed from per-tissue corrected P-1/P-2 verdicts (window kappa on the manual references); stress winner %s" % ws)
cpe = rd(RES + "/S1.4e/S1.4e_cp.csv")
cpE = {}
for _, r in cpe.iterrows():
    A = r["archetype"]
    if str(r["id"]).startswith("CP-1"):
        bz = "sotto zero" if bool(r["below_zero"]) else "non sotto zero"
        j = judge(r["est"], r["lo2"], r["hi2"], -INF, 0.0, "win_se", "fixed", "nuc", A)
        f = lambda v: "sotto zero" if v == "confermata" else "non sotto zero"
        rep = (verdict(r["lo2"], r["hi2"], -INF, 0.0) == "confermata") == bool(r["below_zero"])
        add(step="S1.4e", prediction_id="%s %s" % (r["id"], r["contrast"]), archetype=A, quantity="paired contact difference on Luca's judged pairs <6 um (bench A)", est=r["est"], lo=r["lo2"], hi=r["hi2"], a=-INF, b=0.0,
            est_unit="window", margin_unit="fixed", qtype="nuc", verdict_recorded=bz, reproduces=str(rep), lo_new=j["lo_new"], hi_new=j["hi_new"], kappa_est=j["kappa_est"],
            verdict_new=f(j["verdict_new"]), verdict_new_stress=f(j["verdict_new_stress"]), touches_conclusion="CP-1 passed -> orientation, not seeding, removes contacts (DEC-026 wording)",
            notes="estimate = mean of 4 manual windows +- 2 SE between windows (one-stage) -> window kappa_se %.4f x 2/2 (stress x1.343); margin 'below 0' fixed" % KWS)
        cpE[(r["id"], A, r["contrast"])] = ROWS[-1]
    else:
        cpE[("CP-2", A, r["contrast"])] = pred_row("S1.4e", pd.Series({"id": "P-2", "archetype": A, "arm": r["contrast"], "est": r["est"], "se": r["se"], "lo": r["lo"], "hi": r["hi"], "verdict": r["verdict"], "source": "R-1"}), lambda p, a: "CP-2 -> rule N judged on P-1 and P-2 as written", extra_id=" (CP-2, noisy centres)")
        ROWS[-1]["prediction_id"] = "CP-2 %s (P-2 with centre noise)" % r["contrast"]
for lab in ("CP-1", "CP-1 con le non giudicabili"):
    for key, nm in (("verdict_recorded", "o"), ("verdict_new", "n"), ("verdict_new_stress", "s")):
        pass
    def cnt(key):
        et = sum(cpE[(lab, a_, "ET-K2")][key] == "sotto zero" for a_ in ARCHS[:4]); er = sum(cpE[(lab, a_, "ER-K2")][key] != "sotto zero" for a_ in ARCHS[:4])
        return ("superata" if (et >= 3 and er >= 3) else ("ET part met, ER part not" if et >= 3 else "fallita")), et, er
    o, n, s_ = cnt("verdict_recorded"), cnt("verdict_new"), cnt("verdict_new_stress")
    add(step="S1.4e", prediction_id="%s (esito)" % lab, archetype="A1-A4", quantity="ET-K2 < 0 beyond 2 SE in >=3/4 and ER-K2 not negative beyond 2 SE in >=3/4", verdict_recorded=o[0], verdict_new=n[0], verdict_new_stress=s_[0],
        est_unit="window", margin_unit="fixed", qtype="nuc", kappa_est=KW, reproduces=str(o[0] == ("superata" if lab == "CP-1" else "ET part met, ER part not")), touches_conclusion="CP-1 passed (main reading) -> DEC-026",
        notes="counts ET<0 / ER not<0: old %d/%d, new %d/%d, stress %d/%d" % (o[1], o[2], n[1], n[2], s_[1], s_[2]))
def cp2cnt(key): return all(cpE[("CP-2", a_, "K2")][key] == "smentita" for a_ in ("A1", "A4"))
add(step="S1.4e", prediction_id="CP-2 (esito)", archetype="A1,A4", quantity="K2 stays smentita in intestine and lymph node with centre noise", verdict_recorded="atteso: %s" % cp2cnt("verdict_recorded"),
    verdict_new="atteso: %s" % cp2cnt("verdict_new"), verdict_new_stress="atteso: %s" % cp2cnt("verdict_new_stress"), est_unit="seed", margin_unit="window_two_stage", qtype="nuc", kappa_margin=KW, reproduces="NA",
    touches_conclusion="if K2 re-entered in all 4 tissues, rule N would be judged on P-1 only (not triggered)", notes="recorded text: intestine indeterminata, lymph node smentita -> expectation literally not met in intestine; kept as recorded rows")
for g_ in ("G-1", "G-2", "G-3", "G-4", "G-5"):
    add(step="S1.4e", prediction_id=g_, archetype="A1-A6", quantity="gate (see step file)", verdict_recorded="PASS", est_unit="fixed", margin_unit="fixed", qtype="nuc", notes="gate, from step file table 'Cancelli'")

# --- S1.4f
def tc_f(pid, arm):
    if pid in ("P-1", "P-2"): return "rule N2 -> preset disc K2 kept; remaining contact excess declared limit (DEC-028, DEC-029)"
    if pid in ("P-3", "P-5"): return "EP law = compression V2 (rule O, DEC-028), package option, not preset; BL-175"
    return "none (weak test by construction)"
vf = rd(RES + "/S1.4f/S1.4f_verdicts.csv"); Fd = {}
for _, r in vf.iterrows():
    row = pred_row("S1.4f", r, tc_f); Fd[(r["id"], r["archetype"], r["arm"]) if pd.isna(r["cls"]) else (r["id"], r["archetype"], r["arm"], r["cls"])] = row
wo, so = rule_N(Fd, "verdict_recorded", ["K2", "ET", "EP"], cond_c=True); wn, sn = rule_N(Fd, "verdict_new", ["K2", "ET", "EP"], cond_c=True); ws, ss = rule_N(Fd, "verdict_new_stress", ["K2", "ET", "EP"], cond_c=True)
rn2 = rd(RES + "/S1.4f/S1.4f_ruleN2.csv")
add(step="S1.4f", prediction_id="Regola N2", archetype="A1-A6", quantity="arm proposed for presets (N + P-3 refuted in <=2/5)", verdict_recorded="%s (scores %s)" % (wo, so), verdict_new="%s (scores %s)" % (wn, sn),
    verdict_new_stress="%s (scores %s)" % (ws, ss), transition="winner %s->%s" % (wo, wn), est_unit="area", margin_unit="window", qtype="nuc", kappa_margin=KW,
    reproduces=str(wo == rn2.winner.iloc[0] and so == dict(zip(rn2.arm, rn2.score_P1_P2))), touches_conclusion="preset disc K2 (DEC-028/DEC-029)",
    notes="recomputed from corrected P-1/P-2 (window) and P-3 (area+range) verdicts; stress winner %s" % ws)
ro = rd(RES + "/S1.4f/S1.4f_ruleO_cells.csv"); rO = rd(RES + "/S1.4f/S1.4f_ruleO.csv"); sc_o = {}; sc_n = {}
for _, r in ro.iterrows():
    A = r["archetype"]; kr = K[A]["rnuc"]; mid = (r["lo"] + r["hi"]) / 2; a2, b2 = mid - kr * (mid - r["lo"]), mid + kr * (r["hi"] - mid)
    ins_o = bool(r["lo"] <= r["est"] <= r["hi"]); ins_n = bool(a2 <= r["est"] <= b2)
    sc_o[r["law"]] = sc_o.get(r["law"], 0) + ins_o; sc_n[r["law"]] = sc_n.get(r["law"], 0) + ins_n
    add(step="S1.4f", prediction_id="Regola O cell %s %s" % (r["law"], r["cls"]), archetype=A, quantity="alignment %s, seeds 1-5, mean of 5 ROIs (point) inside R-4" % r["cls"], est=r["est"], lo=r["est"], hi=r["est"],
        a=r["lo"], b=r["hi"], a_new=a2, b_new=b2, est_unit="other", margin_unit="range", qtype="nuc", kappa_margin=kr, verdict_recorded="inside" if bool(r["inside"]) else "outside",
        verdict_new="inside" if ins_n else "outside", reproduces=str(ins_o == bool(r["inside"])), touches_conclusion="rule O -> EP law (compression V2, DEC-028)",
        notes="rule uses the point estimate only (no SE) -> estimate unit 'other'; R-4 = [min,max] of 5 real ROIs stretched by kappa_range_nuc")
def win(sc):
    ab = rO.set_index("law")["sum_absdev"]; best = max(sc.values()); c = [l for l in sc if sc[l] == best]
    return min(c, key=lambda l: ab[l])
add(step="S1.4f", prediction_id="Regola O", archetype="A1-A4", quantity="law with most P-5 cells inside R-4 (tie: smallest sum |est - real mean|)", verdict_recorded="%s (scores %s)" % (win(sc_o), sc_o),
    verdict_new="%s (scores %s)" % (win(sc_n), sc_n), transition="winner %s->%s" % (win(sc_o), win(sc_n)), est_unit="other", margin_unit="range", qtype="nuc",
    kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS[:4]), reproduces=str(win(sc_o) == rO.winner.iloc[0] and all(sc_o[l] == s for l, s in zip(rO.law, rO.score))),
    touches_conclusion="EP law = compression V2 (DEC-028)", notes="per-tissue kappa_range_nuc on each cell; tie-break sums unchanged (do not depend on R-4)")
cpf = rd(RES + "/S1.4f/S1.4f_cp.csv"); cF = {}
for _, r in cpf.iterrows():
    A = r["archetype"]
    if r["id"] in ("CP-1", "S-1"):
        j = judge(r["est"], r["lo2"], r["hi2"], -INF, 0.0, "win_se", "fixed", "nuc", A); f = lambda v: "sotto zero" if v == "confermata" else "non sotto zero"
        add(step="S1.4f", prediction_id="%s %s" % (r["id"], r["contrast"]), archetype=A, quantity="paired contact difference, bench A%s" % (" (duplicate centres removed)" if r["id"] == "S-1" else ""), est=r["est"], lo=r["lo2"], hi=r["hi2"], a=-INF, b=0.0,
            est_unit="window", margin_unit="fixed", qtype="nuc", verdict_recorded="sotto zero" if bool(r["below_zero"]) else "non sotto zero", lo_new=j["lo_new"], hi_new=j["hi_new"], kappa_est=j["kappa_est"],
            verdict_new=f(j["verdict_new"]), verdict_new_stress=f(j["verdict_new_stress"]), reproduces=str((verdict(r["lo2"], r["hi2"], -INF, 0.0) == "confermata") == bool(r["below_zero"])),
            touches_conclusion="CP-1 of S1.4f passed -> EP effect from nucleus drawing" if r["id"] == "CP-1" else "none (sensitivity)", notes="mean of 4 manual windows +- 2 SE between windows (one-stage) -> window kappa_se %.4f x 2/2 (stress x1.343)" % KWS)
        cF[(r["id"], A, r["contrast"])] = ROWS[-1]
    else:
        add(step="S1.4f", prediction_id="cuore %s (duplicate of P-1 A6)" % r["contrast"], archetype=A, quantity="contacts in heart (negative control)", verdict_recorded=r["verdict"],
            verdict_new=Fd[("P-1", "A6", r["contrast"])]["verdict_new"], verdict_new_stress=Fd[("P-1", "A6", r["contrast"])]["verdict_new_stress"], est_unit="seed", margin_unit="window", qtype="nuc", kappa_margin=KW,
            reproduces=str(r["verdict"] == Fd[("P-1", "A6", r["contrast"])]["verdict_recorded"]), notes="same numbers as the P-1 A6 row")
def cpf_out(key):
    ep = sum(cF[("CP-1", a_, "EP-K2")][key] == "sotto zero" for a_ in ARCHS[:4])
    pt = sum(cF[("CP-1", a_, "EP-K2")]["est"] >= cF[("CP-1", a_, "ET-K2")]["est"] for a_ in ARCHS[:4])
    return ("superata" if ep >= 2 and pt >= 3 else "fallita"), ep, pt
o, n, s_ = cpf_out("verdict_recorded"), cpf_out("verdict_new"), cpf_out("verdict_new_stress")
add(step="S1.4f", prediction_id="CP-1 (esito)", archetype="A1-A4", quantity="EP-K2 < 0 beyond 2 SE in >=2/4 and EP-K2 not more negative than ET-K2 in >=3/4", verdict_recorded=o[0], verdict_new=n[0], verdict_new_stress=s_[0],
    est_unit="window", margin_unit="fixed", qtype="nuc", kappa_est=KW, reproduces=str(o[0] == "superata"), touches_conclusion="CP-1 S1.4f passed", notes="EP<0 count old %d, new %d, stress %d" % (o[1], n[1], s_[1]))
s1t = rd(RES + "/S1.4f/S1.4f_S1.csv")
for _, r in s1t.iterrows():
    A = r["archetype"]; m = 2 * r["se_A"]
    add(step="S1.4f", prediction_id="S-1 change %s" % r["contrast"], archetype=A, quantity="change of CP-1 contrast after removing duplicate centres, within 2 SE", est=r["delta"], lo=r["delta"], hi=r["delta"], a=-m, b=m,
        a_new=-KWS * m, b_new=KWS * m, est_unit="other", margin_unit="window", qtype="nuc", kappa_margin=KW, verdict_recorded="entro 2 SE" if bool(r["within_2se"]) else "oltre",
        verdict_new="entro 2 SE" if abs(r["delta"]) <= KWS * m else "oltre", verdict_new_stress="entro 2 SE" if abs(r["delta"]) <= KWS * STRESS * m else "oltre",
        reproduces=str((abs(r["delta"]) <= m) == bool(r["within_2se"])), notes="tolerance = 2 SE between the 4 manual windows (one-stage) -> window kappa_se x 2/2")
def cp2f(key): return sum(Fd[("P-5", a_, "V0", c_)][key] == "smentita" for a_ in ARCHS[:4] for c_ in ("nn5-6", "nn6-8"))
add(step="S1.4f", prediction_id="CP-2 (esito)", archetype="A1-A4", _sde="simulated", quantity="constant-dispersion law V0 refuted in >=4 of 8 P-5 cells", verdict_recorded="superata (%d/8)" % cp2f("verdict_recorded"),
    verdict_new=("superata" if cp2f("verdict_new") >= 4 else "fallita") + " (%d/8)" % cp2f("verdict_new"), verdict_new_stress=("superata" if cp2f("verdict_new_stress") >= 4 else "fallita") + " (%d/8)" % cp2f("verdict_new_stress"),
    est_unit="area", margin_unit="range", qtype="nuc", kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS[:4]), reproduces=str(cp2f("verdict_recorded") == 6), touches_conclusion="rule O reading (F2) / EP law", notes="aggregate of the V0 P-5 rows")
def p5ep(key):
    c = [Fd[("P-5", a_, "EP", c_)][key] for a_ in ARCHS[:4] for c_ in ("nn5-6", "nn6-8")]; n_ = sum(x in ("confermata", "indeterminata") for x in c)
    return ("confermata" if n_ >= 6 else "smentita") + " (%d/8 conf.+indet.)" % n_
add(step="S1.4f", prediction_id="P-5 EP (esito)", _sde="simulated", archetype="A1-A4", quantity="EP confirmed or indeterminate in >=6 of 8 cells", verdict_recorded=p5ep("verdict_recorded"), verdict_new=p5ep("verdict_new"), verdict_new_stress=p5ep("verdict_new_stress"),
    est_unit="area", margin_unit="range", qtype="nuc", kappa_margin=max(K[a_]["rnuc"] for a_ in ARCHS[:4]), reproduces="True", touches_conclusion="BL-175 (steeper law deferred)", notes="aggregate of the EP P-5 rows")
ph2 = rd(RES + "/S1.4f/S1.4f_PH2_verdicts.csv")
for _, r in ph2.iterrows():
    rr = r.copy(); rr["source"] = "B-127...B-132" if r["id"] == "P-1" else "R-1"; rr["cls"] = np.nan
    pred_row("S1.4f", rr, lambda p, a: "PH-2: rule O tie-break does not change conclusions", extra_id=" (PH-2 post hoc)")
for g_ in ("G-1", "G-2", "G-3", "G-4", "G-5"):
    add(step="S1.4f", prediction_id=g_, archetype="A1-A6", quantity="gate (see step file)", verdict_recorded="PASS", est_unit="fixed", margin_unit="fixed", qtype="nuc", notes="gate, from step file table 'Cancelli'")

# =========================== write ===========================
out = pd.DataFrame(ROWS, columns=COLS)
out.to_csv(T + "/R7_M5_S13_S14.csv", index=False)
aff = out[out.transition != "unaffected"]
chg = aff[aff.verdict_recorded.astype(str) != aff.verdict_new.astype(str)]
chs = aff[aff.verdict_recorded.astype(str) != aff.verdict_new_stress.astype(str)]
summ = {"n_rows": len(out), "n_affected": len(aff), "by_step": out.groupby("step").size().to_dict(), "aff_by_step": aff.groupby("step").size().to_dict(),
        "transitions": aff.transition.value_counts().to_dict(), "changed": chg[["step", "prediction_id", "archetype", "verdict_recorded", "verdict_new", "verdict_new_stress", "touches_conclusion"]].to_dict("records"),
        "changed_stress_only": chs[~chs.index.isin(chg.index)][["step", "prediction_id", "archetype", "verdict_recorded", "verdict_new", "verdict_new_stress"]].to_dict("records"),
        "repro": out.reproduces.value_counts().to_dict(), "log": LOG,
        "a5": out[(out.verdict_geom_A5 != "") & (out.verdict_geom_A5.astype(str) != out.verdict_new.astype(str))][["step", "prediction_id", "archetype", "verdict_recorded", "verdict_new", "verdict_geom_A5"]].to_dict("records"),
        "n_a5rows": int((out.verdict_geom_A5 != "").sum()),
        "trans_by_step": aff.assign(t=np.where(aff.transition == "unchanged", "unchanged", aff.verdict_recorded.astype(str).str[:45] + "->" + aff.verdict_new.astype(str).str[:45])).groupby(["step", "t"]).size().reset_index().values.tolist()}
json.dump(summ, open(T + "/summary.json", "w"), indent=1, default=str)
print("rows", len(out), "affected", len(aff), "changed", len(chg), "stress-changed", len(chs))
