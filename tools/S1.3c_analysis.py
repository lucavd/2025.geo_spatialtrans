"""tools/S1.3c_analysis.py — aggregazione, regola pre-registrata, tabelle. Uso: .venv/bin/python tools/S1.3c_analysis.py [outdir]
Legge <outdir>/roi/*.json, <outdir>/ctrl/*.json, <outdir>/bins/*_bg.json; scrive results/S1.3c/S1.3c_*.csv."""
import sys, json, glob
from pathlib import Path
import numpy as np, pandas as pd
import importlib.util
sys_argv = sys.argv; sys.argv = sys.argv[:1]
spec = importlib.util.spec_from_file_location("e", Path(__file__).resolve().parent / "S1.3c_ecm.py"); E = importlib.util.module_from_spec(spec); spec.loader.exec_module(E)
OUT = Path(sys_argv[1]) if len(sys_argv) > 1 else Path("/mnt/micron/geo_spatialtrans/S1.3c")
RES = E.ROOT / "results/S1.3c"; RES.mkdir(exist_ok=True)
ARCH = [f"A{i}" for i in range(1, 7)]

def load_roi():
    rows, extra = [], []
    for f in sorted(glob.glob(str(OUT / "roi" / "*.json"))):
        d = json.load(open(f))
        c2 = d.get("by_k_c2shift", {})
        for k, v in d["by_k"].items():
            vs = c2.get(k, v)     # M-02: con C-2 FAIL le grandezze basate su d_nuc usano lo spostamento ottimo
            rows.append({"C2_shift_applied": bool(d.get("C2_shift_applied", False)), "BA_c2": vs["BA"], "r_star_c2": vs["r_star"],
                         "R_profile_c2": vs["R_profile"], "R_null_c2": vs["R_null"], "verdict_roi_c2": vs["verdict_roi"],"archetype": d["archetype"], "roi_id": d["roi_id"], "variant": d["variant"], "geneset": d["geneset"],
                         "territories": d["territories"], "mask": d["mask"], "n_gen": d["n_gen"], **{kk: v[kk] for kk in
                         ("k", "s_um", "thr", "n_tissue_windows", "n_empty", "phi_obs", "phi_null", "phi_null_min", "phi_null_max",
                          "E", "r_star", "BA", "empty_area", "large_area", "n_comp", "large_frac", "R_profile", "R_null", "verdict_roi")}})
        if d["variant"] == "primary": extra.append(d)
    return pd.DataFrame(rows), extra

def ctrl_status():
    out = {}
    for f in sorted(glob.glob(str(OUT / "ctrl" / "*_C5.json"))):
        d = json.load(open(f)); A = d["archetype"]
        out[A] = {"roi": d["roi_id"], "C5a_exact": d["C5a"]["exact"], "C5a_delta": d["C5a"]["delta_empty"], "C5a_full": d["C5a"]["full_geom"],
                  "C5a_partial": d["C5a"]["partial_fallen"], "C5a_geom_frac": d["C5a"]["geom_frac"],
                  "base_E": d["base"]["E"], "A_E": d["typeA"]["E"], "A_BA": d["typeA"]["BA"], "A_large": d["typeA"]["large_frac"],
                  "A_verdict": d["typeA"]["verdict_roi"], "A_frac_zeroed": d["typeA"]["frac_zeroed"],
                  "B_E": d["typeB"]["E"], "B_BA": d["typeB"]["BA"], "B_large": d["typeB"]["large_frac"], "B_verdict": d["typeB"]["verdict_roi"],
                  "B_frac_zeroed": d["typeB"]["frac_zeroed"], "B_disc_radius_um": d["typeB"]["disc_radius_um"], "B_n_discs": d["typeB"]["n_discs"]}
        out[A]["C5b_pass"] = out[A]["A_verdict"] == "A"; out[A]["C5c_pass"] = out[A]["B_verdict"] == "B"
    return out

def archetype_table(R, c5):
    rows = []
    for (v, k, A), g in R.groupby(["variant", "k", "archetype"]):
        bac = g.BA_c2 if v == "primary" else g.BA          # pre-registrato: spostamento ottimo dove C-2 FAIL (solo primaria)
        Em, BAm = float(g.E.median()), float(bac.median())   # mediana sui ROI con BA definita (M-14: dichiarato)
        lf = float(g.large_area.sum() / g.empty_area.sum()) if g.empty_area.sum() > 0 else float("nan")
        raw = E.verdict(Em, BAm, lf)
        c = c5.get(A, {}); shape_ok = bool(c.get("C5b_pass", False) and c.get("C5c_pass", False))
        final = raw if (raw == "nessuno" or shape_ok) else "non determinabile (solo E)"
        rows.append({"variant": v, "k": k, "s_um": 2 * k, "archetype": A, "n_roi": len(g), "E_median": Em, "E_min": g.E.min(), "E_max": g.E.max(),
                     "phi_obs_median": g.phi_obs.median(), "phi_null_median": g.phi_null.median(), "BA_median": BAm, "n_BA_defined": int(bac.notna().sum()), "BA_median_noshift": float(g.BA.median()),
                     "large_frac_pooled": lf, "r_star_median": g.r_star.median(), "R_profile_median": float((g.R_profile_c2 if v == "primary" else g.R_profile).median()), "R_profile_median_noshift": float(g.R_profile.median()), "R_null_median": float((g.R_null_c2 if v == "primary" else g.R_null).median()),
                     "verdict_rule": raw, "C5_shape_valid": shape_ok, "verdict": final})
    return pd.DataFrame(rows)

def check_c(extra, c5):
    rows = []
    for d in extra:
        tag = f"{d['archetype']}_{d['roi_id']}"
        dev = d["n_bins_roi"] / d["n_bins_expected_C1"] - 1
        rows.append({"check": "C-1", "unit": tag, "value": f"{d['n_bins_roi']} / {d['n_bins_expected_C1']:.0f} ({100*dev:+.2f} %)", "status": "PASS" if abs(dev) <= 0.01 else "FAIL"})
        sh = d["C2_max_abs_shift_bins"]
        rows.append({"check": "C-2", "unit": tag, "value": f"({d['C2_best_dx_um']}, {d['C2_best_dy_um']}) um, r {d['C2_best_corr']:.3f} (a 0: {d['C2_corr_at_0']:.3f})",
                     "status": "PASS" if sh <= 1 else ("WARN" if sh == 2 else "FAIL")})
        rows.append({"check": "C-4", "unit": tag, "value": f"conservati {d['C4_sums_preserved']}, mutante viola {d['C4_mutant_violates']}",
                     "status": "PASS" if d["C4_sums_preserved"] and d["C4_mutant_violates"] else "FAIL"})
    for f in sorted(glob.glob(str(OUT / "bins" / "*_bg.json"))):
        b = json.load(open(f)); ds = Path(f).name.replace("_bg.json", "")
        rows.append({"check": "C-3 (pre-registrato, tutti i geni)", "unit": ds, "value": f"bin diversi {b['C3']['prereg_all_genes']['n_diff']} / {b['C3']['n_in_tissue']}",
                     "status": "PASS" if b["C3"]["prereg_all_genes"]["identical"] else "FAIL"})
        rows.append({"check": "C-3 (geni comuni, primaria)", "unit": ds, "value": f"bin diversi {b['C3']['common_genes']['n_diff']}; UMI dei geni esclusi {100*b['C3']['frac_umi_excluded_genes_in_tissue']:.2f} %",
                     "status": "PASS" if b["C3"]["common_genes"]["identical"] else "FAIL"})
        for gs in ("com", "all"):
            for k, v in b["bg"][gs].items():
                rows.append({"check": f"sfondo ({gs}, k={k})", "unit": ds, "value": f"soglia {v['thr']:.0f} UMI, {v['n_windows']} finestre, media {v['mean']:.2f}",
                             "status": "PASS" if v["n_windows"] >= 1000 else "WARN"})
    for A, c in c5.items():
        rows.append({"check": "C-5a", "unit": f"{A}_{c['roi']}", "value": f"delta {c['C5a_delta']} = {c['C5a_full']} + {c['C5a_partial']}; geometrica {c['C5a_geom_frac']:.3f}", "status": "PASS" if c["C5a_exact"] else "FAIL"})
        rows.append({"check": "C-5b (tipo A)", "unit": f"{A}_{c['roi']}", "value": f"E {c['A_E']:.3f}, BA {c['A_BA']:.3f}, grandi {c['A_large']:.2f} -> {c['A_verdict']}", "status": "PASS" if c["C5b_pass"] else "FAIL"})
        rows.append({"check": "C-5c (tipo B)", "unit": f"{A}_{c['roi']}", "value": f"E {c['B_E']:.3f}, BA {c['B_BA']:.3f}, grandi {c['B_large']:.2f} -> {c['B_verdict']}", "status": "PASS" if c["C5c_pass"] else "FAIL"})
    return pd.DataFrame(rows)

def predictions(AT):
    p = AT[(AT.variant == "primary") & (AT.k == 3)].set_index("archetype")
    e = p.E_median; v = p.verdict; vr = p.verdict_rule
    tie = bool((e.drop("A4") <= e["A4"]).any())             # M-16: «la piu' bassa» con parita' = non valutabile
    rows = [
        ("P1", "A4 E piu' bassa del pannello e E(A4) < 0.05", f"E(A4) {e['A4']:.4f}; min altri {e.drop('A4').min():.4f}" + (" (parita')" if tie else ""),
         None if tie else bool((e["A4"] < e.drop("A4").min()) and e["A4"] < 0.05)),
        ("P2", "E(A3) >= 0.05 con verdetto B", f"E(A3) {e['A3']:.4f}; verdetto {v['A3']} (regola: {vr['A3']})", (e["A3"] >= 0.05) and v["A3"] == "B"),
        ("P3", "E(A1) <= E(A2) < E(A3)", f"{e['A1']:.4f} / {e['A2']:.4f} / {e['A3']:.4f}", (e["A1"] <= e["A2"]) and (e["A2"] < e["A3"])),
        ("P4", "A5: A piu' probabile di B", f"verdetto {v['A5']} (regola: {vr['A5']})", None if v["A5"] not in ("A", "B") else v["A5"] == "A"),
        ("P5", "E(A6) < 0.05, nessun B", f"E(A6) {e['A6']:.4f}; verdetto {v['A6']}", (e["A6"] < 0.05) and v["A6"] != "B"),
    ]
    return pd.DataFrame([{"id": i, "previsione": t, "osservato": o, "esito": ("non valutabile" if ok is None else ("PASS" if ok else "FAIL"))} for i, t, o, ok in rows])

def controprove(extra, R, AT):
    rows = []
    P = R[R.variant == "primary"].set_index(["archetype", "roi_id", "k"])
    for d in extra:
        A, roi = d["archetype"], d["roi_id"]
        key = "CP1_c2shift" if d.get("C2_shift_applied", False) else "CP1"     # M-02
        for k in ("3", "5"):
            ba_r = P.loc[(A, roi, int(k)), "BA_c2"]; ba_p = d[key][k]["BA"]; Rp = d[f"{key}_R_profile"]
            ok = cp1_pass(Rp, ba_p, ba_r)                                      # M-01: booleano Python
            rows.append({"cp": "CP-1", "archetype": A, "roi_id": roi, "k": int(k), "c2_shift": key != "CP1",
                         "R_real": P.loc[(A, roi, int(k)), "R_profile_c2"], "R_placebo": Rp,
                         "BA_real": ba_r, "BA_placebo": ba_p, "shift_um": d[f"{key}_shift_um"], "frac_changed": d[f"{key}_frac_changed"],
                         "R_placebo_in_band": bool(0.9 <= Rp <= 1.1), "BA_defined": bool(np.isfinite(ba_r) and np.isfinite(ba_p)), "pass": ok})
            rows.append({"cp": "CP-2", "archetype": A, "roi_id": roi, "k": int(k), "n_windows": d["CP2"][k]["n_windows"], "n_empty": d["CP2"][k]["n_empty"]})
    cp = pd.DataFrame(rows)
    summ = []
    for A in ARCH:
        for k in (3, 5):
            c1 = cp[(cp.cp == "CP-1") & (cp.archetype == A) & (cp.k == k)]
            c2 = cp[(cp.cp == "CP-2") & (cp.archetype == A) & (cp.k == k)]
            nw, ne = int(c2.n_windows.sum()), int(c2.n_empty.sum())
            npass = int(c1["pass"].astype(bool).sum()); nband = int(c1["R_placebo_in_band"].astype(bool).sum()); nba = int(c1["BA_defined"].astype(bool).sum())
            summ.append({"cp": "CP-1", "archetype": A, "k": k, "value": f"{npass}/{len(c1)} ROI (R placebo in banda {nband}/{len(c1)}; BA definita {nba}/{len(c1)}); R placebo {c1.R_placebo.median():.3f} (reale {c1.R_real.median():.3f}); BA placebo {c1.BA_placebo.median():.3f} (reale {c1.BA_real.median():.3f})",
                         "status": f"PASS {npass}/{len(c1)}" if npass == len(c1) else f"FAIL ({npass}/{len(c1)} ROI)"})
            summ.append({"cp": "CP-2", "archetype": A, "k": k, "value": f"{ne}/{nw} finestre vuote ({(ne/nw if nw else float('nan')):.3f})",
                         "status": "non valutabile" if nw < 20 else ("PASS" if ne / nw >= 0.90 else "FAIL")})
        a6 = AT[(AT.variant == "primary") & (AT.archetype == A)].set_index("k")
        same = ((a6.loc[3, "E_median"] >= 0.05) == (a6.loc[5, "E_median"] >= 0.05)) and (a6.loc[3, "verdict"] == a6.loc[5, "verdict"])
        summ.append({"cp": "CP-3", "archetype": A, "k": None, "value": f"s6: E {a6.loc[3,'E_median']:.4f} {a6.loc[3,'verdict']}; s10: E {a6.loc[5,'E_median']:.4f} {a6.loc[5,'verdict']}",
                     "status": "PASS" if same else "FAIL"})
    return cp, pd.DataFrame(summ)

def cp1_pass(Rp, ba_p, ba_r):
    """CP-1 per ROI (prereg. §8): R placebo in [0.9, 1.1] E BA placebo < BA reale; BA non definita -> False."""
    return bool((0.9 <= Rp <= 1.1) and np.isfinite(ba_p) and np.isfinite(ba_r) and (ba_p < ba_r))

def he_table(extra):
    return pd.DataFrame([{"archetype": d["archetype"], "roi_id": d["roi_id"], **{f"frac_d>{k}um": v for k, v in d["HE_frac_beyond"].items()}} for d in extra])

if __name__ == "__main__" and len(sys_argv) > 1 and sys_argv[1] == "--test":
    NP = NF = 0
    def ck(n, c):
        global NP, NF; print(("PASS " if c else "FAIL ") + n); NP += bool(c); NF += (not c)
    ck("A1 cp1_pass tipo bool", type(cp1_pass(np.float64(1.0), np.float64(0.5), np.float64(0.6))) is bool)
    s = pd.Series([cp1_pass(1.0, 0.5, 0.6), cp1_pass(1.0, 0.5, 0.6), cp1_pass(0.5, 0.5, 0.6)], dtype=object)
    ck("A2 conteggio 2/3 su colonna object (M-01)", int(s.astype(bool).sum()) == 2)
    ck("A3 BA NaN -> False", cp1_pass(1.0, float("nan"), 0.6) is False)
    ck("A4 R fuori banda -> False", cp1_pass(0.8, 0.5, 0.6) is False)
    T = pd.DataFrame({"variant": ["primary"] * 5, "k": [3] * 5, "archetype": ["A1"] * 5, "E": [0.06, 0.07, 0.08, 0.06, 0.07],
                      "BA": [0.5] * 5, "R_null": [1.0] * 5, "R_null_c2": [1.0] * 5, "BA_c2": [0.8, 0.8, 0.8, np.nan, 0.8], "phi_obs": [0.1] * 5, "phi_null": [0.03] * 5,
                      "large_area": [1.0] * 5, "empty_area": [2.0] * 5, "r_star": [3.0] * 5, "R_profile": [0.8] * 5, "R_profile_c2": [0.7] * 5})
    a = archetype_table(T, {"A1": {"C5b_pass": True, "C5c_pass": True}})
    ck("A5 regola usa BA con spostamento (A)", a.verdict.iloc[0] == "A"); ck("A6 n BA definite 4", a.n_BA_defined.iloc[0] == 4)
    a2 = archetype_table(T, {"A1": {"C5b_pass": False, "C5c_pass": True}})
    ck("A7 C-5 fallito -> non determinabile", a2.verdict.iloc[0].startswith("non determinabile"))
    T2 = T.assign(E=0.01); ck("A8 E < 0.05 -> nessuno", archetype_table(T2, {}).verdict.iloc[0] == "nessuno")
    print(f"\n{NP} PASS, {NF} FAIL"); sys.exit(1 if NF else 0)

if __name__ == "__main__":
    R, extra = load_roi(); c5 = ctrl_status()
    AT = archetype_table(R, c5); CC = check_c(extra, c5); PR = predictions(AT); CPr, CPs = controprove(extra, R, AT); HE = he_table(extra)
    R.to_csv(RES / "S1.3c_roi.csv", index=False); AT.to_csv(RES / "S1.3c_archetype.csv", index=False); CC.to_csv(RES / "S1.3c_checkC.csv", index=False)
    PR.to_csv(RES / "S1.3c_predictions.csv", index=False); CPr.to_csv(RES / "S1.3c_cp_roi.csv", index=False); CPs.to_csv(RES / "S1.3c_controprove.csv", index=False)
    HE.to_csv(RES / "S1.3c_he.csv", index=False); pd.DataFrame(c5).T.to_csv(RES / "S1.3c_C5.csv")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(AT[AT.variant == "primary"][["k", "archetype", "E_median", "E_min", "E_max", "phi_obs_median", "phi_null_median", "BA_median", "large_frac_pooled", "R_profile_median", "verdict_rule", "C5_shape_valid", "verdict"]].round(4).to_string(index=False))
    print(PR.to_string(index=False)); print(CPs.to_string(index=False)); print(CC.groupby(["check", "status"]).size().to_string())
