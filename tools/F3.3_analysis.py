"""tools/F3.3_analysis.py - F3.3: statistiche pre-registrate (docs/steps/F3.3.md) dalle tabelle per area.
Ingressi: /mnt/micron/geo_spatialtrans/F3.3/calib/*.json (taratura e simulato di F3.1), F3.3/measures/*.json (aree di verifica e aree di
taratura del tumore del colon con i geni comuni), results/F3.3/F3.3_prereg_numbers.csv (margini h, fissati prima).
Uscite in results/F3.3/: F3.3_area_values.csv, F3.3_bridge_cg.csv, F3.3_cells.csv, F3.3_summary.csv, F3.3_varcomp.csv, F3.3_r4a.csv.
Uso: .venv/bin/python tools/F3.3_analysis.py [sim_cg]"""
import sys, glob, json, itertools, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "tools"); import F3_3_stats as F
M = "/mnt/micron/geo_spatialtrans/F3.3"; OUT = "results/F3.3"; B = 2000; SEED = 20261119
MEAS = ["M1", "M2", "M3", "M4", "M5", "M6"]; MID = {m: f"M-{m[1]}" for m in MEAS + ["M7"]}
SIM_CG = len(sys.argv) > 1 and sys.argv[1] == "sim_cg"
H = pd.read_csv(f"{OUT}/F3.3_prereg_numbers.csv").set_index(["archetype", "id"])
cal = pd.DataFrame([r for f in sorted(glob.glob(f"{M}/calib/*.json")) for r in json.load(open(f))])
hol = pd.DataFrame([r for f in sorted(glob.glob(f"{M}/measures/*.json")) for r in json.load(open(f))])
AR = pd.read_csv(f"{OUT}/F3.3_areas.csv").set_index("area")
# ---------- valori per unita' ----------
rows = []
for _, r in cal[(cal.labels == "tl")].iterrows():
    if r.A == "A4":
        if not (r.part.startswith("q") and "_f" in r.tag): continue
        unit = f"{r.tag}_{r.part}"
    else:
        if r.part != "whole": continue
        unit = r.tag
    if r.side == "real" and r.A in ("A2", "A3"):
        src = "IN_panel"                      # tumore e stroma: il confronto usa i geni comuni (sotto); il pannello intero serve al ponte
    else:
        src = "IN" if r.side == "real" else "SIM"
    rows.append({"unit": unit, "run": r.run, "A": r.A, "src": src, "set": f"{r.A}_IN" if r.side == "real" else f"{r.A}_SIM", "level": "IN" if r.side == "real" else "SIM",
                 **{m: r[m] for m in MEAS + ["M7"]}})
for _, r in hol[hol.labels == "tl"].iterrows():
    a = AR.loc[r.area]; lev = r.level
    if lev == "IN": src, st, lv = "IN", f"{r.A}_IN", "IN"
    else: src, st, lv = "OUT", f"{a.section}_{a.level}_{r.A}", a.level
    rows.append({"unit": r.area, "run": r.run, "A": r.A, "src": src, "set": st, "level": lv, **{m: r[m] for m in MEAS + ["M7"]},
                 **{k: r[k] for k in ("M2r", "M3r", "M3u") if k in r}})
if SIM_CG:
    simcg = pd.DataFrame([r for f in sorted(glob.glob(f"{M}/sim_cg_measures/*.json")) for r in json.load(open(f))])
    rows = [x for x in rows if not (x["src"] == "SIM" and x["A"] in ("A2", "A3"))]
    for _, r in simcg[simcg.labels == "tl"].iterrows():
        rows.append({"unit": r.run.split("_k")[0], "run": r.run, "A": r.A, "src": "SIM", "set": f"{r.A}_SIM", "level": "SIM", **{m: r[m] for m in MEAS + ["M7"]}})
U = pd.DataFrame(rows)
# simulato: media dei semi per unita' (area o quarto)
U["u2"] = np.where(U.src == "SIM", U.run.str.replace(r"_k[01]", "", regex=True) if True else U.unit, U.unit)
U.loc[(U.src == "SIM") & (U.A == "A4"), "u2"] = U.loc[(U.src == "SIM") & (U.A == "A4"), "unit"]
U.to_csv(f"{OUT}/F3.3_area_values.csv", index=False)
SIMU = U[U.src == "SIM"].groupby(["A", "u2"])[MEAS + ["M7"]].mean().reset_index()
# ---------- ponte dei geni comuni (DEC-097 punto 4) ----------
br = []
for A in ("A2", "A3"):
    p = U[(U.A == A) & (U.src == "IN_panel")].set_index("unit")[MEAS]; c = U[(U.A == A) & (U.src == "IN")].set_index("unit")[MEAS]
    for m in MEAS:
        d = float(c[m].mean() - p[m].mean()); h = float(H.loc[(A, MID[m]), "h"])
        br.append({"A": A, "measure": MID[m], "panel": p[m].mean(), "common": c[m].mean(), "delta": d, "h": h, "abs_delta_over_h": abs(d) / h, "exceeds_h10": abs(d) > h / 10})
BR = pd.DataFrame(br); BR.to_csv(f"{OUT}/F3.3_bridge_cg.csv", index=False)
need_sim_cg = bool(BR.exceeds_h10.any())
print("ponte dei geni comuni: rifare il simulato del tumore e dello stroma =", need_sim_cg, "| max |delta|/h =", round(BR.abs_delta_over_h.max(), 3))
# ---------- celle ----------
REAL = U[U.src.isin(["IN", "OUT"])]
cells, crow = [], []
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    R = REAL[REAL.A == A]; IN = R[R.set == f"{A}_IN"]; S = SIMU[SIMU.A == A]
    outs = {s: g for s, g in R[R.src == "OUT"].groupby("set")}
    if not len(outs): continue
    for mi, m in enumerate(MEAS):
        h = float(H.loc[(A, MID[m]), "h"]); rin = IN[m].to_numpy(float); sv = S[m].to_numpy(float)
        sp, df = F.pooled_sd([rin] + [g[m].to_numpy(float) for g in outs.values()]); t = stats.t.ppf(0.975, df)
        Sm = float(sv.mean()); se_s = float(sv.std(ddof=1) / np.sqrt(len(sv)))
        g_in = Sm - rin.mean(); w_in = t * np.sqrt(se_s ** 2 + sp ** 2 / len(rin)); v_in = F.vrule(g_in, w_in, h)
        resid = np.vstack([np.column_stack([gg[x].to_numpy(float) - gg[x].to_numpy(float).mean() for x in MEAS]) for gg in [IN] + list(outs.values())])
        base = {"in": rin, "h": h, "sd": sp, "df": df, "sim": Sm, "se_sim": se_s, "good": v_in != "smentita", "A": A, "m": m, "mi": mi, "resid": resid}
        for lev in ("L1", "L2"):
            o = {s: g[m].to_numpy(float) for s, g in outs.items() if g.level.iloc[0] == lev}
            if o: cells.append({**base, "out": o, "level": lev, "group": f"{A}_{lev}"})
        for s, g in outs.items():
            rv = g[m].to_numpy(float); n = len(rv)
            g_out = Sm - rv.mean(); w_out = t * np.sqrt(se_s ** 2 + sp ** 2 / n)
            D = rv.mean() - rin.mean(); wD = t * sp * np.sqrt(1 / n + 1 / len(rin))
            crow.append({"A": A, "measure": MID[m], "set": s, "level": g.level.iloc[0], "n_out": n, "n_in": len(rin), "h": h, "sd_pooled": sp, "df": df,
                         "real_in": rin.mean(), "real_out": rv.mean(), "sim": Sm, "gap_in": g_in, "w_in": w_in, "verdict_in": v_in,
                         "gap_out": g_out, "w_out": w_out, "verdict_out": F.vrule(g_out, w_out, h), "shift": D, "w_shift": wD, "verdict_shift": F.vrule(D, wD, h),
                         "gap_in_over_h": g_in / h, "gap_out_over_h": g_out / h, "shift_over_h": D / h,
                         "same_sign": bool(np.sign(g_out) == np.sign(g_in))})
C = pd.DataFrame(crow); C.to_csv(f"{OUT}/F3.3_cells.csv", index=False)
summ = []
def band(name, lev, est):
    summ.append({"id": name, "level": lev, **{k: est[k] for k in ("obs", "null_median", "est", "lo", "hi")}, "verdict": F.verdict_band(est["lo"], est["hi"], (-0.5, 0.5))})
def qband(name, lev, q):
    v = "smentita" if q["lo"] > 0.5 else ("smentita per difetto" if q["hi"] < -0.5 else "non rilevato a questa profondita'")
    summ.append({"id": name, "level": lev, "obs": q["est"], "est": q["est"], "lo": q["lo"], "hi": q["hi"], "null_lo": q.get("null_lo"), "null_hi": q.get("null_hi"),
                 "obs_q": q["obs_q"], "verdict": v})
for lev, mode in (("L1", "perm"), ("L2", "param")):
    cl = [c for c in cells if c["level"] == lev]
    qband("P-1" if lev == "L1" else "P-2", lev, F.q_estimate(cl, F.stat_QS, B, SEED, null_mode=mode))          # DEC-098: stimatori principali
    qband("P-3", lev, F.q_estimate(cl, F.stat_QG, B, SEED + 1, null_mode=mode))
    band("P-1/2 come scritta (S*)", lev, F.corrected(cl, F.stat_S, mode, B, SEED))
    band("P-3 come scritta (G*)", lev, F.corrected(cl, F.stat_G, mode, B, SEED + 1))
    # sensibilita' non pre-registrata: residui per area ricampionati insieme (correlazione fra misure)
    band("S* sens. residui congiunti", lev, F.corrected(cl, F.stat_S, "resid", B, SEED + 2))
    band("G* sens. residui congiunti", lev, F.corrected(cl, F.stat_G, "resid", B, SEED + 3))
    if lev == "L2":
        qband("P-3 solo altri pazienti (P1, P5), descrittiva", lev, F.q_estimate([{**c, "out": {k: v for k, v in c["out"].items() if "CRC_P" in k}} for c in cl if any("CRC_P" in k for k in c["out"])], F.stat_QG, B, SEED + 4))
        qband("P-2 solo altri pazienti (P1, P5), descrittiva", lev, F.q_estimate([{**c, "out": {k: v for k, v in c["out"].items() if "CRC_P" in k}} for c in cl if any("CRC_P" in k for k in c["out"])], F.stat_QS, B, SEED + 5))
S_ = pd.DataFrame(summ)
# persistenza (P-4) e condizione 4 di R-3 (P-5)
for lev in ("L1", "L2"):
    c = C[(C.level == lev) & (C.verdict_in == "smentita")]
    k = int(((c.verdict_out == "smentita") & c.same_sign).sum()); n = len(c)
    summ.append({"id": "P-4", "level": lev, "obs": k / n if n else np.nan, "n": n, "k": k, "verdict": "confermata" if n and k / n >= 0.8 else ("smentita" if n else "non valutabile")})
    c = C[(C.level == lev) & (C.verdict_in == "confermata")]
    bad = c[c.verdict_out == "smentita"]; per_m = bad.groupby("measure").set.nunique() if len(bad) else pd.Series(dtype=int)
    summ.append({"id": "P-5", "level": lev, "obs": int(per_m.max()) if len(per_m) else 0, "n": len(c), "k": len(bad),
                 "detail": "; ".join(f"{a}:{s}" for a, s in zip(bad.A + " " + bad.measure, bad.set)), "verdict": "rispettata" if (per_m <= 1).all() else "violata"})
# P-7: condizione (a) di R-4 fuori campione (runs simulati: aree x semi; linfonodo: quarti x semi)
r4 = []
SIMR = U[U.src == "SIM"]
for _, r in REAL[REAL.src == "OUT"].iterrows():
    sv = SIMR[SIMR.A == r.A]
    for m in MEAS:
        lo, hi = np.percentile(sv[m].to_numpy(float), [5, 95]); r4.append({"unit": r.unit, "A": r.A, "set": r.set, "measure": MID[m], "value": r[m], "sim_q05": lo, "sim_q95": hi, "inside": bool(lo <= r[m] <= hi)})
R4 = pd.DataFrame(r4); R4.to_csv(f"{OUT}/F3.3_r4a.csv", index=False)
summ.append({"id": "P-7", "level": "L1+L2", "obs": R4.inside.mean(), "k": int(R4.inside.sum()), "n": len(R4), "verdict": "come atteso" if R4.inside.mean() <= 0.10 else "diversa dall'attesa"})
# P-6: componenti di varianza fra pazienti (P2 = taratura e aree nuove; P1; P5), tumore del colon da 6,5 mm per sensibilita'
vc = []
for A in ("A2", "A3"):
    R = REAL[REAL.A == A]; grp = {"P2": R[R.set.isin([f"{A}_IN", f"A23_L1_{A}"])], "P1": R[R.set == f"CRC_P1_SR300_L2_{A}"], "P5": R[R.set == f"CRC_P5_SR300_L2_{A}"],
                                  "6p5": R[R.set == f"CRC_6p5mm_SR410_L2_{A}"]}
    for m in MEAS:
        for nm, keys in (("P1-P2-P5", ["P2", "P1", "P5"]), ("con 6,5 mm", ["P2", "P1", "P5", "6p5"])):
            gs = [grp[k][m].to_numpy(float) for k in keys if len(grp[k])]
            v = F.var_components(gs, B, SEED + 5); verdict = F.verdict_band(v["lo"], v["hi"], (0.0, 0.5))
            vc.append({"A": A, "measure": MID[m], "groups": nm, "n": ",".join(str(len(g)) for g in gs), **v, "verdict": verdict})
VC = pd.DataFrame(vc); VC.to_csv(f"{OUT}/F3.3_varcomp.csv", index=False)
# CP-1: cuore, M-7 sulle aree nuove (vero - simulato)
r6 = REAL[(REAL.A == "A6") & (REAL.src == "OUT")].M7.to_numpy(float); s6 = SIMU[SIMU.A == "A6"].M7.to_numpy(float)
d7 = r6.mean() - s6.mean(); se7 = np.sqrt(r6.var(ddof=1) / len(r6) + s6.var(ddof=1) / len(s6))
other = []
for A in ("A1", "A2", "A4", "A5"):
    ro = REAL[(REAL.A == A) & (REAL.src == "OUT") & (REAL.level == "L1")].M7.to_numpy(float); so = SIMU[SIMU.A == A].M7.to_numpy(float)
    if len(ro): other.append({"A": A, "diff": ro.mean() - so.mean()})
summ.append({"id": "CP-1", "level": "L1", "obs": d7, "se": se7, "verdict": "passata" if (d7 > 0.032 and d7 - 2 * se7 > 0) else "fallita",
             "detail": "; ".join(f"{o['A']} {o['diff']:.4f}" for o in other)})
# CP-2: falso fuori campione (3 + 2 aree di taratura, 10 modi, media per tessuto e totale)
cp2 = []
for A in ("A1", "A2", "A3", "A5", "A6"):
    IN = REAL[(REAL.A == A) & (REAL.set == f"{A}_IN")].reset_index(drop=True)
    for k, idx in enumerate(itertools.combinations(range(len(IN)), 2)):
        ins = [i for i in range(len(IN)) if i not in idx]; cl = []
        for m in MEAS:
            h = float(H.loc[(A, MID[m]), "h"]); a, b = IN.loc[ins, m].to_numpy(float), IN.loc[list(idx), m].to_numpy(float)
            sp, df = F.pooled_sd([a, b]); cl.append({"in": a, "out": {"p": b}, "h": h, "sd": sp, "df": df, "sim": 0.0, "good": True, "group": "g"})
        e = F.corrected(cl, F.stat_S, "perm", 500, SEED + 10 + k); cp2.append({"A": A, "split": k, **e})
CP2 = pd.DataFrame(cp2); CP2.to_csv(f"{OUT}/F3.3_cp2.csv", index=False)
summ.append({"id": "CP-2", "level": "taratura", "obs": CP2.est.mean(), "detail": "; ".join(f"{A} {v:.3f}" for A, v in CP2.groupby("A").est.mean().items()),
             "verdict": "passata" if abs(CP2.est.mean()) <= 0.25 else "fallita"})
SU = pd.DataFrame(summ); SU["need_sim_cg"] = need_sim_cg; SU["sim_cg_used"] = SIM_CG; SU.to_csv(f"{OUT}/F3.3_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
print(SU.drop(columns=["need_sim_cg"]).round(3).to_string(index=False))
