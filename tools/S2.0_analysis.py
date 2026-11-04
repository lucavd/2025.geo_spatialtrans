"""tools/S2.0_analysis.py — S2.0: tabelle, esiti delle previsioni, regole di decisione R-A, R-B, R-C e controprova CP-1.
Pre-registrazione: docs/steps/S2.0.md (afef06b, a8690b3, f2a23d7, AM-5 in 86410af). Seme 20261103, B = 2000.
Uso: .venv/bin/python tools/S2.0_analysis.py <outdir> <stage>     stage: edges | inner | final
  edges: tabelle delle finestre e parametri mediani per sezione -> results/S2.0/S2.0_edges_*.csv, <outdir>/params_rules.json
  inner: controllo interno (P-3) e regola R-A -> results/S2.0/S2.0_inner_*.csv, S2.0_RA.json
  final: strisce (P-1, P-2, R-B), dispersione (P-4, R-C), lumi (CP-1) -> results/S2.0/S2.0_predictions.csv e tabelle"""
import sys, os, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, "tools")
import S2_0_common as S
OUT, STAGE = sys.argv[1], sys.argv[2]
RES = os.environ.get("S20_RES", "results/S2.0"); os.makedirs(RES, exist_ok=True)   # G-4: la seconda esecuzione scrive altrove
R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
SEC = ["A1", "A23", "A4", "A5", "A6"]
SEC_NAME = {"A1": "intestino", "A23": "tumore del colon (e stroma)", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
VOTE = ["A1", "A23", "A4", "A5"]                # il cuore e' descrittivo (DEC-037)
PRM = {"C0": ["sigma", "w", "ell", "dr", "dc", "a", "b"], "C1": ["q", "p", "sigma_m", "w", "ell", "dr", "dc", "a", "b"],
       "C2": ["sigma_e", "D", "w", "ell", "dr", "dc", "a", "b"], "C3": ["sigma_g", "w", "ell", "dr", "dc", "a", "b"]}

def boot_med(v, seed=S.SEED):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) == 0: return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed); reps = np.median(v[rng.integers(0, len(v), (S.B, len(v)))], axis=1)
    return float(np.median(v)), float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5))

def boot_mean(v, seed=S.SEED):
    v = np.asarray(v, float)
    rng = np.random.default_rng(seed); reps = np.mean(v[rng.integers(0, len(v), (S.B, len(v)))], axis=1)
    return float(np.mean(v)), float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5))

def verdict(lo, hi, a, b):
    if not (np.isfinite(lo) and np.isfinite(hi)): return "non calcolabile"
    return S.C.verdict(lo, hi, a, b)

def stage_edges():
    rows = []
    for f in sorted(glob.glob(f"{OUT}/edges/*.json")):
        d = json.load(open(f))
        if "error" in d:
            rows.append({"key": d["key"], "tile": d["tile"], "error": True}); continue
        r = {"key": d["key"], "tile": d["tile"], "n_bins": d["n_bins"], "error": False,
             "nll_C0_r5": d["C0"]["nll"], "nll_C0_recomputed": d["C0"]["nll_recomputed"]}
        for rl in ("C0", "C1", "C2", "C3"):
            for k in PRM[rl] + ["nll", "aic"]:
                if k in d[rl]: r[f"{rl}_{k}"] = d[rl][k]
        for rl in ("C1", "C2", "C3"):
            r[f"dAIC_{rl}"] = r[f"{rl}_aic"] - r["C0_aic"]
        rows.append(r)
    T = pd.DataFrame(rows); T.to_csv(f"{RES}/S2.0_edges_tiles.csv", index=False)
    ok = T[~T.error]
    print("tiles", len(T), "errors", int(T.error.sum()), "C0 recomputed max |diff|", float((ok.nll_C0_r5 - ok.nll_C0_recomputed).abs().max()))
    secrows, params = [], {}
    for k in SEC:
        t = ok[ok.key == k]; params[k] = {}
        for rl in ("C0", "C1", "C2", "C3"):
            params[k][rl] = {}
            for p in PRM[rl]:
                m, lo, hi = boot_med(t[f"{rl}_{p}"])
                params[k][rl][p] = m
                secrows.append({"key": k, "rule": rl, "param": p, "median": m, "lo": lo, "hi": hi, "n_tiles": len(t)})
        for rl in ("C1", "C2", "C3"):
            m, lo, hi = boot_med(t[f"dAIC_{rl}"])
            secrows.append({"key": k, "rule": rl, "param": "dAIC_vs_C0_per_window", "median": m, "lo": lo, "hi": hi, "n_tiles": len(t)})
            secrows.append({"key": k, "rule": rl, "param": "sum_AIC", "median": float(t[f"{rl}_aic"].sum()), "lo": np.nan, "hi": np.nan, "n_tiles": len(t)})
        secrows.append({"key": k, "rule": "C0", "param": "sum_AIC", "median": float(t["C0_aic"].sum()), "lo": np.nan, "hi": np.nan, "n_tiles": len(t)})
    # cuore: mediana delle altre quattro sezioni (DEC-037); i suoi valori restano descrittivi in params_rules_own
    own = json.loads(json.dumps(params["A6"]))
    for rl in params["A6"]:
        for p in params["A6"][rl]:
            params["A6"][rl][p] = float(np.median([params[k][rl][p] for k in VOTE]))
    # C0 per i controlli: i parametri di R5 (params.json di run3), per riprodurre R5
    p5 = json.load(open(f"{R5RUN}/params.json"))
    for k in SEC:
        params[k]["C0_R5"] = p5[k]
    json.dump({"params": params, "A6_own": own, "note": "A6 = mediana delle mediane di A1, A23, A4, A5 (DEC-037)"}, open(f"{OUT}/params_rules.json", "w"), indent=1)
    flat = {k: {rl: params[k][rl] for rl in ("C0", "C1", "C2", "C3")} for k in SEC}
    for k in SEC: flat[k]["C0"] = {**p5[k]}
    json.dump(flat, open(f"{OUT}/params_flat.json", "w"), indent=1)
    pd.DataFrame(secrows).to_csv(f"{RES}/S2.0_edges_sections.csv", index=False)

def stage_inner():
    rows = [json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/inner/*.json"))]
    T = pd.DataFrame(rows).drop(columns=["sec"], errors="ignore"); T.to_csv(f"{RES}/S2.0_inner_roi.csv", index=False)
    out = []
    for A in ("A1", "A2", "A3", "A4", "A5", "A6"):
        for rl in ("C0", "C0_r5u", "C1", "C2", "C3"):
            t = T[(T.archetype == A) & (T.rule == rl)] if "archetype" in T else T.iloc[0:0]
            v = t.ratio_far.to_numpy(float) if len(t) else np.array([])
            n_nan = int((~np.isfinite(v)).sum()) if len(v) else 0
            if len(v) and n_nan == 0:
                m, lo, hi = boot_mean(v); vd = verdict(lo, hi, 0.909, 1.1)
            else:
                m = lo = hi = np.nan; vd = "non calcolabile come scritto" if len(v) else "mancante"
            vv = v[np.isfinite(v)] if len(v) else v
            ms, los, his = boot_mean(vv) if len(vv) else (np.nan,) * 3
            out.append({"archetype": A, "rule": rl, "n_roi": len(v), "n_not_invertible": n_nan, "ratio_mean": m, "lo": lo, "hi": hi, "verdict": vd,
                        "ratio_mean_ok_only": ms, "lo_ok_only": los, "hi_ok_only": his,
                        "r_median": float(np.nanmedian(t.r)) if len(t) else np.nan})
    I = pd.DataFrame(out); I.to_csv(f"{RES}/S2.0_inner_summary.csv", index=False)
    # R-A
    E = pd.read_csv(f"{RES}/S2.0_edges_tiles.csv"); E = E[~E.error]
    elig = {}
    for rl in ("C1", "C2", "C3"):
        vds = {A: I[(I.archetype == A) & (I.rule == rl)].verdict.iloc[0] for A in ("A2", "A5")}
        elig[rl] = {"eligible": all(v != "smentita" for v in vds.values()), "verdicts": vds}
    cand = [rl for rl in ("C1", "C2", "C3") if elig[rl]["eligible"]]
    ra = {"eligibility": elig, "candidates": cand}
    if cand:
        wins = {rl: 0 for rl in cand}; per = {}
        for k in VOTE:
            t = E[E.key == k]; s = {rl: float(t[f"{rl}_aic"].sum()) for rl in cand}
            w = min(s, key=s.get); wins[w] += 1; per[k] = {"sum_aic": s, "winner": w}
        top = max(wins.values()); lead = [rl for rl in cand if wins[rl] == top]
        tot = {rl: float(E[E.key.isin(VOTE)][f"{rl}_aic"].sum()) for rl in cand}
        win = lead[0] if len(lead) == 1 else min(lead, key=lambda r: tot[r])
        ra.update({"per_section": per, "wins": wins, "total_aic": tot, "winner_before_tie_rule": win})
        order = sorted(cand, key=lambda r: tot[r])
        if len(order) > 1:
            a, b = order[0], order[1]
            near = {k: float(np.median(E[E.key == k][f"{a}_aic"] - E[E.key == k][f"{b}_aic"])) for k in VOTE}
            ra["median_dAIC_first_two"] = {"pair": [a, b], "by_section": near}
            if all(abs(x) < 2 for x in near.values()):
                rank = {"C1": (S.K_AIC["C1"], 0), "C3": (S.K_AIC["C3"], 1), "C2": (S.K_AIC["C2"], 2)}
                win = min([a, b], key=lambda r: rank[r]); ra["tie_rule_applied"] = True
        ra["winner"] = win
    else:
        ra["winner"] = None
    json.dump(ra, open(f"{RES}/S2.0_RA.json", "w"), indent=1)
    print(json.dumps(ra, indent=1)[:3000])


def stage_final():
    P = []
    # M1: strisce (P-1, P-2, R-B)
    rb = {"nonperiodic": {}, "parity_spatial": {}}
    for k in SEC:
        js = json.load(open(f"{OUT}/stripes/{k}.json"))
        for ax in ("rows", "cols"):
            r = js[ax]; vd = verdict(r["sd_lo"], r["sd_hi"], 0.0, 0.05)
            P.append({"id": "P-1", "quantity": f"SD riproducibile dei fattori non periodici ({ax})", "unit": k, "estimate": r["sd_repro"], "lo": r["sd_lo"],
                      "hi": r["sd_hi"], "interval": "[0, 0.05]", "verdict": vd, "extra": f"parita' {r['odd_over_even']:.3f}; r fra meta' {r['corr_halves']:.3f}; passo 1 {r['lag1_repro']:.3f}"})
            rb["nonperiodic"][f"{k}_{ax}"] = {"verdict": vd, "sd": r["sd_repro"], "lag1": r["lag1_repro"]}
        for ax in ("cols", "rows"):
            r = js["tiles"][ax]; vd = verdict(r["sd_lo"], r["sd_hi"], 0.0, 0.05)
            P.append({"id": "P-2", "quantity": f"SD riproducibile fra tessere del log dispari/pari ({ax})", "unit": k, "estimate": r["sd_between_repro"],
                      "lo": r["sd_lo"], "hi": r["sd_hi"], "interval": "[0, 0.05]", "verdict": vd, "extra": f"mediana {r['median_ratio']:.3f} (p10 {r['p10_ratio']:.3f}, p90 {r['p90_ratio']:.3f}); tessere {js['tiles']['n_tiles']}"})
            rb["parity_spatial"][f"{k}_{ax}"] = vd
    v1 = [x["verdict"] for x in rb["nonperiodic"].values()]
    rb["nonperiodic_layer"] = "fattori casuali per riga e colonna" if "smentita" in v1 else ("solo parita'" if all(v == "confermata" for v in v1) else "opzione spenta")
    rb["parity_spatial_sections"] = sorted({k.split("_")[0] for k, v in rb["parity_spatial"].items() if v == "smentita"})
    json.dump(rb, open(f"{RES}/S2.0_RB.json", "w"), indent=1)
    # P-3: controllo interno della regola C1 (due popolazioni)
    I = pd.read_csv(f"{RES}/S2.0_inner_summary.csv")
    for A_ in ("A2", "A5"):
        x = I[(I.archetype == A_) & (I.rule == "C1")].iloc[0]
        P.append({"id": "P-3", "quantity": "previsto/osservato senza i bin entro 2 µm dal bordo del nucleo, regola C1", "unit": A_, "estimate": x.ratio_mean,
                  "lo": x.lo, "hi": x.hi, "interval": "[0.909, 1.1]", "verdict": x.verdict,
                  "extra": f"ROI non invertibili {int(x.n_not_invertible)}/5; sulle ROI riuscite {x.ratio_mean_ok_only:.3f} [{x.lo_ok_only:.3f}, {x.hi_ok_only:.3f}]"})
    # M3: P-4 e R-C
    D = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/disp/*.json"))])
    D.drop(columns=["sec", "psi_sim_reps"], errors="ignore").to_csv(f"{RES}/S2.0_disp_roi.csv", index=False)
    rc = {}
    for A_ in ("A1", "A2", "A3", "A4", "A5", "A6"):
        t = D[D.archetype == A_]
        m, lo, hi = boot_mean(t.ratio.to_numpy(float)); vd = verdict(lo, hi, 0.909, 1.1)
        vm, vlo, vhi = boot_mean(t.v.to_numpy(float))
        P.append({"id": "P-4", "quantity": "indice di dispersione dentro il nucleo, osservato/simulato", "unit": A_, "estimate": m, "lo": lo, "hi": hi,
                  "interval": "[0.909, 1.1]", "verdict": vd, "extra": f"v {vm:.4f} [{vlo:.4f}, {vhi:.4f}]; psi oss {t.psi_obs.mean():.3f}, sim {t.psi_sim.mean():.3f}; r da {', '.join(sorted(set(t.r_source)))}"})
        rc[A_] = {"verdict": vd, "ratio": m, "v": vm, "v_lo": vlo, "v_hi": vhi}
    vs = [x["verdict"] for x in rc.values()]
    over = [a for a, x in rc.items() if x["verdict"] == "smentita" and x["ratio"] > 1.1]
    rc["layer"] = "rumore gamma dentro la cellula, v per tessuto" if over else ("nessuno strato" if all(v == "confermata" for v in vs) else "opzione spenta")
    rc["tissues_over"] = over
    json.dump(rc, open(f"{RES}/S2.0_RC.json", "w"), indent=1)
    # CP-1: lumi con la regola vincente
    ra = json.load(open(f"{RES}/S2.0_RA.json")); win = ra["winner"]
    rows = []; cp_fail = 0
    for k in SEC:
        for rl in ("C0", win):
            f = f"{OUT}/lumens/{k}_{rl}.csv"
            if not os.path.exists(f): continue
            L = pd.read_csv(f)
            if len(L) == 0: continue
            g = L.groupby("tile").sum(numeric_only=True); n = len(g)      # come R5 (tools/R5_analysis.py, modo lumens): bootstrap per tessere
            import R5_common as C5
            for lo_, hi_ in ((2, 4), (4, 8), (8, 16), (16, 32)):
                o, p = g[f"obs_{lo_}_{hi_}"].to_numpy(), g[f"pred_{lo_}_{hi_}"].to_numpy()
                rng = np.random.default_rng([C5.SEED, SEC.index(k), 40 + lo_]); idx = rng.integers(0, n, size=(C5.B, n))
                reps = np.array([o[ix].sum() / max(p[ix].sum(), 1e-12) for ix in idx])
                est = o.sum() / p.sum(); lo, hi = float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5))
                vd = verdict(lo, hi, 0.8, 1.25)
                judged = (rl == win and k != "A5" and hi_ <= 8)
                if judged and vd == "smentita": cp_fail += 1
                rows.append({"key": k, "rule": rl, "band": f"{lo_}-{hi_}", "n_lumens": len(L), "n_tiles": n, "obs_over_pred": est, "lo": lo, "hi": hi, "verdict": vd, "judged_CP1": judged})
    Lm = pd.DataFrame(rows); Lm.to_csv(f"{RES}/S2.0_CP1_lumens.csv", index=False)
    P.append({"id": "CP-1", "quantity": f"lumi interni, osservati/previsti con {win}, fasce 2-4 e 4-8 µm", "unit": "A1, A23, A4, A6", "estimate": float(cp_fail),
              "lo": np.nan, "hi": np.nan, "interval": "0 celle smentite su 8", "verdict": "superata" if cp_fail == 0 else "fallita", "extra": f"celle smentite {cp_fail}"})
    pd.DataFrame(P).to_csv(f"{RES}/S2.0_predictions.csv", index=False)
    print(pd.DataFrame(P)[["id", "unit", "estimate", "lo", "hi", "verdict"]].round(4).to_string(index=False))
    print(json.dumps({"RB": rb["nonperiodic_layer"], "RB_spatial": rb["parity_spatial_sections"], "RC": rc["layer"], "RC_over": over}))

if __name__ == "__main__":
    {"edges": stage_edges, "inner": stage_inner, "final": stage_final}[STAGE]()
