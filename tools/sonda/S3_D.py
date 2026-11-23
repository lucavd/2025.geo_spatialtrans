"""tools/sonda/S3_D.py - S3: misura principale D sul polmone, pavimento, tetto, soglia e controlli del metro (docs/sonda/REGOLA_S3.md).
Usa le funzioni invariate di tools/riesame_D.py (prepare, analyse, compare; B = 2000, seme 20261123; delta = 1,0 h) con T = H1 (un insieme).
Bracci: B (media dei 2 semi per area), P (controllo positivo C-4), B0 e B1 (B con il solo seme 0 o 1, controllo negativo C-5).
C-1: prova sintetica di G-3 (REGOLA_B) con l'insieme H1; C-2: lo stesso percorso sulla tabella del colon ritrova B_D_common.json.
M-4 si toglie per un tessuto se la media dei singoletti di V e' >= 0,99 (decisione sulla sola tabella di V).
Uscite: results/sonda/S3/D.json, per_measure.csv; con --misura anche results/misura.json (A1 = D(B), se = deviazione standard bootstrap), solo se C-1, C-2, C-4, C-5 passano;\naltrimenti results/sonda/S3/misura_provvisoria.json (REGOLA_S3, «Esito e conseguenze»).
Uso: .venv/bin/python -B tools/sonda/S3_D.py [--table T] [--h H] [--out D.json] [--per per_measure.csv] [--misura results/misura.json] [--skip-c2]"""
import os, sys, json, argparse, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); import riesame_D as RD
MEAS = RD.MEAS; MID = RD.MID; TS = ("H1",)

def synth(cells, tsets, reps=200, seed=20261124):
    """G-3 di REGOLA_B (come RD.synth_check) con gli insiemi di T dati."""
    rng = np.random.default_rng(seed); out = {}
    for shift in (0.0, 1.0, 2.0):
        ests = []
        for _ in range(reps):
            cs = []
            for c in cells:
                sV = c["V"].std(ddof=1); V = rng.normal(0, sV, len(c["V"]))
                T = {k: rng.normal(shift * c["h"], c["sp"], len(o)) for k, o in c["T"].items()}
                sp, df = RD.pooled_sd([V] + list(T.values())); cs.append({**c, "V": V, "T": T, "sp": sp, "df": df})
            ests.append(RD.q_arm(cs, "T", tsets))
        out[str(shift)] = {"mean_square": float(np.mean(ests)), "root": RD.ssqrt(np.mean(ests))}
    ok = abs(out["0.0"]["mean_square"]) <= 0.15 and abs(out["1.0"]["root"] - 1) <= 0.25 and abs(out["2.0"]["root"] - 2) <= 0.25
    return out, bool(ok)

def ci_of(Dobs, Dbs, k):
    if np.isnan(Dobs[k]): return None
    sh = np.nanmedian(Dbs[k]) - Dobs[k]; return [float(np.nanpercentile(Dbs[k], 2.5) - sh), float(np.nanpercentile(Dbs[k], 97.5) - sh)]

def colon_identity():
    """C-2: tabella del colon (geni comuni) con lo stesso percorso: D di T, S, B contro results/riesame/B_D_common.json."""
    tab = pd.read_csv("results/riesame/B_area_values.csv"); H = pd.read_csv("results/riesame/B_h_panel.csv"); ref = json.load(open("results/riesame/B_D_common.json"))
    arms = sorted(set(tab[tab.genes == "common"].arm) - {"V", "T"}); out = {}; worst = 0.0
    for g, ts in {"A2+A3": ["A2", "A3"], "A2": ["A2"], "A3": ["A3"]}.items():
        cells = RD.prepare(tab, H, "common", ts, {"A2": False, "A3": False}); Dobs, _ = RD.analyse(cells, arms)
        for k in ("T", "S", "B"):
            r = float(ref["by"][g]["D"][k]); rel = abs(Dobs[k] - r) / abs(r); worst = max(worst, rel)
            out[f"{g}:{k}"] = {"recomputed": float(Dobs[k]), "reference": r, "rel_dev": rel}
    return {"values": out, "max_rel_dev": worst, "pass": bool(worst <= 1e-9)}

def per_measure(tab, H, drop):
    rows = []
    for A in ("tumour", "stroma"):
        for c in RD.prepare(tab, H, "all", [A], drop):
            nV = len(c["V"]); r = {"A": A, "measure": c["m"], "h": c["h"], "V_mean": float(c["V"].mean()), "s_pooled": c["sp"]}
            r["T_q"] = float(np.mean([((o.mean() - c["V"].mean()) ** 2 - c["sp"] ** 2 / len(o) - c["sp"] ** 2 / nV) / c["h"] ** 2 for o in c["T"].values()]))
            r["T_mean"] = float(np.mean([o.mean() for o in c["T"].values()]))
            for a in ("B", "P"):
                if a in c["arms"]:
                    x = c["arms"][a]; r[f"{a}_mean"] = float(x.mean()); r[f"{a}_q"] = float(((x.mean() - c["V"].mean()) ** 2 - c["sp"] ** 2 / nV - x.var(ddof=1) / len(x)) / c["h"] ** 2)
            rows.append(r)
    P = pd.DataFrame(rows)
    for a in ("T", "B", "P"):
        if f"{a}_q" in P: P[f"{a}_d"] = np.sign(P[f"{a}_q"]) * np.sqrt(np.abs(P[f"{a}_q"]))
    P["B_minus_T_q"] = P["B_q"] - P["T_q"]; tot = P.B_minus_T_q.sum(); P["share_of_gap"] = P.B_minus_T_q / tot
    return P.sort_values("B_minus_T_q", ascending=False).reset_index(drop=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--table", default="results/sonda/S3/area_values.csv"); ap.add_argument("--h", default="results/sonda/S3/h.csv")
    ap.add_argument("--out", default="results/sonda/S3/D.json"); ap.add_argument("--per", default="results/sonda/S3/per_measure.csv")
    ap.add_argument("--misura", default=None); ap.add_argument("--provv", default="results/sonda/S3/misura_provvisoria.json"); ap.add_argument("--skip-c2", action="store_true"); a = ap.parse_args()
    tab = pd.read_csv(a.table); H = pd.read_csv(a.h); tab["set"] = tab["set"].astype(str)
    B = tab[tab.arm == "B"]; tab = pd.concat([tab, B[B.k.astype(str) == "0"].assign(arm="B0"), B[B.k.astype(str) == "1"].assign(arm="B1")], ignore_index=True)
    Vs = tab[tab.arm == "V"].groupby("A")["M4"].mean(); drop = {A: bool(Vs.get(A, 0) >= 0.99) for A in ("tumour", "stroma")}
    arms = [x for x in ("B", "P", "B0", "B1") if x in set(tab.arm)]
    res = {"rule": "docs/sonda/REGOLA_S3.md", "delta": RD.DELTA, "B_boot": RD.B, "seed": RD.SEED, "drop_M4": drop, "V_singlets": Vs.to_dict(), "arms": arms,
           "n_rows": tab.groupby("arm").size().to_dict(), "by": {}}
    for g, ts in {"tumour+stroma": ["tumour", "stroma"], "tumour": ["tumour"], "stroma": ["stroma"]}.items():
        cells = RD.prepare(tab, H, "all", ts, drop); Dobs, Dbs = RD.analyse(cells, arms, tsets=TS)
        keys = ["T"] + arms
        comps = [RD.compare(Dobs, Dbs, x, y) for x, y in (("B", "T"), ("P", "T"), ("B0", "B1"), ("P", "B")) if x in Dobs and y in Dobs]
        res["by"][g] = {"D": {k: float(Dobs[k]) for k in keys}, "D_ci": {k: ci_of(Dobs, Dbs, k) for k in keys},
                        "D_se": {k: float(np.nanstd(Dbs[k], ddof=1)) for k in keys}, "comparisons": comps, "n_cells": len(cells)}
        if g == "tumour+stroma":
            g3, ok = synth(cells, TS); res["C1_synthetic"] = {"values": g3, "pass": ok}
    m = res["by"]["tumour+stroma"]; cmp = {(c["a"], c["b"]): c for c in m["comparisons"]}
    res["C4_positive"] = {"comparison": cmp.get(("P", "T")), "pass": bool(cmp.get(("P", "T"), {}).get("outcome") == "different" and cmp[("P", "T")]["better"] == "T")}
    res["C5_negative"] = {"comparison": cmp.get(("B0", "B1")), "pass": bool(cmp.get(("B0", "B1"), {}).get("outcome") in ("equal", "undetermined"))}
    if not a.skip_c2: res["C2_colon_identity"] = colon_identity()
    pav, tet = m["D"]["B"], m["D"]["T"]
    res["floor_ceiling"] = {"pavimento": pav, "tetto": tet, "soglia": tet + 1.0, "soglia_delta_0.5": tet + 0.5, "soglia_delta_2": tet + 2.0,
                            "floor_minus_ceiling": pav - tet, "G1_precondition_gt_delta": bool(pav - tet > 1.0), "B_vs_T": cmp.get(("B", "T")),
                            "A1_today": pav, "A1_se": m["D_se"]["B"], "A1_holds": bool(pav + 2 * m["D_se"]["B"] <= tet + 1.0)}
    P = per_measure(tab, H, drop); P.to_csv(a.per, index=False)
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    allpass = all(res.get(k, {}).get("pass", True) for k in ("C1_synthetic", "C2_colon_identity", "C4_positive", "C5_negative"))
    if a.misura and not allpass:   # REGOLA_S3: con un controllo fallito A1 e' provvisorio e non va in results/misura.json
        if os.path.exists(a.misura): os.remove(a.misura)
        a.misura = a.provv
    if a.misura:
        voci = {f"{r.A} {r.measure}": round(float(r.B_d), 4) for r in P.itertuples()}
        json.dump({"valore": round(pav, 6), "se": round(m["D_se"]["B"], 6), "voci": voci, "affermazioni": {},
                   "fonte": a.out, "tetto_misurato": round(tet, 6), "nota": "A1 = D(B): simulatore = catena banale (S3)", "controlli_passati": allpass}, open(a.misura, "w"), indent=1)
    print("D", {g: {k: round(v, 3) for k, v in r["D"].items()} for g, r in res["by"].items()})
    print("C1", res["C1_synthetic"]["pass"], {k: round(v["root"], 3) for k, v in res["C1_synthetic"]["values"].items()}, "C4", res["C4_positive"]["pass"], "C5", res["C5_negative"]["pass"],
          "C2", res.get("C2_colon_identity", {}).get("pass"), res.get("C2_colon_identity", {}).get("max_rel_dev"))
    for c in m["comparisons"]: print("  ", c["a"], "-", c["b"], round(c["delta_ab"], 3), "w", round(c["w"], 3), c["outcome"], c.get("better", ""))
    print("floor", round(pav, 3), "ceiling", round(tet, 3), "threshold", round(tet + 1.0, 3), "G1 precondition", res["floor_ceiling"]["G1_precondition_gt_delta"])
