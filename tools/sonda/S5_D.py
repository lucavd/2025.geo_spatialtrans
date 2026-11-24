"""tools/sonda/S5_D.py - S5 (docs/sonda/REGOLA_S5.md): D_p sul pannello G_p con lo stimatore appaiato di S3c, funzioni di tools/sonda/S3c_D.py
senza modifiche (load, q_paired, D_of, boot_idx, intervals, compare, synthetic, cp_perm). Tabella results/sonda/S5/area_values.csv (genes = panel),
h_p results/sonda/S5/h.csv. Bracci appaiati a Vp per area k: T (Tp), B (Bp, media dei 2 semi), P (Pp), B0, B1, XB, XV. M-4 si toglie per un
tessuto se la media dei singoletti di Vp e' >= 0,99 (deciso da S3c_D.load sulla tabella di Vp prima di D). Bootstrap per coppie stratificato per
tessuto, stessi indici per tutti i bracci (default_rng(20261123), B = 2000); intervallo (ii) con f = sqrt(n/(n-1)) t(0,975; n-1) / 1,96 = 1,584.
Confronti con la regola di S3c (uguale se |Delta| + w <= 1 h, diverso se |Delta| - w > 1 h). A3 = D_p(XB) - D_p(B), w e se dalla distribuzione
bootstrap della differenza; regge se A3 - 2 se >= 1. Cancelli G-MET: C-1p (S3c_D.synthetic, 1000 repliche, seme 20261127), C-4 (P contro T
«diverso», P peggio), C-5 (B0 contro B1 non «diverso»). Controprova CP-XB: D_p(XB) < D_p(P) e XB contro P «diverso». Attese P1-P5 (come S4:
confermata se l'intervallo sta dentro, smentita se sta tutto fuori, altrimenti indecisa). Descrittivi: XB contro XV, XV contro Vp, D di XB con XV
come riferimento al posto di V, D sui geni completi di S3c. Uscite: results/sonda/S5/D.json, results/sonda/S5/per_measure.csv.
Uso: .venv/bin/python -B tools/sonda/S5_D.py [--nrep 1000]"""
import os, sys, json, argparse, numpy as np, pandas as pd
from scipy import stats
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); sys.path.insert(0, "tools/sonda")
import S3c_D as SC
TAB, HF, OUT, PER = "results/sonda/S5/area_values.csv", "results/sonda/S5/h.csv", "results/sonda/S5/D.json", "results/sonda/S5/per_measure.csv"
NEW = ("XB", "XV"); ARMS = ("T", "B", "P", "B0", "B1", "XB", "XV")
LAB = {"V": "Vp", "T": "Tp", "B": "Bp", "P": "Pp", "B0": "B0p", "B1": "B1p", "XB": "XB", "XV": "XV"}
PAIRS = [("B", "T"), ("P", "T"), ("B0", "B1"), ("P", "B"), ("XB", "B"), ("XB", "T"), ("XB", "P"), ("XB", "XV"), ("XV", "T"), ("XV", "B")]

def load_all():
    data, drop, Vs = SC.load(TAB, HF); U = pd.read_csv(TAB); assert set(U.genes) == {"panel"}
    U["area"] = U.unit.str.extract(r"_(tumour|stroma)_(\d+)$")[1].astype(int)
    for A in SC.TIS:
        keep = [m for m in SC.MEAS if SC.MID[m] in data[A]["meas"]]
        for a in NEW:
            u = U[(U.A == A) & (U.arm == a)]; x = u.groupby("area")[keep].mean()
            assert list(x.index) == data[A]["areas"] and len(u) == 5, (A, a, list(x.index)); data[A]["X"][a] = x.to_numpy(float)
    return data, drop, Vs

def judge(ci, lo, hi):
    if lo <= ci[0] and ci[1] <= hi: return "confermata"
    if ci[1] < lo or ci[0] > hi: return "smentita"
    return "indecisa"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--nrep", type=int, default=1000); a = ap.parse_args()
    data, drop, Vs = load_all(); n = 5; f = float(np.sqrt(n / (n - 1)) * stats.t.ppf(0.975, n - 1) / stats.norm.ppf(0.975)); assert abs(f - 1.584) < 5e-4
    idx = SC.boot_idx(np.random.default_rng(SC.SEED), SC.NB)
    res = {"rule": "docs/sonda/REGOLA_S5.md", "estimator": "paired (S3c_D)", "genes": "panel (G_p, 284)", "method": "ii", "f": f, "B_boot": SC.NB, "seed": SC.SEED,
           "arms": list(ARMS), "labels": LAB, "drop_M4": drop, "V_singlets": Vs, "areas": {A: data[A]["areas"] for A in SC.TIS},
           "measures_kept": {A: data[A]["meas"] for A in SC.TIS}, "by": {}}
    for g, ts in SC.GROUPS.items():
        Dobs = {k: float(v) for k, v in SC.D_of(data, ARMS, ts).items()}; Dbs = SC.D_of(data, ARMS, ts, idx)
        ci = {k: SC.intervals(Dobs[k], Dbs[k], f)[1] for k in ARMS}
        comps = [SC.compare(Dobs, Dbs, x, y, f) for x, y in PAIRS]
        d = Dobs["XB"] - Dobs["B"]; x = Dbs["XB"] - Dbs["B"]; se = float(f * np.std(x, ddof=1)); plo, phi = np.percentile(x, [2.5, 97.5]); w = float(f * (phi - plo) / 2)
        lo, hi = SC.intervals(d, x, f)[1]
        A3 = {"value": float(d), "se": se, "w": w, "lo": lo, "hi": hi, "A3_minus_2se": float(d - 2 * se), "regge": bool(d - 2 * se >= 1.0)}
        # descrittivo: XB con XV come riferimento (stesso stimatore, stessi indici, h_p di Vp)
        qs = [SC.q_paired(data[A]["X"]["XB"], data[A]["X"]["XV"], data[A]["h"]) for A in ts]; dx = float(SC.ssqrt(np.concatenate(qs).mean()))
        qb = [SC.q_paired(data[A]["X"]["XB"][idx[A]], data[A]["X"]["XV"][idx[A]], data[A]["h"]) for A in ts]; dxb = SC.ssqrt(np.concatenate(qb, axis=-1).mean(-1))
        res["by"][g] = {"D": Dobs, "D_ci": ci, "D_se": {k: float(f * np.std(Dbs[k], ddof=1)) for k in ARMS}, "comparisons": comps, "A3": A3,
                        "D_XB_ref_XV": {"D": dx, "ci": SC.intervals(dx, dxb, f)[1], "se": float(f * np.std(dxb, ddof=1))},
                        "XB_minus_T": float(Dobs["XB"] - Dobs["T"])}
    m = res["by"]["tumour+stroma"]; cmp = {(c["a"], c["b"]): c for c in m["comparisons"]}
    syn = SC.synthetic(data, f, a.nrep)
    res["C1p"] = {"values": {k: v["root"] for k, v in syn["scenarios"].items()}, "mean_square_0": syn["scenarios"]["0"]["mean_square"]["T"], "pass": syn["C1p_pass"],
                  "coverage_min": syn["min_coverage"], "coverage_rule_would_choose": syn["method"], "nrep": a.nrep}
    res["C4"] = {"comparison": cmp[("P", "T")], "pass": bool(cmp[("P", "T")]["outcome"] == "different" and cmp[("P", "T")]["better"] == "T")}
    res["C5"] = {"comparison": cmp[("B0", "B1")], "pass": bool(cmp[("B0", "B1")]["outcome"] in ("equal", "undetermined"))}
    res["G_MET_pass"] = bool(res["C1p"]["pass"] and res["C4"]["pass"] and res["C5"]["pass"])
    res["CP_XB"] = {"D_XB": m["D"]["XB"], "D_P": m["D"]["P"], "comparison": cmp[("XB", "P")],
                    "pass": bool(m["D"]["XB"] < m["D"]["P"] and cmp[("XB", "P")]["outcome"] == "different" and cmp[("XB", "P")]["better"] == "XB")}
    res["CP_perm"] = SC.cp_perm(data)
    res["A3"] = m["A3"]
    res["objective_level"] = {"D_XB_minus_D_T": m["XB_minus_T"], "alarm_A3_unreachable": bool(m["XB_minus_T"] <= 1.0),
                              "note": "se D_p(XB) - D_p(T) <= 1 h nessun simulatore batte lo Xenium di 1 h: allarme sul criterio, decide Luca"}
    res["ceiling_floor_panel"] = {"tetto_D_T": m["D"]["T"], "tetto_ci": m["D_ci"]["T"], "pavimento_D_B": m["D"]["B"], "pavimento_ci": m["D_ci"]["B"]}
    th = json.load(open("results/sonda/S5/thinning.json")); fz = {t: v["f"] for t, v in th["by_tissue"].items()}
    P = {}
    P["P1"] = {"observed": fz, "verdict": "confermata" if all(0.1 <= v <= 0.6 for v in fz.values()) else "smentita"}
    P["P2"] = {"observed": {"D": m["D"]["T"], "ci": m["D_ci"]["T"]}, "verdict": judge(m["D_ci"]["T"], 0.0, 1.0)}
    P["P3"] = {"observed": {"D": m["D"]["B"], "ci": m["D_ci"]["B"]}, "verdict": judge(m["D_ci"]["B"], 2.0, 4.0)}
    P["P4"] = {"observed": {"D": m["D"]["XB"], "ci": m["D_ci"]["XB"]}, "verdict": judge(m["D_ci"]["XB"], 1.0, 3.0)}
    P["P5"] = {"observed": {"A3": m["A3"]["value"], "ci": [m["A3"]["lo"], m["A3"]["hi"]], "regge": m["A3"]["regge"]},
               "verdict": "confermata" if (m["A3"]["hi"] < 0 and not m["A3"]["regge"]) else "smentita"}
    res["predictions"] = P
    try:
        s3c = json.load(open("results/sonda/S3c/D.json"))["by"]
        res["full_genes_S3c"] = {g: {"D": s3c[g]["D"], "D_ci": s3c[g]["D_ci"]} for g in s3c}
    except Exception as e: res["full_genes_S3c"] = {"error": repr(e)}
    rows = []
    for A in SC.TIS:
        X = data[A]["X"]; h = data[A]["h"]; q = {k: SC.q_paired(X[k], X["V"], h) for k in ARMS}
        for j, mm in enumerate(data[A]["meas"]):
            r = {"A": A, "measure": mm, "genes": "panel", "h": float(h[j]), "V_mean": float(X["V"][:, j].mean())}
            for k in ARMS: r[f"{k}_mean"] = float(X[k][:, j].mean()); r[f"{k}_q"] = float(q[k][j]); r[f"{k}_d"] = float(SC.ssqrt(q[k][j]))
            r["XB_minus_B_q"] = float(q["XB"][j] - q["B"][j]); rows.append(r)
    Pm = pd.DataFrame(rows); Pm.to_csv(PER, index=False)
    res["per_measure_top"] = {k: Pm.sort_values(f"{k}_q", ascending=False).head(4)[["A", "measure", f"{k}_q"]].to_dict("records") for k in ("XB", "B", "XV", "T")}
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    for g, r in res["by"].items(): print(g, "D", {LAB[k]: round(v, 3) for k, v in r["D"].items()}, "A3", round(r["A3"]["value"], 3), "se", round(r["A3"]["se"], 3))
    print("CI", {LAB[k]: [round(u, 3) for u in v] for k, v in m["D_ci"].items()})
    for c in m["comparisons"]: print("  ", c["a"], "-", c["b"], round(c["delta_ab"], 3), "w", round(c["w"], 3), c["outcome"], c["better"])
    print("C1p", res["C1p"]["pass"], "C4", res["C4"]["pass"], "C5", res["C5"]["pass"], "CP-XB", res["CP_XB"]["pass"], "drop", drop, "P", {k: v["verdict"] for k, v in P.items()})
