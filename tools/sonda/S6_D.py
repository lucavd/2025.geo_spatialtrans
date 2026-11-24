"""tools/sonda/S6_D.py - S6: D dei bracci con i tipi non copiati (docs/sonda/REGOLA_S6.md), stimatore appaiato di S3c invariato
(tools/sonda/S3c_D.py). Bracci di S3c dalla tabella di S3b, di S4 (Y, G, U, A, AE) dalla tabella di S4, nuovi (Yalt, YH1, Yshuf) dalla tabella
di S6; stessi indici bootstrap (default_rng(20261123), B = 2000), stesso metodo d'intervallo di results/sonda/S3c/D.json.
Cancello G-S4: D e intervalli dei bracci di S3c e di S4 identici a results/sonda/S4/D.json (1e-12).
Esito della regola di docs/riesame/RIESAME_C.md par. 6 con delta = 1 (pavimento - delta = 2,151) e con delta = 0,5 (2,651).
Uscite: results/sonda/S6/D.json, results/sonda/S6/per_measure.csv. Uso: .venv/bin/python -B tools/sonda/S6_D.py"""
import os, sys, json, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); sys.path.insert(0, "tools/sonda")
import S3c_D as SC
T3, H3 = "results/sonda/S3b/area_values.csv", "results/sonda/S3b/h.csv"; T4, T6 = "results/sonda/S4/area_values.csv", "results/sonda/S6/area_values.csv"
REF3, REF4 = "results/sonda/S3c/D.json", "results/sonda/S4/D.json"; OUT, PER = "results/sonda/S6/D.json", "results/sonda/S6/per_measure.csv"
OLD, NEW = ("Y", "G", "U", "A", "AE"), ("Yalt", "YH1", "Yshuf")


def add(data, U, names):
    U = U.copy(); U["area"] = U.unit.str.extract(r"_(tumour|stroma)_(\d+)$")[1].astype(int); present = [a for a in names if (U.arm == a).any()]
    for A in SC.TIS:
        keep = [m for m in SC.MEAS if SC.MID[m] in data[A]["meas"]]
        for a in present:
            x = U[(U.A == A) & (U.arm == a)].groupby("area")[keep].mean()
            assert list(x.index) == data[A]["areas"] and (U[(U.A == A) & (U.arm == a)].shape[0] == 5), (A, a, list(x.index))
            data[A]["X"][a] = x.to_numpy(float)
    return present


def outcome(ci, thr):
    """regola di RIESAME_C par. 6: 1 statistiche bastano, 2 serve una sezione vera, 3 guadagno copiato, 4 indeterminato."""
    la, ua = ci["Yalt"]; lh, uh = ci["YH1"]
    if ua <= thr: return 1
    if uh <= thr: return 2
    if la > thr and lh > thr: return 3
    return 4


if __name__ == "__main__":
    data, drop, Vs = SC.load(T3, H3); p4 = add(data, pd.read_csv(T4), OLD); p6 = add(data, pd.read_csv(T6), NEW) if os.path.exists(T6) else []
    r3 = json.load(open(REF3)); r4 = json.load(open(REF4)); meth = r3["COV"]["method"]; f = r3["f"]; fac = 1.0 if meth == "i" else f
    arms = list(SC.ARMS) + p4 + p6; idx = SC.boot_idx(np.random.default_rng(SC.SEED), SC.NB)
    res = {"rule": "docs/sonda/REGOLA_S6.md", "method": meth, "f": f, "arms": arms, "by": {}}; dev = []
    for g, ts in SC.GROUPS.items():
        Dobs = {k: float(v) for k, v in SC.D_of(data, arms, ts).items()}; Dbs = SC.D_of(data, arms, ts, idx)
        ci = {k: SC.intervals(Dobs[k], Dbs[k], f)[0 if meth == "i" else 1] for k in arms}
        for k in list(SC.ARMS) + p4:
            dev.append(abs(Dobs[k] - r4["by"][g]["D"][k])); dev += [abs(u - v) for u, v in zip(ci[k], r4["by"][g]["D_ci"][k])]
        pairs = [(a, b) for a in p6 for b in ("Y", "B", "B0", "T")] + [(a, b) for a, b in (("YH1", "Yalt"), ("Yshuf", "Yalt")) if a in p6 and b in p6]
        comps = [SC.compare(Dobs, Dbs, a, b, fac) for a, b in pairs]
        qm = {k: np.sign(Dbs[k]) * Dbs[k] ** 2 for k in arms}; q0 = {k: np.sign(Dobs[k]) * Dobs[k] ** 2 for k in arms}; closed = {}
        for a in p4 + p6:
            fr = (q0["B"] - q0[a]) / (q0["B"] - q0["T"]); frb = (qm["B"] - qm[a]) / (qm["B"] - qm["T"])
            lo, hi = np.percentile(frb, [2.5, 97.5]); closed[a] = {"fraction": float(fr), "pct_interval": [float(lo), float(hi)]}
        res["by"][g] = {"D": Dobs, "D_ci": ci, "D_se": {k: float(fac * np.std(Dbs[k], ddof=1)) for k in arms}, "comparisons": comps, "gap_closed_fraction": closed}
    res["G_S4_identity"] = {"max_abs_dev": float(max(dev)), "pass": bool(max(dev) <= 1e-12)}
    m = res["by"]["tumour+stroma"]; Dfloor = m["D"]["B"]
    if all(a in p6 for a in ("Yalt", "YH1")):
        res["rule_outcome"] = {f"delta_{dl}": {"threshold": float(Dfloor - dl), "case": outcome(m["D_ci"], Dfloor - dl)} for dl in (1.0, 0.5)}
    rows = []
    for A in SC.TIS:
        X = data[A]["X"]; h = data[A]["h"]; q = {a: SC.q_paired(X[a], X["V"], h) for a in arms}
        for j, mm in enumerate(data[A]["meas"]):
            r = {"A": A, "measure": mm, "h": float(h[j]), "V_mean": float(X["V"][:, j].mean())}
            for a in arms: r[f"{a}_mean"] = float(X[a][:, j].mean()); r[f"{a}_q"] = float(q[a][j])
            r["B_minus_T_q"] = float(q["B"][j] - q["T"][j]); rows.append(r)
    pd.DataFrame(rows).sort_values("B_minus_T_q", ascending=False).to_csv(PER, index=False)
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    print("G-S4", res["G_S4_identity"]); print("D", {k: round(v, 3) for k, v in m["D"].items()})
    print("CI", {k: [round(u, 3) for u in v] for k, v in m["D_ci"].items() if k in p4 + p6 + ["B", "T"]})
    print("closed", {k: round(v["fraction"], 3) for k, v in m["gap_closed_fraction"].items()}); print("rule", res.get("rule_outcome"))
    for c in m["comparisons"]: print("  ", c["a"], "-", c["b"], round(c["delta_ab"], 3), "w", round(c["w"], 3), c["outcome"], c["better"])
