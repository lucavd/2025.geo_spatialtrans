"""tools/sonda/S4_D.py - S4: D dei bracci a componenti veri (docs/sonda/REGOLA_S4.md) con lo stimatore appaiato di S3c, invariato
(tools/sonda/S3c_D.py: load, q_paired, D_of, boot_idx, intervals, compare). Bracci di S3c (T, B, P, B0, B1) dalla tabella di S3b; bracci
nuovi (Y, G, U, A, AE, un seme) dalla tabella di S4; stessi indici bootstrap di S3c (default_rng(20261123), B = 2000), stesso metodo
d'intervallo (quello scelto in results/sonda/S3c/D.json). Cancello G-S3c: D e intervalli di T, B, P, B0, B1 identici a S3c (1e-12).
Uscite: results/sonda/S4/D.json, results/sonda/S4/per_measure.csv. Uso: .venv/bin/python -B tools/sonda/S4_D.py"""
import os, sys, json, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); sys.path.insert(0, "tools/sonda")
import S3c_D as SC
NEW = ("Y", "G", "U", "A", "AE"); T3, H3, T4 = "results/sonda/S3b/area_values.csv", "results/sonda/S3b/h.csv", "results/sonda/S4/area_values.csv"
REF = "results/sonda/S3c/D.json"; OUT, PER = "results/sonda/S4/D.json", "results/sonda/S4/per_measure.csv"

def load_all():
    data, drop, Vs = SC.load(T3, H3); U = pd.read_csv(T4)
    U["area"] = U.unit.str.extract(r"_(tumour|stroma)_(\d+)$")[1].astype(int); present = [a for a in NEW if (U.arm == a).any()]
    for A in SC.TIS:
        keep = [m for m in SC.MEAS if SC.MID[m] in data[A]["meas"]]
        for a in present:
            x = U[(U.A == A) & (U.arm == a)].groupby("area")[keep].mean()
            assert list(x.index) == data[A]["areas"] and (U[(U.A == A) & (U.arm == a)].shape[0] == 5), (A, a, list(x.index))
            data[A]["X"][a] = x.to_numpy(float)
    return data, present

if __name__ == "__main__":
    data, present = load_all(); ref = json.load(open(REF)); meth = ref["COV"]["method"]; f = ref["f"]; fac = 1.0 if meth == "i" else f
    arms = list(SC.ARMS) + present; idx = SC.boot_idx(np.random.default_rng(SC.SEED), SC.NB)
    res = {"rule": "docs/sonda/REGOLA_S4.md", "method": meth, "f": f, "arms": arms, "by": {}}
    gs3c = []
    for g, ts in SC.GROUPS.items():
        Dobs = {k: float(v) for k, v in SC.D_of(data, arms, ts).items()}; Dbs = SC.D_of(data, arms, ts, idx)
        ci = {k: SC.intervals(Dobs[k], Dbs[k], f)[0 if meth == "i" else 1] for k in arms}
        for k in SC.ARMS:
            gs3c.append(abs(Dobs[k] - ref["by"][g]["D"][k])); gs3c += [abs(u - v) for u, v in zip(ci[k], ref["by"][g]["D_ci"][k])]
        pairs = [(a, b) for a in present for b in ("B", "B0", "T")] + [(a, b) for a, b in (("A", "Y"), ("AE", "A")) if a in present and b in present]
        comps = [SC.compare(Dobs, Dbs, a, b, fac) for a, b in pairs]
        qm = {k: np.sign(Dbs[k]) * Dbs[k] ** 2 for k in arms}; q0 = {k: np.sign(Dobs[k]) * Dobs[k] ** 2 for k in arms}
        closed = {}
        for a in present:
            fr = (q0["B"] - q0[a]) / (q0["B"] - q0["T"]); frb = (qm["B"] - qm[a]) / (qm["B"] - qm["T"])
            lo, hi = np.percentile(frb, [2.5, 97.5]); closed[a] = {"fraction": float(fr), "pct_interval": [float(lo), float(hi)]}
        res["by"][g] = {"D": Dobs, "D_ci": ci, "D_se": {k: float(fac * np.std(Dbs[k], ddof=1)) for k in arms}, "comparisons": comps,
                        "gap_closed_fraction": closed}
    res["G_S3c_identity"] = {"max_abs_dev": float(max(gs3c)), "pass": bool(max(gs3c) <= 1e-12)}
    rows = []
    for A in SC.TIS:
        X = data[A]["X"]; h = data[A]["h"]; q = {a: SC.q_paired(X[a], X["V"], h) for a in arms}
        for j, m in enumerate(data[A]["meas"]):
            r = {"A": A, "measure": m, "h": float(h[j]), "V_mean": float(X["V"][:, j].mean())}
            for a in arms: r[f"{a}_mean"] = float(X[a][:, j].mean()); r[f"{a}_q"] = float(q[a][j])
            gap = q["B"][j] - q["T"][j]
            for a in present: r[f"{a}_closed_q"] = float(q["B"][j] - q[a][j]); r[f"{a}_closed_of_gap"] = float((q["B"][j] - q[a][j]) / gap) if abs(gap) > 1e-9 else None
            r["B_minus_T_q"] = float(gap); rows.append(r)
    P = pd.DataFrame(rows).sort_values("B_minus_T_q", ascending=False); P.to_csv(PER, index=False)
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    m = res["by"]["tumour+stroma"]
    print("G-S3c", res["G_S3c_identity"]); print("D", {k: round(v, 3) for k, v in m["D"].items()})
    print("CI", {k: [round(u, 3) for u in v] for k, v in m["D_ci"].items()})
    print("closed", {k: (round(v["fraction"], 3), [round(u, 3) for u in v["pct_interval"]]) for k, v in m["gap_closed_fraction"].items()})
    for c in m["comparisons"]: print("  ", c["a"], "-", c["b"], round(c["delta_ab"], 3), "w", round(c["w"], 3), c["outcome"], c["better"])
