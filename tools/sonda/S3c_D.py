"""tools/sonda/S3c_D.py - S3c: misura principale D con lo stimatore appaiato (docs/sonda/REGOLA_S3c.md).
Bracci T, B (media dei 2 semi per area), P, B0, B1 appaiati a V per area k. Per cella (tessuto x misura): d_k = a_k - V_k,
q = (mean(d)^2 - var(d)/n) / h^2; D(a) = radice con segno della media di q sulle celle.
Incertezza: bootstrap per coppie stratificato per tessuto (stessi indici per V e tutti i bracci), B = 2000, seme 20261123;
candidati (i) bootstrap e (ii) bootstrap con semi-ampiezze x f = sqrt(n/(n-1)) t(0.975, n-1) / 1.96, scelti da COV.
Cancelli: G-ID (media sulle 120 permutazioni = stimatore indipendente di S3_D/riesame_D), C-1p e COV (sintetico appaiato,
Sigma = 0.95 S + 0.05 diag S, 1000 repliche, seme 20261127), C-2 (colon, invariato), C-4, C-5. Controprova CP-perm.
Uscite: --out D.json, --per per_measure.csv, --misura results/misura.json solo se tutti i cancelli passano (altrimenti --provv).
Uso: .venv/bin/python -B tools/sonda/S3c_D.py [--table T] [--h H] [--out D.json] [--per per.csv] [--misura M] [--provv P] [--skip-c2] [--nrep 1000]"""
import os, sys, json, argparse, itertools, numpy as np, pandas as pd
from scipy import stats
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); sys.path.insert(0, "tools/sonda")
import riesame_D as RD
import S3_D
MEAS = ["M1", "M2", "M3", "M4", "M5", "M6"]; MID = {m: f"M-{m[1]}" for m in MEAS}
TIS = ("tumour", "stroma"); ARMS = ("T", "B", "P", "B0", "B1"); SYN_ARMS = ("T", "B", "P")
DELTA, NB, SEED, SYN_SEED = 1.0, 2000, 20261123, 20261127
GROUPS = {"tumour+stroma": TIS, "tumour": ("tumour",), "stroma": ("stroma",)}

def ssqrt(x): return np.sign(x) * np.sqrt(np.abs(x))

def load(table, hfile):
    """per tessuto: matrici (5 aree x misure tenute) per V e per i bracci, ordinate per area k; h; misure tenute."""
    U = pd.read_csv(table); H = pd.read_csv(hfile)
    U["area"] = U.unit.str.extract(r"_(tumour|stroma)_(\d+)$")[1].astype(int)
    Vs = U[U.arm == "V"].groupby("A")["M4"].mean(); drop = {A: bool(Vs.get(A, 0) >= 0.99) for A in TIS}
    out = {}
    for A in TIS:
        keep = [m for m in MEAS if not (drop[A] and m == "M4")]
        h = np.array([float(H[(H.A == A) & (H.measure == MID[m])].h.iloc[0]) for m in keep])
        u = U[U.A == A]; X = {}
        X["V"] = u[u.arm == "V"].groupby("area")[keep].mean()
        X["T"] = u[u.arm == "T"].groupby("area")[keep].mean()
        X["P"] = u[u.arm == "P"].groupby("area")[keep].mean()
        b = u[u.arm == "B"]; X["B"] = b.groupby("area")[keep].mean()
        X["B0"] = b[b.k.astype(str) == "0"].groupby("area")[keep].mean(); X["B1"] = b[b.k.astype(str) == "1"].groupby("area")[keep].mean()
        areas = list(X["V"].index)
        for a, x in X.items():
            assert list(x.index) == areas and len(areas) == 5, (A, a, list(x.index))
            assert u[(u.arm == a)].shape[0] in (5, 10) or a in ("B0", "B1")
        out[A] = {"X": {a: x.to_numpy(float) for a, x in X.items()}, "h": h, "meas": [MID[m] for m in keep], "areas": areas}
    return out, drop, Vs.to_dict()

def q_paired(Xa, Xv, h):
    """Xa, Xv: (..., n, m). Restituisce q (..., m)."""
    d = Xa - Xv; n = d.shape[-2]
    return (d.mean(-2) ** 2 - d.var(-2, ddof=1) / n) / h ** 2

def D_of(data, arms, tissues, idx=None):
    """D per braccio sulle celle dei tessuti dati; idx: {tessuto: (nb, n) indici} per il bootstrap -> vettori (nb,)."""
    res = {}
    for a in arms:
        qs = []
        for A in tissues:
            X = data[A]["X"]; h = data[A]["h"]
            if idx is None: qs.append(q_paired(X[a], X["V"], h))
            else: qs.append(q_paired(X[a][idx[A]], X["V"][idx[A]], h))
        res[a] = ssqrt(np.concatenate(qs, axis=-1).mean(-1))
    return res

def boot_idx(rng, nb, n=5): return {A: rng.integers(0, n, (nb, n)) for A in TIS}

def intervals(Dobs, Dbs, f):
    sh = np.median(Dbs) - Dobs; lo, hi = np.percentile(Dbs, [2.5, 97.5]) - sh
    return [float(lo), float(hi)], [float(Dobs + f * (lo - Dobs)), float(Dobs + f * (hi - Dobs))]

def compare(Dobs, Dbs, a, b, fac):
    d = Dobs[a] - Dobs[b]; lo, hi = np.percentile(Dbs[a] - Dbs[b], [2.5, 97.5]); w = fac * (hi - lo) / 2
    out = "equal" if abs(d) + w <= DELTA else ("different" if abs(d) - w > DELTA else "undetermined")
    return {"a": a, "b": b, "delta_ab": float(d), "w": float(w), "outcome": out, "better": (a if d < 0 else b) if out == "different" else ""}

def g_id(data, table, hfile):
    """G-ID: media di q appaiato di T sulle 120 permutazioni di k = termine di T di S3_D.per_measure (deviazione combinata)."""
    tab = pd.read_csv(table); H = pd.read_csv(hfile); tab["set"] = tab["set"].astype(str)
    drop = {A: len(data[A]["meas"]) < 6 for A in TIS}; worst = 0.0; rows = []
    for A in TIS:
        X = data[A]["X"]; h = data[A]["h"]; perms = list(itertools.permutations(range(5)))
        qbar = np.mean([q_paired(X["T"][list(p)], X["V"], h) for p in perms], axis=0)
        cells = RD.prepare(tab, H, "all", [A], drop)
        for j, cc in enumerate(cells):
            o = list(cc["T"].values())[0]; nV = len(cc["V"])
            qi = ((o.mean() - cc["V"].mean()) ** 2 - cc["sp"] ** 2 / len(o) - cc["sp"] ** 2 / nV) / cc["h"] ** 2
            rel = abs(qbar[j] - qi) / max(abs(qi), 1e-12); worst = max(worst, rel)
            rows.append({"A": A, "measure": cc["m"], "q_perm_mean": float(qbar[j]), "q_independent": float(qi), "rel_dev": float(rel)})
    return {"cells": rows, "max_rel_dev": float(worst), "pass": bool(worst <= 1e-9)}

def cp_perm(data):
    out = {}
    for A in TIS:
        X = data[A]["X"]; h = data[A]["h"]; perms = list(itertools.permutations(range(5)))
        R = np.array([((X["T"][list(p)] - X["V"]).var(0, ddof=1) / h ** 2).sum() for p in perms]); Rid = R[0]
        assert perms[0] == (0, 1, 2, 3, 4)
        rank = int((R < Rid - 1e-12).sum()) + 1; ties = int((np.abs(R - Rid) <= 1e-12).sum()) - 1
        out[A] = {"R_identity": float(Rid), "R_perm_median": float(np.median(R)), "rank": rank, "ties": ties, "n_perm": len(perms), "pass": bool(rank <= 6)}
    out["pass"] = bool(all(out[A]["pass"] for A in TIS)); return out

def synthetic(data, f, nrep, nb=NB):
    """C-1p e COV. Per tessuto: normale multivariata sui 24 valori per area (misure x V, T, B, P)."""
    base = {}
    for A in TIS:
        X = data[A]["X"]; Y = np.hstack([X[a] for a in ("V",) + SYN_ARMS]); S = np.cov(Y, rowvar=False)
        Sig = 0.95 * S + 0.05 * np.diag(np.diag(S)); base[A] = {"Sig": Sig, "m": X["V"].shape[1], "mu_obs": Y.mean(0), "muV": X["V"].mean(0)}
        assert np.linalg.eigvalsh(Sig).min() > 0
    scen = {"0": 0.0, "1": 1.0, "2": 2.0, "observed": None}; res = {}
    rng = np.random.default_rng(SYN_SEED)
    for si, (sname, s) in enumerate(scen.items()):
        mus = {}
        for A in TIS:
            b = base[A]; h = data[A]["h"]
            mus[A] = b["mu_obs"].copy() if s is None else np.concatenate([b["muV"]] + [b["muV"] + s * h for _ in SYN_ARMS])
        truth = {}
        for a_i, a in enumerate(SYN_ARMS):
            qs = []
            for A in TIS:
                m = base[A]["m"]; muV = mus[A][:m]; mua = mus[A][(a_i + 1) * m:(a_i + 2) * m]; qs.append((mua - muV) ** 2 / data[A]["h"] ** 2)
            truth[a] = float(ssqrt(np.concatenate(qs).mean()))
        qsum = {a: [] for a in SYN_ARMS}; cov = {"i": {}, "ii": {}}; keys = list(SYN_ARMS) + ["P-T", "B-T"]
        for k in keys: cov["i"][k] = []; cov["ii"][k] = []
        for r in range(nrep):
            sd = {}
            for A in TIS:
                m = base[A]["m"]; Y = rng.multivariate_normal(mus[A], base[A]["Sig"], size=5)
                sd[A] = {"X": {a: Y[:, j * m:(j + 1) * m] for j, a in enumerate(("V",) + SYN_ARMS)}, "h": data[A]["h"]}
            Dobs = {a: float(v) for a, v in D_of(sd, SYN_ARMS, TIS).items()}
            for a in SYN_ARMS:
                qq = np.concatenate([q_paired(sd[A]["X"][a], sd[A]["X"]["V"], sd[A]["h"]) for A in TIS]).mean(); qsum[a].append(qq)
            if s == 0.0: continue
            Dbs = D_of(sd, SYN_ARMS, TIS, boot_idx(np.random.default_rng([SYN_SEED, si, r]), nb))
            for a in SYN_ARMS:
                ci_i, ci_ii = intervals(Dobs[a], Dbs[a], f)
                cov["i"][a].append(ci_i[0] <= truth[a] <= ci_i[1]); cov["ii"][a].append(ci_ii[0] <= truth[a] <= ci_ii[1])
            if s is None:
                for a, b in (("P", "T"), ("B", "T")):
                    d = Dobs[a] - Dobs[b]; x = Dbs[a] - Dbs[b]; sh = np.median(x) - d; lo, hi = np.percentile(x, [2.5, 97.5]) - sh; tr = truth[a] - truth[b]
                    cov["i"][f"{a}-{b}"].append(lo <= tr <= hi); cov["ii"][f"{a}-{b}"].append(d + f * (lo - d) <= tr <= d + f * (hi - d))
        ms = {a: float(np.mean(qsum[a])) for a in SYN_ARMS}
        res[sname] = {"truth": truth, "mean_square": ms, "root": {a: float(ssqrt(v)) for a, v in ms.items()},
                      "coverage": {meth: {k: (float(np.mean(v)) if v else None) for k, v in cv.items()} for meth, cv in cov.items()}}
    c1 = abs(res["0"]["mean_square"]["T"]) <= 0.15 and abs(res["1"]["root"]["T"] - 1) <= 0.25 and abs(res["2"]["root"]["T"] - 2) <= 0.25
    covs = {meth: [v for sn in ("1", "2", "observed") for v in res[sn]["coverage"][meth].values() if v is not None] for meth in ("i", "ii")}
    mins = {meth: float(min(v)) for meth, v in covs.items()}
    method = "i" if mins["i"] >= 0.90 else "ii"
    return {"scenarios": res, "C1p_pass": bool(c1), "min_coverage": mins, "method": method, "method_flag": bool(mins[method] < 0.90), "nrep": nrep}

def per_measure(data):
    rows = []
    for A in TIS:
        X = data[A]["X"]; h = data[A]["h"]
        q = {a: q_paired(X[a], X["V"], h) for a in ("T", "B", "P")}
        for j, m in enumerate(data[A]["meas"]):
            d = X["T"][:, j] - X["V"][:, j]
            rows.append({"A": A, "measure": m, "h": float(h[j]), "V_mean": float(X["V"][:, j].mean()), "T_mean": float(X["T"][:, j].mean()), "B_mean": float(X["B"][:, j].mean()),
                         "P_mean": float(X["P"][:, j].mean()), "r_TV": float(np.corrcoef(X["T"][:, j], X["V"][:, j])[0, 1]), "r_BV": float(np.corrcoef(X["B"][:, j], X["V"][:, j])[0, 1]),
                         "T_q": float(q["T"][j]), "B_q": float(q["B"][j]), "P_q": float(q["P"][j])})
    P = pd.DataFrame(rows)
    for a in ("T", "B", "P"): P[f"{a}_d"] = ssqrt(P[f"{a}_q"])
    P["B_minus_T_q"] = P.B_q - P.T_q; P["share_of_gap"] = P.B_minus_T_q / P.B_minus_T_q.sum()
    return P.sort_values("B_minus_T_q", ascending=False).reset_index(drop=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--table", default="results/sonda/S3b/area_values.csv"); ap.add_argument("--h", default="results/sonda/S3b/h.csv")
    ap.add_argument("--out", default="results/sonda/S3c/D.json"); ap.add_argument("--per", default="results/sonda/S3c/per_measure.csv")
    ap.add_argument("--misura", default=None); ap.add_argument("--provv", default="results/sonda/S3c/misura_provvisoria.json")
    ap.add_argument("--skip-c2", action="store_true"); ap.add_argument("--nrep", type=int, default=1000); a = ap.parse_args()
    data, drop, Vs = load(a.table, a.h); n = 5; f = float(np.sqrt(n / (n - 1)) * stats.t.ppf(0.975, n - 1) / stats.norm.ppf(0.975))
    res = {"rule": "docs/sonda/REGOLA_S3c.md", "estimator": "paired", "delta": DELTA, "B_boot": NB, "seed": SEED, "syn_seed": SYN_SEED, "f": f,
           "drop_M4": drop, "V_singlets": Vs, "areas": {A: data[A]["areas"] for A in TIS}, "by": {}}
    res["GID"] = g_id(data, a.table, a.h)
    syn = synthetic(data, f, a.nrep); res["COV"] = syn; meth = syn["method"]; fac = 1.0 if meth == "i" else f
    res["C1_synthetic"] = {"values": {k: v["root"] for k, v in syn["scenarios"].items()}, "mean_square_0": syn["scenarios"]["0"]["mean_square"]["T"], "pass": syn["C1p_pass"]}
    idx = boot_idx(np.random.default_rng(SEED), NB)
    for g, ts in GROUPS.items():
        Dobs = {k: float(v) for k, v in D_of(data, ARMS, ts).items()}; Dbs = D_of(data, ARMS, ts, idx)
        ci = {k: intervals(Dobs[k], Dbs[k], f) for k in ARMS}
        comps = {m: [compare(Dobs, Dbs, x, y, 1.0 if m == "i" else f) for x, y in (("B", "T"), ("P", "T"), ("B0", "B1"), ("P", "B"))] for m in ("i", "ii")}
        res["by"][g] = {"D": Dobs, "D_ci": {k: ci[k][0 if meth == "i" else 1] for k in ARMS}, "D_ci_i": {k: ci[k][0] for k in ARMS}, "D_ci_ii": {k: ci[k][1] for k in ARMS},
                        "D_se": {k: float(fac * np.std(Dbs[k], ddof=1)) for k in ARMS}, "D_se_i": {k: float(np.std(Dbs[k], ddof=1)) for k in ARMS},
                        "comparisons": comps[meth], "comparisons_i": comps["i"], "comparisons_ii": comps["ii"], "method": meth}
    m = res["by"]["tumour+stroma"]; cmp = {(c["a"], c["b"]): c for c in m["comparisons"]}
    res["C4_positive"] = {"comparison": cmp[("P", "T")], "pass": bool(cmp[("P", "T")]["outcome"] == "different" and cmp[("P", "T")]["better"] == "T")}
    res["C5_negative"] = {"comparison": cmp[("B0", "B1")], "pass": bool(cmp[("B0", "B1")]["outcome"] in ("equal", "undetermined"))}
    if not a.skip_c2: res["C2_colon_identity"] = S3_D.colon_identity()
    res["CP_perm"] = cp_perm(data)
    pav, tet = m["D"]["B"], m["D"]["T"]
    res["floor_ceiling"] = {"pavimento": pav, "tetto": tet, "soglia": tet + 1.0, "soglia_delta_0.5": tet + 0.5, "soglia_delta_2": tet + 2.0,
                            "floor_minus_ceiling": pav - tet, "G1_precondition_gt_delta": bool(pav - tet > 1.0), "B_vs_T": cmp[("B", "T")],
                            "A1_today": pav, "A1_se": m["D_se"]["B"], "A1_holds": bool(pav + 2 * m["D_se"]["B"] <= tet + 1.0)}
    P = per_measure(data); P.to_csv(a.per, index=False)
    gates = ("GID", "C1_synthetic", "C2_colon_identity", "C4_positive", "C5_negative")
    allpass = all(res.get(k, {}).get("pass", True) for k in gates); res["gates_pass"] = allpass
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    if a.misura and not allpass:
        if os.path.exists(a.misura): os.remove(a.misura)
        a.misura = a.provv
    if a.misura:
        voci = {f"{r.A} {r.measure}": round(float(r.B_d), 4) for r in P.itertuples()}
        json.dump({"valore": round(pav, 6), "se": round(m["D_se"]["B"], 6), "voci": voci, "affermazioni": {}, "fonte": a.out, "tetto_misurato": round(tet, 6),
                   "nota": f"A1 = D(B), stimatore appaiato (S3c), intervallo {meth}", "controlli_passati": allpass}, open(a.misura, "w"), indent=1)
    print("D", {g: {k: round(v, 3) for k, v in r["D"].items()} for g, r in res["by"].items()})
    print("GID", res["GID"]["pass"], f"{res['GID']['max_rel_dev']:.1e}", "| C1p", res["C1_synthetic"]["pass"], {k: round(v["T"], 3) for k, v in res["C1_synthetic"]["values"].items()},
          "ms0", round(res["C1_synthetic"]["mean_square_0"], 4))
    print("COV min", {k: round(v, 3) for k, v in syn["min_coverage"].items()}, "method", meth, "flag", syn["method_flag"], "f", round(f, 4))
    for sn, sv in syn["scenarios"].items():
        if sv["coverage"]["i"]["T"] is not None: print("  cov", sn, {mm: {k: (round(v, 3) if v is not None else None) for k, v in cc.items()} for mm, cc in sv["coverage"].items()})
    print("C4", res["C4_positive"]["pass"], "C5", res["C5_negative"]["pass"], "C2", res.get("C2_colon_identity", {}).get("pass"), "| CP", {A: (res["CP_perm"][A]["rank"], res["CP_perm"][A]["ties"]) for A in TIS}, res["CP_perm"]["pass"])
    for c in m["comparisons"]: print("  ", c["a"], "-", c["b"], round(c["delta_ab"], 3), "w", round(c["w"], 3), c["outcome"], c["better"])
    print("floor", round(pav, 3), "ci", [round(x, 3) for x in m["D_ci"]["B"]], "se", round(m["D_se"]["B"], 3), "| ceiling", round(tet, 3), "ci", [round(x, 3) for x in m["D_ci"]["T"]],
          "| threshold", round(tet + 1.0, 3), "| G1 precondition", res["floor_ceiling"]["G1_precondition_gt_delta"], "| gates", allpass)
