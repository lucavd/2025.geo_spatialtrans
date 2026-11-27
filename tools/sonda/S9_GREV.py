#!/usr/bin/env python
"""S9 gate G-REV: independent recomputation of A4 (chains U, C), A2 (chain L) and A5-gain.

Written from the specifications only (docs/sonda/REGOLA_S9.md sections A4, A2, A5 col guadagno;
docs/sonda/REGOLA_S3c.md "Stimatore appaiato", "Incertezza", method (ii);
docs/sonda/REGOLA_S8.md splits/bootstrap; comparison rule of docs/riesame/FATTI_B.md:
equal if |D|+w <= 1, different if |D|-w > 1, else undetermined).
Does not import tools/sonda/S9_A.py, S3c_D.py or S8_A5.py.
Bootstrap indices for D/A2: default_rng(20261123), integers(0,5,(2000,5)) for tumour, then stroma.
A5 bootstrap: own reading of REGOLA_S8 (draw order not specified in text) -> se compared separately.
"""
import hashlib, json, os, sys, itertools
import numpy as np, pandas as pd
from scipy import stats

REPO = os.path.expanduser("~/2025.geo_spatialtrans")
OUTD = os.path.join(REPO, "results/sonda/S9/GREV")
REVIEW = os.path.join(REPO, "results/sonda/S9/review_GREV.csv")
X1ROOT = "/mnt/micron/geo_spatialtrans/S8_sonda/X1"
TIS = ["tumour", "stroma"]
MEAS = ["M1", "M2", "M3", "M4", "M5", "M6"]
NB, SEED_D, SEED_A5, DELTA = 2000, 20261123, 20261201, 1.0
N = 5
F = np.sqrt(N / (N - 1)) * stats.t.ppf(0.975, N - 1) / stats.norm.ppf(0.975)

INPUTS = {
    "U_tab": "results/sonda/S9/U/area_values_U.csv", "U_h": "results/sonda/S9/U/h_U.csv",
    "C_tab": "results/sonda/S9/C/area_values_C.csv", "C_h": "results/sonda/S9/C/h_C.csv",
    "L_tab": "results/sonda/S3b/area_values.csv", "L_h": "results/sonda/S3b/h.csv",
    "S8_pu": "results/sonda/S8/per_unit.csv", "S8_A5": "results/sonda/S8/A5.json",
    "dec_A4U": "results/sonda/S9/A4_U.json", "dec_A4C": "results/sonda/S9/A4_C.json",
    "dec_A2": "results/sonda/S9/A2_L.json", "dec_A2v": "results/sonda/S9/A2_per_voce_L.csv",
    "dec_A5": "results/sonda/S9/A5_gain.json", "dec_A5b": "results/sonda/S9/A5_gain_boot.csv",
}


def P(rel): return rel if rel.startswith("/") else os.path.join(REPO, rel)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for ch in iter(lambda: fh.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()


def ssqrt(x): return np.sign(x) * np.sqrt(np.abs(x))


# ---------------------------------------------------------------- data for D
def load_chain(tab, hfile):
    """Return {tissue: {arm: array (5 areas, 6 measures)}}, h {tissue: (6,)}, V M-4 mean."""
    t = pd.read_csv(P(tab)); hh = pd.read_csv(P(hfile))
    out, H, vm4 = {}, {}, {}
    for A in TIS:
        s = t[t.A == A].copy()
        s["kk"] = s.unit.str.extract(r"_(\d+)$")[0].astype(int)  # area index k in unit name
        arms = {}
        def mat(df):
            df = df.sort_values("kk"); assert list(df.kk) == [1, 2, 3, 4, 5], df.kk.tolist()
            return df[MEAS].to_numpy(float)
        arms["V"] = mat(s[s.arm == "V"])
        arms["T"] = mat(s[s.arm == "T"])
        arms["P"] = mat(s[s.arm == "P"])
        b = s[s.arm == "B"]
        arms["B0"] = mat(b[b.k.astype(str) == "0"])
        arms["B1"] = mat(b[b.k.astype(str) == "1"])
        arms["B"] = (arms["B0"] + arms["B1"]) / 2.0  # mean of seeds per area first
        out[A] = arms
        H[A] = np.array([float(hh[(hh.A == A) & (hh.measure == "M-" + m[1])].h.iloc[0]) for m in MEAS])
        vm4[A] = float(arms["V"][:, 3].mean())
    return out, H, vm4


def q_cells(Xa, Xv, h, idx=None):
    """Paired q per measure. Xa, Xv (5,6). idx None -> (6,); idx (nb,5) -> (nb,6)."""
    d = Xa - Xv
    if idx is not None: d = d[idx]  # (nb,5,6)
    n = d.shape[-2]
    return (d.mean(-2) ** 2 - d.var(-2, ddof=1) / n) / h ** 2


def boot_indices():
    rng = np.random.default_rng(SEED_D)
    return {A: rng.integers(0, N, (NB, N)) for A in TIS}


def keep_mask(vm4):
    """measures kept per tissue (M-4 dropped only if V mean >= 0.99)."""
    return {A: np.array([not (m == "M4" and vm4[A] >= 0.99) for m in MEAS]) for A in TIS}


def D_arm(data, H, keep, arm, tissues, idx=None):
    qs = []
    for A in tissues:
        q = q_cells(data[A][arm], data[A]["V"], H[A], None if idx is None else idx[A])
        qs.append(q[..., keep[A]])
    return ssqrt(np.concatenate(qs, axis=-1).mean(-1))


def interval(dobs, dbs):
    sh = np.median(dbs) - dobs
    lo, hi = np.percentile(dbs, [2.5, 97.5]) - sh
    return [float(dobs - F * (dobs - lo)), float(dobs + F * (hi - dobs))]


def compare(da, db, xa, xb):
    d = da - db
    lo, hi = np.percentile(xa - xb, [2.5, 97.5])
    w = F * (hi - lo) / 2.0
    if abs(d) - w > DELTA: out = "different"
    elif abs(d) + w <= DELTA: out = "equal"
    else: out = "undetermined"
    return d, w, out


ARMS = ["T", "B", "P", "B0", "B1"]
PAIRS = [("B", "T"), ("P", "T"), ("B0", "B1"), ("P", "B")]


def a4_chain(tab, hfile):
    data, H, vm4 = load_chain(tab, hfile)
    keep = keep_mask(vm4); idx = boot_indices()
    res = {"M4_dropped": {A: bool(not keep[A][3]) for A in TIS}, "V_M4_mean": vm4}
    for name, tissues in [("tumour+stroma", TIS), ("tumour", ["tumour"]), ("stroma", ["stroma"])]:
        Dobs = {a: float(D_arm(data, H, keep, a, tissues)) for a in ARMS}
        Dbs = {a: D_arm(data, H, keep, a, tissues, idx) for a in ARMS}
        r = {"D": Dobs, "D_se": {a: float(F * np.std(Dbs[a], ddof=1)) for a in ARMS},
             "D_se_ddof0": {a: float(F * np.std(Dbs[a], ddof=0)) for a in ARMS},
             "ci": {a: interval(Dobs[a], Dbs[a]) for a in ARMS}, "comp": {}}
        for a, b in PAIRS:
            d, w, out = compare(Dobs[a], Dobs[b], Dbs[a], Dbs[b])
            better = "" if out != "different" else (a if Dobs[a] < Dobs[b] else b)
            r["comp"][f"{a}-{b}"] = {"delta_ab": float(d), "w": float(w), "outcome": out, "better": better}
        if name == "tumour+stroma":
            e = Dobs["B"] - (Dobs["T"] + DELTA)
            ese = float(F * np.std(Dbs["B"] - Dbs["T"], ddof=1))
            esito = "regge" if e + 2 * ese <= 0 else ("non regge" if e - 2 * ese > 0 else "incerta")
            r.update({"e": float(e), "e_se": ese, "esito": esito, "tetto": Dobs["T"], "pavimento": Dobs["B"],
                      "soglia": Dobs["T"] + DELTA,
                      "CP_chain_pass": r["comp"]["B-T"]["outcome"] == "different",
                      "C4_pass": r["comp"]["P-T"]["outcome"] == "different" and r["comp"]["P-T"]["better"] == "T",
                      "C5_pass": r["comp"]["B0-B1"]["outcome"] != "different"})
        res[name] = r
    return res


# ---------------------------------------------------------------- A2
def a2_L():
    data, H, vm4 = load_chain(INPUTS["L_tab"], INPUTS["L_h"])
    keep = keep_mask(vm4); idx = boot_indices()
    rows, Bboot, Tvals = [], [], []
    for A in TIS:
        qT = q_cells(data[A]["T"], data[A]["V"], H[A]); qB = q_cells(data[A]["B"], data[A]["V"], H[A])
        qBb = q_cells(data[A]["B"], data[A]["V"], H[A], idx[A])
        for j, m in enumerate(MEAS):
            if not keep[A][j]: continue
            rows.append({"A": A, "measure": "M-" + m[1], "T_d": float(ssqrt(qT[j])), "B_d": float(ssqrt(qB[j]))})
            Bboot.append(ssqrt(qBb[:, j]))
    v = pd.DataFrame(rows); Bboot = np.stack(Bboot, 1)  # (nb, nvoci)
    Tmax = float(v.T_d.max()); voce = v.loc[v.T_d.idxmax()]
    v["se"] = F * np.std(Bboot, axis=0, ddof=1)
    v["class"] = np.where(v.B_d - 2 * v.se > Tmax, "oltre", np.where(v.B_d + 2 * v.se <= Tmax, "dentro", "indeterminata"))
    cnt = int((v.B_d > Tmax).sum()); cnt_b = (Bboot > Tmax).sum(1)
    se = float(F * np.std(cnt_b, ddof=1))
    esito = "regge" if cnt + 2 * se <= 0 else ("non regge" if cnt - 2 * se > 0 else "incerta")
    return {"T_max": Tmax, "T_max_voce": f"{voce.A} {voce.measure}", "A2": cnt, "A2_se": se,
            "A2_se_ddof0": float(F * np.std(cnt_b, ddof=0)), "esito": esito,
            "n_oltre": int((v["class"] == "oltre").sum()), "n_indeterminata": int((v["class"] == "indeterminata").sum()),
            "n_dentro": int((v["class"] == "dentro").sum())}, v


# ---------------------------------------------------------------- A5 gain
def taub(x, y): return float(stats.kendalltau(x, y, variant="b").statistic)


def a5_gain():
    methods = json.load(open(P(INPUTS["S8_A5"])))["methods"]
    pu = pd.read_csv(P(INPUTS["S8_pu"])); pu["gain"] = pu.score - pu.score_shuf
    def area_mat(arm):
        g = pu[pu.arm == arm].groupby(["i", "method"]).gain.mean().unstack("method")[methods]
        assert list(g.index) == list(range(10)); return g.to_numpy(float)
    GB, GX = area_mat("B"), area_mat("X1")  # (10 areas, 5 methods); B = mean of 2 seeds
    tis = pu[pu.arm == "X1"].groupby("i").tissue.first()
    tum = [i for i in range(10) if tis[i] == "tumour"]; stro = [i for i in range(10) if tis[i] == "stroma"]
    divs = [(list(a) + list(b), sorted(set(tum) - set(a)) + sorted(set(stro) - set(b)))
            for a in itertools.combinations(tum, 3) for b in itertools.combinations(stro, 2)]
    assert len(divs) == 100
    ordered = divs + [(h2, h1) for h1, h2 in divs]  # 200 ordered pairs
    def score(G, half): return G[half].mean(0)
    T5 = float(np.mean([taub(score(GX, h1), score(GX, h2)) for h1, h2 in ordered]))
    S5c = float(np.mean([taub(score(GB, h1), score(GX, h2)) for h1, h2 in ordered]))
    A5 = S5c - T5
    # bootstrap (own reading): one random division, resample areas with replacement within each half and tissue
    rng = np.random.default_rng(SEED_A5); bt = []
    tset = set(tum)
    for _ in range(NB):
        h1, h2 = divs[rng.integers(0, len(divs))]
        def rs(h):
            out = []
            for grp in ([i for i in h if i in tset], [i for i in h if i not in tset]):
                out += list(rng.choice(grp, len(grp), replace=True))
            return out
        r1, r2 = rs(h1), rs(h2)
        t5 = taub(score(GX, r1), score(GX, r2))
        s5 = (taub(score(GB, r1), score(GX, r2)) + taub(score(GB, r2), score(GX, r1))) / 2.0
        bt.append((t5, s5, s5 - t5))
    bt = np.array(bt)
    se = {k: float(np.std(bt[:, j], ddof=1)) for j, k in enumerate(["T5", "S5c", "A5"])}
    if T5 - 2 * se["T5"] <= 0 <= T5 + 2 * se["T5"]: esito = "non misurabile"
    elif A5 - 2 * se["A5"] >= -0.2: esito = "regge"
    elif A5 + 2 * se["A5"] < -0.2: esito = "non regge"
    else: esito = "incerta"
    # CP-gain: tau-b between gain ranking and panel UMI per cell (mean over kept cells, then over 10 areas)
    umi = []
    for u in sorted(pu[pu.arm == "X1"].unit.unique()):
        c = pd.read_parquet(os.path.join(X1ROOT, u, "cells.parquet"))
        n_kept = int(pu[(pu.arm == "X1") & (pu.unit == u)].n_kept.iloc[0]); assert len(c) == n_kept, (u, len(c), n_kept)
        umi.append([c["panel_" + m].mean() for m in methods])
    umi = np.array(umi).mean(0); gx = GX.mean(0); gb = GB.mean(0)
    cp = taub(gx, umi)
    return {"T5": T5, "S5c": S5c, "A5": A5, "T5_se": se["T5"], "S5c_se": se["S5c"], "A5_se": se["A5"], "esito": esito,
            "methods": methods, "gain_B": dict(zip(methods, map(float, gb))), "gain_X1": dict(zip(methods, map(float, gx))),
            "panel_umi_X1": dict(zip(methods, map(float, umi))), "CP_gain_tau": cp, "CP_gain_pass": bool(cp < 0.8)}, bt


# ---------------------------------------------------------------- compare
def flat(prefix, d, out):
    if isinstance(d, dict):
        for k, v in d.items(): flat(f"{prefix}.{k}" if prefix else str(k), v, out)
    elif isinstance(d, (list, tuple)) and all(isinstance(x, (int, float)) for x in d):
        for j, v in enumerate(d): out[f"{prefix}[{j}]"] = v
    else: out[prefix] = d


def main():
    os.makedirs(OUTD, exist_ok=True)
    x1files = [os.path.join(X1ROOT, f"X1_{A}_{k}", "cells.parquet") for A in TIS for k in range(1, 6)]
    files = [P(v) for v in INPUTS.values()] + x1files
    before = {f: md5(f) for f in files}

    res = {"f": float(F), "U": a4_chain(INPUTS["U_tab"], INPUTS["U_h"]), "C": a4_chain(INPUTS["C_tab"], INPUTS["C_h"])}
    a2, a2v = a2_L(); res["A2"] = a2
    a5, bt = a5_gain(); res["A5"] = a5
    a2v.to_csv(os.path.join(OUTD, "A2_per_voce_GREV.csv"), index=False)
    pd.DataFrame(bt, columns=["T5", "S5c", "A5"]).to_csv(os.path.join(OUTD, "A5_gain_boot_GREV.csv"), index=False)
    json.dump(res, open(os.path.join(OUTD, "GREV_recomputed.json"), "w"), indent=1)

    # declared
    dec = {"U": json.load(open(P(INPUTS["dec_A4U"]))), "C": json.load(open(P(INPUTS["dec_A4C"]))),
           "A2": json.load(open(P(INPUTS["dec_A2"]))), "A5": json.load(open(P(INPUTS["dec_A5"])))}
    rows = []
    for key in ["U", "C", "A2", "A5"]:
        fr, fd = {}, {}; flat(key, res[key], fr); flat(key, dec[key], fd)
        for q, dv in fd.items():
            if q not in fr: continue
            rv = fr[q]
            kind = ("se" if ("_se" in q or q.endswith(".w") or ".ci." in q or "ci[" in q) else "point")
            if key == "A5" and q.split(".")[-1] in ("T5_se", "S5c_se", "A5_se"): kind = "se_A5boot"
            if isinstance(dv, bool) or isinstance(dv, str) or isinstance(rv, (bool, str)):
                rows.append({"quantity": q, "kind": "label", "declared": dv, "recomputed": rv, "abs_diff": 0.0 if dv == rv else np.nan, "match": dv == rv})
            elif isinstance(dv, (int, float)):
                ad = abs(float(dv) - float(rv))
                rows.append({"quantity": q, "kind": kind, "declared": dv, "recomputed": rv, "abs_diff": ad, "match": ad < 1e-9})
    dv = pd.read_csv(P(INPUTS["dec_A2v"])).merge(a2v, on=["A", "measure"], suffixes=("_dec", "_rec"))
    for _, r in dv.iterrows():
        for c in ["T_d", "B_d", "se"]:
            ad = abs(r[c + "_dec"] - r[c + "_rec"])
            rows.append({"quantity": f"A2v.{r.A}.{r.measure}.{c}", "kind": "se" if c == "se" else "point",
                         "declared": r[c + "_dec"], "recomputed": r[c + "_rec"], "abs_diff": ad, "match": ad < 1e-9})
        rows.append({"quantity": f"A2v.{r.A}.{r.measure}.class", "kind": "label", "declared": r.class_dec,
                     "recomputed": r.class_rec, "abs_diff": 0.0 if r.class_dec == r.class_rec else np.nan, "match": r.class_dec == r.class_rec})
    assert len(dv) == 12
    rev = pd.DataFrame(rows); rev.to_csv(REVIEW, index=False)

    after = {f: md5(f) for f in files}
    pd.DataFrame({"file": files, "md5_before": [before[f] for f in files], "md5_after": [after[f] for f in files]}).to_csv(
        os.path.join(OUTD, "md5_inputs.csv"), index=False)
    summ = {"n_rows": len(rev),
            "max_abs_diff_point": float(rev[rev.kind == "point"].abs_diff.max()),
            "max_abs_diff_se": float(rev[rev.kind == "se"].abs_diff.max()),
            "max_abs_diff_se_A5boot": float(rev[rev.kind == "se_A5boot"].abs_diff.max()),
            "labels_mismatch": rev[(rev.kind == "label") & (~rev.match.astype(bool))].quantity.tolist(),
            "numeric_mismatch": rev[(rev.kind != "label") & (rev.kind != "se_A5boot") & (~rev.match.astype(bool))].quantity.tolist(),
            "md5_identical": all(before[f] == after[f] for f in files), "n_files_md5": len(files)}
    json.dump(summ, open(os.path.join(OUTD, "GREV_summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
