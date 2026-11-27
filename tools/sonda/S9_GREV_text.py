#!/usr/bin/env python
"""S9 G-REV, text-only re-derivation of the D part (A4 U/C, A2 L).

Written to address a contamination: the first G-REV script was written after a grep of
tools/sonda/S3c_D.py had shown its q / interval / w lines. Here every formula is taken from the
TEXT of REGOLA_S3c.md / REGOLA_S9.md / FATTI_B.md and coded in a different form (explicit loops,
explicit sums), and every quantity whose text admits more than one reading is computed under ALL
readings, so it is visible which declared values/outcomes rest on a reading not fixed by the text.

Text -> code:
 q   = (mean(d)^2 - var(d)/n)/h^2, var ddof=1, n=5                 [S3c, explicit]
 D   = signed root of mean of q over cells                          [S3c, explicit]
 se  = bootstrap standard deviation x f                             [S3c; ddof NOT stated -> ddof0, ddof1]
 w   = 95% half-width of the bootstrap distribution of D(a)-D(b), x f [S3c; form NOT stated ->
       'pct' = (p97.5 - p2.5)/2 ; 'norm' = 1.96 * sd(ddof1)]
 comparison: equal |d|+w<=1, different |d|-w>1                      [FATTI_B.md line 24]
 e = D(B)-(D(T)+1), e_se = f x sd(D_B* - D_T*)                     [S9 A4]
Bootstrap indices: default_rng(20261123).integers(0,5,(2000,5)), tumour then stroma (task-permitted).
"""
import json, math, os
import numpy as np, pandas as pd
from scipy import stats

REPO = os.path.expanduser("~/2025.geo_spatialtrans")
OUTD = os.path.join(REPO, "results/sonda/S9/GREV")
TIS = ["tumour", "stroma"]; MEAS = ["M1", "M2", "M3", "M4", "M5", "M6"]
NB, n = 2000, 5
f = math.sqrt(n / (n - 1)) * stats.t.ppf(0.975, n - 1) / stats.norm.ppf(0.975)


def read(tab, hf):
    t = pd.read_csv(os.path.join(REPO, tab)); h = pd.read_csv(os.path.join(REPO, hf))
    V, A_, H = {}, {}, {}
    for A in TIS:
        for _, r in t[t.A == A].iterrows():
            k = int(r.unit.split("_")[2]) - 1          # H?_<tissue>_<k>[_k<s>]
            arm = r.arm if r.arm != "B" else "B" + str(r.k)
            A_.setdefault((A, arm), [None] * 5)[k] = [float(r[m]) for m in MEAS]
        for j, m in enumerate(MEAS):
            H[(A, j)] = float(h[(h.A == A) & (h.measure == "M-" + m[1])].h.iloc[0])
    arms = {}
    for (A, arm), rows in A_.items(): arms[(A, arm)] = rows
    for A in TIS:
        arms[(A, "B")] = [[(x + y) / 2 for x, y in zip(a, b)] for a, b in zip(arms[(A, "B0")], arms[(A, "B1")])]
    return arms, H


def q_one(a_rows, v_rows, j, h, ks):
    d = [a_rows[k][j] - v_rows[k][j] for k in ks]
    m = sum(d) / n
    s2 = sum((x - m) ** 2 for x in d) / (n - 1)
    return (m * m - s2 / n) / (h * h)


def D(arms, H, arm, tissues, ks_by_t, keep):
    qs = [q_one(arms[(A, arm)], arms[(A, "V")], j, H[(A, j)], ks_by_t[A])
          for A in tissues for j in range(6) if keep[(A, j)]]
    x = sum(qs) / len(qs)
    return math.copysign(math.sqrt(abs(x)), x)


def sd(x, ddof): x = np.asarray(x); return float(np.sqrt(((x - x.mean()) ** 2).sum() / (len(x) - ddof)))


def run_chain(tab, hf, label, rows_out):
    arms, H = read(tab, hf)
    keep = {(A, j): not (MEAS[j] == "M4" and np.mean([r[3] for r in arms[(A, "V")]]) >= 0.99) for A in TIS for j in range(6)}
    rng = np.random.default_rng(20261123); idx = {A: rng.integers(0, 5, (NB, 5)) for A in TIS}
    ident = {A: list(range(5)) for A in TIS}
    out = {}
    for name, tissues in [("tumour+stroma", TIS), ("tumour", ["tumour"]), ("stroma", ["stroma"])]:
        Do = {a: D(arms, H, a, tissues, ident, keep) for a in ["T", "B", "P", "B0", "B1"]}
        Db = {a: [D(arms, H, a, tissues, {A: list(idx[A][b]) for A in TIS}, keep) for b in range(NB)]
              for a in Do}
        rec = {"D": Do, "se_ddof1": {a: f * sd(Db[a], 1) for a in Do}, "se_ddof0": {a: f * sd(Db[a], 0) for a in Do}, "comp": {}}
        for a, b in [("B", "T"), ("P", "T"), ("B0", "B1"), ("P", "B")]:
            dd = Do[a] - Do[b]; diff = np.array(Db[a]) - np.array(Db[b])
            ws = {"pct": f * (np.percentile(diff, 97.5) - np.percentile(diff, 2.5)) / 2, "norm": f * 1.96 * sd(diff, 1)}
            oc = {}
            for wn, w in ws.items():
                oc[wn] = "different" if abs(dd) - w > 1 else ("equal" if abs(dd) + w <= 1 else "undetermined")
            rec["comp"][f"{a}-{b}"] = {"delta": dd, "w": ws, "outcome": oc}
        if name == "tumour+stroma":
            e = Do["B"] - (Do["T"] + 1); diff = np.array(Db["B"]) - np.array(Db["T"])
            rec["e"] = e; rec["e_se"] = {"ddof1": f * sd(diff, 1), "ddof0": f * sd(diff, 0)}
            rec["esito"] = {k: ("regge" if e + 2 * s <= 0 else "non regge" if e - 2 * s > 0 else "incerta") for k, s in rec["e_se"].items()}
        out[name] = rec
    dec = json.load(open(os.path.join(REPO, f"results/sonda/S9/A4_{label}.json")))
    for name in out:
        r, d = out[name], dec[name]
        for a in r["D"]:
            rows_out.append((f"{label}.{name}.D.{a}", "point", "text", d["D"][a], r["D"][a]))
            for dd in (1, 0):
                rows_out.append((f"{label}.{name}.D_se.{a}", "se", f"ddof{dd}", d["D_se"][a], r[f"se_ddof{dd}"][a]))
        for c, v in r["comp"].items():
            rows_out.append((f"{label}.{name}.comp.{c}.delta_ab", "point", "text", d["comp"][c]["delta_ab"], v["delta"]))
            for wn in ("pct", "norm"):
                rows_out.append((f"{label}.{name}.comp.{c}.w", "se", wn, d["comp"][c]["w"], v["w"][wn]))
                rows_out.append((f"{label}.{name}.comp.{c}.outcome", "label", wn, d["comp"][c]["outcome"], v["outcome"][wn]))
        if "e" in r:
            rows_out.append((f"{label}.{name}.e", "point", "text", d["e"], r["e"]))
            for k in ("ddof1", "ddof0"):
                rows_out.append((f"{label}.{name}.e_se", "se", k, d["e_se"], r["e_se"][k]))
                rows_out.append((f"{label}.{name}.esito", "label", k, d["esito"], r["esito"][k]))


def run_A2(rows_out):
    arms, H = read("results/sonda/S3b/area_values.csv", "results/sonda/S3b/h.csv")
    rng = np.random.default_rng(20261123); idx = {A: rng.integers(0, 5, (NB, 5)) for A in TIS}
    rs = lambda x: math.copysign(math.sqrt(abs(x)), x)
    voci = []
    for A in TIS:
        for j in range(6):
            Td = rs(q_one(arms[(A, "T")], arms[(A, "V")], j, H[(A, j)], range(5)))
            Bd = rs(q_one(arms[(A, "B")], arms[(A, "V")], j, H[(A, j)], range(5)))
            Bb = [rs(q_one(arms[(A, "B")], arms[(A, "V")], j, H[(A, j)], list(idx[A][b]))) for b in range(NB)]
            voci.append((A, "M-" + str(j + 1), Td, Bd, np.array(Bb)))
    Tmax = max(v[2] for v in voci); cnt = sum(v[3] > Tmax for v in voci)
    cb = np.sum([v[4] > Tmax for v in voci], axis=0)
    dec = json.load(open(os.path.join(REPO, "results/sonda/S9/A2_L.json")))
    decv = pd.read_csv(os.path.join(REPO, "results/sonda/S9/A2_per_voce_L.csv"))
    rows_out.append(("A2.T_max", "point", "text", dec["T_max"], Tmax))
    rows_out.append(("A2.A2", "point", "text", dec["A2"], int(cnt)))
    for dd in (1, 0):
        s = f * sd(cb, dd)
        rows_out.append(("A2.A2_se", "se", f"ddof{dd}", dec["A2_se"], s))
        rows_out.append(("A2.esito", "label", f"ddof{dd}", dec["esito"], "regge" if cnt + 2 * s <= 0 else "non regge" if cnt - 2 * s > 0 else "incerta"))
    for A, m, Td, Bd, Bb in voci:
        r = decv[(decv.A == A) & (decv.measure == m)].iloc[0]
        rows_out.append((f"A2v.{A}.{m}.T_d", "point", "text", r.T_d, Td))
        rows_out.append((f"A2v.{A}.{m}.B_d", "point", "text", r.B_d, Bd))
        for dd in (1, 0):
            s = f * sd(Bb, dd)
            cl = "oltre" if Bd - 2 * s > Tmax else ("dentro" if Bd + 2 * s <= Tmax else "indeterminata")
            rows_out.append((f"A2v.{A}.{m}.se", "se", f"ddof{dd}", r.se, s))
            rows_out.append((f"A2v.{A}.{m}.class", "label", f"ddof{dd}", r["class"], cl))


def main():
    rows = []
    run_chain("results/sonda/S9/U/area_values_U.csv", "results/sonda/S9/U/h_U.csv", "U", rows)
    run_chain("results/sonda/S9/C/area_values_C.csv", "results/sonda/S9/C/h_C.csv", "C", rows)
    run_A2(rows)
    df = pd.DataFrame(rows, columns=["quantity", "kind", "reading", "declared", "recomputed"])
    num = df.kind != "label"
    df["abs_diff"] = np.where(num, (pd.to_numeric(df.declared, errors="coerce") - pd.to_numeric(df.recomputed, errors="coerce")).abs(), np.nan)
    df["match"] = np.where(num, df.abs_diff < 1e-9, df.declared.astype(str) == df.recomputed.astype(str))
    df.to_csv(os.path.join(OUTD, "readings_text_only.csv"), index=False)
    summ = {}
    for (k, rd), g in df.groupby(["kind", "reading"]):
        summ[f"{k}|{rd}"] = {"n": int(len(g)), "n_match": int(g.match.sum()),
                             "max_abs_diff": (None if k == "label" else float(g.abs_diff.max()))}
    # outcomes that change between the two w readings / ddof readings (robustness of decisions)
    lab = df[df.kind == "label"]
    flips = []
    for q, g in lab.groupby("quantity"):
        if g.recomputed.nunique() > 1: flips.append({q: dict(zip(g.reading, g.recomputed))})
    summ["label_changes_across_readings"] = flips
    json.dump(summ, open(os.path.join(OUTD, "readings_text_only_summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
