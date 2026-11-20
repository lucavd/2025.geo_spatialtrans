"""Riesame, Parte B: D per braccio, confronti, cancelli G-1, G-2, G-3, G-6 ed esiti, come in docs/riesame/REGOLA_B.md (45fdbbe).

Ingresso: una tabella lunga con colonne arm, genes, A, unit, run, set, M1..M6 (una riga per area o per area x seme).
  arm in {V, T, S, B, X, X1, XC, XC2}; per T la colonna set distingue P1, P5, 6p5mm; genes in {panel, common}.
  margini h: tabella con A, measure (M-1..M-6), h, genes.
D(a) = radice con segno della media sulle celle tessuto x misura di (Delta^2 - Var Delta) / h^2, Delta = media del braccio - media di V
(stimatore quadratico corretto di DEC-098). Bracci con piu' semi: prima la media dei semi per area.
Var Delta = s_V^2 / n_V + s_a^2 / n_a; s_V = deviazione standard combinata dei bracci reali (V e insiemi di T) come in F3.3;
s_a = deviazione standard fra aree del braccio (per T la combinata). Bootstrap parametrico congiunto (V ridisegnata una volta per replica,
cosi' le differenze fra bracci tengono conto della V comune): t con i gradi di liberta' della combinata per i reali, t con n_a - 1 per gli altri.
B = 2000, seme 20261123. w = semi-ampiezza al 95 % della distribuzione bootstrap di D(a) - D(b).
Esito: uguale se |Delta_ab| + w <= delta; diverso se |Delta_ab| - w > delta; altrimenti indeterminato. delta = 1.0.
"""
import argparse, json
from itertools import combinations
import numpy as np, pandas as pd
MEAS = ["M1", "M2", "M3", "M4", "M5", "M6"]; MID = {m: f"M-{m[1]}" for m in MEAS}
DELTA, B, SEED = 1.0, 2000, 20261123

def ssqrt(x): return float(np.sign(x) * np.sqrt(abs(x)))

def pooled_sd(groups):
    num = sum(((g - g.mean()) ** 2).sum() for g in groups if len(g) >= 2); df = sum(len(g) - 1 for g in groups if len(g) >= 2)
    return float(np.sqrt(num / df)), int(df)

def prepare(tab, H, genes, tissues, drop_m4):
    """per tessuto e misura: aree di V, insiemi di T, aree degli altri bracci (media dei semi), h, s combinata."""
    t = tab[tab.genes == genes]; cells = []
    for A in tissues:
        for m in MEAS:
            if drop_m4.get(A) and m == "M4": continue
            h = float(H[(H.A == A) & (H.measure == MID[m]) & (H.genes == genes)].h.iloc[0])
            V = t[(t.arm == "V") & (t.A == A)].groupby("unit")[m].mean().to_numpy(float)
            T = {s: g.groupby("unit")[m].mean().to_numpy(float) for s, g in t[(t.arm == "T") & (t.A == A)].groupby("set")}
            sp, df = pooled_sd([V] + list(T.values()))
            arms = {a: g.groupby("unit")[m].mean().to_numpy(float) for a, g in t[(~t.arm.isin(["V", "T"])) & (t.A == A)].groupby("arm")}
            cells.append({"A": A, "m": MID[m], "h": h, "V": V, "T": T, "arms": arms, "sp": sp, "df": df})
    return cells

def q_arm(cells, arm, tsets=None):
    v = []
    for c in cells:
        nV = len(c["V"])
        if arm == "T":
            for s, o in c["T"].items():
                if tsets and not any(k in s for k in tsets): continue
                v.append(((o.mean() - c["V"].mean()) ** 2 - c["sp"] ** 2 / len(o) - c["sp"] ** 2 / nV) / c["h"] ** 2)
        else:
            if arm not in c["arms"]: return np.nan
            a = c["arms"][arm]; sa2 = a.var(ddof=1) if len(a) > 1 else 0.0
            v.append(((a.mean() - c["V"].mean()) ** 2 - c["sp"] ** 2 / nV - sa2 / len(a)) / c["h"] ** 2)
    return float(np.mean(v)) if v else np.nan

def draw(cells, rng):
    out = []
    for c in cells:
        s, df = c["sp"], c["df"]
        V = np.full(len(c["V"]), c["V"].mean() + s / np.sqrt(len(c["V"])) * rng.standard_t(df))
        T = {k: np.full(len(o), o.mean() + s / np.sqrt(len(o)) * rng.standard_t(df)) for k, o in c["T"].items()}
        arms = {}
        for k, a in c["arms"].items():
            sa = a.std(ddof=1) if len(a) > 1 else 0.0; n = len(a)
            arms[k] = a - a.mean() + a.mean() + sa / np.sqrt(n) * (rng.standard_t(n - 1) if n > 1 else 0.0)
        out.append({**c, "V": V, "T": T, "arms": arms})
    return out
# nota: nel ridisegno V e T diventano costanti (varianza fra aree nulla), quindi Var Delta del ridisegno usa ancora sp (fisso) e,
# per i bracci, la loro varianza originale (le aree si spostano tutte insieme): la correzione resta quella osservata, come in F3.3.

def D_all(cells, arms, tsets):
    r = {a: q_arm(cells, a) for a in arms}; r["T"] = q_arm(cells, "T", tsets)
    for s in ("P1", "P5", "6p5"): r[f"T_{s}"] = q_arm(cells, "T", [s])
    return r

def analyse(cells, arms, tsets=("P1", "P5")):
    obs = D_all(cells, arms, tsets); rng = np.random.default_rng(SEED)
    bs = [D_all(draw(cells, rng), arms, tsets) for _ in range(B)]
    keys = list(obs)
    Dobs = {k: ssqrt(obs[k]) if not np.isnan(obs[k]) else np.nan for k in keys}
    Dbs = {k: np.array([ssqrt(b[k]) if not np.isnan(b[k]) else np.nan for b in bs]) for k in keys}
    return Dobs, Dbs

def compare(Dobs, Dbs, a, b):
    if np.isnan(Dobs.get(a, np.nan)) or np.isnan(Dobs.get(b, np.nan)): return {"a": a, "b": b, "delta_ab": np.nan, "w": np.nan, "outcome": "not available"}
    d = Dobs[a] - Dobs[b]; x = Dbs[a] - Dbs[b]; lo, hi = np.nanpercentile(x, [2.5, 97.5]); w = (hi - lo) / 2
    out = "equal" if abs(d) + w <= DELTA else ("different" if abs(d) - w > DELTA else "undetermined")
    return {"a": a, "b": b, "delta_ab": d, "w": w, "outcome": out, "better": (a if d < 0 else b) if out == "different" else ""}

def split_check(cells):
    """G-3: le 5 aree di V divise 3 + 2 in tutti i 10 modi; D della parte da 2 contro quella da 3 (stimatore corretto)."""
    vals = []
    for idx in combinations(range(5), 2):
        v = []
        for c in cells:
            if len(c["V"]) != 5: continue
            two = c["V"][list(idx)]; three = np.delete(c["V"], list(idx))
            v.append(((two.mean() - three.mean()) ** 2 - c["sp"] ** 2 / 2 - c["sp"] ** 2 / 3) / c["h"] ** 2)
        if v: vals.append(ssqrt(np.mean(v)))
    return float(np.mean(vals)) if vals else np.nan

def synth_check(cells, reps=200, seed=20261124):
    """G-3 (corretto il 2026-11-20): spostamenti costruiti 0, 1, 2 h con le numerosita' e le deviazioni osservate
    (V con la sua deviazione fra aree, gli insiemi di T con la combinata); stima con lo stesso stimatore di D(T)."""
    rng = np.random.default_rng(seed); out = {}
    for shift in (0.0, 1.0, 2.0):
        ests = []
        for _ in range(reps):
            cs = []
            for c in cells:
                sV = c["V"].std(ddof=1); V = rng.normal(0, sV, len(c["V"]))
                T = {k: rng.normal(shift * c["h"], c["sp"], len(o)) for k, o in c["T"].items()}
                sp, df = pooled_sd([V] + list(T.values())); cs.append({**c, "V": V, "T": T, "sp": sp, "df": df})
            ests.append(q_arm(cs, "T", ("P1", "P5")))
        out[str(shift)] = {"mean_square": float(np.mean(ests)), "root": ssqrt(np.mean(ests))}
    ok = abs(out["0.0"]["mean_square"]) <= 0.15 and abs(out["1.0"]["root"] - 1) <= 0.25 and abs(out["2.0"]["root"] - 2) <= 0.25
    return out, bool(ok)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--table"); ap.add_argument("--h"); ap.add_argument("--genes", default="panel")
    ap.add_argument("--out"); ap.add_argument("--tissues", default="A2,A3")
    a = ap.parse_args()
    tab = pd.read_csv(a.table); H = pd.read_csv(a.h)
    Vs = tab[(tab.arm == "V") & (tab.genes == a.genes)].groupby("A")["M4"].mean()
    drop = {A: bool(a.genes == "panel" and Vs.get(A, 0) >= 0.99) for A in a.tissues.split(",")}
    res = {"genes": a.genes, "drop_M4": drop, "V_singlets": Vs.to_dict(), "by": {}}
    groups = {"A2+A3": a.tissues.split(",")} | {A: [A] for A in a.tissues.split(",")}
    arms = sorted(set(tab[(tab.genes == a.genes)].arm) - {"V", "T"})
    for g, ts in groups.items():
        cells = prepare(tab, H, a.genes, ts, drop)
        Dobs, Dbs = analyse(cells, arms)
        ci = {k: [float(np.nanpercentile(Dbs[k], 2.5) - (np.nanmedian(Dbs[k]) - Dobs[k])), float(np.nanpercentile(Dbs[k], 97.5) - (np.nanmedian(Dbs[k]) - Dobs[k]))] if not np.isnan(Dobs[k]) else None for k in Dobs}
        comps = [compare(Dobs, Dbs, x, y) for x, y in [("S", "T"), ("T_P1", "T_P5"), ("XC2", "X1"), ("X1", "T"), ("X", "T"), ("XC", "X1"), ("S", "XC"), ("S", "X1"), ("S", "B"), ("B", "T"), ("XC", "T")]]
        g3, g3ok = synth_check(cells)
        res["by"][g] = {"D": Dobs, "D_ci": ci, "comparisons": comps, "G3_synthetic": g3, "G3_pass": g3ok, "split_3v2_descriptive": split_check(cells), "n_cells": len(cells)}
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    for g, r in res["by"].items():
        print(g, {k: round(v, 2) for k, v in r["D"].items() if not np.isnan(v)}, "G-3", r["G3_pass"], {k: round(v["root"], 3) for k, v in r["G3_synthetic"].items()}, "split", round(r["split_3v2_descriptive"], 3))
        for cmp in r["comparisons"]:
            if cmp["outcome"] != "not available": print("  ", cmp["a"], "-", cmp["b"], round(cmp["delta_ab"], 2), "w", round(cmp["w"], 2), cmp["outcome"], cmp["better"])
