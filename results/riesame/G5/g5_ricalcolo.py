"""Revisione G-5: ricalcolo indipendente di D e dei confronti (geni 'common').
Scritto da zero dalla definizione del compito; non usa lo script originale."""
import sys, json
import numpy as np, pandas as pd

ROOT = "/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/"
F_AREA = ROOT + "c6bf1c49-787d-44c9-900f-d7a1cdef4644/v8a322235_B_area_values.csv"
F_H = ROOT + "fc857421-7019-4f5b-b75d-8afdfb6a0cd2/va977e54a_B_h_panel.csv"
MEAS = ["M1", "M2", "M3", "M4", "M5", "M6"]
TISS = ["A2", "A3"]

def tset_label(s):
    if "_P1_" in s: return "T_P1"
    if "_P5_" in s: return "T_P5"
    if "6p5" in s: return "T_6p5"
    raise ValueError(s)

def load():
    a = pd.read_csv(F_AREA)
    a = a[a.genes == "common"].copy()
    h = pd.read_csv(F_H)
    h = h[h.genes == "common"].copy()
    h["m"] = h.measure.str.replace("-", "")
    H = {(r.A, r.m): r.h for r in h.itertuples()}
    HV = {(r.A, r.m): r.V_mean for r in h.itertuples()}
    a["grp"] = a.arm
    t = a.arm == "T"
    a.loc[t, "grp"] = a.loc[t, "set"].map(tset_label)
    # media dei semi dentro ogni area (S, B); per V, T una riga per area
    area = a.groupby(["A", "grp", "unit"])[MEAS].mean().reset_index()
    return a, area, H, HV

def cell_stats(area, sp_6p5_A3=True):
    """per (tessuto, misura): media, n, sd per gruppo; s_p e df per i bracci reali"""
    st = {}
    for A in TISS:
        for m in MEAS:
            d = {}
            for g, sub in area[area.A == A].groupby("grp"):
                x = sub[m].to_numpy()
                d[g] = dict(mean=x.mean(), n=len(x), sd=x.std(ddof=1) if len(x) > 1 else np.nan,
                            ss=((x - x.mean())**2).sum())
            real = ["V", "T_P1", "T_P5"] + (["T_6p5"] if (A == "A2" or sp_6p5_A3) else [])
            real = [g for g in real if g in d]
            ss = sum(d[g]["ss"] for g in real); df = sum(d[g]["n"] - 1 for g in real)
            st[(A, m)] = dict(groups=d, sp=np.sqrt(ss / df), df=df)
    return st

def q_terms(st, H, means=None):
    """q per cella e gruppo; means = dict (A,m)->{g: media} per il bootstrap"""
    out = {}
    for (A, m), c in st.items():
        d = c["groups"]; sp = c["sp"]; h = H[(A, m)]
        mu = {g: d[g]["mean"] for g in d} if means is None else means[(A, m)]
        nV = d["V"]["n"]
        for g in d:
            if g == "V": continue
            diff2 = (mu[g] - mu["V"])**2
            if g.startswith("T_"):
                corr = sp**2 / d[g]["n"] + sp**2 / nV
            else:
                corr = sp**2 / nV + d[g]["sd"]**2 / d[g]["n"]
            out[(A, m, g)] = (diff2 - corr) / h**2
    return out

def sroot(x): return np.sign(x) * np.sqrt(np.abs(x))

GROUPS = {"A2+A3": TISS, "A2": ["A2"], "A3": ["A3"]}
ARMS = {"T": ["T_P1", "T_P5"], "T_P1": ["T_P1"], "T_P5": ["T_P5"], "T_6p5": ["T_6p5"], "S": ["S"], "B": ["B"]}

def Dvals(q):
    res = {}
    for G, tl in GROUPS.items():
        for arm, sets in ARMS.items():
            v = [q[(A, m, g)] for A in tl for m in MEAS for g in sets if (A, m, g) in q]
            res[(G, arm)] = sroot(np.mean(np.array(v), axis=0)) if v else np.nan  # axis 0 = celle
    return res

def bootstrap(st, H, B=2000, seed=20261123):
    """bootstrap parametrico congiunto, vettorializzato sulle B repliche"""
    rng = np.random.default_rng(seed)
    means = {}
    for (A, m), c in st.items():
        d = c["groups"]; means[(A, m)] = {}
        for g in d:
            if g in ("S", "B"):
                se, df = d[g]["sd"] / np.sqrt(d[g]["n"]), d[g]["n"] - 1
            else:
                se, df = c["sp"] / np.sqrt(d[g]["n"]), c["df"]
            means[(A, m)][g] = d[g]["mean"] + se * rng.standard_t(df, size=B)
    q = q_terms(st, H, means)          # ogni valore e' un vettore di lunghezza B
    return Dvals(q)

COMPS = [("S", "T"), ("S", "B"), ("B", "T"), ("T_P1", "T_P5")]

def outcome(dab, w, delta=1.0):
    if abs(dab) + w <= delta: return "equal"
    if abs(dab) - w > delta: return "different"
    return "undetermined"

def run(sp_6p5_A3=True, seeds=(20261123,), B=2000):
    a, area, H, HV = load()
    st = cell_stats(area, sp_6p5_A3)
    D = Dvals(q_terms(st, H))
    boots = [bootstrap(st, H, B, s) for s in seeds]
    rows = []
    for G in GROUPS:
        for x, y in COMPS:
            ws = []
            for bt in boots:
                dd = bt[(G, x)] - bt[(G, y)]
                lo, hi = np.percentile(dd, [2.5, 97.5]); ws.append((hi - lo) / 2)
            dab = D[(G, x)] - D[(G, y)]
            rows.append(dict(group=G, a=x, b=y, Dab=dab, w=ws[0], w_seeds=ws,
                             outcome=outcome(dab, ws[0])))
    return a, area, H, HV, st, D, rows

F_JSON = ROOT + "ed7d0a17-29bc-4534-a059-3ec609470a80/v6a56e1de_B_D_common.json"
SEEDS = [20261123] + list(range(1, 20))   # 20 semi per l'errore Monte Carlo di w

if __name__ == "__main__":
    J = json.load(open(F_JSON))
    decl = {(G, c["a"], c["b"]): c for G in J["by"] for c in J["by"][G]["comparisons"]
            if c["outcome"] != "not available"}
    a, area, H, HV, st, D, rows = run(True, seeds=SEEDS, B=2000)
    _, _, _, _, st0, D0, _ = run(False, seeds=(1,), B=10)   # variante: s_p stroma senza 6,5 mm
    out = []
    for G in GROUPS:
        for arm in ["T", "T_P1", "T_P5", "T_6p5", "S", "B"]:
            out.append(dict(group=G, quantity=f"D({arm})", mine=D[(G, arm)],
                            declared=J["by"][G]["D"][arm], diff=D[(G, arm)] - J["by"][G]["D"][arm],
                            mine_sp_no6p5_A3=D0[(G, arm)]))
    for r in rows:
        c = decl[(r["group"], r["a"], r["b"])]; ws = np.array(r["w_seeds"])
        lab = f"{r['a']}-{r['b']}"
        out.append(dict(group=r["group"], quantity=f"Dab {lab}", mine=r["Dab"], declared=c["delta_ab"],
                        diff=r["Dab"] - c["delta_ab"]))
        out.append(dict(group=r["group"], quantity=f"w {lab}", mine=r["w"], declared=c["w"],
                        diff=r["w"] - c["w"], w_mean_20seeds=ws.mean(), w_sd_20seeds=ws.std(ddof=1),
                        w_min_20seeds=ws.min(), w_max_20seeds=ws.max(),
                        z_declared=(c["w"] - ws.mean()) / ws.std(ddof=1)))
        out.append(dict(group=r["group"], quantity=f"esito {lab}", mine_text=r["outcome"],
                        declared_text=c["outcome"],
                        esiti_20seeds="|".join(sorted({outcome(r["Dab"], w) for w in ws}))))
    pd.DataFrame(out).to_csv("g5_ricalcolo_numeri.csv", index=False)
    print(pd.DataFrame(out).to_string())
