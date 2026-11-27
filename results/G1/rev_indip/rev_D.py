"""Independent recomputation of the paired D metric (REGOLA_S3c) for G1. Own code; reads only result tables."""
import json, sys, numpy as np, pandas as pd
from scipy import stats
REPO = "/home/user/2025.geo_spatialtrans/"
OUT = "/mnt/micron/geo_spatialtrans/G1_review/rev_indip/"
MEAS = ["M1","M2","M3","M4","M5","M6"]; TIS = ["tumour","stroma"]
def area_k(u): return u.str.rsplit("_", n=1).str[-1].astype(int)
t3 = pd.read_csv(REPO+"results/sonda/S3b/area_values.csv")
t4 = pd.read_csv(REPO+"results/sonda/S4/area_values.csv")
t6 = pd.read_csv(REPO+"results/sonda/S6/area_values.csv")
H = pd.read_csv(REPO+"results/sonda/S3b/h.csv")
allt = pd.concat([t3.assign(src="S3b"), t4.assign(src="S4"), t6.assign(src="S6")], ignore_index=True)
allt["kk"] = area_k(allt.unit); allt["seed"] = allt["k"].astype(str)
# arm -> per tissue 5x6 matrix ordered by area 1..5
def mat(arm, A, seed=None):
    u = allt[(allt.arm == arm) & (allt.A == A)]
    if seed is not None: u = u[u.seed == seed]
    g = u.groupby("kk")[MEAS].mean().sort_index()
    assert list(g.index) == [1,2,3,4,5], (arm, A, list(g.index))
    return g
X = {}; hvec = {}; keep = {}
for A in TIS:
    vm = mat("V", A)
    keep[A] = [m for m in MEAS if not (m == "M4" and vm["M4"].mean() >= 0.99)]
    hvec[A] = np.array([H[(H.A == A) & (H.measure == "M-"+m[1])].h.values[0] for m in keep[A]])
    X[A] = {"V": vm[keep[A]].to_numpy()}
    for arm in ["T","B","P","Y","G","U","A","AE","Yalt","YH1","Yshuf"]: X[A][arm] = mat(arm, A)[keep[A]].to_numpy()
    X[A]["B0"] = mat("B", A, "0")[keep[A]].to_numpy(); X[A]["B1"] = mat("B", A, "1")[keep[A]].to_numpy()
ARMS = ["T","B","P","B0","B1","Y","G","U","A","AE","Yalt","YH1","Yshuf"]
# row counts per arm (duplicates check)
cnt = allt.groupby(["src","arm","A"]).size().to_dict()
def qcells(Xa, Xv, h):
    d = Xa - Xv; n = d.shape[-2]
    return (d.mean(-2)**2 - d.var(-2, ddof=1)/n) / h**2
def Dval(arm, tis, idx=None):
    qs = []
    for A in tis:
        xa, xv = X[A][arm], X[A]["V"]
        if idx is not None: xa, xv = xa[idx[A]], xv[idx[A]]
        qs.append(qcells(xa, xv, hvec[A]))
    m = np.concatenate(qs, -1).mean(-1)
    return np.sign(m)*np.sqrt(np.abs(m))
n = 5; f = np.sqrt(n/(n-1))*stats.t.ppf(0.975, n-1)/stats.norm.ppf(0.975)
def boot(seed, NB=2000):
    rng = np.random.default_rng(seed); return {"tumour": rng.integers(0,5,(NB,5)), "stroma": rng.integers(0,5,(NB,5))}
def ci_ii(D0, Db):
    shift = np.median(Db) - D0; lo, hi = np.percentile(Db, [2.5, 97.5]) - shift
    return D0 + f*(lo-D0), D0 + f*(hi-D0)
GR = {"tumour+stroma": TIS, "tumour": ["tumour"], "stroma": ["stroma"]}
res = {"f": f, "keep": keep, "row_counts": {"|".join(k): v for k, v in cnt.items()}, "by": {}}
idx = boot(20261123)
for g, ts in GR.items():
    r = {}
    for a in ARMS:
        D0 = float(Dval(a, ts)); Db = Dval(a, ts, idx); lo, hi = ci_ii(D0, Db)
        r[a] = {"D": D0, "lo": float(lo), "hi": float(hi), "se": float(f*np.std(Db, ddof=1)), "se_i": float(np.std(Db, ddof=1))}
    res["by"][g] = r
# Monte Carlo sensitivity of interval endpoints over 20 other seeds (tumour+stroma)
mc = {a: [] for a in ARMS}
for s in range(20):
    ix = boot(1000+s)
    for a in ARMS:
        D0 = float(Dval(a, TIS)); mc[a].append(ci_ii(D0, Dval(a, TIS, ix)))
res["mc_seed_sd"] = {a: [float(np.std([x[0] for x in v], ddof=1)), float(np.std([x[1] for x in v], ddof=1))] for a, v in mc.items()}
# exact bootstrap (all 5^5 resamples per tissue, joint by outer sum) for G1-relevant arms
import itertools
allidx = np.array(list(itertools.product(range(5), repeat=5)))
ex = {}
for a in ["T","B","Y","Yalt","YH1","Yshuf"]:
    sums = {A: qcells(X[A][a][allidx], X[A]["V"][allidx], hvec[A]).sum(-1) for A in TIS}
    ncell = sum(len(keep[A]) for A in TIS)
    m = (sums["tumour"][:, None] + sums["stroma"][None, :]).ravel()/ncell
    Db = np.sign(m)*np.sqrt(np.abs(m)); D0 = res["by"]["tumour+stroma"][a]["D"]
    lo, hi = ci_ii(D0, Db); ex[a] = [float(lo), float(hi)]
res["exact_boot_ci_ii"] = ex
json.dump(res, open(OUT+"rev_D.json", "w"), indent=1, default=float)
dec = json.load(open(REPO+"results/sonda/S6/D.json"))
rows = []
for g in GR:
    db = dec["by"][g]
    for a in ARMS:
        m = res["by"][g][a]
        rows.append({"group": g, "arm": a, "D_dec": db["D"][a], "D_rev": m["D"], "lo_dec": db["D_ci"][a][0], "lo_rev": m["lo"], "hi_dec": db["D_ci"][a][1], "hi_rev": m["hi"], "se_dec": db["D_se"][a], "se_rev": m["se"]})
P = pd.DataFrame(rows)
for c in ["D","lo","hi","se"]: P[c+"_absdiff"] = (P[c+"_dec"] - P[c+"_rev"]).abs()
P.to_csv(OUT+"rev_D_compare.csv", index=False)
print("f", f, "keep", keep)
print("max absdiff", {c: P[c+"_absdiff"].max() for c in ["D","lo","hi","se"]})
print(P[P.group=="tumour+stroma"][["arm","D_rev","lo_rev","hi_rev","se_rev"]].round(4).to_string(index=False))
print("MC sd", {a: [round(x,4) for x in v] for a, v in res["mc_seed_sd"].items() if a in ["T","B","Y","Yalt","YH1","Yshuf"]})
print("exact", {a: [round(x,4) for x in v] for a, v in ex.items()})
print("counts", {k: v for k, v in res["row_counts"].items()})
# ---- paired comparisons (rule of S3c: equal if |d|+w<=1, different if |d|-w>1), method (ii)
cmpr = []
for a, b in (("B","T"), ("P","T"), ("B0","B1"), ("Y","T"), ("Yalt","Y"), ("YH1","Y"), ("Yshuf","Y"), ("Yalt","YH1"), ("Yalt","T"), ("YH1","T"), ("Yalt","B"), ("YH1","B")):
    d = float(Dval(a, TIS) - Dval(b, TIS)); bd = Dval(a, TIS, idx) - Dval(b, TIS, idx)
    lo, hi = np.percentile(bd, [2.5, 97.5]); w = f*(hi-lo)/2
    for dl in (1.0, 0.5):
        out = "equal" if abs(d)+w <= dl else ("different" if abs(d)-w > dl else "undetermined")
        cmpr.append({"a": a, "b": b, "delta_margin": dl, "d": d, "w": float(w), "outcome": out})
C = pd.DataFrame(cmpr); C.to_csv(OUT+"rev_D_comparisons.csv", index=False)
dc = {(x["a"], x["b"]): x for x in dec["by"]["tumour+stroma"]["comparisons"]}
C1 = C[C.delta_margin == 1.0].copy()
C1["d_dec"] = [dc.get((r.a, r.b), {}).get("delta_ab") for r in C1.itertuples()]; C1["w_dec"] = [dc.get((r.a, r.b), {}).get("w") for r in C1.itertuples()]
C1["out_dec"] = [dc.get((r.a, r.b), {}).get("outcome") for r in C1.itertuples()]
print(C1.round(6).to_string(index=False)); print(C[C.delta_margin == 0.5][["a","b","d","w","outcome"]].round(4).to_string(index=False))
