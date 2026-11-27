"""Independent recomputation of A5 (gain), CP-gain, A4 (chains U and C) and of the h margins. Own code."""
import json, itertools, numpy as np, pandas as pd
from scipy import stats
REPO = "/home/user/2025.geo_spatialtrans/"; G = "/mnt/micron/geo_spatialtrans/"; OUT = G + "G1_review/rev_indip/"
res = {}
# ---------- A5 on gain
pu = pd.read_csv(REPO + "results/sonda/S8/per_unit.csv"); pu["gain"] = pu.score - pu.score_shuf
meth = ["E-naive", "E-area", "E-gene", "E-cluster", "b2c"]
res["per_unit_rows"] = len(pu); res["per_unit_dups"] = int(pu.duplicated(["arm", "unit", "method"]).sum())
res["X1_i"] = sorted(pu[pu.arm == "X1"].i.unique().tolist()); res["B_i"] = sorted(pu[pu.arm == "B"].i.unique().tolist())
def tab(col):
    out = {}
    for arm in ("B", "X1"):
        w = pu[pu.arm == arm].pivot_table(index="i", columns="method", values=col, aggfunc="mean")[meth]
        assert list(w.index) == list(range(10)); out[arm] = w.to_numpy()
    return out
def tau(a, b):
    t = stats.kendalltau(a, b).statistic; return 0.0 if not np.isfinite(t) else float(t)
tumour, stroma = [0, 1, 2, 3, 4], [5, 6, 7, 8, 9]
splits = []
for a in itertools.combinations(tumour, 3):
    for b in itertools.combinations(stroma, 2):
        h1 = list(a) + list(b); splits.append((h1, [x for x in tumour + stroma if x not in h1]))
def triple(S, h1, h2):
    X, B = S["X1"], S["B"]; mx1, mx2, mb1, mb2 = X[h1].mean(0), X[h2].mean(0), B[h1].mean(0), B[h2].mean(0)
    t5 = tau(mx1, mx2); s5 = (tau(mb1, mx2) + tau(mb2, mx1)) / 2; return t5, s5, s5 - t5
def a5(col):
    S = tab(col); est = np.mean([triple(S, *sp) for sp in splits], 0)
    rng = np.random.default_rng(20261201); bs = []
    for _ in range(2000):
        h1, h2 = splits[rng.integers(len(splits))]
        def rs(h):
            o = []
            for grp in (tumour, stroma):
                part = [x for x in h if x in grp]; o += [int(v) for v in rng.choice(part, size=len(part), replace=True)]
            return o
        r1 = rs(h1); r2 = rs(h2); bs.append(triple(S, r1, r2))
    bs = np.array(bs); se = bs.std(0, ddof=1)
    return {"T5": est[0], "S5c": est[1], "A5": est[2], "T5_se": se[0], "S5c_se": se[1], "A5_se": se[2], "n_splits": len(splits)}, S, bs
res["A5_score"], _, _ = a5("score")
r, S, bs = a5("gain"); res["A5_gain"] = r
T5, T5s, A5, A5s = r["T5"], r["T5_se"], r["A5"], r["A5_se"]
res["A5_gain"]["outcome"] = ("not measurable" if T5 - 2*T5s <= 0 <= T5 + 2*T5s else "holds" if A5 - 2*A5s >= -0.2 else "does not hold" if A5 + 2*A5s < -0.2 else "uncertain")
dec_boot = pd.read_csv(REPO + "results/sonda/S9/A5_gain_boot.csv"); res["A5_boot_maxdiff_vs_declared"] = float(np.abs(dec_boot.to_numpy() - bs).max())
# per-split T5 values for diagnosis
S = tab("gain"); pt = np.array([triple(S, *sp) for sp in splits]); res["T5_gain_per_split_unique"] = sorted(set(np.round(pt[:, 0], 12).tolist()))
res["gain_X1_means"] = dict(zip(meth, S["X1"].mean(0).tolist())); res["gain_B_means"] = dict(zip(meth, S["B"].mean(0).tolist()))
# CP-gain
units = sorted(pu[pu.arm == "X1"].unit.unique()); pan = []
for u in units:
    c = pd.read_parquet(f"{G}S8_sonda/X1/{u}/cells.parquet"); pan.append([c[f"panel_{m}"].mean() for m in meth])
pan = np.array(pan).mean(0); res["panel_umi_X1"] = dict(zip(meth, pan.tolist())); res["CP_gain_tau"] = tau(S["X1"].mean(0), pan)
# ---------- A4 with own paired-D code
f = np.sqrt(5/4) * stats.t.ppf(0.975, 4) / stats.norm.ppf(0.975)
def paired(table, hfile):
    U = pd.read_csv(table); H = pd.read_csv(hfile); U["kk"] = U.unit.str.rsplit("_", n=1).str[-1].astype(int)
    X, h = {}, {}
    for A in ("tumour", "stroma"):
        u = U[U.A == A]; vm = u[u.arm == "V"].groupby("kk")[["M1","M2","M3","M4","M5","M6"]].mean()
        keep = [m for m in vm.columns if not (m == "M4" and vm.M4.mean() >= 0.99)]
        h[A] = np.array([H[(H.A == A) & (H.measure == "M-" + m[1])].h.values[0] for m in keep])
        X[A] = {a: u[u.arm == a].groupby("kk")[keep].mean().sort_index().to_numpy() for a in ("V", "T", "B")}
    return X, h
def Dfun(X, h, arm, idx=None):
    q = []
    for A in ("tumour", "stroma"):
        a, v = X[A][arm], X[A]["V"]
        if idx is not None: a, v = a[idx[A]], v[idx[A]]
        d = a - v; q.append((d.mean(-2)**2 - d.var(-2, ddof=1)/5) / h[A]**2)
    m = np.concatenate(q, -1).mean(-1); return np.sign(m) * np.sqrt(np.abs(m))
rng = np.random.default_rng(20261123); idx = {"tumour": rng.integers(0, 5, (2000, 5)), "stroma": rng.integers(0, 5, (2000, 5))}
for ch, tb, hf in (("U", "results/sonda/S9/U/area_values_U.csv", "results/sonda/S9/U/h_U.csv"), ("C", "results/sonda/S9/C/area_values_C.csv", "results/sonda/S9/C/h_C.csv")):
    X, h = paired(REPO + tb, REPO + hf); DT, DB = float(Dfun(X, h, "T")), float(Dfun(X, h, "B"))
    diff = Dfun(X, h, "B", idx) - Dfun(X, h, "T", idx); e = DB - (DT + 1); se = f * np.std(diff, ddof=1)
    res["A4_" + ch] = {"tetto": DT, "D_B": DB, "e": e, "e_se": se, "outcome": "holds" if e + 2*se <= 0 else "does not hold" if e - 2*se > 0 else "uncertain"}
# ---------- h recomputation (S3b h.csv) by trying seed-index conventions
H = pd.read_csv(REPO + "results/sonda/S3b/h.csv"); V = pd.read_csv(REPO + "results/sonda/S3b/area_values.csv"); V = V[V.arm == "V"].copy()
V["kk"] = V.unit.str.rsplit("_", n=1).str[-1].astype(int); best = None
for itum, istr in ((1, 2), (0, 1), (1, 0), (0, 0)):
    for j0 in (0, 1):
        err = 0.0
        for A, i in (("tumour", itum), ("stroma", istr)):
            v = V[V.A == A].sort_values("kk")
            for jj, m in enumerate(["M1","M2","M3","M4","M5","M6"]):
                x = v[m].to_numpy(); g = np.random.default_rng([20261115, i, jj + j0]); mb = x[g.integers(0, 5, (2000, 5))].mean(1)
                lo, hi = np.percentile(mb, [2.5, 97.5]); row = H[(H.A == A) & (H.measure == "M-" + m[1])].iloc[0]
                hh = (hi - max(0, lo)) / 2 if m == "M2" else (min(1, hi) - max(0, lo)) / 2
                if m == "M4": hh = max(hh, 0.01)
                err = max(err, abs(lo - row.lo), abs(hi - row.hi), abs(hh - row.h))
        if best is None or err < best[0]: best = (err, itum, istr, j0)
res["h_recompute"] = {"max_abs_err": best[0], "i_tumour": best[1], "i_stroma": best[2], "j_offset": best[3]}
json.dump(res, open(OUT + "rev_A.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items() if k not in ("T5_gain_per_split_unique",)}, indent=0, default=float))
print("T5 per split unique values", res["T5_gain_per_split_unique"])
