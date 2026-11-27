"""tools/sonda/S9_A.py - S9 (docs/sonda/REGOLA_S9.md): A2 sulla catena di oggi con incertezza per voce; A5 col guadagno
(score - score_shuf) sui file di S8, con il codice di tools/sonda/S8_A5.py summary() copiato (identita' sul coseno: G-A5ID);
A4 per catena X (U, C) con le funzioni di tools/sonda/S3c_D.py invariate (stessi indici bootstrap, seme 20261123, metodo (ii), f).
Uso: .venv/bin/python -B tools/sonda/S9_A.py a2 [table h tag] | a5 | a4 <X> <table> <h>"""
import os, sys, json, glob, itertools, numpy as np, pandas as pd
from scipy import stats
from scipy.stats import kendalltau
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools/sonda")
import S3c_D as M
RES = "results/sonda/S9"; os.makedirs(RES, exist_ok=True)
N = 5; F = float(np.sqrt(N / (N - 1)) * stats.t.ppf(0.975, N - 1) / stats.norm.ppf(0.975))
S8 = "results/sonda/S8"; S8_OUT = "/mnt/micron/geo_spatialtrans/S8_sonda"; SEED8 = 20261201

def esito_basso(v, se, serve):
    return "regge" if v + 2 * se <= serve else ("non regge" if v - 2 * se > serve else "incerta")

def a2(table="results/sonda/S3b/area_values.csv", h="results/sonda/S3b/h.csv", tag="L"):
    data, drop, Vs = M.load(table, h); idx = M.boot_idx(np.random.default_rng(M.SEED), M.NB)
    rows, bs = [], []
    for A in M.TIS:
        X = data[A]["X"]; hh = data[A]["h"]
        qT = M.q_paired(X["T"], X["V"], hh); qB = M.q_paired(X["B"], X["V"], hh)
        bs.append(M.ssqrt(M.q_paired(X["B"][idx[A]], X["V"][idx[A]], hh)))
        for j, m in enumerate(data[A]["meas"]):
            rows.append(dict(A=A, measure=m, T_d=float(M.ssqrt(qT[j])), B_d=float(M.ssqrt(qB[j]))))
    P = pd.DataFrame(rows); bs = np.concatenate(bs, axis=1); Tmax = float(P.T_d.max())
    P["se"] = F * bs.std(axis=0, ddof=1)
    P["class"] = np.where(P.B_d - 2 * P.se > Tmax, "oltre", np.where(P.B_d + 2 * P.se <= Tmax, "dentro", "indeterminata"))
    cnt = int((P.B_d > Tmax).sum()); cb = (bs > Tmax).sum(axis=1); se = float(F * cb.std(ddof=1))
    out = dict(chain=tag, T_max=Tmax, T_max_voce=" ".join(P.loc[P.T_d.idxmax(), ["A", "measure"]]), A2=cnt, A2_se=se,
               esito=esito_basso(cnt, se, 0), n_oltre=int((P["class"] == "oltre").sum()), n_indeterminata=int((P["class"] == "indeterminata").sum()),
               n_dentro=int((P["class"] == "dentro").sum()), f=F, seed=M.SEED, B=M.NB, table=table, h=h)
    P.to_csv(f"{RES}/A2_per_voce_{tag}.csv", index=False); json.dump(out, open(f"{RES}/A2_{tag}.json", "w"), indent=1)
    print(json.dumps(out)); print(P.to_string())
    return out

def taub(a, b):
    t = kendalltau(a, b).statistic; return float(t) if np.isfinite(t) else 0.0

def a5_core(pa_, col, METHODS):
    S = {arm: pa_[pa_.arm == arm].groupby(["i", "method"])[col].mean().unstack()[METHODS].sort_index().to_numpy() for arm in ("B", "X1")}
    TU, ST = list(range(5)), list(range(5, 10)); splits = [(list(a) + list(b), [x for x in TU if x not in a] + [x for x in ST if x not in b]) for a in itertools.combinations(TU, 3) for b in itertools.combinations(ST, 2)]
    def st(h1, h2):
        m = lambda Z, h: Z[h].mean(axis=0)
        t5 = taub(m(S["X1"], h1), m(S["X1"], h2)); s5 = 0.5 * (taub(m(S["B"], h1), m(S["X1"], h2)) + taub(m(S["B"], h2), m(S["X1"], h1)))
        return t5, s5, s5 - t5
    P = np.array([st(h1, h2) for h1, h2 in splits]); est = P.mean(axis=0)
    rng = np.random.default_rng(SEED8); BS = []
    for _ in range(2000):
        h1, h2 = splits[rng.integers(len(splits))]
        r = lambda h: [int(v) for t in (TU, ST) for v in rng.choice([x for x in h if x in t], size=len([x for x in h if x in t]), replace=True)]
        BS.append(st(r(h1), r(h2)))
    BS = np.array(BS); se = BS.std(axis=0, ddof=1)
    return dict(T5=float(est[0]), S5c=float(est[1]), A5=float(est[2]), T5_se=float(se[0]), S5c_se=float(se[1]), A5_se=float(se[2])), S, BS

def a5():
    ref = json.load(open(f"{S8}/A5.json")); METHODS = ref["methods"]; pa_ = pd.read_csv(f"{S8}/per_unit.csv")
    chk, _, _ = a5_core(pa_, "score", METHODS)
    gid = {k: abs(chk[k] - ref[k]) for k in ("T5", "S5c", "A5", "T5_se", "S5c_se", "A5_se")}
    pa_["gain"] = pa_.score - pa_.score_shuf
    out, S, BS = a5_core(pa_, "gain", METHODS)
    T5, T5s, A5, A5s = out["T5"], out["T5_se"], out["A5"], out["A5_se"]
    if T5 - 2 * T5s <= 0 <= T5 + 2 * T5s: es = "non misurabile a questa profondita"
    elif A5 - 2 * A5s >= -0.2: es = "regge"
    elif A5 + 2 * A5s < -0.2: es = "non regge"
    else: es = "incerta"
    units = sorted(pa_[pa_.arm == "X1"].unit.unique()); pan = []
    for u in units:
        c = pd.read_parquet(f"{S8_OUT}/X1/{u}/cells.parquet"); pan.append([c[f"panel_{m}"].mean() for m in METHODS])
    pan = np.array(pan).mean(axis=0); gainX1 = S["X1"].mean(axis=0); cp = taub(gainX1, pan)
    out.update(esito=es, methods=METHODS, gain_B=dict(zip(METHODS, S["B"].mean(axis=0).tolist())), gain_X1=dict(zip(METHODS, gainX1.tolist())),
               panel_umi_X1=dict(zip(METHODS, pan.tolist())), CP_gain_tau=cp, CP_gain_pass=bool(cp < 0.8), n_X1_units=len(units),
               G_A5ID_maxdiff=max(gid.values()), G_A5ID_pass=bool(max(gid.values()) < 1e-12), seed=SEED8, B=2000)
    json.dump(out, open(f"{RES}/A5_gain.json", "w"), indent=1); pd.DataFrame(BS, columns=["T5", "S5c", "A5"]).to_csv(f"{RES}/A5_gain_boot.csv", index=False)
    print(json.dumps(out)); return out

def a4(X, table, h):
    data, drop, Vs = M.load(table, h); idx = M.boot_idx(np.random.default_rng(M.SEED), M.NB)
    res = {"chain": X, "table": table, "h": h, "f": F, "M4_dropped": drop, "V_M4_mean": Vs}
    for g, ts in M.GROUPS.items():
        Dobs = {k: float(v) for k, v in M.D_of(data, M.ARMS, ts).items()}; Dbs = M.D_of(data, M.ARMS, ts, idx)
        r = {"D": Dobs, "D_se": {k: float(F * np.std(Dbs[k], ddof=1)) for k in M.ARMS}, "ci": {k: M.intervals(Dobs[k], Dbs[k], F)[1] for k in M.ARMS},
             "comp": {f"{a}-{b}": M.compare(Dobs, Dbs, a, b, F) for a, b in (("B", "T"), ("P", "T"), ("B0", "B1"), ("P", "B"))}}
        if g == "tumour+stroma":
            e = Dobs["B"] - (Dobs["T"] + 1.0); se = float(F * np.std(Dbs["B"] - Dbs["T"], ddof=1))
            c = r["comp"]
            r.update(e=e, e_se=se, esito=esito_basso(e, se, 0.0), tetto=Dobs["T"], pavimento=Dobs["B"], soglia=Dobs["T"] + 1.0,
                     CP_chain_pass=c["B-T"]["outcome"] == "different", C4_pass=c["P-T"]["outcome"] == "different" and c["P-T"]["better"] == "T",
                     C5_pass=c["B0-B1"]["outcome"] != "different")
        res[g] = r
    M.per_measure(data).to_csv(f"{RES}/per_measure_{X}.csv", index=False)
    json.dump(res, open(f"{RES}/A4_{X}.json", "w"), indent=1, default=float); print(json.dumps(res["tumour+stroma"], default=float)); return res

if __name__ == "__main__":
    c = sys.argv[1]
    if c == "a2": a2(*sys.argv[2:5]) if len(sys.argv) > 2 else a2()
    elif c == "a5": a5()
    elif c == "a4": a4(*sys.argv[2:5])
