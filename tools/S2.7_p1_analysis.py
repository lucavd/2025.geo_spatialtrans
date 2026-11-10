"""tools/S2.7_p1_analysis.py - S2.7, P-1 e CP-1: stessa legge del motore di S2.5 (docs/steps/S2.7.md). Misure per esecuzione da
tools/S2.5_stats.py (definizione reg del pannello di S2.5) sulle cartelle <PD>/old e <PD>/new. Rapporto pezzi / S2.5 delle medie sui semi,
errore standard col metodo delta dalle deviazioni standard fra semi; esito con stima +- 2 errori standard. CP-1: UMI per quadratino delle
uscite dei pezzi assottigliate con probabilita' 0.95 (binomiale per quadratino, seme 20261221).
Uso (root del repo): .venv/bin/python tools/S2.7_p1_analysis.py <PD> <outdir>"""
import sys, os, json, glob
import numpy as np, pandas as pd, shapely
sys.path.insert(0, "tools")
import importlib.util
PD, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
TIS = {"A1": "intestino", "A2": "tumore", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
def verdict(lo, hi, rlo, rhi):
    if lo >= rlo and hi <= rhi: return "confermata"
    if hi < rlo or lo > rhi: return "smentita"
    return "indeterminata"
def load(v):
    R = [json.load(open(f)) for f in sorted(glob.glob(f"{PD}/stats_{v}/*_k*.json"))]
    bad = [r for r in R if "error" in r]
    assert not bad, bad[0]
    return pd.DataFrame(R)
O = load("old"); N = load("new")
V5 = pd.read_csv("results/S2.5/S2.5_verdicts_main_r2.csv")
rows = []
for A in sorted(TIS):
    o = O[O.archetype == A]; n = N[N.archetype == A]
    for m, kind in [("reg_mean_umi2", "umi"), ("reg_mean_genes2", "genes"), ("reg_mean_genes8", "genes"), ("reg_frac_empty2", "empty")]:
        a = n[m].to_numpy(float); b = o[m].to_numpy(float)
        ma, mb = a.mean(), b.mean(); r = ma / mb if mb > 0 else np.nan
        se = r * np.sqrt(a.var(ddof=1) / len(a) / ma**2 + b.var(ddof=1) / len(b) / mb**2) if ma > 0 and mb > 0 else np.nan
        if kind == "empty":
            v = V5[(V5.archetype == A) & (V5.metric == "frac_empty2")].iloc[0]
            if float(v.rel_hi) == 0 or float(v.real_est) == 0: rlo, rhi, src = np.nan, np.nan, "vero 0: rapporto non definito"
            else: rlo, rhi, src = float(v.rel_lo) / float(v.real_est), float(v.rel_hi) / float(v.real_est), "intervallo corretto del vero (DEC-061), relativo"
        else: rlo, rhi, src = 0.98, 1.02, "margine del 2 % (S2.6)"
        lo, hi = r - 2 * se, r + 2 * se
        es = verdict(lo, hi, rlo, rhi) if np.isfinite(rlo) and np.isfinite(lo) else "non valutabile"
        rows.append(dict(id="P-1", archetype=A, tessuto=TIS[A], metric=m[4:], n_new=len(a), n_old=len(b), mean_new=ma, mean_old=mb, sd_new=a.std(ddof=1), sd_old=b.std(ddof=1),
                         ratio=r, se=se, lo=lo, hi=hi, rel_lo=rlo, rel_hi=rhi, rel_source=src, esito=es))
    # CP-1: assottigliamento 0.95 delle uscite dei pezzi, UMI per quadratino nella stessa regione
    rng = np.random.default_rng(20261221)
    th = []
    for _, rr in n.iterrows():
        base = f"{PD}/new/{rr.tag}_k{int(rr.k)}"
        u2 = np.fromfile(base + ".bin", dtype=np.int32)
        th.append(rr.reg_mean_umi2 * (rng.binomial(u2, 0.95).sum() / max(u2.sum(), 1)))
    a = np.array(th); b = o.reg_mean_umi2.to_numpy(float); ma, mb = a.mean(), b.mean(); r = ma / mb
    se = r * np.sqrt(a.var(ddof=1) / len(a) / ma**2 + b.var(ddof=1) / len(b) / mb**2)
    rows.append(dict(id="CP-1", archetype=A, tessuto=TIS[A], metric="mean_umi2_thinned_0.95", n_new=len(a), n_old=len(b), mean_new=ma, mean_old=mb, sd_new=a.std(ddof=1),
                     sd_old=b.std(ddof=1), ratio=r, se=se, lo=r - 2 * se, hi=r + 2 * se, rel_lo=0.98, rel_hi=1.02, rel_source="margine del 2 % (S2.6)",
                     esito=verdict(r - 2 * se, r + 2 * se, 0.98, 1.02) + (" per difetto" if r + 2 * se < 0.98 else "")))
P = pd.DataFrame(rows); P.to_csv(f"{OUT}/S2.7_P1.csv", index=False)
pd.set_option("display.width", 250); print(P.round(4).to_string(index=False))
print(P.groupby("id").esito.value_counts())
