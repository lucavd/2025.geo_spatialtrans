"""Riesame, Parte A, passo 1: il tetto (scarto fra sezioni vere dello stesso tessuto) e lo scarto del simulatore, in unita' del margine h.

Legge solo tabelle di F3.3 gia' committate (results/F3.3/F3.3_area_values.csv, F3.3_cells.csv); nessun file sigillato,
nessuna simulazione nuova, nessun cambiamento giudicato. Stimatore quadratico corretto di DEC-098 (D^2 - Var D) / h^2,
riportato come radice con segno; intervallo: bootstrap parametrico (modello t con la deviazione standard combinata,
simulato ridisegnato con il suo errore standard fra aree), percentile corretto per lo scarto della mediana. B = 2000, seme 20261121.
Uscite: results/riesame/tetto_pooled.csv, tetto_by_measure.csv, single_section_tissues.csv.
"""
import importlib.util, os
from itertools import combinations
import numpy as np, pandas as pd
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("F", f"{ROOT}/tools/F3_3_stats.py"); F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
OUT = f"{ROOT}/results/riesame"; os.makedirs(OUT, exist_ok=True)
B, SEED = 2000, 20261121
MEAS = ["M1", "M2", "M3", "M4", "M5", "M6"]; MID = {m: f"M-{m[1]}" for m in MEAS}
NAMES = {"M-1": "assigned UMI fraction", "M-2": "UMI per cell", "M-3": "foreign markers", "M-4": "singlets",
         "M-5": "RCTD-nucleus agreement", "M-6": "same-type neighbours"}
cells = pd.read_csv(f"{ROOT}/results/F3.3/F3.3_cells.csv"); av = pd.read_csv(f"{ROOT}/results/F3.3/F3.3_area_values.csv")
H = cells.drop_duplicates(["A", "measure"]).set_index(["A", "measure"])[["h"]]
REAL = av[av.src.isin(["IN", "OUT"])]; SIMU = av[av.src == "SIM"].groupby(["A", "u2"])[MEAS].mean().reset_index()

def build():
    """celle come in tools/F3.3_analysis.py (stessa deviazione standard combinata, stesso giudizio dentro campione)."""
    out = []
    for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
        R = REAL[REAL.A == A]; IN = R[R.set == f"{A}_IN"]; S = SIMU[SIMU.A == A]
        outs = {s: g for s, g in R[R.src == "OUT"].groupby("set")}
        if not outs: continue
        for m in MEAS:
            h = float(H.loc[(A, MID[m]), "h"]); rin = IN[m].to_numpy(float); sv = S[m].to_numpy(float)
            sp, df = F.pooled_sd([rin] + [g[m].to_numpy(float) for g in outs.values()])
            Sm = float(sv.mean()); se = float(sv.std(ddof=1) / np.sqrt(len(sv))); t = stats.t.ppf(0.975, df)
            v_in = F.vrule(Sm - rin.mean(), t * np.sqrt(se ** 2 + sp ** 2 / len(rin)), h)
            for lev in ("L1", "L2"):
                o = {s: g[m].to_numpy(float) for s, g in outs.items() if g.level.iloc[0] == lev}
                if o: out.append({"in": rin, "h": h, "sd": sp, "df": df, "sim": Sm, "se_sim": se, "good": v_in != "smentita",
                                  "A": A, "m": MID[m], "level": lev, "out": o})
    return out

vm = lambda c, n: c["sd"] ** 2 / n
def d_TV(cs, keep=None):
    v = [((o.mean() - c["in"].mean()) ** 2 - vm(c, len(o)) - vm(c, len(c["in"]))) / c["h"] ** 2
         for c in cs for k, o in c["out"].items() if keep is None or keep(k)]
    return np.mean(v) if v else np.nan
def d_TT(cs, keep=None):
    v = []
    for c in cs:
        ks = [k for k in c["out"] if keep is None or keep(k)]
        for a, b in combinations(ks, 2):
            oa, ob = c["out"][a], c["out"][b]
            v.append(((oa.mean() - ob.mean()) ** 2 - vm(c, len(oa)) - vm(c, len(ob))) / c["h"] ** 2)
    return np.mean(v) if v else np.nan
def d_SV(cs, keep=None):
    return np.mean([((c["sim"] - c["in"].mean()) ** 2 - vm(c, len(c["in"])) - c["se_sim"] ** 2) / c["h"] ** 2 for c in cs])
def d_ST(cs, keep=None):
    v = [((c["sim"] - o.mean()) ** 2 - vm(c, len(o)) - c["se_sim"] ** 2) / c["h"] ** 2
         for c in cs for k, o in c["out"].items() if keep is None or keep(k)]
    return np.mean(v) if v else np.nan
def bdraw(cs, rng):
    out = []
    for c in cs:
        s, df = c["sd"], c["df"]
        ni = np.full(len(c["in"]), c["in"].mean() + s / np.sqrt(len(c["in"])) * rng.standard_t(df))
        o2 = {k: np.full(len(o), o.mean() + s / np.sqrt(len(o)) * rng.standard_t(df)) for k, o in c["out"].items()}
        out.append({**c, "in": ni, "out": o2, "sim": c["sim"] + c["se_sim"] * rng.standard_normal()})
    return out
def est(cs, fun, keep=None, seed=SEED):
    obs = fun(cs, keep)
    if np.isnan(obs): return (np.nan, np.nan, np.nan)
    rng = np.random.default_rng(seed); bs = np.array([fun(bdraw(cs, rng), keep) for _ in range(B)])
    lo, hi, med = np.nanpercentile(bs, [2.5, 97.5, 50]); b = med - obs
    return F.ssqrt(obs), F.ssqrt(lo - b), F.ssqrt(hi - b)

CL = build(); L1 = [c for c in CL if c["level"] == "L1"]; L2 = [c for c in CL if c["level"] == "L2"]
# controllo: si ritrovano P-1 e P-2 di F3.3 (stime puntuali)
chk = {"P-1_L1": F.ssqrt(F.stat_QS(L1)), "P-2_L2": F.ssqrt(F.stat_QS(L2)), "P-3_L2": F.ssqrt(F.stat_QG(L2))}
print("check vs F3.3_summary:", {k: round(v, 4) for k, v in chk.items()})
pat = lambda k: "CRC_P" in k
groups = {"A2": [c for c in L2 if c["A"] == "A2"], "A3": [c for c in L2 if c["A"] == "A3"], "A5": [c for c in L2 if c["A"] == "A5"],
          "A2+A3": [c for c in L2 if c["A"] in ("A2", "A3")], "all_L2": L2}
rows = []
for g, cs in groups.items():
    for kn, keep in (("all_other_sections", None), ("P1_P5_only", pat)):
        if g == "A5" and kn != "all_other_sections": continue
        r = {"group": g, "comparison_set": kn, "n_cells": len(cs)}
        for nm, fun in (("D_T_vs_V", d_TV), ("D_T_vs_T", d_TT), ("D_S_vs_V", d_SV), ("D_S_vs_T", d_ST)):
            e = est(cs, fun, keep); r[nm], r[nm + "_lo"], r[nm + "_hi"] = e
        rows.append(r)
pd.DataFrame(rows).to_csv(f"{OUT}/tetto_pooled.csv", index=False)
pm = []
for i, c in enumerate(L2):
    tv, tt, sv, st = d_TV([c]), d_TT([c]), d_SV([c]), d_ST([c])
    ceil = np.nanmean([tv, tt]) if not np.isnan(tt) else tv
    diff = lambda cs, keep=None: d_SV(cs) - (np.nanmean([d_TV(cs), d_TT(cs)]) if not np.isnan(d_TT(cs)) else d_TV(cs))
    rng = np.random.default_rng(SEED + i); bs = np.array([diff(bdraw([c], rng)) for _ in range(B)])
    lo, hi, med = np.percentile(bs, [2.5, 97.5, 50]); b = med - diff([c]); lo, hi = lo - b, hi - b
    pm.append({"tissue": c["A"], "measure": c["m"], "name": NAMES[c["m"]], "h": c["h"], "V": c["in"].mean(), "sim": c["sim"],
               **{f"T_{k}": o.mean() for k, o in c["out"].items()},
               "D_TV": F.ssqrt(tv), "D_TT": F.ssqrt(tt) if not np.isnan(tt) else np.nan, "ceiling": F.ssqrt(ceil),
               "D_SV": F.ssqrt(sv), "D_ST": F.ssqrt(st), "diff_sq_h2": diff([c]), "diff_lo": lo, "diff_hi": hi,
               "class": "beyond ceiling" if lo > 0 else ("within ceiling" if hi < 0 else "undetermined"),
               "verdict_in_F33": "not rejected" if c["good"] else "rejected"})
PM = pd.DataFrame(pm); PM.to_csv(f"{OUT}/tetto_by_measure.csv", index=False)
print(PM.groupby(["tissue", "class"]).size().unstack(fill_value=0).to_string())
mx = float(pd.DataFrame(rows).query("comparison_set == 'all_other_sections' and group in ['A2','A3','A5']").D_T_vs_V.max())
pooled = float(pd.DataFrame(rows).query("group == 'all_L2' and comparison_set == 'all_other_sections'").D_T_vs_V.iloc[0])
one = cells[cells.A.isin(["A1", "A4", "A6"])].drop_duplicates(["A", "measure"])[["A", "measure", "gap_in_over_h", "verdict_in"]].copy()
one["name"] = one.measure.map(NAMES); one["beyond_max_ceiling"] = one.gap_in_over_h.abs() > mx; one["within_pooled_ceiling"] = one.gap_in_over_h.abs() <= pooled
one["max_ceiling_used"] = mx; one["pooled_ceiling_used"] = pooled
one.to_csv(f"{OUT}/single_section_tissues.csv", index=False)
print("max ceiling", round(mx, 3), "pooled", round(pooled, 3)); print(one.groupby("A")[["beyond_max_ceiling", "within_pooled_ceiling"]].sum().to_string())
