#!/usr/bin/env python
"""tools/S1.4e_analysis.py - S1.4e: verdetti P-1...P-4, regola N, controprove CP-1 e CP-2, descrittive (definizioni: docs/steps/S1.4e.md).
Unita': P-1, P-2, CP-2 seme (finestre fisse; stima = media dei 10 semi +- 2 ES, ES = sd fra semi / sqrt(10)); P-3, P-4 ROI (valore della ROI = media
dei 10 semi; 5 ROI); CP-1 finestra (media dei 10 semi e delle coppie; 4 finestre). Esiti: confermata (stima +- 2 ES dentro l'intervallo),
smentita (tutta fuori), indeterminata (a cavallo). Uscite: results/S1.4e/S1.4e_verdicts.csv, S1.4e_ruleN.csv, S1.4e_cp.csv, S1.4e_desc.csv.
Uso (root del repo): .venv/bin/python tools/S1.4e_analysis.py"""
import json, numpy as np, pandas as pd
RES = "results/S1.4e"; ARMS = ["K2", "ER", "ET"]; T6 = ["A1", "A2", "A3", "A4", "A5", "A6"]; DENSE = ["A1", "A2", "A3", "A4"]
CONTACT = {"A1": (0.366, 0.650), "A2": (0.217, 0.567), "A3": (0.083, 0.283), "A4": (0.350, 0.750), "A5": (0.000, 0.217), "A6": (0.000, 0.070)}  # B-127...B-132
REF = pd.read_csv(f"{RES}/S1.4e_refs.csv")
def iv(ref, A, cls): r = REF[(REF.ref == ref) & (REF.archetype == A) & (REF.cls == cls)].iloc[0]; return float(r.lo), float(r.hi)
def est(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; return (v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v)) if len(v) > 1 else (np.nan, np.nan, len(v))
def verdict(m, se, lo, hi):
    if not np.isfinite(m): return "non valutabile"
    a, b = m - 2 * se, m + 2 * se
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"
rows = []
def add(pid, A, arm, v, lo, hi, unit, src):
    m, se, n = est(v); rows.append(dict(id=pid, archetype=A, arm=arm, unit=unit, est=m, se=se, n_units=n, lo=lo, hi=hi, verdict=verdict(m, se, lo, hi), source=src))
W = pd.read_csv(f"{RES}/S1.4e_windows.csv"); PB = pd.read_csv(f"{RES}/S1.4e_pairs_B.csv"); SH = pd.read_csv(f"{RES}/S1.4e_shape_B.csv"); PA = pd.read_csv(f"{RES}/S1.4e_pairs_A.csv")
assert (W.groupby(["arm", "archetype"]).seed.nunique() == 10).all() and (SH.groupby(["arm", "archetype", "roi_id"]).seed.nunique() == 10).all()
ps = W.groupby(["arm", "archetype", "seed"]).s1.mean().reset_index()
pb = PB[PB.noise == "none"].groupby(["arm", "archetype", "seed", "cls"])[["n", "n_contact"]].sum().reset_index(); pb["f"] = pb.n_contact / pb.n
pc = PB[PB.noise == "click"].groupby(["arm", "archetype", "seed", "cls", "draw"])[["n", "n_contact"]].sum().reset_index(); pc["f"] = pc.n_contact / pc.n
pc = pc.groupby(["arm", "archetype", "seed", "cls"]).f.mean().reset_index()
sr = SH.groupby(["arm", "archetype", "roi_id"])[["p3", "p4"]].mean().reset_index()
for arm in ARMS:
    for A in T6:
        lo, hi = CONTACT[A]; add("P-1", A, arm, ps[(ps.arm == arm) & (ps.archetype == A)].s1, lo, hi, "seme", "B-127...B-132")
        lo, hi = iv("R-3", A, "all"); add("P-4", A, arm, sr[(sr.arm == arm) & (sr.archetype == A)].p4, lo, hi, "ROI", "R-3")
        if A in DENSE:
            lo, hi = iv("R-1", A, "lt6"); add("P-2", A, arm, pb[(pb.arm == arm) & (pb.archetype == A) & (pb.cls == "lt6")].f, lo, hi, "seme", "R-1")
        if A in ["A1", "A2", "A3", "A4", "A6"]:
            lo, hi = iv("R-2", A, "nn<5"); add("P-3", A, arm, sr[(sr.arm == arm) & (sr.archetype == A)].p3, lo, hi, "ROI", "R-2")
V = pd.DataFrame(rows); V.to_csv(f"{RES}/S1.4e_verdicts.csv", index=False)
# identita' con S1.5e: P-1 di K2 = P-1 di RV
s5 = pd.read_csv("results/S1.5e/S1.5e_verdicts.csv"); s5 = s5[(s5.id == "P-1") & (s5.arm == "RV")].set_index("archetype")
k2 = V[(V.id == "P-1") & (V.arm == "K2")].set_index("archetype"); dmax = float((k2.est - s5.est.reindex(k2.index)).abs().max())
print("P-1 K2 contro RV di S1.5e, scarto massimo:", dmax)
# ---- regola N ----
sc = {a: int(((V.arm == a) & V.id.isin(["P-1", "P-2"]) & (V.verdict == "confermata")).sum()) for a in ARMS}
p3c = {a: int(((V.arm == a) & (V.id == "P-3") & (V.verdict == "confermata")).sum()) for a in ARMS}
k2c = set(V[(V.arm == "K2") & (V.id == "P-1") & (V.verdict == "confermata")].archetype)
def ok_b(a): return not any(V[(V.arm == a) & (V.id == "P-1") & V.archetype.isin(k2c)].verdict == "smentita")
elig = [a for a in ["ER", "ET"] if sc[a] >= sc["K2"] + 2 and ok_b(a)]
if not elig: win, why = "K2", "nessuna ellisse supera K2 di almeno 2 voci senza peggiorare un tessuto confermato"
elif len(elig) == 1: win, why = elig[0], "unica ammessa"
elif sc["ER"] != sc["ET"]: win = max(elig, key=lambda a: sc[a]); why = "punteggio piu' alto"
elif p3c["ER"] != p3c["ET"]: win = max(elig, key=lambda a: p3c[a]); why = "pari: piu' P-3 confermate"
else: win, why = "K2", "pari anche su P-3"
RN = pd.DataFrame([dict(arm=a, score_P1_P2=sc[a], P3_confirmed=p3c[a], cond_a=(a == "K2") or sc[a] >= sc["K2"] + 2, cond_b=(a == "K2") or ok_b(a)) for a in ARMS])
RN["winner"] = win; RN["why"] = why; RN.to_csv(f"{RES}/S1.4e_ruleN.csv", index=False); print(RN.to_string())
# ---- controprove ----
cp = []
pa = PA[PA.d_um < 6].groupby(["archetype", "win_id", "arm"]).contact.mean().unstack("arm").reset_index()
for A in DENSE:
    g = pa[pa.archetype == A]
    for b in ["ET", "ER"]:
        m, se, n = est(g[b] - g["K2"])
        cp.append(dict(id="CP-1", archetype=A, contrast=f"{b}-K2", est=m, se=se, n_units=n, lo2=m - 2 * se, hi2=m + 2 * se,
                       below_zero=bool(m + 2 * se < 0), n_pairs=int((PA[(PA.archetype == A) & (PA.arm == "K2")].d_um < 6).sum() / 10)))
for arm in ["K2", "ET"]:
    for A in DENSE:
        lo, hi = iv("R-1", A, "lt6"); m, se, n = est(pc[(pc.arm == arm) & (pc.archetype == A) & (pc.cls == "lt6")].f)
        m0, _, _ = est(pb[(pb.arm == arm) & (pb.archetype == A) & (pb.cls == "lt6")].f)
        cp.append(dict(id="CP-2", archetype=A, contrast=arm, est=m, se=se, n_units=n, lo2=m - 2 * se, hi2=m + 2 * se, verdict=verdict(m, se, lo, hi), lo=lo, hi=hi, est_no_noise=m0))
CP = pd.DataFrame(cp); CP.to_csv(f"{RES}/S1.4e_cp.csv", index=False)
n_et = int(CP[(CP.id == "CP-1") & (CP.contrast == "ET-K2")].below_zero.sum()); n_er = int((CP[(CP.id == "CP-1") & (CP.contrast == "ER-K2")].below_zero.astype(bool) == False).sum())   # D-1: ~ su colonna object
print(f"CP-1: ET-K2 < 0 oltre 2 ES in {n_et}/4 (attesa >= 3); ER-K2 non negativa oltre 2 ES in {n_er}/4 (attesa >= 3) ->",
      "superata" if (n_et >= 3 and n_er >= 3) else "FALLITA")
k2cp2 = CP[(CP.id == "CP-2") & (CP.contrast == "K2")]
print("CP-2: K2 con il rumore dei punti:", dict(zip(k2cp2.archetype, k2cp2.verdict)))
# ---- descrittive ----
D = []
for arm in ARMS:                                          # differenze appaiate per seme di P-1 contro K2
    for A in T6:
        if arm == "K2": continue
        a_ = ps[(ps.arm == arm) & (ps.archetype == A)].set_index("seed").s1; b_ = ps[(ps.arm == "K2") & (ps.archetype == A)].set_index("seed").s1
        m, se, n = est((a_ - b_.reindex(a_.index)).to_numpy()); D.append(dict(desc="P-1 diff vs K2 (per seme)", archetype=A, arm=arm, est=m, se=se, n=n))
    for A in DENSE:
        for c in ["lt6", "6to10"]:
            m, se, n = est(pb[(pb.arm == arm) & (pb.archetype == A) & (pb.cls == c)].f); D.append(dict(desc=f"P-2 classe {c}", archetype=A, arm=arm, est=m, se=se, n=n))
it = pd.read_csv("results/R2c/R2c_items_Luca.csv"); an = pd.DataFrame(json.load(open("results/R2c/annotations/R2c_answers_Luca.json"))["answers"])
J = it.merge(an, on="item"); J = J[(J.kind == "pair") & (J.answer != "unsure")]; J["A"] = J.win_id.str[:2]
BINS = [0, 5, 6, 7, 8, 10, np.inf]; LAB = ["<5", "5-6", "6-7", "7-8", "8-10", ">=10"]
PA["bin"] = pd.cut(PA.d_um, BINS, right=False, labels=LAB); J["bin"] = pd.cut(J.d_um, BINS, right=False, labels=LAB)
for (A, b), g in J.groupby(["A", "bin"], observed=True): D.append(dict(desc="banco A curva", archetype=A, arm="manuale", bin=b, est=float((g.answer == "touch").mean()), n=len(g)))
for (A, arm, b), g in PA.groupby(["archetype", "arm", "bin"], observed=True): D.append(dict(desc="banco A curva", archetype=A, arm=arm, bin=b, est=float(g.contact.mean()), n=int(len(g) / 10)))
pd.DataFrame(D).to_csv(f"{RES}/S1.4e_desc.csv", index=False)
cnt = V.groupby(["id", "arm"]).verdict.value_counts().unstack(fill_value=0); print(cnt.to_string())
print(V[V.id.isin(["P-1", "P-2"])].pivot_table(index=["id", "archetype"], columns="arm", values="est").round(3).to_string())
