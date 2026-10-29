#!/usr/bin/env python
"""tools/S1.4f_analysis.py - S1.4f: verdetti P-1...P-5, regola N2, controprove CP-1 e CP-2, sensibilita' S-1, descrittive (docs/steps/S1.4f.md).
Unita': P-1, P-2 seme (finestre fisse; stima = media dei 10 semi +- 2 ES, ES = sd fra semi / sqrt(10)); P-3, P-4, P-5 ROI (valore del ROI =
media dei 10 semi; 5 ROI); CP-1 e S-1 finestra (media dei 10 semi e delle coppie; 4 finestre). K2 ed ET: P-1 e P-2 dai file di S1.4e (nuclei
identici, cancello G-3). EP = legge vincente della regola O (S1.4f_ruleO.csv). Uscite: results/S1.4f/S1.4f_verdicts.csv, S1.4f_ruleN2.csv,
S1.4f_cp.csv, S1.4f_desc.csv. Uso (root del repo): .venv/bin/python tools/S1.4f_analysis.py"""
import json, numpy as np, pandas as pd
RES = "results/S1.4f"; E = "results/S1.4e"; T6 = ["A1", "A2", "A3", "A4", "A5", "A6"]; DENSE = ["A1", "A2", "A3", "A4"]; P3T = ["A1", "A2", "A3", "A4", "A6"]
CONTACT = {"A1": (0.366, 0.650), "A2": (0.217, 0.567), "A3": (0.083, 0.283), "A4": (0.350, 0.750), "A5": (0.000, 0.217), "A6": (0.000, 0.070)}  # B-127...B-132
def est(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; return (v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v)) if len(v) > 1 else (np.nan, np.nan, len(v))
def verdict(m, se, lo, hi):
    if not np.isfinite(m): return "non valutabile"
    a, b = m - 2 * se, m + 2 * se
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"
def rule_n2(V, arms=("K2", "ET", "EP")):
    """Regola N2: punteggio = confermate in P-1 e P-2; ET o EP solo se (a) >= K2 + 2, (b) nessun tessuto confermato da K2 in P-1 diventa smentito,
    (c) P-3 smentita in al piu' 2 dei 5 tessuti; pari: piu' P-3 confermate; ancora pari: K2."""
    sc = {a: int(((V.arm == a) & V.id.isin(["P-1", "P-2"]) & (V.verdict == "confermata")).sum()) for a in arms}
    p3c = {a: int(((V.arm == a) & (V.id == "P-3") & (V.verdict == "confermata")).sum()) for a in arms}
    p3s = {a: int(((V.arm == a) & (V.id == "P-3") & (V.verdict == "smentita")).sum()) for a in arms}
    k2c = set(V[(V.arm == "K2") & (V.id == "P-1") & (V.verdict == "confermata")].archetype)
    ok_b = lambda a: not any(V[(V.arm == a) & (V.id == "P-1") & V.archetype.isin(k2c)].verdict == "smentita")
    cond = {a: (sc[a] >= sc["K2"] + 2, ok_b(a), p3s[a] <= 2) for a in arms if a != "K2"}
    elig = [a for a in cond if all(cond[a])]
    if not elig: win, why = "K2", "nessun braccio ammesso dalle condizioni (a), (b), (c)"
    elif len(elig) == 1: win, why = elig[0], "unico ammesso"
    elif len({sc[a] for a in elig}) > 1: win = max(elig, key=lambda a: sc[a]); why = "punteggio piu' alto"
    elif len({p3c[a] for a in elig}) > 1: win = max(elig, key=lambda a: p3c[a]); why = "pari: piu' P-3 confermate"
    else: win, why = "K2", "pari anche su P-3"
    T = pd.DataFrame([dict(arm=a, score_P1_P2=sc[a], P3_confirmed=p3c[a], P3_refuted=p3s[a], cond_a=cond.get(a, (True,))[0] if a != "K2" else True,
                           cond_b=cond[a][1] if a != "K2" else True, cond_c=cond[a][2] if a != "K2" else True) for a in arms])
    T["winner"] = win; T["why"] = why
    return win, why, T

def label_ep(SH, win):
    """Etichetta EP = la legge vincente della regola O; la legge resta anche con il suo nome (righe duplicate con arm diverso)."""
    SH = SH.copy(); SH["arm"] = SH.tag.replace({win: "EP"}); return pd.concat([SH, SH[SH.tag == win].assign(arm=win)], ignore_index=True)

def p5_iv(RF, A, c):
    """Intervallo R-4 (o R-2 per nn<5) della classe c del tessuto A."""
    r = RF[(RF.archetype == A) & (RF.cls == c) & RF.ref.isin(["R-2", "R-4"])].iloc[0]; return float(r.lo), float(r.hi)

def cp_rows(PA, bench, cid):
    """Contrasti EP-K2 ed ET-K2 per tessuto denso (coppie giudicate, 'unsure' escluse, d < 6 um; unita' = finestra)."""
    PA = PA[PA.answer != "unsure"]; out = []
    pa = PA[(PA.bench == bench) & (PA.d_um < 6)].groupby(["archetype", "win_id", "arm"]).contact.mean().unstack("arm").reset_index()
    for A in DENSE:
        g = pa[pa.archetype == A]
        for b in ["EP", "ET"]:
            m, se, n = est(g[b] - g["K2"])
            out.append(dict(id=cid, bench=bench, archetype=A, contrast=f"{b}-K2", est=m, se=se, n_units=n, lo2=m - 2 * se, hi2=m + 2 * se,
                            below_zero=bool(m + 2 * se < 0), n_pairs=int(len(PA[(PA.bench == bench) & (PA.archetype == A) & (PA.arm == "K2") & (PA.d_um < 6)]) / 10)))
    return out

def main():
    rows = []
    def add(pid, A, arm, v, lo, hi, unit, src, cls=""):
        m, se, n = est(v); rows.append(dict(id=pid, archetype=A, arm=arm, cls=cls, unit=unit, est=m, se=se, n_units=n, lo=lo, hi=hi, verdict=verdict(m, se, lo, hi), source=src))
    win = pd.read_csv(f"{RES}/S1.4f_ruleO.csv").winner.iloc[0]
    W = pd.concat([pd.read_csv(f"{E}/S1.4e_windows.csv").query("arm in ['K2', 'ET']"), pd.read_csv(f"{RES}/S1.4f_windows.csv")], ignore_index=True)
    PB = pd.concat([pd.read_csv(f"{E}/S1.4e_pairs_B.csv").query("arm in ['K2', 'ET'] and noise == 'none'"), pd.read_csv(f"{RES}/S1.4f_pairs_B.csv")], ignore_index=True)
    SH = label_ep(pd.read_csv(f"{RES}/S1.4f_shape_B.csv"), win)
    R1 = pd.read_csv(f"{E}/S1.4e_refs.csv"); RF = pd.read_csv(f"{RES}/S1.4f_refs.csv")
    iv = lambda D, ref, A, c: tuple(float(x) for x in D[(D.ref == ref) & (D.archetype == A) & (D.cls == c)].iloc[0][["lo", "hi"]])
    assert (W.groupby(["arm", "archetype"]).seed.nunique() == 10).all() and (SH.groupby(["arm", "archetype", "roi_id", "cls"]).seed.nunique() == 10).all()
    ps = W.groupby(["arm", "archetype", "seed"]).s1.mean().reset_index()
    pb = PB.groupby(["arm", "archetype", "seed", "cls"])[["n", "n_contact"]].sum().reset_index(); pb["f"] = pb.n_contact / pb.n
    sr = SH.groupby(["arm", "archetype", "roi_id", "cls"])[["p3", "p4"]].mean().reset_index()
    for arm in ["K2", "ET", "EP"]:
        for A in T6:
            lo, hi = CONTACT[A]; add("P-1", A, arm, ps[(ps.arm == arm) & (ps.archetype == A)].s1, lo, hi, "seme", "B-127...B-132")
            lo, hi = iv(R1, "R-3", A, "all"); add("P-4", A, arm, sr[(sr.arm == arm) & (sr.archetype == A) & (sr.cls == "nn<5")].p4, lo, hi, "ROI", "R-3")
            if A in DENSE:
                lo, hi = iv(R1, "R-1", A, "lt6"); add("P-2", A, arm, pb[(pb.arm == arm) & (pb.archetype == A) & (pb.cls == "lt6")].f, lo, hi, "seme", "R-1")
            if A in P3T:
                lo, hi = iv(RF, "R-2", A, "nn<5"); add("P-3", A, arm, sr[(sr.arm == arm) & (sr.archetype == A) & (sr.cls == "nn<5")].p3, lo, hi, "ROI", "R-2")
    for arm in ["K2", "ET", "EP", "V0", "V1", "V2"]:
        for A in DENSE:
            for c in ["nn5-6", "nn6-8"]:
                lo, hi = p5_iv(RF, A, c); add("P-5", A, arm, sr[(sr.arm == arm) & (sr.archetype == A) & (sr.cls == c)].p3, lo, hi, "ROI", "R-4", cls=c)
    V = pd.DataFrame(rows); V["ep_law"] = win; V.to_csv(f"{RES}/S1.4f_verdicts.csv", index=False)
    s4 = pd.read_csv(f"{E}/S1.4e_verdicts.csv"); chk = []
    for pid in ["P-1", "P-2", "P-3", "P-4"]:
        for arm in ["K2", "ET"]:
            a = V[(V.id == pid) & (V.arm == arm)].set_index("archetype").est; b = s4[(s4.id == pid) & (s4.arm == arm)].set_index("archetype").est.reindex(a.index)
            chk.append(dict(id=pid, arm=arm, max_abs_diff_vs_S14e=float((a - b).abs().max())))
    print(pd.DataFrame(chk).to_string())
    w2, why, RN = rule_n2(V[V.id.isin(["P-1", "P-2", "P-3"])]); RN.to_csv(f"{RES}/S1.4f_ruleN2.csv", index=False); print(RN.to_string())
    # ---- controprove e sensibilita' ----
    PA = pd.read_csv(f"{RES}/S1.4f_pairs_A.csv"); cp = cp_rows(PA, "A", "CP-1") + cp_rows(PA, "A_dd", "S-1"); PA = PA[PA.answer != "unsure"]
    CP = pd.DataFrame(cp)
    a = CP[CP.bench == "A"].set_index(["archetype", "contrast"]); d = CP[CP.bench == "A_dd"].set_index(["archetype", "contrast"])
    S1 = (d.est - a.est).rename("delta").to_frame(); S1["se_A"] = a.se; S1["within_2se"] = S1.delta.abs() < 2 * S1.se_A
    v0 = V[(V.id == "P-5") & (V.arm == "V0")]
    for A in ["A6"]:
        for arm in ["K2", "ET", "EP"]:
            r = V[(V.id == "P-1") & (V.arm == arm) & (V.archetype == A)].iloc[0]; cp.append(dict(id="cuore", archetype=A, contrast=arm, est=r.est, se=r.se, verdict=r.verdict))
    CP = pd.DataFrame(cp); CP.to_csv(f"{RES}/S1.4f_cp.csv", index=False); S1.reset_index().to_csv(f"{RES}/S1.4f_S1.csv", index=False)
    q = CP[CP.id == "CP-1"]; ep = q[q.contrast == "EP-K2"].set_index("archetype"); et = q[q.contrast == "ET-K2"].set_index("archetype")
    n_below = int(ep.below_zero.astype(bool).sum()); n_le = int((ep.est >= et.est.reindex(ep.index)).sum())
    print(f"CP-1: EP-K2 < 0 oltre 2 ES in {n_below}/4 (attesa >= 2); EP-K2 >= ET-K2 in {n_le}/4 (attesa >= 3) ->", "superata" if (n_below >= 2 and n_le >= 3) else "FALLITA")
    print("CP-2: V0 smentita in", int((v0.verdict == "smentita").sum()), "/ 8 caselle di P-5 (attesa >= 4)")
    print("S-1:", int(S1.within_2se.sum()), "/", len(S1), "contrasti entro 2 ES"); print(S1.round(4).to_string())
    # ---- descrittive ----
    D = []
    for A in DENSE:
        g = V[(V.id == "P-1") & (V.archetype == A)].set_index("arm").est; D.append(dict(desc="P-1 quota dell'effetto di ET conservata da EP", archetype=A, est=float((g["EP"] - g["K2"]) / (g["ET"] - g["K2"]))))
        h = SH[(SH.cls == "nn<5") & (SH.archetype == A)].groupby(["arm", "roi_id"]).p3.mean().groupby("arm").mean()
        e5 = pd.read_csv(f"{E}/S1.4e_verdicts.csv"); er = e5[(e5.id == "P-3") & (e5.arm == "ER") & (e5.archetype == A)].est.iloc[0]
        D.append(dict(desc="P-3 posizione di EP fra ER (0) ed ET (1)", archetype=A, est=float((h["EP"] - er) / (h["ET"] - er))))
        for arm in ["EP"]:
            a_ = ps[(ps.arm == arm) & (ps.archetype == A)].set_index("seed").s1; b_ = ps[(ps.arm == "K2") & (ps.archetype == A)].set_index("seed").s1
            m, se, n = est((a_ - b_.reindex(a_.index)).to_numpy()); D.append(dict(desc="P-1 differenza EP - K2 appaiata per seme", archetype=A, est=m, se=se, n=n))
    for A in DENSE:
        for c in ["lt6", "6to10"]:
            m, se, n = est(pb[(pb.arm == "EP") & (pb.archetype == A) & (pb.cls == c)].f); D.append(dict(desc=f"P-2 classe {c}", archetype=A, est=m, se=se, n=n))
    for A in DENSE:
        g = PA[(PA.bench == "A") & (PA.archetype == A) & (PA.d_um < 6)]
        for arm in ["K2", "ET", "EP"]: D.append(dict(desc="banco A < 6 um", archetype=A, arm=arm, est=float(g[g.arm == arm].contact.mean()), n=int(len(g[g.arm == arm]) / 10)))
    pd.DataFrame(D).to_csv(f"{RES}/S1.4f_desc.csv", index=False)
    cnt = V.groupby(["id", "arm"]).verdict.value_counts().unstack(fill_value=0); print(cnt.to_string())
    print(V[V.id.isin(["P-1", "P-2", "P-3"])].pivot_table(index=["id", "archetype"], columns="arm", values="est").round(3).to_string())

if __name__ == "__main__":
    main()
