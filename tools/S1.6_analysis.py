"""tools/S1.6_analysis.py - S1.6: previsione P-1, controprove CP-1 e CP-2, descrittive D-1...D-3 (definizioni: docs/steps/S1.6.md).
Valori per ROI (per seme nei contatti) da tools/S1.6_lib.py, stesse definizioni di S1.5e. P-1: differenza appaiata Q6 (S1.6) meno Q (S1.5e) con
gli intervalli di S1.5e; margine +- un quarto della larghezza dell'intervallo vero (scelta di Luca, 2026-10-29). CP-1: stesso confronto fra XA e RV
di S1.5e. CP-2: quota di proposte accettate dal diradamento, Q6 su Q, per tessuto. D-2: verdetti P-1...P-17 di U6 e Q6 con i riferimenti per tipo
di tools/S1.6_real_refs.py. Uscite: results/S1.6/S1.6_P1.csv, S1.6_CP1.csv, S1.6_CP2.csv, S1.6_verdicts.csv, S1.6_counts.csv, S1.6_desc.csv."""
import importlib.util, numpy as np, pandas as pd
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
L = mod("s16lib", "tools/S1.6_lib.py"); RES = "results/S1.6"
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
Rold = L.load_refs(); Rnew = L.load_refs("results/S1.6/S1.6_real_q90.csv", "results/S1.6/S1.6_real_refs_roi_type.csv")
B6 = L.load_bench(RES); B5 = L.load_bench("results/S1.5e")
assert (B6["U"].groupby(["arm", "archetype", "roi_id"]).seed.count() == 10).all() and set(B6["U"].arm) == {"U6", "Q6"}
# ---- P-1 ----
E3 = pd.read_csv(f"{RES}/S1.6_explore_se.csv"); E3 = E3[(E3.arm == "Q") & (E3.ratio_2se_margin > 1)][["id", "archetype", "unit"]].assign(by_construction=True)
P1 = L.paired(L.roi_values(B6, Rold, "Q6"), L.roi_values(B5, Rold, "Q")).merge(E3, on=["id", "archetype", "unit"], how="left")
P1["by_construction"] = P1.by_construction.fillna(False).astype(bool)
P1.to_csv(f"{RES}/S1.6_P1.csv", index=False)
print("P-1 righe", len(P1), "|", P1.verdict.value_counts().to_dict(), "| indeterminate fra le 12 di E3:", int((P1.by_construction & (P1.verdict == "indeterminata")).sum()),
      "| indeterminate fuori dalle 12:", int((~P1.by_construction & (P1.verdict == "indeterminata")).sum()))
print(P1[P1.verdict != "confermata"][["id", "archetype", "unit", "n_pairs", "diff", "se", "margin", "verdict", "by_construction"]].round(4).to_string(index=False))
# ---- CP-1 ----
C1 = L.paired(L.roi_values(B5, Rold, "XA"), L.roi_values(B5, Rold, "RV"))
need = pd.DataFrame([("P-15", "A1"), ("P-15", "A4"), ("P-15", "A5"), ("P-17", "A1"), ("P-3", "A4")], columns=["id", "archetype"])
c1 = need.merge(C1[C1.unit == "tessuto"], on=["id", "archetype"], how="left")
C1.to_csv(f"{RES}/S1.6_CP1.csv", index=False)
print("CP-1 attese smentite:"); print(c1[["id", "archetype", "diff", "se", "margin", "verdict"]].round(4).to_string(index=False))
print("CP-1", "superata" if (c1.verdict == "smentita").all() else "fallita", "| smentite in tutto il confronto:", int((C1.verdict == "smentita").sum()), "su", len(C1))
# ---- CP-2 ----
U5 = B5["U"]; U6 = B6["U"]
acc = lambda d: d.groupby("archetype").apply(lambda g: g.n_thin_accept.sum() / g.n_attempts.sum(), include_groups=False)
C2 = pd.DataFrame(dict(acc_Q6=acc(U6[U6.arm == "Q6"]), acc_Q=acc(U5[U5.arm == "Q"]))); C2["ratio"] = C2.acc_Q6 / C2.acc_Q
B4 = pd.read_csv(f"{RES}/S1.6_explore_bound.csv"); B4w = B4.groupby("archetype").apply(lambda g: (g.ratio * g.area_um2).sum() / g.area_um2.sum(), include_groups=False)
C2["bound_ratio_E4"] = B4w; C2["expected"] = np.where(C2.index == "A5", "< 0.95", "> 1.05")
C2["ok"] = np.where(C2.index == "A5", C2.ratio < 0.95, C2.ratio > 1.05)
C2.reset_index().to_csv(f"{RES}/S1.6_CP2.csv", index=False); print("CP-2"); print(C2.round(4).to_string()); print("CP-2", "superata" if C2.ok.all() else "fallita")
# ---- D-2: verdetti di U6 e Q6 con i riferimenti nuovi; S1.5e U e Q con i riferimenti vecchi e nuovi ----
V = pd.concat([L.verdicts(L.roi_values(B6, Rnew, "U6"), "U6"), L.verdicts(L.roi_values(B6, Rnew, "Q6"), "Q6"),
               L.verdicts(L.roi_values(B5, Rnew, "U"), "U (S1.5e, riferimenti nuovi)"), L.verdicts(L.roi_values(B5, Rnew, "Q"), "Q (S1.5e, riferimenti nuovi)"),
               L.verdicts(L.roi_values(B5, Rold, "U"), "U (S1.5e)"), L.verdicts(L.roi_values(B5, Rold, "Q"), "Q (S1.5e)")], ignore_index=True)
V.to_csv(f"{RES}/S1.6_verdicts.csv", index=False)
CNT = V.groupby(["arm", "verdict"]).size().unstack(fill_value=0); CNT.to_csv(f"{RES}/S1.6_counts.csv"); print("D-2 conteggi"); print(CNT.to_string())
SPATIAL = ["P-1", "P-2", "P-3", "P-4", "P-5", "P-6", "P-7"]; TYPES = ["P-4", "P-5", "P-6"]
for a in ["U6", "Q6", "Q (S1.5e, riferimenti nuovi)", "Q (S1.5e)"]:
    v = V[V.arm == a]; print(f"D-2 {a}: confermate spaziali {int(((v.verdict == 'confermata') & v.id.isin(SPATIAL)).sum())}, per tipo {int(((v.verdict == 'confermata') & v.id.isin(TYPES)).sum())}")
u = V[V.arm == "U6"].merge(V[V.arm == "U (S1.5e, riferimenti nuovi)"], on=["id", "archetype", "unit"], suffixes=("", "_s15e"))
print("D-2 U6 contro U di S1.5e con gli stessi riferimenti: verdetti diversi", int((u.verdict != u.verdict_s15e).sum()), "max |stima|", float(np.nanmax(np.abs(u.est - u.est_s15e))))
q = V[V.arm == "Q6"].merge(V[V.arm == "Q (S1.5e, riferimenti nuovi)"], on=["id", "archetype", "unit"], suffixes=("", "_s15e"))
dq = q[q.verdict != q.verdict_s15e]; print("D-2 Q6 contro Q di S1.5e con gli stessi riferimenti: verdetti diversi", len(dq))
print(dq[["id", "archetype", "unit", "est_s15e", "est", "lo", "hi", "verdict_s15e", "verdict"]].round(4).to_string(index=False))
# ---- D-1, D-3 ----
D = []
for A, g in U6.groupby("archetype"):
    for arm6, arm5 in [("U6", "U"), ("Q6", "Q")]:
        D.append(dict(desc="D-1 tempo mediano per unita' (s)", archetype=A, arm=arm6, value=float(g[g.arm == arm6].elapsed_s.median()),
                      ref=float(U5[(U5.arm == arm5) & (U5.archetype == A)].elapsed_s.median())))
T6 = B6["T"]; T5 = B5["T"]
def frac(T, arm):
    t = T[T.arm == arm].copy(); t["f"] = t.n_type / t.groupby(["archetype", "roi_id", "seed"]).n_type.transform("sum")
    return t.groupby(["archetype", "cell_type"]).f.mean()
F = pd.DataFrame(dict(Q6=frac(T6, "Q6"), Q=frac(T5, "Q"), RV=frac(T5, "RV")))
for (A, t), r in F.iterrows():
    D.append(dict(desc="D-3 frazione realizzata del tipo (media delle unita')", archetype=A, arm=t, value=r.Q6, ref=r.RV, ref2=r.Q))
D = pd.DataFrame(D); D.to_csv(f"{RES}/S1.6_desc.csv", index=False)
print("D-1"); print(D[D.desc.str.startswith("D-1")].round(2).to_string(index=False))
f3 = F.assign(d_Q6_RV=F.Q6 - F.RV, d_Q_RV=F.Q - F.RV); print("D-3 scarto massimo |Q6 - RV|", round(float(f3.d_Q6_RV.abs().max()), 4), "|Q - RV|", round(float(f3.d_Q_RV.abs().max()), 4), "|Q6 - Q|", round(float((F.Q6 - F.Q).abs().max()), 4))
