#!/usr/bin/env python
"""tools/S1.4f_test.py - S1.4f: test della taratura, della misura del profilo, delle regole O e N2 (cancello G-1, T7, T9, T11, T12; mutanti G-2).
Uso (root del repo): [S14F_PYSRC=<cartella con tools/>] .venv/bin/python tools/S1.4f_test.py"""
import importlib.util, os
import numpy as np, pandas as pd, shapely
from shapely import affinity
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; SRC = Path(os.environ.get("S14F_PYSRC", ROOT))
def load(name, f):
    sp = importlib.util.spec_from_file_location(name, SRC / f); m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m
CA = load("s14f_calib", "tools/S1.4f_calib.py"); MS = load("s14f_measure", "tools/S1.4f_measure.py"); AN = load("s14f_analysis", "tools/S1.4f_analysis.py")
res = []
def ck(name, ok, val=""):
    res.append(bool(ok)); print(f"{'PASS' if ok else 'FAIL'} {name} {val}")
# ---- T7: taratura ----
t = np.array([0, .2, .4, .6, .8, 1.0]); c = 0.6 + 0.3 * t
ts, cl = CA.solve_t(t, c, 0.6 + 0.3 * 0.53)
ck("T7 inversa esatta su una curva lineare monotona", abs(ts - 0.53) < 1e-12 and not cl, f"t = {ts:.12f}")
ck("T7 fuori dalla curva: estremo piu' vicino, segnalato", CA.solve_t(t, c, 0.5) == (0.0, True) and CA.solve_t(t, c, 0.95) == (1.0, True))
ck("T7 curva non monotona: primo incrocio", abs(CA.solve_t(t, np.array([.6, .7, .8, .7, .75, .9]), 0.65)[0] - 0.1) < 1e-12)
ck("T7 lambda da t: 0 -> 0, 1 -> inf, 1 - e^-2 -> 2", CA.lam_of(0.0) == 0 and np.isinf(CA.lam_of(1.0)) and abs(CA.lam_of(1 - np.exp(-2)) - 2) < 1e-12)
rois = ["r1", "r2", "r3", "r4", "r5"]; slope = {"r1": .2, "r2": .25, "r3": .3, "r4": .35, "r5": .4}
cur = pd.DataFrame([dict(roi_id=r, t=x, p3=0.6 + slope[r] * x) for r in rois for x in t]); real = {r: 0.6 + slope[r] * 0.5 for r in rois}
L1 = pd.DataFrame(CA.loo(cur, real, "A1")); cur2 = cur.copy(); cur2.loc[cur2.roi_id == "r3", "p3"] += 0.2; real2 = dict(real); real2["r3"] += 0.1
L2 = pd.DataFrame(CA.loo(cur2, real2, "A1"))
ck("T7 lascia-fuori-un-ROI: lambda del ROI di prova non dipende dai dati di quel ROI", abs(L1.set_index("roi_id").t["r3"] - L2.set_index("roi_id").t["r3"]) < 1e-12
   and abs(L1.set_index("roi_id").t["r1"] - L2.set_index("roi_id").t["r1"]) > 1e-6 and all("r3" not in s.split(";") for s in L1[L1.roi_id == "r3"].train))
# ---- T9: profilo dell'allineamento per classe di distanza ----
def ell(cx, cy, a, b, rot):
    return affinity.translate(affinity.rotate(affinity.scale(shapely.Point(0, 0).buffer(1, 64), a, b), rot, origin=(0, 0)), cx, cy)
G = [ell(0, 0, 1, 1.6, 0), ell(4, 0, 1, 1.6, 0),            # coppia a 4 um, di fianco (asse lungo y, vicino lungo x): psi = 90
     ell(100, 0, 1.6, 1, 0), ell(105.5, 0, 1.6, 1, 0),        # coppia a 5.5 um, in fila: psi = 0
     ell(200, 0, 1, 1.6, 0), ell(206.5, 0, 1, 1.6, 0)]        # coppia a 6.5 um, di fianco: psi = 90 (D-1: a 7 um il mutante M6 non era rilevato)
S = np.array([MS.ME.poly_shape(g) for g in G]); pr, p4 = MS.profile(S, np.ones(len(G), bool))
ck("T9 profilo: < 5 um 1, 5-6 um 0, 6-8 um 1 (2 nuclei per classe)", pr == {"nn<5": (2, 1.0), "nn5-6": (2, 0.0), "nn6-8": (2, 1.0)}, f"{pr}")
# ---- T11: verdetti e regola N2 su tabelle costruite a mano ----
ck("T11 verdetti: dentro, fuori, a cavallo", [AN.verdict(0.5, 0.01, 0.4, 0.6), AN.verdict(0.8, 0.01, 0.4, 0.6), AN.verdict(0.59, 0.01, 0.4, 0.6)] == ["confermata", "smentita", "indeterminata"])
def tab(conf, p3s):
    r = []
    for arm, n in conf.items():
        for k in range(10): r.append(dict(id="P-1" if k < 6 else "P-2", archetype=f"A{k % 6 + 1}", arm=arm, verdict="confermata" if k < n else "smentita"))
        for k in range(5): r.append(dict(id="P-3", archetype=f"A{k + 1}", arm=arm, verdict="smentita" if k < p3s[arm] else "confermata"))
    return pd.DataFrame(r)
w1 = AN.rule_n2(tab({"K2": 3, "ET": 5, "EP": 4}, {"K2": 0, "ET": 5, "EP": 0}))[0]
w2 = AN.rule_n2(tab({"K2": 3, "ET": 5, "EP": 5}, {"K2": 0, "ET": 5, "EP": 2}))[0]
w3 = AN.rule_n2(tab({"K2": 3, "ET": 6, "EP": 5}, {"K2": 0, "ET": 3, "EP": 3}))[0]
ck("T11 regola N2: +1 resta K2; ET escluso da P-3 smentita 5/5, EP a +2 vince; P-3 smentita in 3 tessuti esclude tutti", (w1, w2, w3) == ("K2", "EP", "K2"), f"{w1} {w2} {w3}")
# ---- T12: regola O su una tabella costruita a mano ----
refs = pd.DataFrame([dict(ref="R-4", archetype=A, cls=c, lo=0.5, hi=0.6, est=0.55) for A in CA.DENSE for c in ["nn5-6", "nn6-8"]])
vals = {"V0": 0.70, "V1": 0.56, "V2": 0.58}
sh = pd.DataFrame([dict(tag=law, archetype=A, roi_id=f"r{k}", seed=f"s{s:02d}", cls=c, p3=v) for law, v in vals.items() for A in CA.DENSE for c in ["nn5-6", "nn6-8"]
                   for k in range(1, 6) for s in range(1, 11)])
w, T, Cc = CA.rule_o(sh, refs)
ck("T12 regola O: V1 e V2 dentro 8/8, vince la deviazione minore (V1); V0 0/8", w == "V1" and T.set_index("law").score.to_dict() == {"V0": 0, "V1": 8, "V2": 8}, f"{T.to_dict('records')}")
sh2 = sh.copy(); sh2.loc[(sh2.tag == "V1") & sh2.seed.isin([f"s{s:02d}" for s in range(6, 11)]), "p3"] = 0.9
ck("T12 regola O usa solo i semi 1-5", CA.rule_o(sh2, refs)[0] == "V1")
# ---- T13-T20: aggiunti dopo la revisione avversariale (mutanti X1-X12 del revisore del codice non rilevati) ----
m, se, n = AN.est([1, 2, 3, 4])
ck("T13 est(): media e errore standard con deviazione standard campionaria (ddof = 1)", abs(m - 2.5) < 1e-12 and abs(se - np.std([1, 2, 3, 4], ddof=1) / 2) < 1e-12 and n == 4)
tb = tab({"K2": 3, "ET": 3, "EP": 6}, {"K2": 0, "ET": 5, "EP": 0}); tb.loc[(tb.arm == "EP") & (tb.id == "P-1") & (tb.archetype == "A1"), "verdict"] = "smentita"
T14 = AN.rule_n2(tb)[2].set_index("arm")
ck("T14 regola N2, condizione (b): EP a K2 + 2 ma smentisce un tessuto confermato da K2 in P-1 -> resta K2", AN.rule_n2(tb)[0] == "K2" and bool(T14.loc["EP", "cond_a"])
   and not bool(T14.loc["EP", "cond_b"]) and bool(T14.loc["EP", "cond_c"]), AN.rule_n2(tb)[1])
sh0 = pd.DataFrame([dict(tag=t_, archetype="A1", roi_id="r1", seed="s01", cls="nn<5", n=1, p3=v, p4=0.5) for t_, v in [("V0", 0.1), ("V1", 0.2), ("V2", 0.3), ("K2", 0.4), ("ET", 0.5)]])
le = AN.label_ep(sh0, "V2"); ep = le[le.arm == "EP"]
ck("T15 etichetta EP = legge vincente, e la legge resta anche con il suo nome", len(ep) == 1 and float(ep.p3.iloc[0]) == 0.3 and len(le[le.arm == "V2"]) == 1 and len(le) == 6)
PA0 = pd.DataFrame([dict(bench="A", archetype=A, win_id=f"{A}_w{w}", arm=arm, seed="s01", d_um=3.0, answer=ans, contact=c_) for A in CA.DENSE for w in (1, 2)
                    for arm, c_ in [("K2", 1), ("ET", 0), ("EP", 1)] for ans in ["touch"]] +
                   [dict(bench="A", archetype=A, win_id=f"{A}_w{w}", arm=arm, seed="s01", d_um=3.0, answer="unsure", contact=c_) for A in CA.DENSE for w in (1, 2)
                    for arm, c_ in [("K2", 0), ("ET", 1), ("EP", 0)]])
cr = pd.DataFrame(AN.cp_rows(PA0, "A", "CP-1")).set_index(["archetype", "contrast"])
ck("T16 CP-1 esclude le coppie 'non giudicabile' (ET - K2 = -1, EP - K2 = 0)", all(cr.loc[(A, "ET-K2")].est == -1 and cr.loc[(A, "EP-K2")].est == 0 for A in CA.DENSE))
G2 = [ell(0, 0, 1, 1.6, 0), ell(4, 0, 1, 1.6, 0), ell(100, 0, 1, 1.6, 0), ell(104, 0, 1, 1.6, 0), ell(200, 0, 1, 1.6, 50), ell(204, 0, 1, 1.6, 50)]
S2 = np.array([MS.ME.poly_shape(g) for g in G2]); pr2, _ = MS.profile(S2, np.array([True, True, False, False, True, True]))
ck("T17 profilo: solo i nuclei interni; psi = 40 gradi non e' di traverso (soglia 45)", pr2["nn<5"] == (4, 0.5), f"{pr2['nn<5']}")
GR = np.array(CA.GRID); kk = np.arange(1, 13); rows_ = []
for law in CA.LAWS:
    for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
        for r in rois:
            for k in kk:
                for s in range(1, 6): rows_.append(dict(tag=f"{law}_k{k:02d}", archetype=A, roi_id=r, seed=f"s{s:02d}", cls="nn<5", n=1, p3=0.6 + 0.3 * (1 - np.exp(-GR[k - 1]))))
e0 = pd.DataFrame([dict(arm="ET", archetype=A, roi_id=r, seed=f"s{s:02d}", p3=0.9) for A in ["A1", "A2", "A3", "A4", "A5", "A6"] for r in rois for s in range(1, 6)])
lam_true = {"A1": 0.5, "A2": 1.0, "A3": 2.0, "A4": 3.0, "A6": 0.75}
real0 = {A: {r: 0.6 + 0.3 * (1 - np.exp(-lam_true.get(A, 1.0))) for r in rois} for A in ["A1", "A2", "A3", "A4", "A5", "A6"]}
_, L0 = CA.calibrate_tables(pd.DataFrame(rows_), e0, real0); L0 = L0[L0.law == "V0"].groupby("archetype")["lambda"].agg(["min", "max"])
ok_l = all(abs(L0.loc[A, "min"] - v) < 1e-9 and abs(L0.loc[A, "max"] - v) < 1e-9 for A, v in lam_true.items())
ck("T18 taratura su tabelle sintetiche: lambda veri ritrovati; cervello = mediana dei 25 lambda (= 1.0; la media sarebbe 1.45)", ok_l and abs(L0.loc["A5", "min"] - 1.0) < 1e-9, f"{L0.round(6).to_dict()}")
gr_r = [float(x) for x in __import__("re").search(r"S14F_GRID <- c\(([^)]*)\)", (SRC / "tools/S1.4f_bench_lib.R").read_text()).group(1).split(",")]
ck("T19 griglia di lambda uguale in R (banco) e in Python (taratura)", gr_r == [float(x) for x in CA.GRID], f"{gr_r}")
RF0 = pd.DataFrame([dict(ref="R-4", archetype="A1", cls="nn5-6", lo=0.1, hi=0.2), dict(ref="R-4", archetype="A1", cls="nn6-8", lo=0.3, hi=0.4), dict(ref="R-2", archetype="A1", cls="nn<5", lo=0.5, hi=0.6)])
ck("T20 intervallo di P-5 per classe giusta", AN.p5_iv(RF0, "A1", "nn6-8") == (0.3, 0.4) and AN.p5_iv(RF0, "A1", "nn5-6") == (0.1, 0.2))
sh3 = sh.copy(); sh3.loc[sh3.tag == "V1", "p3"] = 0.605; sh3.loc[sh3.tag == "V2", "p3"] = 0.58
ck("T21 regola O: una stima appena fuori da R-4 (0.605 > 0.6) non conta", CA.rule_o(sh3, refs)[1].set_index("law").score.to_dict() == {"V0": 0, "V1": 0, "V2": 8})
print(f"RESULT {sum(res)}/{len(res)} PASS")
