#!/usr/bin/env python
"""tools/S1.4f_calib.py - S1.4f: taratura lascia-fuori-un-ROI del lambda delle leggi V0, V1, V2 e regola O (docs/steps/S1.4f.md).
calibrate: curva = quota P-3 sotto i 5 um (media dei semi 1-5) in funzione di t = 1 - exp(-lambda) sulla griglia S14F_GRID, piu' t = 1 (ET di
  S1.4e, stessi semi); per il ROI di prova r: curva media dei 4 altri ROI, bersaglio = media dei valori veri R-2 degli stessi 4 ROI; t_r per
  interpolazione lineare al primo incrocio; fuori dalla curva l'estremo (t = 0 o t = 1) piu' vicino al bersaglio. Cervello: mediana dei 25
  lambda_r degli altri cinque tessuti. Uscite: results/S1.4f/S1.4f_lambda.csv, S1.4f_cal_curves.csv.
ruleO: quota P-3 a 5-6 e 6-8 um nei 4 tessuti densi, semi 1-5, lambda_r; ROI = media dei semi, stima = media dei 5 ROI; punteggio = caselle
  dentro R-4; pari: somma |stima - media vera| piu' piccola. Uscita: results/S1.4f/S1.4f_ruleO.csv (colonna winner).
Uso (root del repo): .venv/bin/python tools/S1.4f_calib.py calibrate|ruleO"""
import sys, numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/S1.4f"
GRID = [0, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10]; LAWS = ["V0", "V1", "V2"]; T6 = ["A1", "A2", "A3", "A4", "A5", "A6"]
CAL_T = ["A1", "A2", "A3", "A4", "A6"]; DENSE = ["A1", "A2", "A3", "A4"]; CAL_SEEDS = [f"s{k:02d}" for k in range(1, 6)]

def roi_values(refs, ref, cls):
    """Valori veri per ROI (colonna roi_values di S1.4f_refs.csv) -> {archetipo: {roi: valore}}."""
    out = {}
    for r in refs[(refs.ref == ref) & (refs.cls == cls)].itertuples():
        out[r.archetype] = {k: float(v) for k, v in (x.split(":") for x in r.roi_values.split(";"))}
    return out

def solve_t(t, c, target):
    """Primo incrocio della curva lineare a tratti c(t) con il bersaglio; fuori dalla curva, l'estremo piu' vicino. -> (t*, clamped)."""
    t = np.asarray(t, float); c = np.asarray(c, float); s = c - target
    for j in range(len(t) - 1):
        if s[j] == 0: return float(t[j]), False
        if s[j] * s[j + 1] < 0: return float(t[j] + (t[j + 1] - t[j]) * (target - c[j]) / (c[j + 1] - c[j])), False
    if s[-1] == 0: return float(t[-1]), False
    return (float(t[0]) if abs(s[0]) <= abs(s[-1]) else float(t[-1])), True

def lam_of(t): return np.inf if t >= 1 else float(-np.log(1 - t))

def loo(curves, real, A):
    """curves: DataFrame (roi_id, t, p3) di un tessuto e una legge; real: {roi: valore vero}. -> righe con lambda_r per ROI di prova."""
    rois = sorted(real); rows = []
    for r in rois:
        train = [x for x in rois if x != r]
        cv = curves[curves.roi_id.isin(train)].groupby("t").p3.mean().sort_index(); target = float(np.mean([real[x] for x in train]))
        ts, cl = solve_t(cv.index.to_numpy(), cv.to_numpy(), target)
        rows.append(dict(archetype=A, roi_id=r, t=ts, **{"lambda": lam_of(ts)}, clamped=cl, target=target, train=";".join(train)))
    return rows

def calibrate():
    refs = pd.read_csv(RES / "S1.4f_refs.csv"); real = roi_values(refs, "R-2", "nn<5")
    sh = pd.read_csv(RES / "S1.4f_shape_cal.csv"); sh = sh[(sh.cls == "nn<5") & sh.seed.isin(CAL_SEEDS)].copy()
    sh["law"] = sh.tag.str.split("_").str[0]; sh["k"] = sh.tag.str.split("_k").str[1].astype(int); sh["t"] = 1 - np.exp(-np.array(GRID)[sh.k - 1])
    e = pd.read_csv(ROOT / "results/S1.4e/S1.4e_shape_B.csv"); e = e[(e.arm == "ET") & e.seed.isin(CAL_SEEDS)]
    assert sh.groupby(["law", "archetype", "roi_id", "k"]).seed.nunique().eq(5).all() and e.groupby(["archetype", "roi_id"]).seed.nunique().eq(5).all()
    C = []
    for law in LAWS:
        a = sh[sh.law == law].groupby(["archetype", "roi_id", "t"]).p3.mean().reset_index()
        b = e.groupby(["archetype", "roi_id"]).p3.mean().reset_index(); b["t"] = 1.0
        cc = pd.concat([a, b], ignore_index=True); cc["law"] = law; C.append(cc)
    C = pd.concat(C, ignore_index=True); C.to_csv(RES / "S1.4f_cal_curves.csv", index=False); rows = []
    for law in LAWS:
        lr = []
        for A in CAL_T: lr += [dict(law=law, **x) for x in loo(C[(C.law == law) & (C.archetype == A)], real[A], A)]
        med = float(np.median([x["lambda"] for x in lr]))
        for r in sorted(C[C.archetype == "A5"].roi_id.unique()):
            lr.append(dict(law=law, archetype="A5", roi_id=r, t=1 - np.exp(-med), **{"lambda": med}, clamped=False, target=np.nan, train="median of other tissues"))
        rows += lr
    L = pd.DataFrame(rows); L.to_csv(RES / "S1.4f_lambda.csv", index=False)
    print(L.pivot_table(index=["archetype", "roi_id"], columns="law", values="lambda").round(3).to_string()); print("clamped:", int(L.clamped.sum()))

def rule_o(sh, refs):
    """sh: righe di S1.4f_shape_B.csv; restituisce (vincente, tabella per legge, tabella per casella)."""
    r4 = refs[refs.ref == "R-4"].set_index(["archetype", "cls"]); cells = []
    for law in LAWS:
        g = sh[(sh.tag == law) & sh.seed.isin(CAL_SEEDS) & sh.archetype.isin(DENSE) & sh.cls.isin(["nn5-6", "nn6-8"])]
        v = g.groupby(["archetype", "cls", "roi_id"]).p3.mean().groupby(["archetype", "cls"]).mean()
        for (A, c), est in v.items():
            lo, hi, mu = r4.loc[(A, c), ["lo", "hi", "est"]]
            cells.append(dict(law=law, archetype=A, cls=c, est=float(est), lo=float(lo), hi=float(hi), real=float(mu), inside=bool(lo <= est <= hi), absdev=float(abs(est - mu))))
    Cc = pd.DataFrame(cells); assert len(Cc) == 3 * 8
    T = Cc.groupby("law").agg(score=("inside", "sum"), sum_absdev=("absdev", "sum")).reset_index()
    T = T.sort_values(["score", "sum_absdev"], ascending=[False, True]).reset_index(drop=True); T["winner"] = T.law.iloc[0]
    return T.law.iloc[0], T, Cc

def ruleO():
    sh = pd.read_csv(RES / "S1.4f_shape_B.csv"); refs = pd.read_csv(RES / "S1.4f_refs.csv")
    w, T, Cc = rule_o(sh, refs); T.to_csv(RES / "S1.4f_ruleO.csv", index=False); Cc.to_csv(RES / "S1.4f_ruleO_cells.csv", index=False)
    print(Cc.round(3).to_string()); print(T.to_string()); print("legge vincente:", w)

if __name__ == "__main__":
    {"calibrate": calibrate, "ruleO": ruleO}[sys.argv[1]]()
