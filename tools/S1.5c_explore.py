"""tools/S1.5c_explore.py - S1.5c: descrittive di esplorazione PRIMA della pre-registrazione (nessun verdetto).
Confronta i centri dei nuclei veri (etichette R6 PH-1, territori reali esportati da R6) con la semina di S1.5b
(braccio ricentrato, seme 1, cells_seed1.parquet) sulle stesse 30 ROI. Serve a scegliere la famiglia di modelli di semina.
Uso (root del repo): .venv/bin/python tools/S1.5c_explore.py    Uscita: results/S1.5c/S1.5c_explore.csv"""
import os, importlib.util, numpy as np, pandas as pd
from scipy.spatial import cKDTree
spec = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(spec); spec.loader.exec_module(cal)
RES = "results/S1.5c"; os.makedirs(RES, exist_ok=True)
S = pd.read_parquet("/mnt/micron/geo_spatialtrans/S1.5b/cells_seed1.parquet"); S = S[S.arm == "recentred"]
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
BINS = [0, 10, 20, 40, 80, 160]; RLOC = [15, 30, 60]
rng = np.random.default_rng(20261023)
rows = []
def add(A, r, src, metric, value, cell_type="", bin_=""):
    rows.append(dict(archetype=A, roi_id=r, source=src, metric=metric, cell_type=cell_type, bin=bin_, value=float(value)))
def describe(A, r, src, x, y, an, at, interior, typ):
    P = np.c_[x, y]; tr = cKDTree(P); n = len(P)
    dd, ii = tr.query(P, k=7); nn = dd[:, 1]; j = ii[:, 1]
    rad = np.sqrt(an / np.pi); gap = nn / (rad + rad[j])
    I = interior
    for q in (0.05, 0.25, 0.5): add(A, r, src, f"nn_q{int(q*100):02d}", np.quantile(nn[I], q))
    for q in (0.05, 0.5): add(A, r, src, f"gapratio_q{int(q*100):02d}", np.quantile(gap[I], q))
    add(A, r, src, "gapratio_lt1", np.mean(gap[I] < 1)); add(A, r, src, "gapratio_lt1.2", np.mean(gap[I] < 1.2))
    add(A, r, src, "cv_terr", at[I].std(ddof=1) / at[I].mean())
    for R in RLOC:
        cnt = np.array([len(v) - 1 for v in tr.query_ball_point(P[I], R)]); add(A, r, src, "cv_localcount", cnt.std(ddof=1) / cnt.mean(), bin_=str(R))
    typed = ~pd.isna(typ) & ~pd.Series(typ).isin(list(cal.UNTYPED)).values
    tt = np.asarray(typ, dtype=object)
    for t in MAIN[A]:
        f_t = np.mean(tt[typed] == t)
        foc = np.where(I & (tt == t))[0]
        if len(foc) < 20: continue
        add(A, r, src, "terr_median", np.median(at[foc]), t); add(A, r, src, "nn_median_type", np.median(nn[foc]), t)
        add(A, r, src, "nt_median_type", np.median(an[foc] / at[foc]), t)
        nb = ii[foc, 1:]; ok = typed[nb]; same = (tt[nb] == t) & ok
        add(A, r, src, "same6_ratio", same.sum() / ok.sum() / f_t, t)
        sub = foc if len(foc) <= 1500 else rng.choice(foc, 1500, replace=False)
        nbl = tr.query_ball_point(P[sub], BINS[-1])
        num = np.zeros(len(BINS) - 1); den = np.zeros(len(BINS) - 1)
        for k, lst in zip(sub, nbl):
            lst = np.array([m for m in lst if m != k and typed[m]], dtype=int)
            if not len(lst): continue
            d = np.hypot(*(P[lst] - P[k]).T); b = np.digitize(d, BINS) - 1; s = tt[lst] == t
            np.add.at(num, b, s); np.add.at(den, b, 1)
        for kb in range(len(BINS) - 1):
            if den[kb] > 0: add(A, r, src, "markconn_ratio", num[kb] / den[kb] / f_t, t, f"{BINS[kb]}-{BINS[kb+1]}")
for A in cal.ARCH:
    lab, d = cal.load(A)
    for r, e in d.groupby("roi_id"):
        describe(A, r, "real", e.x.values, e.y.values, e.a_n.values, e.a_t.values, e.interior.values, e.type.values)
    for r, e in S[S.archetype == A].groupby("roi_id"):
        describe(A, r, "sim_S1.5b", e.x.values, e.y.values, e.a_n.values, e.a_t.values, e.interior.values, e.cell_type.values)
    print(A, "done", flush=True)
pd.DataFrame(rows).to_csv(f"{RES}/S1.5c_explore.csv", index=False)
