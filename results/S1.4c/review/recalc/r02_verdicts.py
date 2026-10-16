#!/usr/bin/env python
"""results/S1.4c/review/recalc/r02_verdicts.py — revisione avversariale S1.4c (recalc).
Ricalcolo indipendente (codice nuovo, Python; nessuna funzione di tools/S1.4c_analysis.R o S1.4c_decision.R) di:
riferimenti reali per ROI, metri inter-ROI, verdetti di fold e d'archetipo (K-c, K-1, K-2 x B1,B2,B7,B8,S1,S2,S3),
CP-2 (declassazione), punteggi, bootstrap sui semi, decisione, CP-1, B-4. Confronto con results/S1.4c/*.csv.
Input: dump di r01_dump.R (solo I/O) + parquet grezzi. Uso: .venv/bin/python r02_verdicts.py [obs_csv] [tag]
"""
import sys, json, numpy as np, pandas as pd
from scipy.stats import ks_2samp
B = "/mnt/micron/geo_spatialtrans/S1.4c/"; D = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/dump/"
OUTD = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/"; RES = "results/S1.4c/"
SEG = dict(A1="cellpose_rgb", A2="cellpose_rgb", A3="cellpose_rgb", A4="stardist_he", A5="cellpose_rgb", A6="stardist_he")
ARCHS = list(SEG); VARS = ["K-c", "K-1", "K-2"]; M7 = ["B1", "B2", "B7", "B8", "S1", "S2", "S3"]
obs_csv = sys.argv[1] if len(sys.argv) > 1 else D + "obs_all.csv"; TAG = sys.argv[2] if len(sys.argv) > 2 else "obs"
ks = lambda a, b: float(ks_2samp(np.asarray(a, float), np.asarray(b, float), method="asymp").statistic)
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
ps = pd.read_csv(D + "prep_scalars.csv")
log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)

# ---------------- riferimenti reali per ROI (insieme ev, maschera di forma) ----------------
real = {}
for _, rw in rois.iterrows():
    A, roi = rw.archetype, rw.roi_id
    cc = pd.read_parquet(D + f"cells_{A}_{roi}.parquet")
    ev = (cc.interior.astype(bool) & cc.matched.astype(bool) & np.isfinite(cc.a_n) & (cc.a_n > 0)).to_numpy()
    assert (ev == cc.ev_prep.to_numpy()).all(), (A, roi, "ev diverso da prep")
    rs = pd.read_parquet(B + f"real_shape/{A}_{roi}_{SEG[A]}_shape.parquet").set_index("label")
    lab = cc.shape_label.to_numpy()[ev]
    assert len(np.unique(lab)) == len(lab), (A, roi, "etichette ev duplicate")
    z = rs.loc[lab]
    pr = ps[(ps.archetype == A) & (ps.roi == roi)].iloc[0]
    real[(A, roi)] = dict(s1=float((z.n_c2 > 0).mean()), s2=(z.n_c2 / z.n_bnd).to_numpy(), s3=z.sol.to_numpy(),
                          b7=float(pr.b7_real), b8=float(pr.b8_real), n_ev=int(ev.sum()))
refs = pd.read_csv(RES + "S1.4c_real_refs.csv")
dmax = {}
for k in ["s1", "s2_median", "s3_median", "b7", "b8", "n_ev"]:
    mine = np.array([(real[(a, r)][k] if k in ("s1", "b7", "b8", "n_ev") else np.median(real[(a, r)][k[:2]])) for a, r in zip(refs.archetype, refs.roi)])
    dmax[k] = float(np.max(np.abs(mine - refs[k].to_numpy())))
P("REFS max |mio - riportato|:", json.dumps(dmax))

# ---------------- metri inter-ROI per fold ----------------
yard = pd.read_csv(D + "loo_yard.csv")
yards = {}; fold_roi = {}
for (A, k), z in yard.groupby(["archetype", "fold"]):
    kr = z.roi_test.iloc[0]; fold_roi[(A, k)] = kr
    js = list(z.j); allr = list(rois.roi_id[rois.archetype == A])
    assert sorted(js) == sorted([r for r in allr if r != kr])
    R = real[(A, kr)]; Rj = [real[(A, j)] for j in js]
    yards[(A, k)] = dict(B1=z.ks_nc.to_numpy(), B2=z.ks_an.to_numpy(),
                         B7=np.array([abs(R["b7"] - q["b7"]) for q in Rj]), B8=np.array([abs(R["b8"] - q["b8"]) for q in Rj]),
                         S1=np.array([abs(R["s1"] - q["s1"]) for q in Rj]), S2=np.array([ks(R["s2"], q["s2"]) for q in Rj]),
                         S3=np.array([ks(R["s3"], q["s3"]) for q in Rj]))

# ---------------- distanze per seme ----------------
obs = pd.read_csv(obs_csv)
assert len(obs) == 900 and obs.groupby(["archetype", "fold", "variant"]).size().eq(10).all()
for (A, k), z in obs.groupby(["archetype", "fold"]):
    assert (z.roi == fold_roi[(A, k)]).all(), (A, k)
def fv(d, y):
    return "PASS" if d <= np.median(y) else ("WARN" if d <= np.max(y) else "FAIL")
def av(v):
    v = list(v); return "PASS" if v.count("PASS") >= 4 else ("FAIL" if v.count("FAIL") >= 2 else "WARN")
arr = {}   # (A, V, k, m) -> vettore dei 10 semi ordinati per seed_idx
for (A, V, k), z in obs.groupby(["archetype", "variant", "fold"]):
    z = z.sort_values("seed_idx")
    for m in M7: arr[(A, V, k, m)] = z[m].to_numpy(float)
def verdicts(draw=None):
    out = {}; fold = {}
    for A in ARCHS:
        for V in VARS:
            for m in M7:
                vs = []
                for k in range(1, 6):
                    x = arr[(A, V, k, m)]
                    if draw is not None: x = x[draw[(A, k)]]
                    d = float(np.median(x)); v = fv(d, yards[(A, k)][m]); vs.append(v)
                    fold[(A, V, k, m)] = (d, v)
                out[(A, V, m)] = (av(vs), vs.count("PASS"), vs.count("WARN"), vs.count("FAIL"))
    return out, fold
vt, fold = verdicts()

# confronto verdetti di fold
rf = pd.read_csv(RES + "S1.4c_fold_verdicts.csv"); rf = rf[rf.metric.isin(M7)]
nd = 0; dd = 0.0; dy = 0.0
for _, r in rf.iterrows():
    d, v = fold[(r.archetype, r.variant, int(r.fold), r.metric)]
    y = yards[(r.archetype, int(r.fold))][r.metric]
    dd = max(dd, abs(d - r.dist_median)); dy = max(dy, abs(np.median(y) - r.yard_median), abs(np.max(y) - r.yard_max))
    nd += v != r.verdict
P(f"FOLD verdetti: {len(rf)} righe, discordanti {nd}; max|d_med| {dd:.2e}; max|metro| {dy:.2e}")
rv = pd.read_csv(RES + "S1.4c_verdicts.csv"); rv = rv[rv.metric.isin(M7)]
na = 0
for _, r in rv.iterrows():
    v = vt[(r.archetype, r.variant, r.metric)]
    na += (v[0] != r.verdict) or (v[1:] != (r.n_pass, r.n_warn, r.n_fail))
P(f"ARCHETIPO verdetti: {len(rv)} righe, discordanti {na}")

# ---------------- CP-2 indipendente ----------------
def cp2v(A, roi, seg):
    z = pd.read_parquet(B + f"cp2/{A}_{roi}_{seg}_cp2.parquet"); z = z[z.ev]
    return dict(S1=float((z.n_c2 > 0).mean()), S2=(z.n_c2 / z.n_bnd).to_numpy(), S3=z.sol.to_numpy())
def dS(m, a, b): return abs(a["S1"] - b["S1"]) if m == "S1" else ks(a[m], b[m])
cp2rows = []; demoted = []
for m in ["S1", "S2", "S3"]:
    n_ex = 0
    for A in ARCHS:
        Pm = SEG[A]; O = "stardist_he" if Pm == "cellpose_rgb" else "cellpose_rgb"; rr = list(rois.roi_id[rois.archetype == A])
        vP = [cp2v(A, r, Pm) for r in rr]; vO = [cp2v(A, r, O) for r in rr]
        dseg = np.array([dS(m, vP[i], vO[i]) for i in range(5)])
        ymax = np.array([max(dS(m, vP[i], vP[j]) for j in range(5) if j != i) for i in range(5)])
        above = int((dseg > ymax).sum()); ex = above >= 3; n_ex += ex
        cp2rows.append(dict(metric=m, archetype=A, n_roi_above=above, exceeds=ex, d_seg=";".join(f"{x:.3g}" for x in dseg), yard_max=";".join(f"{x:.3g}" for x in ymax)))
    if n_ex >= 3: demoted.append(m)
    P(f"CP-2 {m}: archetipi che superano {n_ex}/6 -> {'DECLASSATA' if n_ex >= 3 else 'resta'}")
cp2 = pd.DataFrame(cp2rows); rc = pd.read_csv(RES + "S1.4c_cp2.csv")
mm = cp2.merge(rc[["metric", "archetype", "n_roi_above", "exceeds", "d_seg", "yard_max"]], on=["metric", "archetype"], suffixes=("", "_rep"))
P("CP-2 n_roi_above discordanti:", int((mm.n_roi_above != mm.n_roi_above_rep).sum()), "| exceeds discordanti:", int((mm.exceeds != mm.exceeds_rep).sum()),
  "| d_seg stringhe diverse:", int((mm.d_seg != mm.d_seg_rep).sum()), "| yard_max stringhe diverse:", int((mm.yard_max != mm.yard_max_rep).sum()))
# margine dei ROI rispetto alla soglia (quanto e' robusta la declassazione)
panel = [m for m in M7 if m not in demoted]
P("PANNELLO:", panel, "declassate:", demoted)

# ---------------- punteggi, bootstrap, decisione ----------------
def score(vt_, ms): return {V: sum(vt_[(A, V, m)][0] == "PASS" for A in ARCHS for m in ms) for V in VARS}
sc = score(vt, panel); sc7 = score(vt, M7)
P("PUNTEGGI pannello:", sc, "| su 7 metriche:", sc7)
def boot(seed, Bn=200):
    rng = np.random.default_rng(seed); res = []
    for b in range(Bn):
        draw = {(A, k): rng.integers(0, 10, 10) for A in ARCHS for k in range(1, 6)}
        v_, _ = verdicts(draw); res.append(score(v_, panel))
    return pd.DataFrame(res)
bs = boot(777); bs2 = boot(31415)
pb = {V: float((bs[V] > bs["K-c"]).mean()) for V in ["K-1", "K-2"]}; pb2 = {V: float((bs2[V] > bs2["K-c"]).mean()) for V in ["K-1", "K-2"]}
p21 = float((bs["K-2"] > bs["K-1"]).mean()); p21e = float((bs["K-2"] >= bs["K-1"]).mean())
P("BOOT (seme 777) P(V > K-c):", pb, "| (seme 31415):", pb2, "| P(K-2 > K-1):", p21, "P(K-2 >= K-1):", p21e)
P("BOOT media punteggi:", bs.mean().round(2).to_dict(), "| quantili 5-95%:", {V: (float(bs[V].quantile(.05)), float(bs[V].quantile(.95))) for V in VARS})
rb = pd.read_csv(RES + "S1.4c_bootstrap.csv")
P("BOOT riportato media:", rb.mean().round(2).to_dict(), "P(K-1>K-c)", float((rb["K-1"] > rb["K-c"]).mean()), "P(K-2>K-c)", float((rb["K-2"] > rb["K-c"]).mean()), "P(K-2>K-1)", float((rb["K-2"] > rb["K-1"]).mean()))
win = [V for V in ["K-1", "K-2"] if sc[V] >= sc["K-c"] + 3 and pb[V] >= 0.9]
choice = "K-c" if not win else sorted(win, key=lambda V: (-sc[V], ["K-c", "K-1", "K-2"].index(V)))[0]
P("DECISIONE (eleggibilita' assunta dai check C):", choice, "vincenti:", win)

# ---------------- CP-1 ----------------
def npass(V, As, ms): return sum(vt[(A, V, m)][0] == "PASS" for A in As for m in ms)
g14 = npass("K-2", ARCHS[:4], M7) - npass("K-c", ARCHS[:4], M7); g56 = npass("K-2", ARCHS[4:], M7) - npass("K-c", ARCHS[4:], M7)
g14p = npass("K-2", ARCHS[:4], panel) - npass("K-c", ARCHS[:4], panel); g56p = npass("K-2", ARCHS[4:], panel) - npass("K-c", ARCHS[4:], panel)
P(f"CP-1: guadagno A1-A4 {g14} (per coppia {g14/2}), A5+A6 {g56} -> smentita {g56 >= g14/2} | pannello: {g14p}, {g56p}")

# ---------------- B-4 ----------------
b4 = pd.read_csv(D + "gen_b4.csv")
assert b4.groupby(["archetype", "fold", "variant"]).size().eq(10).all()
b4a = b4.groupby(["archetype", "variant", "fold"]).trunc_all.median().groupby(["archetype", "variant"]).mean()
rb4 = pd.read_csv(RES + "S1.4c_b4.csv").set_index(["archetype", "variant"])
d4 = float(np.max(np.abs(b4a.loc[rb4.index].to_numpy() - rb4.trunc_all.to_numpy())))
ver = np.where(b4a.loc[rb4.index] <= 0.05, "PASS", "FAIL")
P(f"B-4: max|mio - riportato| {d4:.2e}; verdetti discordanti {int((ver != rb4.verdict.to_numpy()).sum())}; K-2 PASS {int((b4a.xs('K-2', level='variant') <= .05).sum())}/6")
P("B-4 K-2:", b4a.xs("K-2", level="variant").round(4).to_dict(), "| K-c:", b4a.xs("K-c", level="variant").round(4).to_dict())

# ---------------- salvataggi ----------------
rows = [dict(archetype=A, variant=V, metric=m, verdict=v[0], n_pass=v[1], n_warn=v[2], n_fail=v[3], in_panel=m in panel) for (A, V, m), v in vt.items()]
pd.DataFrame(rows).to_csv(OUTD + f"recalc_verdicts_{TAG}.csv", index=False)
pd.DataFrame([dict(archetype=A, variant=V, fold=k, metric=m, dist_median=d, verdict=v, yard_median=float(np.median(yards[(A, k)][m])), yard_max=float(np.max(yards[(A, k)][m])))
              for (A, V, k, m), (d, v) in fold.items()]).to_csv(OUTD + f"recalc_fold_verdicts_{TAG}.csv", index=False)
cp2.to_csv(OUTD + f"recalc_cp2.csv", index=False); bs.to_csv(OUTD + f"recalc_bootstrap_{TAG}.csv", index=False)
json.dump(dict(score=sc, score7=sc7, pboot=pb, pboot2=pb2, p21=p21, choice=choice, panel=panel, demoted=demoted, cp1=[g14, g56, g14p, g56p], b4_maxdiff=d4, refs=dmax),
          open(OUTD + f"recalc_summary_{TAG}.json", "w"), indent=1)
open(OUTD + f"r02_{TAG}.log", "w").write("\n".join(log) + "\n")
