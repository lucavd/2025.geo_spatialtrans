#!/usr/bin/env python
"""review/code/03_decision_recompute.py — ricalcolo indipendente di metri, verdetti, punteggi, bootstrap, CP-0/1/2, phi, B-4
e analisi di robustezza (descrittive) della decisione di S1.4c. Codice proprio; nessuna funzione di S1.4c_analysis/decision.
Uso (root del repo): .venv/bin/python results/S1.4c/review/code/03_decision_recompute.py [obs_csv]
"""
import os, sys, json, itertools
import numpy as np, pandas as pd
from scipy.spatial import cKDTree, Delaunay
from scipy.stats import rankdata

S = "/mnt/micron/geo_spatialtrans/S1.4c"; SC = "/mnt/micron/geo_spatialtrans/S1.4c_review/code"; EXP = SC + "/exp"
OUTD = "results/S1.4c/review/code"; TAG = os.environ.get("RV_TAG", "")
SEG = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "stardist_he", "A5": "cellpose_rgb", "A6": "stardist_he"}
ARCH = list(SEG); PANEL7 = ["B1", "B2", "B7", "B8", "S1", "S2", "S3"]; VAR = ["K-c", "K-1", "K-2"]
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
RO = {A: list(rois.roi_id[rois.archetype == A]) for A in ARCH}
out = {}; rec = []
def R(q, rep, val, match): rec.append(dict(quantity=q, reported=str(rep), recomputed=str(val), match=bool(match)))

def ks(a, b):
    a = np.sort(np.asarray(a, float)); b = np.sort(np.asarray(b, float)); v = np.unique(np.concatenate([a, b]))
    return float(np.max(np.abs(np.searchsorted(a, v, "right") / len(a) - np.searchsorted(b, v, "right") / len(b))))
def hits(x, y, r):
    t = cKDTree(np.column_stack([x, y])); h = np.zeros(len(x), bool); p = t.query_pairs(2 * r.max(), output_type="ndarray")
    if len(p):
        i, j = p[:, 0], p[:, 1]; bad = r[i] + r[j] > np.hypot(x[i] - x[j], y[i] - y[j]); h[i[bad]] = True; h[j[bad]] = True
    return h

# ---------------------------------------------------------------- riferimenti reali per ROI (ricalcolati)
prep = pd.read_csv(EXP + "/prep.csv"); real = {}; chk_ev = []
for A in ARCH:
    for roi in RO[A]:
        cc = pd.read_parquet(f"{EXP}/cells_{A}_{roi}.parquet")
        ev = (cc.interior & cc.matched & np.isfinite(cc.a_n) & (cc.a_n > 0)).to_numpy()
        chk_ev.append(bool(np.array_equal(ev, cc.ev_prep.to_numpy())))
        rs = pd.read_parquet(f"{S}/real_shape/{A}_{roi}_{SEG[A]}_shape.parquet")
        z = rs.set_index("label").loc[cc.shape_label.to_numpy()[ev]]
        ra = np.where(cc.matched & np.isfinite(cc.a_n) & (cc.a_n > 0), cc.a_n, cc.a_own)
        b7 = float(np.mean(np.sqrt(cc.a_n.to_numpy()[ev] / np.pi) > cc.d_gen.to_numpy()[ev]))
        b8 = float(np.mean(hits(cc.x.to_numpy(), cc.y.to_numpy(), np.sqrt(ra / np.pi))[ev]))
        real[(A, roi)] = dict(nc=(cc.a_n / cc.a_t).to_numpy()[ev], an=cc.a_n.to_numpy()[ev], b7=b7, b8=b8,
            s1=float(np.mean(z.n_c2 > 0)), s2=(z.n_c2 / z.n_bnd).to_numpy(), s3=z.sol.to_numpy(), rs=rs,
            ev_lab=set(cc.shape_label.to_numpy()[ev].tolist()), n=len(cc))
p2 = prep.set_index(["archetype", "roi"])
db7 = max(abs(real[k]["b7"] - p2.loc[k, "b7_real"]) for k in real); db8 = max(abs(real[k]["b8"] - p2.loc[k, "b8_real"]) for k in real)
R("insieme ev (interior & matched & a_n>0) = prep$ev, 30 ROI", "30/30", f"{sum(chk_ev)}/30", all(chk_ev))
R("b7_real (S1.4b prep) max |diff| su 30 ROI", "prep", f"{db7:.1e}", db7 < 1e-12)
R("b8_real (S1.4b prep) max |diff| su 30 ROI", "prep", f"{db8:.1e}", db8 < 1e-12)
refs = pd.read_csv("results/S1.4c/S1.4c_real_refs.csv").set_index(["archetype", "roi"])
dmax = max(max(abs(real[k]["s1"] - refs.loc[k, "s1"]), abs(np.median(real[k]["s2"]) - refs.loc[k, "s2_median"]),
               abs(np.median(real[k]["s3"]) - refs.loc[k, "s3_median"])) for k in real)
R("riferimenti reali S-1, mediane S-2, S-3 (real_refs.csv), 30 ROI", "real_refs.csv", f"max |diff| {dmax:.1e}", dmax < 1e-12)

# ---------------------------------------------------------------- metri inter-ROI, fold e ROI di test
loo = pd.read_csv(EXP + "/loo.csv"); y14 = pd.read_csv(EXP + "/yard14.csv")
yards = {}; okj = []; dy = 0.0
for A in ARCH:
    for k in range(1, 6):
        kr = loo.roi_test[(loo.archetype == A) & (loo.fold == k)].iloc[0]
        okj.append(kr == RO[A][k - 1])
        js = [j for j in RO[A] if j != kr]; T = real[(A, kr)]
        yy = y14[(y14.archetype == A) & (y14.fold == k)].sort_values("pos")
        okj.append(list(yy.j) == js)
        b1 = np.array([ks(T["nc"], real[(A, j)]["nc"]) for j in js]); b2 = np.array([ks(T["an"], real[(A, j)]["an"]) for j in js])
        dy = max(dy, np.max(np.abs(b1 - yy.ks_nc.to_numpy())), np.max(np.abs(b2 - yy.ks_an.to_numpy())))
        yards[(A, k)] = dict(B1=b1, B2=b2, B7=np.array([abs(T["b7"] - real[(A, j)]["b7"]) for j in js]),
            B8=np.array([abs(T["b8"] - real[(A, j)]["b8"]) for j in js]), S1=np.array([abs(T["s1"] - real[(A, j)]["s1"]) for j in js]),
            S2=np.array([ks(T["s2"], real[(A, j)]["s2"]) for j in js]), S3=np.array([ks(T["s3"], real[(A, j)]["s3"]) for j in js]))
R("fold k -> ROI di test = k-esimo ROI (R2_rois_checked) e yard$j = altri ROI nello stesso ordine (30 fold)", "stopifnot in analysis", f"{sum(okj)}/60", all(okj))
R("metro B-1/B-2 (yard$ks_nc, ks_an di S1.4) ricalcolato dai reali", "S1.4 loo", f"max |diff| {dy:.1e}", dy < 1e-12)
# addestramento: phi, a_ref, n del fit M2 solo dai ROI di addestramento
gen = pd.read_csv(EXP + "/gen.csv"); dphi = 0.0; dref = 0.0; okn = []
for A in ARCH:
    for k in range(1, 6):
        kr = loo.roi_test[(loo.archetype == A) & (loo.fold == k)].iloc[0]; tr = [j for j in RO[A] if j != kr]
        q = np.quantile(np.concatenate([real[(A, j)]["nc"] for j in tr]), 0.99)
        g = gen[(gen.archetype == A) & (gen.fold == k)].iloc[0]; dphi = max(dphi, abs(min(1, q) - g.phi), abs(q - g.q99))
        cells = [pd.read_parquet(f"{EXP}/cells_{A}_{j}.parquet") for j in tr]
        a_ref = sum(c.a_t.sum() for c in cells) / sum(len(c) for c in cells)
        l = loo[(loo.archetype == A) & (loo.fold == k)].iloc[0]; dref = max(dref, abs(a_ref - l.a_ref) / l.a_ref)
        okn.append(int(l.n) == sum(len(real[(A, j)]["an"]) for j in tr))
        # contro-prova di fuga: phi con il ROI di test incluso
        out.setdefault("phi_with_test", []).append(dict(archetype=A, fold=k, phi_train=min(1, q),
            phi_all=min(1, np.quantile(np.concatenate([real[(A, j)]["nc"] for j in RO[A]]), 0.99))))
R("phi = min(1, q99 tipo 7 di a_n/a_t su ev dei 4 ROI di addestramento), 30 fold", "S1.4c_phi.csv", f"max |diff| {dphi:.1e}", dphi < 1e-12)
R("a_ref di M2 = sum a_t / n sui soli ROI di addestramento (S1.4 loo)", "fits$M2['a_ref']", f"max rel diff {dref:.1e}", dref < 1e-9)
R("n del fit M2 = celle ev dei soli ROI di addestramento", "fits$M2['n']", f"{sum(okn)}/30", all(okn))

# ---------------------------------------------------------------- verdetti e punteggi
def fold_v(d, y):
    return "PASS" if d <= np.median(y) else ("WARN" if d <= np.max(y) else "FAIL")
def arch_v(v):
    return "PASS" if v.count("PASS") >= 4 else ("FAIL" if v.count("FAIL") >= 2 else "WARN")
def verdicts(obs, metrics=PANEL7, seeds=None, archs=ARCH, strict=False, yd=None):
    yd = yd or yards; M = {}
    for A in archs:
        for V in VAR:
            for m in metrics:
                fv = []
                for k in range(1, 6):
                    z = obs[(A, V, k)][m]
                    sel = z if seeds is None else z[seeds[(A, k)]]
                    d = np.median(sel); y = yd[(A, k)][m]
                    fv.append(("PASS" if d < np.median(y) else ("WARN" if d <= np.max(y) else "FAIL")) if strict else fold_v(d, y))
                M[(A, V, m)] = arch_v(fv)
    return M
def score(M, metrics, archs=ARCH):
    return {V: sum(M[(A, V, m)] == "PASS" for A in archs for m in metrics) for V in VAR}
def to_obs(df):
    o = {}
    for (A, V, k), z in df.groupby(["archetype", "variant", "fold"]):
        z = z.sort_values("seed_idx"); assert list(z.seed_idx) == list(range(1, 11))
        o[(A, V, int(k))] = {m: z[m].to_numpy() for m in PANEL7}
    return o
obs_csv = sys.argv[1] if len(sys.argv) > 1 else "results/S1.4c/S1.4c_obs_raw.csv"
obs_df = pd.read_csv(obs_csv); OBS = to_obs(obs_df)
M = verdicts(OBS)
vt = pd.read_csv("results/S1.4c/S1.4c_verdicts.csv")
vt_m = sum(M[(r.archetype, r.variant, r.metric)] == r.verdict for r in vt.itertuples() if r.metric in PANEL7)
R("verdetti d'archetipo (6 x 3 varianti x 7 metriche)", "S1.4c_verdicts.csv", f"{vt_m}/126 uguali", vt_m == 126)

# ---------------------------------------------------------------- CP-2
def cp2v(A, roi, seg):
    z = pd.read_parquet(f"{S}/cp2/{A}_{roi}_{seg}_cp2.parquet"); z = z[z.ev]
    return dict(S1=float(np.mean(z.n_c2 > 0)), S2=(z.n_c2 / z.n_bnd).to_numpy(), S3=z.sol.to_numpy())
def dS(m, a, b): return abs(a[m] - b[m]) if m == "S1" else ks(a[m], b[m])
cp2 = {}; demoted = []; demoted_b = []; demoted_c = []; cp2rows = []
for m in ["S1", "S2", "S3"]:
    ex = exb = exc = 0
    for A in ARCH:
        P = SEG[A]; O = "stardist_he" if P == "cellpose_rgb" else "cellpose_rgb"
        vP = [cp2v(A, r, P) for r in RO[A]]; vO = [cp2v(A, r, O) for r in RO[A]]
        dseg = np.array([dS(m, vP[i], vO[i]) for i in range(5)])
        ymax = np.array([max(dS(m, vP[i], vP[j]) for j in range(5) if j != i) for i in range(5)])
        ybench = np.array([yards[(A, i + 1)][m].max() for i in range(5)])
        above = int(np.sum(dseg > ymax)); ex += above >= 3; exb += int(np.sum(dseg > ybench)) >= 3; exc += np.median(dseg) > np.median(ymax)
        cp2rows.append(dict(metric=m, archetype=A, n_above=above, n_above_bench_yard=int(np.sum(dseg > ybench)),
                            dseg_med=float(np.median(dseg)), ymax_med=float(np.median(ymax)), ybench_med=float(np.median(ybench))))
    if ex >= 3: demoted.append(m)
    if exb >= 3: demoted_b.append(m)
    if exc >= 3: demoted_c.append(m)
cp2rep = pd.read_csv("results/S1.4c/S1.4c_cp2.csv")
cmatch = all(int(r.n_above) == int(cp2rep[(cp2rep.metric == r["metric"]) & (cp2rep.archetype == r["archetype"])].n_roi_above.iloc[0]) for _, r in pd.DataFrame(cp2rows).iterrows())
R("CP-2: ROI sopra il massimo del metro, 18 celle metrica x archetipo; metriche declassate", "S1.4c_cp2.csv; S1,S2,S3", f"{'uguali' if cmatch else 'DIVERSI'}; {','.join(demoted)}", cmatch and demoted == ["S1", "S2", "S3"])
pd.DataFrame(cp2rows).to_csv(f"{OUTD}/03_cp2_alternatives{TAG}.csv", index=False)
out["cp2_demoted"] = dict(as_implemented=demoted, bench_yard=demoted_b, median_rule=demoted_c)
panel = [m for m in PANEL7 if m not in demoted]
sc = score(M, panel)
R("punteggi sul pannello {B1,B2,B7,B8}", "K-c 3, K-1 7, K-2 8", ", ".join(f"{v} {sc[v]}" for v in VAR), (sc["K-c"], sc["K-1"], sc["K-2"]) == (3, 7, 8))
sc7 = score(M, PANEL7); out["score_full7"] = sc7

# ---------------------------------------------------------------- CP-1
def npass(V, As, ms): return sum(M[(A, V, m)] == "PASS" for A in As for m in ms)
g14 = npass("K-2", ARCH[:4], PANEL7) - npass("K-c", ARCH[:4], PANEL7); g56 = npass("K-2", ARCH[4:], PANEL7) - npass("K-c", ARCH[4:], PANEL7)
R("CP-1 guadagno K-2 vs K-c (7 metriche): A1-A4 | A5+A6 | smentita (g56 >= g14/2)", "6 | 3 | TRUE", f"{g14} | {g56} | {g56 >= g14 / 2}", (g14, g56) == (6, 3))
g14p = npass("K-2", ARCH[:4], panel) - npass("K-c", ARCH[:4], panel); g56p = npass("K-2", ARCH[4:], panel) - npass("K-c", ARCH[4:], panel)
out["cp1"] = dict(g14=g14, g56=g56, refuted=g56 >= g14 / 2, g14_panel=g14p, g56_panel=g56p, refuted_panel=g56p >= g14p / 2,
                  alt_per_archetype=dict(gain_per_arch_A1A4=g14 / 4, gain_per_arch_A5A6=g56 / 2))

# ---------------------------------------------------------------- bootstrap (RNG proprio) e decisione
def boot(obs, metrics, B, seed, archs=ARCH, nseed=10):
    rng = np.random.default_rng(seed); res = []
    for b in range(B):
        dr = {(A, k): rng.integers(0, nseed, nseed) for A in archs for k in range(1, 6)}
        s = score(verdicts(obs, metrics, seeds=dr, archs=archs), metrics, archs); res.append([s[v] for v in VAR])
    return np.array(res)
def decide(s, bs, margin=3, pmin=0.9):
    win = []
    for i, V in enumerate(VAR[1:], start=1):
        p = float(np.mean(bs[:, i] > bs[:, 0]))
        if s[V] >= s["K-c"] + margin and p >= pmin: win.append((-s[V], i, V))
    return sorted(win)[0][2] if win else "K-c"
bs = boot(OBS, panel, 2000, 7)
pK1 = float(np.mean(bs[:, 1] > bs[:, 0])); pK2 = float(np.mean(bs[:, 2] > bs[:, 0]))
R("P(bootstrap) punteggio > K-c (RNG proprio, B = 2000)", "K-1 1.0, K-2 1.0", f"K-1 {pK1:.4f}, K-2 {pK2:.4f}", pK1 >= 0.9 and pK2 >= 0.9)
choice = decide(sc, bs)
R("scelta della regola (eleggibili tutte, margine 3, p 0.9)", "K-2", choice, choice == "K-2")
out["boot"] = dict(p_K2_gt_K1=float(np.mean(bs[:, 2] > bs[:, 1])), p_K2_eq_K1=float(np.mean(bs[:, 2] == bs[:, 1])),
                   p_K1_gt_K2=float(np.mean(bs[:, 1] > bs[:, 2])), score_quantiles={v: np.percentile(bs[:, i], [2.5, 50, 97.5]).tolist() for i, v in enumerate(VAR)},
                   p_K2_margin3=float(np.mean(bs[:, 2] >= bs[:, 0] + 3)), p_K1_margin3=float(np.mean(bs[:, 1] >= bs[:, 0] + 3)))
# winner per replica (regola sui punteggi della replica, margine 3, senza secondo livello di bootstrap)
wr = []
for row in bs:
    s = dict(zip(VAR, row)); c = [(-s[V], i, V) for i, V in enumerate(VAR[1:], 1) if s[V] >= s["K-c"] + 3]; wr.append(sorted(c)[0][2] if c else "K-c")
out["boot"]["winner_freq"] = pd.Series(wr).value_counts(normalize=True).to_dict()
el_strict = {"K-c": True, "K-1": True, "K-2": False}
out["decision_if_Cperf_blocks"] = "K-1" if (sc["K-1"] >= sc["K-c"] + 3 and pK1 >= 0.9) else "K-c"

# ---------------------------------------------------------------- robustezza (descrittiva)
rob = []
for mg in range(1, 7): rob.append(dict(scenario=f"margine {mg}", **sc, choice=decide(sc, bs, margin=mg)))
for pm in (0.8, 0.95, 0.99): rob.append(dict(scenario=f"p_min {pm}", **sc, choice=decide(sc, bs, pmin=pm)))
for s_out in range(10):
    keep = [i for i in range(10) if i != s_out]
    o9 = {k: {m: v[m][keep] for m in PANEL7} for k, v in OBS.items()}
    M9 = verdicts(o9); s9 = score(M9, panel); b9 = boot(o9, panel, 200, 100 + s_out, nseed=9)
    rob.append(dict(scenario=f"senza seme {s_out + 1}", **s9, choice=decide(s9, b9)))
for A_out in ARCH:
    aa = [A for A in ARCH if A != A_out]; s5 = score(M, panel, aa); b5 = boot(OBS, panel, 200, 200 + ARCH.index(A_out), archs=aa)
    rob.append(dict(scenario=f"senza {A_out}", **s5, choice=decide(s5, b5)))
b7f = boot(OBS, PANEL7, 200, 300); rob.append(dict(scenario="pannello completo 7 metriche (nessun declassamento)", **sc7, choice=decide(sc7, b7f)))
Mst = verdicts(OBS, strict=True); sst = score(Mst, panel); rob.append(dict(scenario="PASS solo se d < mediana (stretto)", **sst, choice=decide(sst, bs)))
rob = pd.DataFrame(rob); rob.to_csv(f"{OUTD}/03_robustness{TAG}.csv", index=False)
def verdicts_meanyard(obs, metrics):
    Mx = {}
    for A in ARCH:
        for V in VAR:
            for m in metrics:
                fv = []
                for k in range(1, 6):
                    d = np.median(obs[(A, V, k)][m]); y = yards[(A, k)][m]
                    fv.append("PASS" if d <= np.mean(y) else ("WARN" if d <= np.max(y) else "FAIL"))
                Mx[(A, V, m)] = arch_v(fv)
    return Mx
smy = score(verdicts_meanyard(OBS, panel), panel); rob = pd.concat([rob, pd.DataFrame([dict(scenario="PASS se d <= media del metro (mutante DM3)", **smy, choice=decide(smy, bs))])])
rng = np.random.default_rng(11); ub = []
for b in range(500):
    s = {}
    for V in VAR:
        dr = {(A, k): rng.integers(0, 10, 10) for A in ARCH for k in range(1, 6)}
        s[V] = score(verdicts(OBS, panel, seeds=dr), panel)[V]
    ub.append([s[v] for v in VAR])
ub = np.array(ub); out["boot_unpaired"] = dict(B=500, p_K1_gt_Kc=float(np.mean(ub[:, 1] > ub[:, 0])), p_K2_gt_Kc=float(np.mean(ub[:, 2] > ub[:, 0])), p_K2_gt_K1=float(np.mean(ub[:, 2] > ub[:, 1])))
rob.to_csv(f"{OUTD}/03_robustness{TAG}.csv", index=False)


# ---------------------------------------------------------------- CP-0 (rapporto, come nel codice; differenza, come nella pre-registrazione)
def auc(s, y):
    y = np.asarray(y, bool); r = rankdata(s); n1 = float(y.sum()); n0 = float((~y).sum())
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
cp0rep = pd.read_csv("results/S1.4c/S1.4c_cp0.csv").set_index("archetype"); cp0 = []
for A in ARCH:
    SR, SD, YY = [], [], []
    for roi in RO[A]:
        Rr = real[(A, roi)]; rs = Rr["rs"]; P = rs[["cx", "cy"]].to_numpy(); lab = rs.label.to_numpy(); ar = rs.area.to_numpy()
        tri = Delaunay(P).simplices
        e = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]]), axis=1), axis=0)
        evl = np.isin(lab, list(Rr["ev_lab"])); e = e[evl[e[:, 0]] | evl[e[:, 1]]]
        d = np.hypot(*(P[e[:, 0]] - P[e[:, 1]]).T); ri = np.sqrt(ar[e[:, 0]] / np.pi); rj = np.sqrt(ar[e[:, 1]] / np.pi)
        pr = pd.read_parquet(f"{S}/real_shape/{A}_{roi}_{SEG[A]}_pairs.parquet")
        M_ = int(lab.max()) + 1; key = set((pr.label_i.astype(np.int64) * M_ + pr.label_j).tolist())
        lo = np.minimum(lab[e[:, 0]], lab[e[:, 1]]).astype(np.int64); hi = np.maximum(lab[e[:, 0]], lab[e[:, 1]]).astype(np.int64)
        y = np.fromiter(((a * M_ + b) in key for a, b in zip(lo, hi)), bool, len(lo))
        SR.append((ri + rj) / d); SD.append(ri + rj - d); YY.append(y)
    sr, sd, yy = map(np.concatenate, (SR, SD, YY))
    cp0.append(dict(archetype=A, n_pairs=len(yy), auc_ratio=auc(sr, yy), auc_diff=auc(sd, yy), reported=cp0rep.loc[A, "auc"], n_rep=int(cp0rep.loc[A, "n_pairs"])))
cp0 = pd.DataFrame(cp0); cp0.to_csv(f"{OUTD}/03_cp0{TAG}.csv", index=False)
dauc = float(np.max(np.abs(cp0.auc_ratio - cp0.reported)))
R("CP-0 AUC (score (r_i + r_j)/d come nel codice), 6 archetipi; n coppie", "0.947-0.997", f"max |diff| {dauc:.1e}; n coppie {'uguali' if (cp0.n_pairs == cp0.n_rep).all() else 'DIVERSE'}", dauc < 1e-9 and (cp0.n_pairs == cp0.n_rep).all())
R("CP-0 AUC con lo score pre-registrato r_i + r_j - d_ij (A1-A4)", "ratio 0.947-0.974", "; ".join(f"{r.archetype} {r.auc_diff:.3f}" for r in cp0.itertuples() if r.archetype in ARCH[:4]), bool((cp0.auc_diff[:4] >= 0.75).all()))

# ---------------------------------------------------------------- B-4
b4 = pd.read_csv(EXP + "/gen_b4.csv"); b4r = pd.read_csv("results/S1.4c/S1.4c_b4.csv")
a = b4.groupby(["archetype", "variant", "fold"]).trunc_all.median().groupby(["archetype", "variant"]).mean()
db4 = max(abs(a[(r.archetype, r.variant)] - r.trunc_all) for r in b4r.itertuples())
k2pass = int(sum(a[(A, "K-2")] <= 0.05 for A in ARCH))
R("B-4 (media sui fold della mediana sui semi), K-2 PASS", "6/6", f"{k2pass}/6, max |diff| {db4:.1e}", k2pass == 6 and db4 < 1e-12)

# ---------------------------------------------------------------- C-0..C-area dai check del banco
gc = pd.read_csv(EXP + "/gen_checks.csv")
c0 = gc[gc.variant == "K-c"].c0; k12 = gc[gc.variant != "K-c"]
R("C-0 identico (fold x semi 0..10)", "330/330", f"{int(c0.sum())}/{len(c0)}", int(c0.sum()) == 330)
R("C3 | C-NO | vuoti (K-1, K-2, 300 configurazioni ciascuna)", "0 | 0 | 0", f"{int(k12.c3_viol.sum() + k12.c3_sf_viol.sum())} | {int(k12.cno_pairs.sum())} | {int(k12.n_empty.sum())}", (k12.c3_viol.sum() + k12.c3_sf_viol.sum() + k12.cno_pairs.sum() + k12.n_empty.sum()) == 0)
R("C-area max", "5.84e-07", f"{gc[gc.variant == 'K-2'].c_area_rel_max.max():.2e}", gc[gc.variant == "K-2"].c_area_rel_max.max() <= 1e-6)

pd.DataFrame(rec).to_csv(f"{OUTD}/03_recomputed{TAG}.csv", index=False)
json.dump(out, open(f"{OUTD}/03_extra{TAG}.json", "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
print(pd.DataFrame(rec).to_string()); print(rob.to_string()); print(json.dumps(out, indent=1, default=str)[:4000]); print(cp0.to_string())
