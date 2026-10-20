#!/usr/bin/env python
"""tools/R2c_posthoc.py — R2c: analisi di sensibilita' POST HOC, scritte dopo i risultati e la revisione avversariale
(2026-10-20). NON cambiano i verdetti pre-registrati: servono a dire quanto sono robusti. Ogni blocco cita il rilievo.
 ci        (F-1, F-4, RV-1)  IC alternativi della S-1 manuale: esatto a grappolo (enumerazione G^G), Wilson sul totale, bootstrap
                             a due stadi (finestre, poi focali); B-dec con ciascuno.
 unsure    (RV-4)            B-dec sulla variante 'non giudicabile = contatto'.
 paired    (F-3, RV-2)       differenza segmentatore - manuale sugli stessi focali appaiati (IC a grappolo, B = 2000).
 lowo      (RV-3)            B-dec lasciando fuori una finestra per volta (IC esatto a grappolo sulle 3 rimaste).
 partners  (RV-7)            S-1 dei segmentatori ristretta ai contatti con i partner di Luca (coppie giudicate).
 k2        (F-2, RV-8)       K-2 con peso uguale per finestra; numero di semi dentro l'IC.
 agreement (RV-5, F-9)       kappa a 3 categorie, accordo positivo/negativo su 'touch', tabelle 2x2 per archetipo.
 fast      (RV-14)           coppie vicine giudicate 'separate' in < 0.8 s (lista per un eventuale riesame).
Uscite: results/R2c/R2c_posthoc_<blocco>.csv.  Uso: .venv/bin/python tools/R2c_posthoc.py"""
import itertools, importlib.util
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import norm
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/R2c"
spec = importlib.util.spec_from_file_location("ana", ROOT / "tools/R2c_analysis.py"); A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
ARCH, A14 = A.ARCH, A.A14
rng0 = 20261020


def exact_cluster(win, y):
    W = list(dict.fromkeys(win)); g = [y[win == w] for w in W]
    ms = [np.nanmean(np.concatenate([g[k] for k in t])) for t in itertools.product(range(len(W)), repeat=len(W))]
    return tuple(np.percentile(ms, [2.5, 97.5]))


def wilson(y):
    y = y[np.isfinite(y)]; n = len(y); p = y.mean(); z = norm.ppf(0.975)
    d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def two_stage(win, y, B=2000, seed=rng0):
    rng = np.random.default_rng(seed); W = list(dict.fromkeys(win)); g = [y[win == w] for w in W]; bs = []
    for _ in range(B):
        pick = rng.choice(len(W), len(W), replace=True)
        bs.append(np.nanmean(np.concatenate([rng.choice(g[k], len(g[k]), replace=True) for k in pick])))
    return tuple(np.percentile(bs, [2.5, 97.5]))


def main():
    focal = pd.read_csv(RES / "R2c_focal.csv"); focal = focal[focal.focal][["win_id", "archetype", "idx"]].reset_index(drop=True)
    sf = pd.read_csv(RES / "R2c_seg_focal.csv"); s1 = pd.read_csv(RES / "R2c_s1.csv").set_index("archetype")
    k2 = pd.read_csv(RES / "R2c_k2_windows.csv"); pairs = pd.read_csv(RES / "R2c_pairs.csv")
    L = A.load_answers("Luca"); Fd = A.load_answers("Federica")
    P = L[L.kind == "pair"][["win_id", "i", "j", "answer", "pair_id", "d_um", "t_ms"]]
    ci, unsure, paired, lowo, partners, k2r = [], [], [], [], [], []
    for a in ARCH:
        f = focal[focal.archetype == a].reset_index(drop=True); w = f.win_id.to_numpy()
        y = A.focal_s1(P[P.win_id.str.startswith(a + "_")], f, "primary"); est = float(np.nanmean(y)); r = s1.loc[a]
        scp, ssd = r.seg_cellpose_rgb, r.seg_stardist_he
        for nm, (lo, hi) in [("cluster percentile (pre-reg)", (r.lo_primary, r.hi_primary)), ("cluster esatto G^G", exact_cluster(w, y)),
                             ("Wilson sul totale", wilson(y)), ("due stadi", two_stage(w, y))]:
            ci.append(dict(archetype=a, ci=nm, est=est, lo=lo, hi=hi, b_dec=A.b_dec(lo, hi, est, scp, ssd),
                           k2_verdict=A.k2_verdict(r.k2_s1, lo, hi, r.widen if np.isfinite(r.widen) else 0.0)))
        unsure.append(dict(archetype=a, est=r.s1_unsure_touch, lo=r.lo_unsure_touch, hi=r.hi_unsure_touch,
                           b_dec=A.b_dec(r.lo_unsure_touch, r.hi_unsure_touch, r.s1_unsure_touch, scp, ssd), b_dec_primary=r.b_dec))
        yy = pd.Series(y, index=pd.MultiIndex.from_frame(f[["win_id", "idx"]]))
        for seg in ["cellpose_rgb", "stardist_he"]:
            m = sf[(sf.archetype == a) & (sf.seg == seg) & sf.matched].merge(f[["win_id", "idx"]], on=["win_id", "idx"])
            ym = yy.reindex(pd.MultiIndex.from_frame(m[["win_id", "idx"]])).to_numpy(); dlt = m.seg_contact.astype(float).to_numpy() - ym
            e, lo, hi = A.cluster_boot(m.win_id.to_numpy(), dlt, seed=rng0)
            paired.append(dict(archetype=a, seg=seg, n=len(m), seg_s1=float(np.nanmean(m.seg_contact.astype(float))),
                               manual_same=float(np.nanmean(ym)), diff=e, lo=lo, hi=hi, compatible=bool(lo <= 0 <= hi)))
            # RV-7: contatto del segmentatore solo con partner di Luca fra le coppie giudicate
            J = pairs[pairs.judged & pairs.win_id.str.startswith(a + "_")]
            con = set()
            for q in J[J[f"contact_{seg}"] == True].itertuples():
                con.add((q.win_id, q.i)); con.add((q.win_id, q.j))
            v = np.array([(x.win_id, x.idx) in con for x in m.itertuples()], float)
            partners.append(dict(archetype=a, seg=seg, n=len(m), seg_s1_any_label=float(m.seg_contact.astype(float).mean()),
                                 seg_s1_luca_partners=float(v.mean()) if len(v) else np.nan, manual_same=float(np.nanmean(ym))))
        for out_w in sorted(set(w)):
            keep = w != out_w; lo, hi = exact_cluster(w[keep], y[keep]); e = float(np.nanmean(y[keep]))
            fk = f[keep]; s_ = {}
            for seg in ["cellpose_rgb", "stardist_he"]:
                m = sf[(sf.archetype == a) & (sf.seg == seg) & sf.matched].merge(fk[["win_id", "idx"]], on=["win_id", "idx"])
                s_[seg] = float(m.seg_contact.astype(float).mean())
            lowo.append(dict(archetype=a, left_out=out_w, est=e, lo=lo, hi=hi, cellpose=s_["cellpose_rgb"], stardist=s_["stardist_he"],
                             b_dec=A.b_dec(lo, hi, e, s_["cellpose_rgb"], s_["stardist_he"])))
        kk = k2[k2.archetype == a]
        agg = kk.groupby("seed")[["n", "n_contact"]].sum(); agg = agg.n_contact / agg.n
        eqw = kk.groupby("seed").s1.mean()
        k2r.append(dict(archetype=a, k2_aggregate=float(np.median(agg)), k2_equal_windows=float(np.median(eqw)),
                        seed_min=float(agg.min()), seed_max=float(agg.max()), lo=r.lo_primary, hi=r.hi_primary, widen=r.widen,
                        seeds_in_ci=int(((agg >= r.lo_primary) & (agg <= r.hi_primary)).sum()),
                        verdict_aggregate=r.k2_verdict, verdict_equal_windows=A.k2_verdict(float(np.median(eqw)), r.lo_primary, r.hi_primary, r.widen)))
    # agreement
    ag = []
    a_ = L[L.kind == "pair"].set_index("pair_id").answer; b_ = Fd[Fd.kind == "pair"].set_index("pair_id").answer
    m = pd.concat([a_.rename("L"), b_.rename("F")], axis=1, join="inner").dropna()
    for scope, mm in [("tutti", m)] + [(x, m[m.index.str.startswith(x + "_")]) for x in ARCH]:
        tL, tF = mm.L == "touch", mm.F == "touch"; both = int((tL & tF).sum()); onlyL = int((tL & ~tF).sum()); onlyF = int((~tL & tF).sum())
        nn = mm[(mm.L != "unsure") & (mm.F != "unsure")]
        ag.append(dict(scope=scope, n=len(mm), kappa_3cat=A.kappa(mm.L, mm.F) if len(mm) else np.nan,
                       kappa_2cat=A.kappa(nn.L, nn.F) if len(nn) else np.nan, touch_both=both, touch_only_Luca=onlyL, touch_only_Federica=onlyF,
                       pos_agreement=2 * both / (2 * both + onlyL + onlyF) if (both + onlyL + onlyF) else np.nan,
                       LsepFtouch=int(((nn.L == "separate") & (nn.F == "touch")).sum()), LtouchFsep=int(((nn.L == "touch") & (nn.F == "separate")).sum())))
    for nm, X in [("Luca", L), ("Federica", Fd)]:
        o = X[X.kind == "pair"].set_index("pair_id").answer; r_ = X[X.kind == "repeat"].set_index("pair_id").answer
        mm = pd.concat([o.rename("a"), r_.rename("b")], axis=1, join="inner").dropna()
        ag.append(dict(scope=f"intra {nm}", n=len(mm), kappa_3cat=A.kappa(mm.a, mm.b)))
    # fast separate near pairs (RV-14)
    dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype").d_max_um
    Q = P.copy(); Q["arch"] = Q.win_id.str[:2]; Q["rel"] = Q.d_um / Q.arch.map(dm)
    q75 = float(Q[Q.answer == "touch"].rel.quantile(0.75))
    fast = Q[(Q.answer == "separate") & (Q.t_ms < 800) & (Q.rel <= q75)].copy(); fast["q75_rel_touch"] = q75
    outs = dict(ci=ci, unsure=unsure, paired=paired, lowo=lowo, partners=partners, k2=k2r, agreement=ag)
    for k, v in outs.items():
        pd.DataFrame(v).to_csv(RES / f"R2c_posthoc_{k}.csv", index=False)
    fast.to_csv(RES / "R2c_posthoc_fast_separate.csv", index=False)
    pd.set_option("display.width", 220)
    for k, v in outs.items():
        print("==", k); print(pd.DataFrame(v).round(3).to_string(index=False))
    print("== fast", len(fast), fast.arch.value_counts().to_dict(), "q75 rel", round(q75, 3))


if __name__ == "__main__":
    main()
