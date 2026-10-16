#!/usr/bin/env python
"""tools/R2c_analysis.py — R2c: analisi delle risposte manuali (pre-registrazione + addendum D-1..D-4).

Funzioni pure (testate da tools/R2c_analysis_test.py prima di qualsiasi risposta, BL-125) + pipeline `main()`.
Mutanti (variabile R2C_MUTANT, C-4): MU1 bootstrap senza reinserimento; MU2 bootstrap non a grappolo (ricampiona i nuclei);
MU3 kappa con le marginali scambiate (p_e dalle marginali del primo annotatore per entrambi).
Ingressi: results/R2c/annotations/R2c_answers_<rater>.json (export dell'annotatore), R2c_items_<rater>.csv, R2c_pairs.csv,
R2c_focal.csv, R2c_seg_focal.csv, R2c_k2_windows.csv.  Uscite: results/R2c/R2c_*.csv (vedi main).
Uso: .venv/bin/python tools/R2c_analysis.py [--partial]
"""
import os, sys, json
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import beta

MUT = os.environ.get("R2C_MUTANT", "none")
ROOT = Path(__file__).resolve().parents[1]; RES = Path(os.environ.get("R2C_RES", ROOT / "results/R2c"))
B_BOOT, SEED_BOOT = 2000, 20261016
A14 = ["A1", "A2", "A3", "A4"]; ARCH = A14 + ["A5", "A6"]


# ----------------------------------------------------------------------------------------------- funzioni pure
def kappa(a, b):
    """Kappa di Cohen per due vettori di etichette (stesse categorie)."""
    a = np.asarray(a); b = np.asarray(b); assert len(a) == len(b) and len(a) > 0
    cats = np.union1d(a, b)
    po = float(np.mean(a == b))
    pa = np.array([np.mean(a == k) for k in cats]); pb = np.array([np.mean(b == k) for k in cats])
    pe = float(np.sum(pa * pa)) if MUT == "MU3" else float(np.sum(pa * pb))
    return np.nan if pe == 1 else (po - pe) / (1 - pe)


def focal_s1(pair_ans, focal, variant="primary"):
    """S-1 per nucleo focale. pair_ans: DataFrame win_id, i, j, answer in {touch, separate, unsure} (una riga per coppia).
    focal: DataFrame win_id, idx. variant: primary (non giudicabili escluse; nucleo NaN se nessun 'touch' e solo 'unsure'),
    unsure_touch, unsure_sep. Un focale senza coppie candidate vale 0."""
    lk = {}
    for r in pair_ans.itertuples():
        lk.setdefault((r.win_id, int(r.i)), []).append(r.answer); lk.setdefault((r.win_id, int(r.j)), []).append(r.answer)
    out = []
    for r in focal.itertuples():
        a = lk.get((r.win_id, int(r.idx)), [])
        if "touch" in a:
            y = 1.0
        elif variant == "unsure_touch":
            y = 1.0 if "unsure" in a else 0.0
        elif variant == "unsure_sep":
            y = 0.0
        else:
            y = np.nan if (len(a) and all(x == "unsure" for x in a)) else 0.0
        out.append(y)
    return np.array(out, float)


def cluster_boot(win, y, B=B_BOOT, seed=SEED_BOOT):
    """IC 95 % percentile della media (nanmean) con ricampionamento delle finestre con reinserimento."""
    win = np.asarray(win); y = np.asarray(y, float); rng = np.random.default_rng(seed)
    W = np.unique(win); groups = [y[win == w] for w in W]
    est = float(np.nanmean(y)); bs = np.empty(B)
    for b in range(B):
        if MUT == "MU2":
            yy = rng.choice(y, size=len(y), replace=True)
        else:
            pick = rng.choice(len(W), size=len(W), replace=(MUT != "MU1"))
            yy = np.concatenate([groups[k] for k in pick])
        bs[b] = np.nanmean(yy)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return est, float(lo), float(hi)


def cp_upper(k, n, alpha=0.05):
    """Limite superiore dell'IC (1 - alpha) bilaterale di Clopper-Pearson."""
    return 1.0 if k >= n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))


def auc(score, y):
    """AUC (Mann-Whitney) di score per y = 1; pareggi a 1/2."""
    s = np.asarray(score, float); y = np.asarray(y).astype(bool)
    p, q = s[y], s[~y]
    if not len(p) or not len(q):
        return np.nan
    return float(((p[:, None] > q[None, :]).sum() + 0.5 * (p[:, None] == q[None, :]).sum()) / (len(p) * len(q)))


def b_dec(lo, hi, est, s_cp, s_sd, tie=0.02):
    """Regola B-dec pre-registrata. Restituisce cellpose / StarDist / nessuno / indeciso."""
    in_cp = lo <= s_cp <= hi; in_sd = lo <= s_sd <= hi
    if in_cp and in_sd:
        dcp, dsd = abs(s_cp - est), abs(s_sd - est)
        return "indeciso" if abs(dcp - dsd) < tie else ("cellpose" if dcp < dsd else "StarDist")
    return "cellpose" if in_cp else "StarDist" if in_sd else "nessuno"


def band_verdict(k, n):
    """D-3: PASS se cp_upper <= 0.05; WARN se k == 1 e limite > 0.05; FAIL se k >= 2.
    D-6 (lacuna della regola trovata dai test prima dei dati): k == 0 con limite > 0.05 (n < 72) = 'non valutabile'."""
    if n == 0:
        return "NA"
    if cp_upper(k, n) <= 0.05:
        return "PASS"
    if k == 0:
        return "non valutabile"
    return "WARN" if k == 1 else "FAIL"


def k2_verdict(s_k2, lo, hi, widen):
    if lo <= s_k2 <= hi:
        return "PASS"
    return "WARN" if (lo - widen) <= s_k2 <= (hi + widen) else "FAIL"


def kappa_verdict(k, pass_=0.60, warn=0.40):
    return "NA" if np.isnan(k) else "PASS" if k >= pass_ else "WARN" if k >= warn else "FAIL"


# ----------------------------------------------------------------------------------------------- pipeline
def load_answers(rater):
    f = RES / f"annotations/R2c_answers_{rater}.json"
    if not f.exists():
        return None
    d = json.load(open(f)); it = pd.read_csv(RES / f"R2c_items_{rater}.csv")
    assert d["doc_id"] == f"R2c_{rater}_v10" and d["n_items"] == len(it), (d["doc_id"], d["n_items"], len(it))
    a = pd.DataFrame(d["answers"]); assert (a["item"].to_numpy() == it["item"].to_numpy()).all()
    return it.merge(a, on="item", validate="1:1")


def s1_table(ans, focal, seg_focal, k2=None, label=""):
    """Per archetipo: S-1 manuale (tre varianti) con IC a grappolo, S-1 dei segmentatori sui focali appaiati, manuale sui
    soli focali appaiati, B-dec, K-2."""
    P = ans[ans.kind == "pair"][["win_id", "i", "j", "answer"]]
    rows = []
    for A in ARCH:
        f = focal[focal.archetype == A]
        if not len(f):
            continue
        r = dict(set=label, archetype=A, n_focal=len(f))
        for v in ["primary", "unsure_touch", "unsure_sep"]:
            y = focal_s1(P[P.win_id.str.startswith(A + "_")], f, v)
            est, lo, hi = cluster_boot(f.win_id.to_numpy(), y)
            r.update({f"s1_{v}": est, f"lo_{v}": lo, f"hi_{v}": hi, f"n_{v}": int(np.isfinite(y).sum())})
            if v == "primary":
                ycache = pd.Series(y, index=pd.MultiIndex.from_frame(f[["win_id", "idx"]]))
        for seg in ["cellpose_rgb", "stardist_he", "spaceranger"]:
            sf = seg_focal[(seg_focal.archetype == A) & (seg_focal.seg == seg) & seg_focal.matched]
            sf = sf.merge(f[["win_id", "idx"]], on=["win_id", "idx"])
            if not len(sf):
                continue
            r[f"seg_{seg}"] = float(sf.seg_contact.astype(float).mean()); r[f"n_{seg}"] = len(sf)
            ym = ycache.reindex(pd.MultiIndex.from_frame(sf[["win_id", "idx"]])).to_numpy()
            r[f"manual_on_{seg}"] = float(np.nanmean(ym))
        r["b_dec"] = b_dec(r["lo_primary"], r["hi_primary"], r["s1_primary"], r.get("seg_cellpose_rgb", np.nan),
                           r.get("seg_stardist_he", np.nan))
        if k2 is not None:
            kk = k2[(k2.archetype == A) & k2.win_id.isin(f.win_id.unique())].groupby("seed")[["n", "n_contact"]].sum()
            r["k2_s1"] = float(np.median(kk.n_contact / kk.n))
        rows.append(r)
    return pd.DataFrame(rows)


def main(partial=False):
    focal = pd.read_csv(RES / "R2c_focal.csv"); focal = focal[focal.focal][["win_id", "archetype", "idx"]].reset_index(drop=True)
    pairs = pd.read_csv(RES / "R2c_pairs.csv"); seg_focal = pd.read_csv(RES / "R2c_seg_focal.csv")
    k2 = pd.read_csv(RES / "R2c_k2_windows.csv")
    L = load_answers("Luca"); Fd = load_answers("Federica")
    assert L is not None, "manca results/R2c/annotations/R2c_answers_Luca.json"
    for nm, X in [("Luca", L), ("Federica", Fd)]:
        if X is not None and not partial:
            assert X.answer.notna().all(), f"{nm}: {int(X.answer.isna().sum())} item senza risposta"
    out = {}
    # --- tempi e completezza
    tm = []
    for nm, X in [("Luca", L), ("Federica", Fd)]:
        if X is None:
            continue
        for k, g in X.groupby("kind"):
            tm.append(dict(rater=nm, kind=k, n=len(g), n_done=int(g.answer.notna().sum()), t_med_s=g.t_ms.median() / 1000,
                           touch=int((g.answer == "touch").sum()), separate=int((g.answer == "separate").sum()),
                           unsure=int((g.answer == "unsure").sum())))
    out["answers_summary"] = pd.DataFrame(tm)
    # --- CP-2 intra, CP-3 trappole, D-3 fascia
    cp = []
    for nm, X in [("Luca", L), ("Federica", Fd)]:
        if X is None:
            continue
        o = X[X.kind == "pair"].set_index("pair_id").answer; r_ = X[X.kind == "repeat"].set_index("pair_id").answer
        m = pd.concat([o.rename("a"), r_.rename("b")], axis=1, join="inner").dropna()
        m = m[(m.a != "unsure") & (m.b != "unsure")]
        k = kappa(m.a, m.b) if len(m) else np.nan
        cp.append(dict(check="CP-2", rater=nm, n=len(m), value=k, agree=float((m.a == m.b).mean()) if len(m) else np.nan,
                       verdict="NA" if np.isnan(k) else ("PASS" if k >= 0.70 else "FAIL")))
        c_ = X[X.kind == "catch"].dropna(subset=["answer"])
        fsep = float((c_.answer == "separate").mean()) if len(c_) else np.nan
        cp.append(dict(check="CP-3", rater=nm, n=len(c_), value=fsep, verdict="PASS" if fsep >= 0.95 else "FAIL"))
        b_ = X[(X.kind == "band") & X.answer.isin(["touch", "separate"])]
        kk, nn = int((b_.answer == "touch").sum()), len(b_)
        cp.append(dict(check="D-3", rater=nm, n=nn, value=kk, upper=cp_upper(kk, nn) if nn else np.nan,
                       verdict=band_verdict(kk, nn) if nm == "Luca" else "descrittivo"))
    # --- CP-1 inter
    if Fd is not None:
        a = L[L.kind == "pair"].set_index("pair_id").answer; b = Fd[Fd.kind == "pair"].set_index("pair_id").answer
        m = pd.concat([a.rename("L"), b.rename("F")], axis=1, join="inner").dropna()
        n_all = len(m); m = m[(m.L != "unsure") & (m.F != "unsure")]
        k = kappa(m.L, m.F) if len(m) else np.nan
        cp.append(dict(check="CP-1", rater="Luca-Federica", n=len(m), n_with_unsure=n_all, value=k,
                       agree=float((m.L == m.F).mean()) if len(m) else np.nan, verdict=kappa_verdict(k)))
        for A in ARCH:
            mm = m[m.index.str.startswith(A + "_")]
            if len(mm):
                cp.append(dict(check="CP-1 per archetipo", rater=A, n=len(mm), value=kappa(mm.L, mm.F),
                               agree=float((mm.L == mm.F).mean()), verdict="descrittivo"))
    # --- CP-4 AUC
    P = L[L.kind == "pair"].dropna(subset=["answer"]); P = P[P.answer != "unsure"]
    for A in ARCH:
        g = P[P.win_id.str.startswith(A + "_")]
        v = auc(-g.d_um.to_numpy(), (g.answer == "touch").to_numpy())
        cp.append(dict(check="CP-4", rater=A, n=len(g), n_touch=int((g.answer == "touch").sum()), value=v,
                       verdict=("PASS" if v >= 0.80 else "FAIL") if A in A14 and np.isfinite(v) else "descrittivo"))
    out["counterproofs"] = pd.DataFrame(cp)
    # --- S-1, B-dec, K-2
    s1L = s1_table(L, focal, seg_focal, k2, "Luca")
    if Fd is not None:
        fw = Fd.win_id.unique(); f2 = focal[focal.win_id.isin(fw)]
        s1F = s1_table(Fd, f2, seg_focal, None, "Federica")
        s1Ls = s1_table(L, f2, seg_focal, None, "Luca_sottoinsieme_Federica")
        widen = (s1Ls.set_index("archetype").s1_primary - s1F.set_index("archetype").s1_primary).abs()
        s1L["widen"] = s1L.archetype.map(widen)
        out["s1_federica"] = pd.concat([s1F, s1Ls])
        cp5 = s1L.set_index("archetype").b_dec.reindex(A14) == s1F.set_index("archetype").b_dec.reindex(A14)
        out["cp5"] = pd.DataFrame(dict(archetype=A14, b_dec_Luca=s1L.set_index("archetype").b_dec.reindex(A14).to_numpy(),
                                       b_dec_Federica=s1F.set_index("archetype").b_dec.reindex(A14).to_numpy(), same=cp5.to_numpy()))
    else:
        s1L["widen"] = 0.0
    s1L["k2_verdict"] = [k2_verdict(r.k2_s1, r.lo_primary, r.hi_primary, r.widen if np.isfinite(r.widen) else 0.0)
                         for r in s1L.itertuples()]
    out["s1"] = s1L
    # --- P-3 sensibilita' / specificita' per coppia
    pp = P.merge(pairs[["pair_id", "both_cellpose_rgb", "both_stardist_he", "contact_cellpose_rgb", "contact_stardist_he"]],
                 on="pair_id", how="left", validate="1:1")
    pp = pp[(pp.both_cellpose_rgb == True) & (pp.both_stardist_he == True)]
    ss = []
    for A in ARCH:
        g = pp[pp.win_id.str.startswith(A + "_")]; t = g.answer == "touch"
        for seg in ["cellpose_rgb", "stardist_he"]:
            c_ = g[f"contact_{seg}"] == True
            ss.append(dict(archetype=A, seg=seg, n=len(g), n_touch=int(t.sum()),
                           sens=float(c_[t].mean()) if t.any() else np.nan, spec=float((~c_[~t]).mean()) if (~t).any() else np.nan))
    out["sens_spec"] = pd.DataFrame(ss)
    for k, v in out.items():
        v.to_csv(RES / f"R2c_{k}.csv", index=False)
    return out


if __name__ == "__main__":
    o = main(partial="--partial" in sys.argv)
    pd.set_option("display.width", 220)
    for k, v in o.items():
        print("==", k); print(v.round(3).to_string(index=False))
