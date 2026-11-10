"""tools/S2.6_analysis.py - S2.6: tabelle ed esiti di P-1...P-8 e CP-1, CP-2; misure J-2, J-3. Definizioni: docs/steps/S2.6.md.
Uso (root del repo): .venv/bin/python tools/S2.6_analysis.py [outdir]   -> <outdir>/S2.6_*.csv, S2.6_verdicts.csv, S2.6_verdict_counts.csv"""
import sys, os, json, glob, numpy as np, pandas as pd
M = "/mnt/micron/geo_spatialtrans/S2.6"; OUT = sys.argv[1] if len(sys.argv) > 1 else "results/S2.6"; os.makedirs(OUT, exist_ok=True)
A6 = ["A1", "A2", "A3", "A4", "A5", "A6"]; V = []
def verdict(est, se2, lo, hi):
    a, b = est - se2, est + se2
    if not np.isfinite(est): return "non valutabile"
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"
def verdict_iv(a, b, lo, hi):
    if not (np.isfinite(a) and np.isfinite(b)): return "non valutabile"
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"
def add(pid, A, what, est, se2, lo, hi, note="", iv=None):
    a, b = (est - se2, est + se2) if iv is None else iv     # per le quantita' derivate delle sezioni: intervallo percentile del bootstrap
    V.append(dict(id=pid, archetype=A, quantity=what, estimate=est, two_se=se2, low=a, high=b, rel_lo=lo, rel_hi=hi,
                  outcome=verdict_iv(a, b, lo, hi), note=note))
J = lambda f: json.load(open(f))
# ---------------- domanda 1
R = []
for A in A6:
    for side in (1000, 2000, 3000):
        for k in (0, 1):
            o = glob.glob(f"{M}/q1/one_{A}_{side}_k{k}/run.json")
            if not o: continue
            j = J(o[0]); lj = J(glob.glob(f"{M}/layers/{j['layer'].replace('.rds', '.json')}")[0])
            R.append(dict(archetype=A, side_um=side, mm2=(side / 1000) ** 2, k=k, step1_s=lj["elapsed_s"], step1_peak_gb=lj["peak_mb"] / 1024,
                          n_cells=lj["n_cells"], sim_s=j["elapsed_s"], sim_peak_gb=j["peak_mb"] / 1024, n_generated=j["n_generated"], M=j["M"],
                          **{f"stage_{s}": v for s, v in j["elapsed"].items()}, step1_tess_s=lj["rprof_total_s"].get("tessellate_voronoi", np.nan)))
Q1 = pd.DataFrame(R); Q1.to_csv(f"{OUT}/S2.6_q1_runs.csv", index=False)
def seedpair(A, f):
    vals = []
    for k in (0, 1):
        d = {s: Q1[(Q1.archetype == A) & (Q1.side_um == s) & (Q1.k == k)] for s in (1000, 2000, 3000)}
        if any(len(x) == 0 for x in d.values()): continue
        vals.append(f(d))
    v = np.array(vals); return (v.mean(), 2 * v.std(ddof=1) / np.sqrt(len(v))) if len(v) >= 2 else (np.nan, np.nan)
for A in A6:
    if not len(Q1[Q1.archetype == A]): continue
    e, s = seedpair(A, lambda d: (d[3000].step1_s.iloc[0] / 9) / d[1000].step1_s.iloc[0]); add("P-1", A, "Step 1, tempo per mm2 9/1", e, s, 0.8, 1.25)
    e, s = seedpair(A, lambda d: (d[3000].sim_s.iloc[0] / 9) / d[1000].sim_s.iloc[0]); add("P-2a", A, "simulazione, tempo per mm2 9/1", e, s, 0.8, 1.25)
    e, s = seedpair(A, lambda d: ((d[3000].sim_peak_gb.iloc[0] - d[2000].sim_peak_gb.iloc[0]) / 5) / ((d[2000].sim_peak_gb.iloc[0] - d[1000].sim_peak_gb.iloc[0]) / 3))
    add("P-2b", A, "memoria aggiunta per mm2, 4-9 / 1-4", e, s, 0.8, 1.25)
    e, s = seedpair(A, lambda d: d[3000].sim_s.iloc[0] / 9); add("P-4b", A, "un pezzo solo, s per mm2 utile (9 mm2)", e, s, 0, 300)
D = []
for A in A6:
    for f in sorted(glob.glob(f"{M}/q1/design_{A}/p*/run.json")):
        j = J(f); cut = J(f"{M}/q1/design_{A}/cut.json"); p = int(os.path.basename(os.path.dirname(f))[1:])
        ci = [c for c in cut["pieces"] if c["p"] == p][0]
        D.append(dict(archetype=A, piece=p, useful_mm2=ci["useful_mm2"], box_mm2=ci["box_mm2"], n_cells=ci["n_cells"], sim_s=j["elapsed_s"], peak_gb=j["peak_mb"] / 1024,
                      s_per_useful_mm2=j["elapsed_s"] / ci["useful_mm2"], n_generated=j["n_generated"]))
DS = pd.DataFrame(D); DS.to_csv(f"{OUT}/S2.6_q1_design_pieces.csv", index=False)
for A in A6:
    x = DS[DS.archetype == A]
    if len(x) < 2: continue
    add("P-3", A, "memoria di un pezzo del disegno (GB)", x.peak_gb.mean(), 2 * x.peak_gb.std(ddof=1) / np.sqrt(len(x)), 0, 8)
    add("P-4a", A, "disegno, s per mm2 utile", x.s_per_useful_mm2.mean(), 2 * x.s_per_useful_mm2.std(ddof=1) / np.sqrt(len(x)), 0, 300)
# ---------------- domanda 2
MS, PR = [], []
for f in glob.glob(f"{M}/q2/*/measures_T*.csv"):
    t = os.path.basename(os.path.dirname(f)); arm, A, k = t.rsplit("_", 2)[0], t.rsplit("_", 2)[1], int(t.rsplit("_", 2)[2][1:])
    m = pd.read_csv(f); m["arm"], m["archetype"], m["k"] = arm, A, k; MS.append(m)
    p = pd.read_csv(f.replace("measures_", "pairs_")); p["arm"], p["archetype"], p["k"] = arm, A, k; PR.append(p)
if MS:
    MS = pd.concat(MS); PR = pd.concat(PR); MS.to_csv(f"{OUT}/S2.6_q2_measures.csv", index=False); PR.to_csv(f"{OUT}/S2.6_q2_pairs.csv", index=False)
    J1 = []
    for A in A6:
        for arm, Tu in [("P1000", 1000), ("P504", 504), ("P1000H0", 1000)]:
            for b in ["0-8", "8-16", "16-40", "40-100", "100-250", "all"]:
                p = MS[(MS.archetype == A) & (MS.arm == arm) & (MS.T_um == Tu) & (MS.band == b)].mean_umi.values
                o = MS[(MS.archetype == A) & (MS.arm == "ONE") & (MS.T_um == Tu) & (MS.band == b)].mean_umi.values
                if len(p) < 2 or len(o) < 2: continue
                r = p.mean() / o.mean(); se = r * np.sqrt(p.var(ddof=1) / len(p) / p.mean() ** 2 + o.var(ddof=1) / len(o) / o.mean() ** 2)
                dp = MS[(MS.archetype == A) & (MS.arm == arm) & (MS.T_um == Tu) & (MS.band == b)].dup_frac.mean()
                J1.append(dict(archetype=A, arm=arm, T_um=Tu, band=b, ratio=r, two_se=2 * se, n_seeds_pieces=len(p), n_seeds_one=len(o), dup_frac=dp))
                if b != "all": add("CP-1" if arm == "P1000H0" else "P-5", A, f"J-1 {arm} fascia {b} um", r, 2 * se, 0.98, 1.02) if (arm != "P1000H0" or b == "0-8") else None
    pd.DataFrame(J1).to_csv(f"{OUT}/S2.6_q2_J1_J2.csv", index=False)
    j3 = PR.groupby(["archetype", "arm", "T_um"]).agg(n_pairs=("n_pairs", "mean"), pearson=("pearson_umi", "mean"), pearson_sd=("pearson_umi", "std"),
                                                     cosine=("mean_cosine", "mean"), cosine_sd=("mean_cosine", "std"), n_seeds=("k", "count")).reset_index()
    j3.to_csv(f"{OUT}/S2.6_q2_J3.csv", index=False)
# ---------------- domanda 3
SL = []
for f in glob.glob(f"{M}/sections/sec_*.json"):
    j = J(f); SL.append(dict(archetype=j["archetype"], n_cells=j["n_cells"], n_regions=j["n_regions"], tissue_mm2=j["tissue_mm2_layer"], step1_s=j["elapsed_s"],
                             step1_peak_gb=j["peak_mb"] / 1024, s_per_mm2=j["elapsed_s"] / j["tissue_mm2_layer"], tess_s=j["rprof_total_s"].get("tessellate_voronoi")))
SL = pd.DataFrame(SL)
if len(SL):
    SL.to_csv(f"{OUT}/S2.6_q3_step1.csv", index=False)
    for _, r in SL.iterrows():
        add("P-8a", r.archetype, "Step 1 sezione, s per mm2 di tessuto", r.s_per_mm2, 0, 0, 300, "un seme")
        add("P-8b", r.archetype, "Step 1 sezione, picco (GB)", r.step1_peak_gb, 0, 0, 8, "un seme")
SR = []
for A in ["A1", "A2", "A4", "A5", "A6"]:
    row = dict(archetype=A)
    for arm in ("one", "pieces"):
        f = f"{M}/q3/{arm}_{A}/run.json"
        if not os.path.exists(f): continue
        j = J(f)
        if arm == "one": row.update(one_s=j["elapsed_s"], one_peak_gb=j["peak_mb"] / 1024, one_generated=j["n_generated"])
        else:
            P = [p for p in j["pieces"] if not p.get("empty")]
            row.update(pieces_n=j["n_pieces"], pieces_empty=j["n_pieces"] - len(P), pieces_wall_s=j["elapsed_wall_s"], pieces_process_s=j["elapsed_process_s"],
                       piece_peak_max_gb=max(p["peak_mb"] for p in P) / 1024, piece_peak_increment_max_gb=max(p["peak_mb"] - p["rss_start_mb"] for p in P) / 1024)
    SR.append(row)
pd.DataFrame(SR).to_csv(f"{OUT}/S2.6_q3_runs.csv", index=False)
for A in ["A1", "A4", "A6"]:
    f = f"{OUT}/S2.6_section_{A}.csv"
    if not os.path.exists(f): continue
    s = pd.read_csv(f)
    for _, r in s.iterrows():
        if "ratio_pieces_one" not in r or not np.isfinite(r.get("ratio_pieces_one", np.nan)): 
            if r.band != "tissue": add("CP-2" if A == "A6" else "P-6", A, f"vetro {r.band}", np.nan, np.nan, 0.8, 1.25, "nessun blocco o media nulla")
            continue
        est = r.ratio_pieces_one; lo, hi = r.lo_pieces_one, r.hi_pieces_one; se2 = max(est - lo, hi - est)
        nt = f"bootstrap percentile delle tessere ({int(r.n_tiles)} tessere, {int(r.n_blocks)} blocchi)"
        if r.band == "tissue": add("P-7", A, "tessuto della sezione, a pezzi / pezzo solo", est, se2, 0.98, 1.02, nt, iv=(lo, hi))
        else: add("CP-2" if A == "A6" else "P-6", A, f"vetro {r.band}", est, se2, 0.8, 1.25, nt, iv=(lo, hi))
VV = pd.DataFrame(V); VV.to_csv(f"{OUT}/S2.6_verdicts.csv", index=False)
cnt = VV.groupby(["id", "outcome"]).size().unstack(fill_value=0); cnt.to_csv(f"{OUT}/S2.6_verdict_counts.csv")
print(cnt.to_string())
