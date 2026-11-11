#!/usr/bin/env python3
"""tools/S2.8_analysis.py - S2.8: previsioni P-1...P-3, regola R-1, misure e descrittiva di BL-237 dai json di
/mnt/micron/geo_spatialtrans/S2.8/ (docs/steps/S2.8.md). Uso (root del repo): .venv/bin/python tools/S2.8_analysis.py [results/S2.8]
Tempo = tempo di orologio dell'esecuzione (elapsed_total_s); picco = massimo fra processo principale e pezzi (MB / 1024 = GB, come S2.7).
Stima: media dei 3 semi; 2 errori standard = 2 sd / sqrt(3). Per il picco la stima e' il massimo sui semi, con lo stesso errore standard.
Esito: stima +- 2 es tutta dentro [0, soglia] confermata, tutta fuori smentita, altrimenti indeterminata."""
import json, os, sys, math
import numpy as np, pandas as pd
O = "/mnt/micron/geo_spatialtrans/S2.8"; RES = sys.argv[1] if len(sys.argv) > 1 else "results/S2.8"; os.makedirs(RES, exist_ok=True)
TAGS = ["A1_r1", "A2_r1", "A3_r1", "A4_f1", "A5_r1", "A6_r1"]
TIS = {"A1": "intestino", "A2": "tumore del colon", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
CONF = ["ref", "e504", "e1008", "einf"]
rows = []
for cf in CONF:
    for t in TAGS:
        for k in (0, 1, 2):
            f = f"{O}/{cf}/{t}_k{k}.json"
            if not os.path.exists(f): continue
            j = json.load(open(f))
            pk = j.get("peak_mb") if cf == "ref" else max(j["main_peak_mb"], j["piece_peak_max_mb"])
            rows.append(dict(conf=cf, tag=t, A=t[:2], k=k, time_s=j["elapsed_total_s"], peak_gb=pk / 1024,
                             main_peak_gb=(j.get("main_peak_mb") or pk) / 1024, piece_peak_gb=(j.get("piece_peak_max_mb") or float("nan")) / 1024,
                             n_pieces=j.get("n_pieces", 1), load_end=j.get("load_end", "")))
R = pd.DataFrame(rows); R.to_csv(f"{RES}/S2.8_runs.csv", index=False)
def summ(g):
    n = len(g); es_t = g.time_s.std(ddof=1) / math.sqrt(n) if n > 1 else float("nan"); es_p = g.peak_gb.std(ddof=1) / math.sqrt(n) if n > 1 else float("nan")
    return pd.Series(dict(n=n, time_mean=g.time_s.mean(), time_es2=2 * es_t, peak_max=g.peak_gb.max(), peak_es2=2 * es_p,
                          main_peak_max=g.main_peak_gb.max(), piece_peak_max=g.piece_peak_gb.max()))
S = R.groupby(["conf", "A"]).apply(summ).reset_index(); S.to_csv(f"{RES}/S2.8_summary.csv", index=False)
def verdict(est, e2, hi):
    lo_, hi_ = est - e2, est + e2
    if hi_ <= hi: return "confermata"
    if lo_ > hi: return "smentita"
    return "indeterminata"
P = []
for A in TIS:
    s = S[(S.conf == "e504") & (S.A == A)]
    if len(s): s = s.iloc[0]; P.append(dict(id="P-1", A=A, tessuto=TIS[A], quantity="time_s E504", n=s.n, est=s.time_mean, es2=s.time_es2, lo=s.time_mean - s.time_es2, hi=s.time_mean + s.time_es2, limit=300, esito=verdict(s.time_mean, s.time_es2, 300)))
for A in TIS:
    s = S[(S.conf == "e504") & (S.A == A)]
    if len(s): s = s.iloc[0]; P.append(dict(id="P-2", A=A, tessuto=TIS[A], quantity="peak_gb E504", n=s.n, est=s.peak_max, es2=s.peak_es2, lo=s.peak_max - s.peak_es2, hi=s.peak_max + s.peak_es2, limit=8, esito=verdict(s.peak_max, s.peak_es2, 8)))
for A in TIS:
    s = S[(S.conf == "einf") & (S.A == A)]
    if len(s): s = s.iloc[0]; P.append(dict(id="P-3", A=A, tessuto=TIS[A], quantity="peak_gb EInf", n=s.n, est=s.peak_max, es2=s.peak_es2, lo=s.peak_max - s.peak_es2, hi=s.peak_max + s.peak_es2, limit=8, esito=verdict(s.peak_max, s.peak_es2, 8)))
P = pd.DataFrame(P); P.to_csv(f"{RES}/S2.8_predictions.csv", index=False)
# regola R-1
Q = []
for A in TIS:
    r0 = S[(S.conf == "ref") & (S.A == A)]
    if not len(r0): continue
    t0 = r0.iloc[0].time_mean; cand = []
    for cf in ("e504", "e1008", "einf"):
        s = S[(S.conf == cf) & (S.A == A)]
        if len(s) and s.iloc[0].peak_max <= 8: cand.append((s.iloc[0].time_mean, cf))
    best = min(cand) if cand else (float("nan"), "none")
    Q.append(dict(A=A, tessuto=TIS[A], t_ref=t0, best=best[1], t_best=best[0], ratio=best[0] / t0, within_1_5=best[0] / t0 <= 1.5,
                  t_e504=S[(S.conf == "e504") & (S.A == A)].time_mean.iloc[0], ratio_e504=S[(S.conf == "e504") & (S.A == A)].time_mean.iloc[0] / t0))
Q = pd.DataFrame(Q); Q.to_csv(f"{RES}/S2.8_rule_R1.csv", index=False)
n15 = int(Q.within_1_5.sum()) if len(Q) else 0; nbest504 = int((Q.best == "e504").sum()) if len(Q) else 0
rule = dict(n_tissues=len(Q), n_ratio_le_1_5=n15, single_entry_point=bool(n15 >= 5), n_best_e504=nbest504, backlog_auto_piece=bool(nbest504 < 5))
json.dump(rule, open(f"{RES}/S2.8_rule_R1.json", "w"), indent=1)
print(P.to_string(index=False)); print(Q.to_string(index=False)); print(rule)
# descrittiva di BL-237: A1_r1, finestra contro predefinito, stessi semi
def load(cf, k):
    b = f"{O}/{cf}/A1_r1_k{k}"; j = json.load(open(b + ".json")); g = j["grid"]
    u = np.fromfile(b + ".bin", dtype="<i4"); tis = np.fromfile(b + ".tis.bin", dtype="i1").astype(bool)
    u8 = np.fromfile(b + ".u8.bin", dtype="<i4"); return g, u, tis, u8
D = []
for k in (0, 1, 2):
    if not (os.path.exists(f"{O}/bl237/A1_r1_k{k}.json") and os.path.exists(f"{O}/e504/A1_r1_k{k}.json")): continue
    g, ud, tis, _ = load("e504", k); _, uw, tis_w, _ = load("bl237", k)
    nr, nc = g["n_row"], g["n_col"]; x = g["x0"] + (np.arange(nc) + 0.5) * g["bin_um"]; X = np.tile(x, nr)
    edge = 5.476725  # lato sinistro del tessuto in A1_r1 (riquadro delle regioni)
    zones = {"tissue_all": tis, "tissue_0_8": tis & (X - edge >= 0) & (X - edge < 8), "tissue_8_16": tis & (X - edge >= 8) & (X - edge < 16),
             "tissue_16_32": tis & (X - edge >= 16) & (X - edge < 32), "tissue_32_64": tis & (X - edge >= 32) & (X - edge < 64),
             "glass_window_strip_0_5.48": (~tis) & (X >= 0) & (X < edge), "glass_outside_window_left": (~tis) & (X < 0)}
    for z, msk in zones.items():
        D.append(dict(k=k, zone=z, n_bins=int(msk.sum()), umi_default=float(ud[msk].mean()) if msk.any() else float("nan"),
                      umi_window=float(uw[msk].mean()) if msk.any() else float("nan"), same_tissue_mask=bool((tis == tis_w).all())))
D = pd.DataFrame(D)
if len(D):
    D["ratio_window_default"] = D.umi_window / D.umi_default; D.to_csv(f"{RES}/S2.8_BL237.csv", index=False)
    print(D.groupby("zone")[["n_bins", "umi_default", "umi_window", "ratio_window_default"]].mean().to_string())
