"""tools/F3.1_PH4b_blocks.py - F3.1, PH-4 domanda 2 a livello di blocco: da quali blocchi vengono le chiamate 'plasmacellule' di RCTD nel linfonodo.
Per i blocchi non scartati: primo tipo di RCTD contro tipo del nucleo dominante (vero e simulato) e contro il tipo che domina gli UMI del blocco (solo simulato),
e contro la purezza del blocco (quota degli UMI di cellule del tipo dominante, solo simulato). Uso: .venv/bin/python tools/F3.1_PH4b_blocks.py <main_dir> <out_dir>"""
import sys, glob, os, numpy as np, pandas as pd
main, out = sys.argv[1:3]; rows = []
for side, pat in [("real", f"{main}/real/A4_*"), ("sim", f"{main}/sim/A4_*_k*"), ("simq6", f"{main}/simq6/A4_*_k0")]:
    for d in sorted(glob.glob(pat)):
        m = pd.read_csv(f"{d}/rctd_blocks.csv"); m = m[m.spot_class.notna() & (m.spot_class != "reject")].copy(); m["side"] = side; m["area"] = os.path.basename(d)
        rows.append(m)
M = pd.concat(rows, ignore_index=True)
pd.set_option("display.width", 220)
t1 = pd.crosstab([M.side, M.nuc_type], M["first"], normalize="index").round(3); print("primo tipo dato il tipo del nucleo dominante\n", t1.to_string())
S = M[M.side != "real"].copy()
t2 = pd.crosstab([S.side, S.true_type], S["first"], normalize="index").round(3); print("primo tipo dato il tipo che domina gli UMI (simulato)\n", t2.to_string())
S["purity_bin"] = pd.cut(S.true_share, [0, 0.5, 0.7, 0.9, 1.0001], labels=["<0.5", "0.5-0.7", "0.7-0.9", ">=0.9"])
B = S[S.true_type == "linfociti B"]
t3 = pd.crosstab([B.side, B.purity_bin], B["first"], normalize="index").round(3); print("blocchi dominati da B, per purezza\n", t3.to_string())
print("quota di blocchi dominati da B per purezza:", B.groupby("side").purity_bin.value_counts(normalize=True).round(3).to_dict())
pc = M.groupby(["side"]).apply(lambda x: (x["first"] == "plasmacellule").mean()).round(3); print("quota di blocchi chiamati plasmacellule:", pc.to_dict())
# scomposizione dell'eccesso di chiamate 'plasmacellule' (simulato k = 0, 1 contro vero) in composizione dei blocchi e tasso per tipo del nucleo
M["nt"] = M.nuc_type.fillna("nessun nucleo tipizzato"); M["pl"] = M["first"] == "plasmacellule"
g = M.groupby(["side", "nt"]).agg(n=("pl", "size"), rate=("pl", "mean")); g["w"] = g.n / g.groupby("side").n.transform("sum")
r_, s_ = g.loc["real"], g.loc["sim"]; idx = r_.index.intersection(s_.index); tot = (g.w * g.rate).groupby("side").sum()
dec = {"real": tot["real"], "sim": tot["sim"], "excess": tot["sim"] - tot["real"], "composition": ((s_.w - r_.w) * r_.rate).reindex(idx).sum(),
       "rate": (s_.w * (s_.rate - r_.rate)).reindex(idx).sum(), "rate_nonB": (s_.w * (s_.rate - r_.rate)).reindex(idx).drop("linfociti B").sum()}
print("scomposizione:", {k: round(v, 3) for k, v in dec.items()}); pd.Series(dec).to_csv(f"{out}/F3.1_PH4b_decomposition.csv"); g.to_csv(f"{out}/F3.1_PH4b_rates.csv")
for nm, t in [("given_nucleus", t1), ("given_truth", t2), ("B_by_purity", t3)]: t.to_csv(f"{out}/F3.1_PH4b_{nm}.csv")
