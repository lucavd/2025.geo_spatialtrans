"""tools/F3.1_PH1_mixing.py - F3.1, analisi post hoc PH-1 (dichiarata dopo gli esiti; non cambia gli esiti): i vicini di una cellula sono di
un altro tipo piu' spesso nel simulato che nel vero, e i marcatori estranei (M-3) crescono con questa mescolanza?
Per ogni nucleo con tipo (centro = media dei quadratini del nucleo; vero: tipo delle etichette di R6 nel vocabolario dei preset; simulato:
tipo vero): quota dei K = 6 nuclei con tipo piu' vicini che hanno un altro tipo. Per area: media della quota (cellule dell'area) e
pendenza di M-3 sulla quota (minimi quadrati sulle cellule con M-3). Uso: .venv/bin/python tools/F3.1_PH1_mixing.py <main_dir> <out_csv>"""
import sys, os, glob, numpy as np, pandas as pd
from scipy.spatial import cKDTree
main, out = sys.argv[1:3]; K = 6; rows = []; LV = []
for side, pat in [("real", f"{main}/real/*"), ("sim", f"{main}/sim/*_k*"), ("simnuc", f"{main}/simnuc/*_k*"), ("simq6", f"{main}/simq6/*_k*")]:
    for d in sorted(glob.glob(pat)):
        b = pd.read_parquet(f"{d}/bins.parquet", columns=["x_um", "y_um", "label"]); b = b[b.label > 0]
        cen = b.groupby("label")[["x_um", "y_um"]].mean()
        ct = pd.read_parquet(f"{d}/cells.parquet").drop_duplicates("label").set_index("label").cell_type
        A = os.path.basename(d)[:2]; types = set(pd.read_csv(f"{main}/ref/markers_{A}.csv").cell_type)
        cen["t"] = ct.reindex(cen.index).to_numpy(); cen = cen[cen.t.isin(types)]
        if len(cen) <= K: continue
        tr = cKDTree(cen[["x_um", "y_um"]].to_numpy()); _, nb = tr.query(cen[["x_um", "y_um"]].to_numpy(), k=K + 1)
        t = cen.t.to_numpy(); mix = (t[nb[:, 1:]] != t[:, None]).mean(axis=1)
        cen["mix"] = mix
        cl = pd.read_parquet(f"{d}/cells_label.parquet").set_index("label")
        m = cen.join(cl[["foreign_exp"]], how="inner").dropna(subset=["foreign_exp"])
        slope = np.polyfit(m.mix, m.foreign_exp, 1)[0] if len(m) > 20 and m.mix.std() > 0 else np.nan
        # M-3 nelle cellule senza vicini di altro tipo (mescolanza 0) e con almeno la meta' dei vicini di altro tipo
        f0 = m.foreign_exp[m.mix == 0].mean(); f5 = m.foreign_exp[m.mix >= 0.5].mean()
        lev = np.rint(m.mix.to_numpy() * K).astype(int)
        for l_ in range(K + 1):
            sel = lev == l_; LV.append({"side": side, "run": os.path.basename(d), "A": A, "level": l_, "n": int(sel.sum()), "sum_M3": float(m.foreign_exp.to_numpy()[sel].sum())})
        rows.append({"side": side, "run": os.path.basename(d), "A": A, "n_cells": len(cen), "mix_mean": float(mix.mean()), "slope_M3_on_mix": float(slope),
                     "M3_mix0": float(f0), "M3_mix_ge_half": float(f5), "share_mix0": float((mix == 0).mean())})
T = pd.DataFrame(rows); T.to_csv(out, index=False)
T["tag"] = T.run.str.replace(r"_k\d$", "", regex=True)
S = T.groupby(["side", "A", "tag"])[["mix_mean", "slope_M3_on_mix", "M3_mix0", "M3_mix_ge_half", "share_mix0"]].mean().groupby(["side", "A"]).agg(["mean", "std"])
pd.set_option("display.width", 250); print(S.round(3).to_string())

# riponderazione: M-3 medio del simulato per livello di mescolanza (0..6 vicini di altro tipo), pesato con la distribuzione dei livelli del vero
L = pd.DataFrame(LV); L.to_csv(out.replace(".csv", "_levels.csv"), index=False)
g = L.groupby(["side", "A", "level"])[["n", "sum_M3"]].sum().reset_index()
rw = []
for A in sorted(g.A.unique()):
    re = g[(g.side == "real") & (g.A == A)].set_index("level"); p_real = re.n / re.n.sum()
    for side in ["sim", "simnuc", "simq6"]:
        si = g[(g.side == side) & (g.A == A)].set_index("level")
        if si.empty: continue
        mu = si.sum_M3 / si.n.replace(0, np.nan); p_sim = si.n / si.n.sum()
        ok = mu.notna() & p_real.reindex(mu.index).notna()
        exp_rw = float((mu[ok] * p_real.reindex(mu.index)[ok]).sum() / p_real.reindex(mu.index)[ok].sum())
        obs = float(si.sum_M3.sum() / si.n.sum()); real_obs = float(re.sum_M3.sum() / re.n.sum())
        rw.append({"A": A, "side": side, "M3_real": real_obs, "M3_sim": obs, "M3_sim_reweighted_to_real_mixing": exp_rw,
                   "fraction_of_gap_explained": (obs - exp_rw) / (obs - real_obs) if obs != real_obs else np.nan})
RW = pd.DataFrame(rw); RW.to_csv(out.replace(".csv", "_reweight.csv"), index=False); print(RW.round(3).to_string(index=False))
# PH-2: composizione dei primi tipi di RCTD (vero contro simulato), dalle tabelle per blocco
import json
comp = []
for side, pat in [("real", f"{main}/real/*"), ("sim", f"{main}/sim/*_k*")]:
    for d in sorted(glob.glob(pat)):
        j = json.load(open(f"{d}/rctd.json")); fs = j.get("first_type_share") or {}
        for t_, v_ in fs.items(): comp.append({"side": side, "A": os.path.basename(d)[:2], "type": t_, "share": v_})
C = pd.DataFrame(comp).groupby(["A", "type", "side"]).share.mean().unstack("side").round(3); C.to_csv(out.replace(".csv", "_rctd_composition.csv")); print(C.to_string())
