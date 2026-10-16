#!/usr/bin/env python
"""results/S1.4c/review/recalc/r08_phi_b4.py — revisione avversariale S1.4c (recalc): phi di K-2 ricalcolato
(q99 tipo 7 di a_n/a_t sull'insieme ev dei 4 ROI di addestramento del fold) e B-4 per fold (mediana sui semi)."""
import numpy as np, pandas as pd
D = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/dump/"; O = D + "../"
rois = pd.read_csv("results/R2/R2_rois_checked.csv"); ph = pd.read_csv(D + "gen_phi.csv"); yd = pd.read_csv(D + "loo_yard.csv")
nc = {}
for A, r in zip(rois.archetype, rois.roi_id):
    cc = pd.read_parquet(D + f"cells_{A}_{r}.parquet"); ev = cc.ev_prep.to_numpy().astype(bool)
    nc[(A, r)] = (cc.a_n / cc.a_t).to_numpy()[ev]
rows = []
for _, z in ph.iterrows():
    kr = yd[(yd.archetype == z.archetype) & (yd.fold == z.fold)].roi_test.iloc[0]
    tr = [r for r in rois.roi_id[rois.archetype == z.archetype] if r != kr]
    q = float(np.quantile(np.concatenate([nc[(z.archetype, r)] for r in tr]), 0.99))
    rows.append(dict(archetype=z.archetype, fold=z.fold, q99_mine=q, q99_rep=z.q99, phi_mine=min(1, q), phi_rep=z.phi))
p = pd.DataFrame(rows); p.to_csv(O + "r08_phi.csv", index=False)
print("phi: max|mio - gen|", float(np.abs(p.phi_mine - p.phi_rep).max()), "| phi per archetipo (min-max):", p.groupby("archetype").phi_mine.agg(["min", "max"]).round(3).to_dict("index"))
b4 = pd.read_csv(D + "gen_b4.csv")
f = b4.groupby(["archetype", "variant", "fold"]).trunc_all.median().unstack("fold")
f["mean"] = f.mean(axis=1); f.to_csv(O + "r08_b4_by_fold.csv")
print(f.xs("K-2", level="variant").round(4).to_string())
