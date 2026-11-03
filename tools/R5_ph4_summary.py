"""tools/R5_ph4_summary.py — riassunto della sensibilita' R-4 (R5_ph4_r4shift.py): media per tessuto del rapporto pesato per cellula
(diviso maschera spostata di 20 µm) con bootstrap sulle ROI come R4 PH-1 (stesso seme e stesso generatore: sign = 0 riproduce R4_PH1.csv).
Uso: .venv/bin/python tools/R5_ph4_summary.py <outdir> -> <outdir>/tables/R5_PH4_r4shift_summary.csv"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R4_common as C4
OUT = sys.argv[1]
T = pd.read_csv(f"{OUT}/tables/R5_PH4_r4shift.csv")
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
r4 = pd.read_csv("results/R4/R4_PH1.csv").query("quantity == 'cellweighted_star'").set_index("archetype")
rows = []
for a in ARCH:
    for sign, g in T[T.archetype == a].groupby("sign"):
        e, lo, hi, _ = C4.boot_one_stage(g.star.values, rng=np.random.default_rng([C4.SEED, ARCH.index(a), 1000 + 2]))
        rows.append({"archetype": a, "sign": sign, "mask": {0: "R4 (no shift)", 1: "shifted by section delta (aligned)", -1: "shifted by -delta (control)"}[sign],
                     "estimate": e, "lo": lo, "hi": hi, "class_vs_0.909_1.1": C4.verdict(lo, hi, 0.909, 1.1),
                     "frac_nuc_bins_same_label_mean": g.frac_nuc_bins_same_label.mean(), "spearman_nuc_umi_mean": g.spearman_nuc_umi.mean(),
                     "R4_PH1_estimate": r4.loc[a, "estimate"], "R4_PH1_lo": r4.loc[a, "lo"], "R4_PH1_hi": r4.loc[a, "hi"]})
S = pd.DataFrame(rows); S.to_csv(f"{OUT}/tables/R5_PH4_r4shift_summary.csv", index=False)
pd.set_option("display.width", 250); print(S.round(4).to_string(index=False))
