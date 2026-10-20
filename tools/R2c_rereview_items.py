#!/usr/bin/env python
"""tools/R2c_rereview_items.py — R2c, CP-6 POST HOC (RV-14; richiesta di Luca 2026-10-20 «Riesame ora, come controprova post hoc»).
Gruppo 'fast': coppie 'pair' di Luca giudicate 'separate' in < 0.8 s con d/d_max <= q75 dei 'touch' (results/R2c/R2c_posthoc_fast_separate.csv).
Controlli: 15 'ctrl_sep' (stesso intervallo di distanza relativa, 'separate', t >= 1.5 s) e 15 'ctrl_touch' (stesso intervallo, 'touch').
Ordine casuale (seme 20261020). Uscita: results/R2c/R2c_items_Luca_rev.csv (stesso formato degli item; annotatore con
tools/R2c_annotator.py Luca_rev)."""
import importlib.util, numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/R2c"
spec = importlib.util.spec_from_file_location("ana", ROOT / "tools/R2c_analysis.py"); A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
L = A.load_answers("Luca"); dm = pd.read_csv(RES / "R2c_dmax.csv").set_index("archetype").d_max_um
P = L[L.kind == "pair"].copy(); P["rel"] = P.d_um / P.win_id.str[:2].map(dm)
fast = pd.read_csv(RES / "R2c_posthoc_fast_separate.csv"); q75 = float(fast.q75_rel_touch.iloc[0])
near = P[P.rel <= q75]; rng = np.random.default_rng(20261020)
cs = near[(near.answer == "separate") & (near.t_ms >= 1500)]; ct = near[near.answer == "touch"]
pick = lambda d, n: d.iloc[rng.choice(len(d), size=min(n, len(d)), replace=False)]
G = pd.concat([P[P.pair_id.isin(fast.pair_id)].assign(kind="fast"), pick(cs, 15).assign(kind="ctrl_sep"), pick(ct, 15).assign(kind="ctrl_touch")])
G = G[["pair_id", "win_id", "i", "j", "d_um", "kind", "answer", "t_ms"]].rename(columns={"answer": "answer_first", "t_ms": "t_first_ms"})
G = G.iloc[rng.permutation(len(G))].reset_index(drop=True); G.insert(0, "item", np.arange(1, len(G) + 1))
G.to_csv(RES / "R2c_items_Luca_rev.csv", index=False)
print(len(G), G.kind.value_counts().to_dict(), "pool ctrl_sep", len(cs), "ctrl_touch", len(ct), "q75", round(q75, 3))
