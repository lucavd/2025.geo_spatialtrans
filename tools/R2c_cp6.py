#!/usr/bin/env python
"""tools/R2c_cp6.py — R2c, CP-6 POST HOC (addendum 2026-10-20 in R2c_preregistration.md, scritto e committato prima del riesame).
f = frazione 'touch' al riesame nel gruppo fast (prima 'separate' in < 0.8 s); c = idem in ctrl_sep; coerenza = 'touch' in ctrl_touch.
Regola: abitudine sospetta se f - c >= 0.10 e Fisher esatto unilaterale (fast > ctrl_sep) p < 0.05; altrimenti non sostenuta.
Sensibilita': S-1 manuale ricalcolata sostituendo le risposte del riesame alle originali (IC a grappolo come pre-registrato e
bootstrap a due stadi come in R2c_posthoc). Uscite: results/R2c/R2c_cp6.csv, R2c_cp6_s1.csv."""
import importlib.util, json, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/R2c"
spec = importlib.util.spec_from_file_location("ana", ROOT / "tools/R2c_analysis.py"); A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
spec = importlib.util.spec_from_file_location("ph", ROOT / "tools/R2c_posthoc.py"); H = importlib.util.module_from_spec(spec); spec.loader.exec_module(H)
it = pd.read_csv(RES / "R2c_items_Luca_rev.csv"); d = json.load(open(RES / "annotations/R2c_answers_Luca_rev.json"))
assert d["doc_id"] == "R2c_Luca_rev_v10" and d["n_items"] == len(it)
a = pd.DataFrame(d["answers"]); assert (a["item"].to_numpy() == it["item"].to_numpy()).all()
R = it.merge(a, on="item"); R["flip_touch"] = R.answer == "touch"
rows = []
for k, g in R.groupby("kind"):
    rows.append(dict(group=k, n=len(g), touch=int((g.answer == "touch").sum()), separate=int((g.answer == "separate").sum()),
                     unsure=int((g.answer == "unsure").sum()), frac_touch=float((g.answer == "touch").mean()), t_med_s=g.t_ms.median() / 1000))
T = pd.DataFrame(rows).set_index("group")
f, cc = T.frac_touch["fast"], T.frac_touch["ctrl_sep"]
tab = [[T.touch["fast"], T.n["fast"] - T.touch["fast"]], [T.touch["ctrl_sep"], T.n["ctrl_sep"] - T.touch["ctrl_sep"]]]
p = fisher_exact(tab, alternative="greater")[1]
verdict = "sospetta" if (f - cc >= 0.10 and p < 0.05) else "non sostenuta"
T = T.reset_index(); T["f_minus_c"] = f - cc; T["fisher_p_one_sided"] = p; T["verdict"] = verdict
T.to_csv(RES / "R2c_cp6.csv", index=False)
# sensibilita': S-1 con le risposte del riesame
L = A.load_answers("Luca"); new = R.set_index("pair_id").answer
L2 = L.copy(); m = (L2.kind == "pair") & L2.pair_id.isin(new.index)
L2.loc[m, "answer"] = L2.loc[m, "pair_id"].map(new)
focal = pd.read_csv(RES / "R2c_focal.csv"); focal = focal[focal.focal][["win_id", "archetype", "idx"]].reset_index(drop=True)
s1 = pd.read_csv(RES / "R2c_s1.csv").set_index("archetype"); out = []
for X, lab in [(L, "originale"), (L2, "con riesame")]:
    P = X[X.kind == "pair"][["win_id", "i", "j", "answer"]]
    for ar in A.ARCH:
        fa = focal[focal.archetype == ar].reset_index(drop=True); w = fa.win_id.to_numpy()
        for v in ["primary", "unsure_touch"]:
            y = A.focal_s1(P[P.win_id.str.startswith(ar + "_")], fa, v)
            e, lo, hi = A.cluster_boot(w, y); l2, h2 = H.two_stage(w, y)
            out.append(dict(answers=lab, archetype=ar, variant=v, s1=e, lo_cluster=lo, hi_cluster=hi, lo_two_stage=l2, hi_two_stage=h2,
                            b_dec_cluster=A.b_dec(lo, hi, e, s1.seg_cellpose_rgb[ar], s1.seg_stardist_he[ar]),
                            b_dec_two_stage=A.b_dec(l2, h2, e, s1.seg_cellpose_rgb[ar], s1.seg_stardist_he[ar]),
                            k2=s1.k2_s1[ar], k2_verdict_cluster=A.k2_verdict(s1.k2_s1[ar], lo, hi, s1.widen[ar])))
S = pd.DataFrame(out); S.to_csv(RES / "R2c_cp6_s1.csv", index=False)
pd.set_option("display.width", 220)
print(T.round(3).to_string(index=False)); print(S.round(3).to_string(index=False))
print("fast flipped by archetype:", R[(R.kind == "fast") & R.flip_touch].win_id.str[:2].value_counts().to_dict())
