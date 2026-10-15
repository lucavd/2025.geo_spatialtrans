
import pandas as pd, numpy as np
me = pd.read_csv("results/S1.4b/review/code/cp1_indep_raw.csv"); ref = pd.read_csv("results/S1.4b/S1.4b_cp1_raw.csv")
j = me.merge(ref[["archetype", "roi", "cond", "m1", "m2"]], on=["archetype", "roi", "cond"], suffixes=("", "_ref"))
j["d1"] = j.m1 - j.m1_ref; j["d2"] = j.m2 - j.m2_ref
j["kind"] = np.where(j.cond == "vor", "vor", np.where(j.cond == "lag", "lag", "perm"))
print(j.groupby("kind").agg(n=("d1", "size"), n_diff=("d1", lambda v: int((v != 0).sum())), max_abs_d1=("d1", lambda v: v.abs().max()), mean_d1=("d1", "mean"), max_abs_d2=("d2", lambda v: v.abs().max())).to_string())
r = pd.read_csv("results/S1.4b/review/code/cp1_indep_roi.csv"); print(r[~r.lag_beats_vor | ~r.vor_beats_permmin].to_string())
