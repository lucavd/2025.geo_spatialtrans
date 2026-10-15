
import numpy as np, pandas as pd
M = pd.read_csv("results/S1.4b/review/code/nonoverlap_sens_raw.csv"); Y = pd.read_csv("/tmp/code_tgt/yards.csv")
M = M.merge(Y[["A", "k", "b7_real"]], on=["A", "k"]); M["b7d"] = (M.b7 - M.b7_real).abs()
Fd = M.groupby(["alg", "A", "k"]).median(numeric_only=True).reset_index().drop(columns="b7_real").merge(Y, on=["A", "k"])
cl = lambda v, a, b: np.where(v <= a, "P", np.where(v <= b, "W", "F"))
Fd["v1"] = cl(Fd.ks_nc, Fd.y_nc_med, Fd.y_nc_max); Fd["v2"] = cl(Fd.ks_an, Fd.y_an_med, Fd.y_an_max); Fd["v3"] = cl(Fd.d_sp, Fd.y_sp_med, Fd.y_sp_max); Fd["v7"] = cl(Fd.b7d, Fd.y_b7_med, Fd.y_b7_max)
av = lambda o: "PASS" if (o == "P").sum() >= 4 else ("FAIL" if (o == "F").sum() >= 2 else "WARN")
rows = []
for (a, A), z in Fd.sort_values("k").groupby(["alg", "A"]):
    rows.append(dict(alg=a, A=A, B1=av(z.v1), B2=av(z.v2), B3=av(z.v3), B7=av(z.v7), f1="".join(z.v1), f2="".join(z.v2), f3="".join(z.v3), f7="".join(z.v7),
                     B4=round(z.trunc_all.mean(), 4), loss=round(z.loss.mean(), 4), b7_sim=round(z.b7.mean(), 3)))
At = pd.DataFrame(rows); At.to_csv("results/S1.4b/review/code/nonoverlap_sens_archetypes.csv", index=False); Fd.to_csv("results/S1.4b/review/code/nonoverlap_sens_folds.csv", index=False)
pd.set_option("display.width", 200); print(At.to_string())
for a in At.alg.unique(): print(a, "score", int((At[At.alg == a][["B1", "B2", "B3", "B7"]] == "PASS").sum().sum()))
