
import pandas as pd, numpy as np
d = pd.read_csv("results/S1.4b/review/code/nonoverlap_opt.csv"); y = pd.read_csv("/tmp/code_vc/yards_A2A4.csv"); d = d.merge(y, on="file")
cl = lambda v, a, b: np.where(v <= a, "P", np.where(v <= b, "W", "F"))
rows = []
for alg in ["onepass", "greedy", "lp_r", "lp_a"]:
    for A in ["A2", "A4"]:
        z = d[d.file.str.startswith(A)]
        f1 = "".join(cl(z[f"{alg}_ks_nc"], z.y_nc_med, z.y_nc_max)); f2 = "".join(cl(z[f"{alg}_ks_an"], z.y_an_med, z.y_an_max))
        rows.append(dict(alg=alg, A=A, B4=round(z[f"{alg}_trunc"].mean(), 4), loss=round(z[f"{alg}_loss"].mean(), 4), f_B1=f1, f_B2=f2))
    rows.append(dict(alg="mvc(min)", A="A2", B4=round(d[d.file.str.startswith("A2")].mvc_trunc.mean(), 4)))
    rows.append(dict(alg="mvc(min)", A="A4", B4=round(d[d.file.str.startswith("A4")].mvc_trunc.mean(), 4)))
r = pd.DataFrame(rows).drop_duplicates(); print(r.to_string()); r.to_csv("results/S1.4b/review/code/nonoverlap_alt_verdicts_s1.csv", index=False)
print(d[["file","y_an_med","onepass_ks_an","ab_ks_an_s1","ab_ks_an_med10"]].to_string())
