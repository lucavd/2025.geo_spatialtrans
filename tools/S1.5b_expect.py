# tools/S1.5b_expect.py - S1.5b: stime a mano degli attesi (non simulazioni del banco).
# I modelli di taglia per tipo (results/S1.5b/S1.5b_type_calibration.csv, bracci "group" e "recentred") e quello a un tipo di S1.5
# si applicano ai territori delle celle interne simulate in S1.5 (seme 1), con tipi estratti secondo la composizione e area del
# nucleo = min(bersaglio, phi * territorio) (approssimazione del ritaglio). 10 estrazioni per ROI, seme 1.
# Uscite: results/S1.5b/S1.5b_expected_roi.csv, results/S1.5b/S1.5b_expected_type.csv
import numpy as np, pandas as pd
S = pd.read_parquet("/mnt/micron/geo_spatialtrans/S1.5b/s15_seed1_terr.parquet")
C = pd.read_csv("results/S1.5b/S1.5b_type_calibration.csv")
R5 = pd.read_csv("results/S1.5b/S1.5b_real_refs_roi.csv"); THR = float(R5.small_threshold_um2.dropna().iloc[0])
PHI = {"A1": 0.757, "A2": 1, "A3": 0.947, "A4": 0.791, "A5": 0.474, "A6": 0.149}
PRE = {"A1": (12419, 10591.0114668242, 0.380, 0.114, 0.475), "A2": (8924, 7723.28633246204, 0.510, 0.320, 0.547),
       "A3": (3100, 4278.87001473923, 0.374, 0.147, 0.553), "A4": (28096, 27897.5, 0.679, 0.455, 0.210),
       "A5": (1185, 1039.7, 0.282, 0.621, 0.813), "A6": (988, 1062.4, 0.154, 0.122, 0.490)}
rng = np.random.default_rng(1); res, tres = [], []
for (A, r), g in S.groupby(["archetype", "roi_id"]):
    at = g.a_t.values; it = g.interior.values; c = C[C.archetype == A].reset_index(drop=True)
    D, Dm, rho, b, s = PRE[A]; rc = rho * (D / Dm) ** ((1 - b) / 2); aref = 1e6 / D
    for rep in range(10):
        z = rng.standard_normal(len(at)); ti = rng.choice(len(c), size=len(at), p=c.fraction_tissue / c.fraction_tissue.sum())
        arms = {"single": np.exp(np.log(rc ** 2 * aref) + b * (np.log(at) - np.log(aref)) + s * z)}
        for arm in ("group", "recentred"):
            a_t = 1e6 / c.density_type.values[ti]; rr = c[f"rho_catalog_{arm}"].values[ti]
            arms[arm] = np.exp(np.log(rr ** 2 * a_t) + c.beta_group.values[ti] * (np.log(at) - np.log(a_t)) + c.sigma_group.values[ti] * z)
        for arm, tgt in arms.items():
            an = np.minimum(tgt, PHI[A] * at); nt = an / at
            rec = dict(archetype=A, roi_id=r, rep=rep, arm=arm, median_area_nuc=np.median(an[it]), median_nt=np.median(nt[it]),
                       frac_nt095=np.mean(nt[it] >= 0.95), frac_trunc=np.mean(tgt > PHI[A] * at))
            if A == "A5":
                rec["frac_band"] = np.mean((an[it] > THR / np.sqrt(2)) & (an[it] < THR * np.sqrt(2)))
                if arm != "single":
                    ng = c.size_group[c.cell_type == "neuroni"].iloc[0]; grp = c.size_group.values[ti]
                    rec["ratio_neuron_glia"] = np.median(an[it & (grp == ng)]) / np.median(an[it & (grp != ng)])
            res.append(rec)
            for k, t in enumerate(c.cell_type):
                m = it & (ti == k)
                if m.any(): tres.append(dict(archetype=A, roi_id=r, rep=rep, arm=arm, cell_type=t, median_area_nuc=np.median(an[m]), median_nt=np.median(nt[m])))
E = pd.DataFrame(res).groupby(["archetype", "roi_id", "arm"]).mean(numeric_only=True).reset_index().drop(columns="rep")
T = pd.DataFrame(tres).groupby(["archetype", "roi_id", "arm", "cell_type"]).mean(numeric_only=True).reset_index().drop(columns="rep")
E.to_csv("results/S1.5b/S1.5b_expected_roi.csv", index=False); T.to_csv("results/S1.5b/S1.5b_expected_type.csv", index=False)
pd.set_option("display.width", 250)
print(E.groupby(["archetype", "arm"]).mean(numeric_only=True).round(4).to_string())
print(T.groupby(["archetype", "cell_type", "arm"]).agg(area=("median_area_nuc", "mean"), nt=("median_nt", "mean")).unstack("arm").round(3).to_string())
