"""S2.2, analisi post hoc (dichiarate dopo i risultati; non cambiano gli esiti pre-registrati).
PH-1: braccio oracolo (sottotipo = gruppo vero di R6 della cellula, disposizione spaziale vera) contro il vincitore e contro il vero:
      scarti relativi di m1...m4 e indice degli zeri in eccesso dell'intestino. Ipotesi: lo scarto di m3 (intestino) e degli zeri viene
      dalla disposizione dei sottotipi nello spazio e dall'assegnazione per cellula, non dai profili.
PH-2: scomposizione del fattore di profondità per area: f ~ (UMI dell'area / UMI nelle impronte dei nuclei veri) x (nuclei veri / cellule simulate)
      x (totale medio vero per nucleo / totale medio simulato per cellula) - fondo.
Uso: .venv/bin/python tools/S2.2_posthoc_summary.py -> results/S2.2/S2.2_PH1.csv, S2.2_PH1_P5.csv, S2.2_PH2.csv"""
import glob, json, numpy as np, pandas as pd
M = "/mnt/micron/geo_spatialtrans"; RES = "results/S2.2"
MET = {"m1": "m1", "m2": "m2", "m3": "m3", "m4": "m4_ratio"}
real = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{M}/S2.2/race/*_real.json"))]).set_index(["archetype", "roi_id"])
ph = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{M}/S2.2/posthoc/*_ph1.csv"))], ignore_index=True)
per = pd.read_csv(f"{RES}/S2.2_race_per_roi.csv"); win = pd.read_csv(f"{RES}/S2.2_race_winners.csv").set_index("archetype").winner
rows = []
for (a, arm), g in ph.groupby(["archetype", "arm"]):
    gm = g.groupby("roi_id")[list(MET.values())].mean()
    for k, col in MET.items():
        e = (gm[col] - real.loc[a][col].reindex(gm.index)) / real.loc[a][col].reindex(gm.index)
        rows.append(dict(archetype=a, arm=arm, metric=k, e_mean=e.mean(), e_se=e.std(ddof=1) / np.sqrt(len(e))))
    for base in ("exp_M", "exp_G", win[a]):
        gb = per[(per.archetype == a) & (per.arm == base)].set_index("roi_id")
        for k in MET:
            rows.append(dict(archetype=a, arm=base, metric=k, e_mean=gb[f"e_{k}"].mean(), e_se=gb[f"e_{k}"].std(ddof=1) / np.sqrt(len(gb))))
out = pd.DataFrame(rows).drop_duplicates(["archetype", "arm", "metric"]); out.to_csv(f"{RES}/S2.2_PH1.csv", index=False)
print(out.pivot_table(index=["archetype", "arm"], columns="metric", values="e_mean").round(3).to_string())
# zeri in eccesso dell'intestino per i bracci oracolo
files = sorted(glob.glob(f"{M}/S2.2/race/A1_*_p5.npz")); pfiles = sorted(glob.glob(f"{M}/S2.2/posthoc/A1_*_ph1_p5.npz"))
real_p = [np.concatenate([np.load(f)[k] for k in np.load(f).files if k.startswith("real|")]) for f in files]
r5 = []
for arm in ("exp_M_oracle", "exp_G_oracle"):
    sims = []
    for r in range(3):
        sims.append(np.median(np.concatenate([np.concatenate([np.load(f)[k] for k in np.load(f).files if k.startswith(f"{arm}|{r}|")]) for f in pfiles])))
    r5.append(dict(arm=arm, sim=np.mean(sims), real=np.median(np.concatenate(real_p))))
pd.DataFrame(r5).to_csv(f"{RES}/S2.2_PH1_P5.csv", index=False); print(pd.DataFrame(r5).round(4).to_string(index=False))
# PH-2
dp = pd.read_csv(f"{RES}/S2.2_depth_per_roi.csv"); rows = []
for _, r in dp.iterrows():
    info = json.load(open(f"{M}/R6/counts/{r.archetype}_{r.roi_id}_info.json"))
    cov = info["umi_assigned_exp"] / info["umi_roi"]
    n_real = info["n_nuclei_keep"]; mean_real = info["umi_assigned_exp"] / n_real
    mean_sim = r.exp_tot / r.n_cells
    area_ratio = r.area_ratio_sim_real
    f_dec = (1 / cov) * (n_real / r.n_cells) * (mean_real / mean_sim) * area_ratio * (r.mean_umi2 - r.b) / r.mean_umi2
    rows.append(dict(archetype=r.archetype, roi_id=r.roi_id, f=r.f, coverage=cov, density_ratio_real_sim=n_real / r.n_cells * area_ratio,
                     mean_ratio_real_sim=mean_real / mean_sim, background_term=(r.mean_umi2 - r.b) / r.mean_umi2, f_decomposed=f_dec))
d2 = pd.DataFrame(rows); d2.to_csv(f"{RES}/S2.2_PH2.csv", index=False)
print(d2.groupby("archetype")[["f", "f_decomposed", "coverage", "density_ratio_real_sim", "mean_ratio_real_sim"]].mean().round(3).to_string())
