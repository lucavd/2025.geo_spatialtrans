"""tools/S1.5d_posthoc_size.py - S1.5d, analisi post hoc (dopo i verdetti, senza verdetto) sul modello di taglia, solo dati veri.
PH-2 (CP-3 smentita): sui territori veri del gruppo, frazione attesa di nuclei >= 0.95 del territorio con tre forme dei residui intorno alla media del
preset: gaussiana (sigma del gruppo), empirica per decimo di log A_T (M3), empirica di tutti i decimi insieme (W3p). Probabilita' esatte per cella
(nessuna estrazione): P(r > log(0.95 A_T) - media(A_T)). Riportata anche la dispersione dei residui per decimo.
PH-3 (P-14 nel linfonodo): 90deg percentile dei residui rispetto alla media del preset: gaussiana (1.2816 sigma), tabella di R3b (decimo 0), nuclei di R6
(insieme degli intervalli di P-8b e P-14). Uscite: results/S1.5d/S1.5d_posthoc_size_tail.csv, S1.5d_posthoc_size_q90.csv"""
import importlib.util, numpy as np, pandas as pd
from scipy.stats import norm
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py")
tab = pd.read_csv("inst/extdata/size_resid.csv"); PT = pd.read_csv("results/S1.5d/S1.5d_preset_types.csv")
def cdf_table(q, x):            # P(r <= x) per la distribuzione con quantili q (interpolazione lineare in p)
    return np.interp(x, q, np.linspace(0, 1, 101), left=0.0, right=1.0)
rows, rows2 = [], []
for key, tk in tab.groupby("key"):
    pre, g = key.split(":"); ty = PT[(PT.preset == pre) & (PT.size_group == int(g))].iloc[0]
    cc = pd.read_parquet(f"/mnt/micron/geo_spatialtrans/S1.5d/size_resid_cells/{pre}_{g}.parquet"); x = np.log(cc.a_t.values)
    mu = 2 * np.log(ty.rho_g) + (1 - ty.beta_g) * np.log(1e6 / ty.density_model) + ty.beta_g * x; thr = np.log(0.95) + x - mu
    r = np.log(cc.a_n.values) - mu
    hi = tk[tk.decile > 0].groupby("decile").x_hi.first().values; dc = np.searchsorted(hi[:-1], x, side="left") + 1
    Q = {d: tk[tk.decile == d].sort_values("p").r.values for d in range(0, 11)}
    p_dec = np.array([1 - cdf_table(Q[d], t) for d, t in zip(dc, thr)]); p_pool = 1 - cdf_table(Q[0], thr); p_gau = norm.sf(thr / ty.sigma_g)
    rows.append(dict(key=key, archetype=ty.archetype, n=len(x), real=np.mean(r >= thr), gaussian=p_gau.mean(), empirical_decile=p_dec.mean(), empirical_pooled=p_pool.mean(),
                     sd_decile1=r[dc == 1].std(ddof=1), sd_decile5=r[dc == 5].std(ddof=1), sd_decile10=r[dc == 10].std(ddof=1),
                     share_tail_from_decile1_3=float(p_pool[dc <= 3].sum() / p_pool.sum())))
    rows2.append(dict(key=key, archetype=ty.archetype, gaussian_q90=1.2816 * ty.sigma_g, table_q90=float(np.interp(0.9, np.linspace(0, 1, 101), Q[0]))))
A4 = PT[PT.archetype == "A4"].iloc[0]; _, d = cal.load("A4"); e = d[d.interior.astype(bool)]
r6 = np.log(e.a_n.values) - (2 * np.log(A4.rho_g) + (1 - A4.beta_g) * np.log(1e6 / A4.density_model) + A4.beta_g * np.log(e.a_t.values))
for R in rows2:
    if R["archetype"] == "A4": R["R6_q90"] = float(np.quantile(r6, 0.9)); R["R6_sd"] = float(r6.std(ddof=1))
pd.DataFrame(rows).to_csv("results/S1.5d/S1.5d_posthoc_size_tail.csv", index=False); pd.DataFrame(rows2).to_csv("results/S1.5d/S1.5d_posthoc_size_q90.csv", index=False)
pd.set_option("display.width", 250); print(pd.DataFrame(rows).round(4).to_string(index=False)); print(pd.DataFrame(rows2).round(4).to_string(index=False))
