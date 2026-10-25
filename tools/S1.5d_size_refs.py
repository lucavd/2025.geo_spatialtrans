"""tools/S1.5d_size_refs.py - S1.5d: riferimenti per l'aggiunta sul modello di taglia, calcolati prima di ogni simulazione.
(R1) per gruppo di taglia del preset (S1.5b, .TP_TYPES): pendenza c1 del logaritmo della dispersione su log A_T (massima verosimiglianza, media lineare
libera), xbar = media di log A_T delle celle interne del gruppo; modello M1 del pacchetto: sigma(A_T) = sigma_g * exp(c1 * (log A_T - xbar)).
Tutte le ROI (come il modello di taglia di S1.5b) e pieghe lascia-fuori-un-ROI (per la stabilita'). (R2) controllo sui territori veri, lascia-fuori:
media del gruppo dal modello di S1.5b ricalcolato sulle 4 ROI, frazione N/C >= 0.95 con M0 e M1 contro il vero. (R3) intervalli fra ROI vere del
90deg percentile dell'area del nucleo per tipo principale (celle interne). Uscite: results/S1.5d/S1.5d_size_groups.csv, S1.5d_size_check.csv, S1.5d_real_q90.csv"""
import importlib.util, numpy as np, pandas as pd
from scipy.optimize import minimize
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py")
TC = pd.read_csv("results/S1.5b/S1.5b_type_calibration.csv")
PRE = {"A1": "epithelial", "A2": "tumor_microenv", "A3": "stromal_rich", "A4": "lymphoid", "A5": "brain_cortex", "A6": "columnar_muscle"}
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
def fit(x, y):
    xm = x.mean(); b = np.polyfit(x, y, 1)
    def nll(p):
        s = np.exp(p[2] + p[3] * (x - xm)); return np.sum(np.log(s) + 0.5 * ((y - p[0] - p[1] * x) / s) ** 2)
    p = minimize(nll, [b[1], b[0], np.log((y - np.polyval(b, x)).std()), 0.0], method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-7, fatol=1e-9)).x
    return dict(a=p[0], b=p[1], c0=p[2], c1=p[3], xbar=xm, sigma0=(y - np.polyval(b, x)).std(ddof=2), b_ols=b[0], a_ols=b[1])
G, C, Q = [], [], []
for A in cal.ARCH:
    _, d = cal.load(A); e = d[d.interior.astype(bool)].copy()
    tc = TC[TC.archetype == A] if "archetype" in TC else TC[TC.preset == PRE[A]]
    gmap = dict(zip(tc.cell_type, tc.size_group)); e["g"] = e.type.map(gmap); e = e[e.g.notna()]
    for gi, gg in e.groupby("g"):
        for fold in ["pooled"] + [f"loo_{r}" for r in sorted(gg.roi_id.unique())]:
            tr = gg if fold == "pooled" else gg[gg.roi_id != fold[4:]]
            f = fit(np.log(tr.a_t.values), np.log(tr.a_n.values))
            G.append(dict(archetype=A, size_group=int(gi), fold=fold, n=len(tr), types="|".join(sorted(tr.type.unique())), **f))
            if fold != "pooled":
                te = gg[gg.roi_id == fold[4:]]; xt = np.log(te.a_t.values); at = te.a_t.values
                rng = np.random.default_rng(20261025); E = rng.standard_normal((40, len(xt)))
                mu = f["a_ols"] + f["b_ols"] * xt
                an0 = np.exp(mu + f["sigma0"] * E); an1 = np.exp(mu + f["sigma0"] * np.exp(f["c1"] * (xt - f["xbar"])) * E)
                for t in sorted(te.type.unique()):
                    k = (te.type == t).values
                    if k.sum() < 30: continue
                    C.append(dict(archetype=A, size_group=int(gi), roi_id=fold[4:], cell_type=t, n=int(k.sum()), real_nc095=np.mean(te.a_n.values[k] / at[k] >= 0.95),
                                  M0_nc095=np.mean(an0[:, k] / at[k] >= 0.95), M1_nc095=np.mean(an1[:, k] / at[k] >= 0.95),
                                  real_q90=np.quantile(te.a_n.values[k], 0.9), M0_q90=np.quantile(an0[:, k], 0.9), M1_q90=np.quantile(an1[:, k], 0.9),
                                  real_med=np.median(te.a_n.values[k]), M0_med=np.median(an0[:, k]), M1_med=np.median(an1[:, k])))
    for t in MAIN[A]:
        for r, gr in e[e.type == t].groupby("roi_id"):
            Q.append(dict(archetype=A, cell_type=t, roi_id=r, n=len(gr), q90_an=np.quantile(gr.a_n, 0.9), q10_an=np.quantile(gr.a_n, 0.1)))
pd.DataFrame(G).to_csv("results/S1.5d/S1.5d_size_groups.csv", index=False)
pd.DataFrame(C).to_csv("results/S1.5d/S1.5d_size_check.csv", index=False)
Qd = pd.DataFrame(Q); Qd.to_csv("results/S1.5d/S1.5d_real_q90_roi.csv", index=False)
Qd.groupby(["archetype", "cell_type"]).agg(q90_min=("q90_an", "min"), q90_max=("q90_an", "max"), n_roi=("q90_an", "size")).reset_index().to_csv("results/S1.5d/S1.5d_real_q90.csv", index=False)
print("ok")
