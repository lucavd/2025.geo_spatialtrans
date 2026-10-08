# rv_task2_recompute.py — revisione avversariale S1.2 / checkB, compito 2.
# Ricalcolo da zero (numpy/pandas) di D, CP-1, CP-2, CP-3, B1b, B2 a partire da:
#   data/real_pcf_export.csv, data/sim_export.csv (esportazioni 1:1 di real_pcf.rds e sim/*.rds, rv_export.R),
#   results/S1.2/S1.2_pcf_curves.csv, S1.2_sim_models.csv, S1.2_B1b.csv; confronto riga per riga con S1.2_checkB_summary.csv.
# Uso: cd ~/2025.geo_spatialtrans && .venv/bin/python results/S1.2/review/checkB/rv_task2_recompute.py
import numpy as np, pandas as pd
RV = "results/S1.2/review/checkB"; RES = "results/S1.2"
R = np.arange(0, 30.0001, 0.5); I = R > 0
ARCH = pd.DataFrame(dict(archetype=[f"A{i}" for i in range(1, 7)], tot=[12419, 8924, 3100, 28096, 1185, 988],
                         evid=[10291, 7376, 2311, 26844, 953, 906], d_nuc=[4.1, 6.3, 6.3, 4.2, 11.3, 4.8]))
ARCH["eq_r"] = np.sqrt(1e6 / (np.pi * ARCH.tot)); ARCH["d_rule"] = 2 / 3 * ARCH.eq_r
real = pd.read_csv(f"{RV}/data/real_pcf_export.csv")
sim = pd.read_csv(f"{RV}/data/sim_export.csv")
gcols = [c for c in sim.columns if c.startswith("g_")]; assert len(gcols) == 61
def D(a, b): return float(np.trapezoid((np.asarray(a)[I] - np.asarray(b)[I]) ** 2, R[I]))
greal = {(a, r_, role): grp.sort_values("r").g.to_numpy() for (a, r_, role), grp in real.groupby(["archetype", "roi_id", "role"])}
checks = []   # (id, oggetto, riportato, ricalcolato, uguale)
# --- coerenza curve: S1.2_pcf_curves.csv vs esportazione rds
cur = pd.read_csv(f"{RES}/S1.2_pcf_curves.csv")
mx = 0.0
for (a, r_, src), grp in cur.groupby(["archetype", "roi_id", "source"]):
    g = grp.sort_values("r").g.to_numpy()
    if src.startswith("reale"):
        ref = greal[(a, r_, "primary" if src == "reale primario" else "secondary")]
    else:
        ref = sim[(sim.archetype == a) & (sim.roi_id == r_) & (sim.model == src)][gcols].to_numpy()[0]
    mx = max(mx, np.max(np.abs(g - ref)))
checks.append(("curve_csv_vs_rds_maxabs", mx))
# --- D per riga
sim["D_rv"] = [D(sim.loc[i, gcols].to_numpy(float), greal[(sim.archetype[i], sim.roi_id[i], "primary")]) if sim.reps[i] > 0 else np.nan for i in sim.index]
sm = pd.read_csv(f"{RES}/S1.2_sim_models.csv")
mm = sim.merge(sm[["archetype", "roi_id", "model", "D"]], on=["archetype", "roi_id", "model"], how="outer", indicator=True)
checks.append(("sim_models_rows_matched", (mm._merge == "both").sum(), len(sm), len(sim)))
checks.append(("D_maxabs_diff", np.nanmax(np.abs(mm.D_rv - mm.D))))
checks.append(("D_NA_mismatch", int((mm.D_rv.isna() != mm.D.isna()).sum())))
# --- integrale da r=0 (pre-registrazione letterale) per confronto
g0 = lambda a: np.asarray(a)
sim["D_from0"] = [float(np.trapezoid((sim.loc[i, gcols].to_numpy(float) - greal[(sim.archetype[i], sim.roi_id[i], "primary")]) ** 2, R)) if sim.reps[i] > 0 else np.nan for i in sim.index]
rows = []
def add(A, check, value, status, extra=""):
    rows.append(dict(archetype=A, check=check, value_rv=value, status_rv=status, extra=extra))
# --- CP-1
just = {}; just0 = {}
cp1rows = []
for A in ARCH.archetype:
    s = sim[sim.archetype == A]
    pd_ = s[s.model == "PD"].set_index("roi_id"); cs = s[s.model == "CSR"].set_index("roi_id")
    red = 1 - pd_.D_rv / cs.D_rv; win = int((pd_.D_rv < cs.D_rv).sum())
    j = win >= 4 and np.median(red) >= 0.25; just[A] = j
    red0 = 1 - pd_.D_from0 / cs.D_from0; just0[A] = int((pd_.D_from0 < cs.D_from0).sum()) >= 4 and np.median(red0) >= 0.25
    add(A, "CP-1", f"{win}/5; {np.median(red):.2f}", "PASS" if j else "FAIL", f"da r=0: {int((pd_.D_from0 < cs.D_from0).sum())}/5; {np.median(red0):.2f}")
    for roi in pd_.index: cp1rows.append(dict(archetype=A, roi_id=roi, D_PD=pd_.D_rv[roi], D_CSR=cs.D_rv[roi], red=red[roi]))
add("tutti", "CP-1", f"{sum(just.values())}/6", "PASS" if sum(just.values()) >= 4 else "FAIL", f"da r=0: {sum(just0.values())}/6")
# --- CP-2
cp2rows = []
for A in ARCH.archetype:
    a = ARCH[ARCH.archetype == A].iloc[0]
    gr = sim[(sim.archetype == A) & sim.grid]
    tab = []
    for dval, s in gr.groupby("d"):
        ok = bool(((s.feasible) & (s.reps == 20)).all()) and len(s) == 5
        tab.append(dict(d=dval, n_ok=int(((s.feasible) & (s.reps == 20)).sum()), Dmed=np.median(s.D_rv) if ok else np.nan))
    tab = pd.DataFrame(tab).sort_values("d")
    t2 = tab.dropna()
    dstar = t2.d.iloc[int(np.argmin(t2.Dmed))]; Dstar = t2.Dmed.min()
    Drule = np.median(sim[(sim.archetype == A) & (sim.model == "PD")].D_rv)
    Dnuc = np.median(sim[(sim.archetype == A) & (sim.model == "NUC")].D_rv)
    # regola valutata sulla griglia (d piu' vicino) come controllo
    dgrid_rule = t2.d.iloc[int(np.argmin(np.abs(t2.d - a.d_rule)))]
    ratio = Drule / Dstar
    # non monotonia/forma: secondo minimo, d massimo fattibile
    add(A, "CP-2", f"{ratio:.2f}; d* {dstar:.1f} (regola {a.d_rule:.2f}, nucleo {a.d_nuc:.1f}: {Dnuc / Dstar:.2f})", "PASS" if ratio <= 1.25 else "FAIL",
        f"d_max_fattibile={t2.d.max():.1f}; primo d non fattibile={tab[tab.Dmed.isna()].d.min() if tab.Dmed.isna().any() else 'nessuno'}; "
        f"D(griglia piu' vicina alla regola, d={dgrid_rule})/D*={t2.Dmed[t2.d == dgrid_rule].iloc[0] / Dstar:.2f}")
    cp2rows.append(dict(archetype=A, d_star=dstar, D_star=Dstar, D_rule=Drule, ratio_rule=ratio, D_nuc=Dnuc, d_max_feasible=t2.d.max(),
                        n_d_feasible=len(t2)))
# --- CP-3
p3 = {}
for A in ARCH.archetype:
    rr = sorted(real[real.archetype == A].roi_id.unique())
    Dseg = [D(greal[(A, r_, "primary")], greal[(A, r_, "secondary")]) for r_ in rr]
    Dcsr = sim[(sim.archetype == A) & (sim.model == "CSR")].D_rv; Dpd = sim[(sim.archetype == A) & (sim.model == "PD")].D_rv
    v = np.median(Dseg) / np.median(Dcsr); p3[A] = v < 0.5
    add(A, "CP-3", f"{v:.2f}; {np.median(Dseg) / np.median(Dpd):.2f}", "PASS" if v < 0.5 else "FAIL")
add("tutti", "CP-3", f"{sum(p3.values())}/6", "PASS" if sum(p3.values()) >= 4 else "FAIL")
# --- B1b (dal CSV del progetto: solo la regola di decisione)
b1 = pd.read_csv(f"{RES}/S1.2_B1b.csv")
for A in ARCH.archetype:
    s = b1[b1.archetype == A]; a = ARCH[ARCH.archetype == A].iloc[0]
    k = int(((s.dens_seed42 >= a.evid) & (s.dens_seed42 <= a.tot)).sum())
    add(A, "B1b", f"{k}/5 (esclusa mediana {np.median(s.frac_excluded):.3f})", "PASS" if k >= 4 else "FAIL",
        f"seed42 > atteso in {int((s.dens_seed42 > s.expected).sum())}/5 ROI")
# --- B1a (coerenza)
for A in ARCH.archetype:
    s = sim[(sim.archetype == A) & (sim.model == "PD")]; a = ARCH[ARCH.archetype == A].iloc[0]
    area = s.n_mean / s.dens_mean
    k = int(((s.dens_min >= a.evid - 1 / area) & (s.dens_max <= a.tot + 1 / area)).sum())
    kexcess = (s.dens_max - a.tot) * area   # cellule oltre il totale nella replica peggiore
    add(A, "B1a", f"{k}/5", "PASS" if k == 5 else "FAIL", f"max cellule oltre tot*area: {kexcess.max():.2f}; |media-tot|/tot max {np.max(np.abs(s.dens_mean - a.tot) / a.tot):.1e}")
# --- B2: inviluppo min-max dei 5 ROI reali vs media dei 5 ROI di g_PD
b2rows = []
for A in ARCH.archetype:
    rr = sorted(real[real.archetype == A].roi_id.unique())
    Rm = np.vstack([greal[(A, r_, "primary")] for r_ in rr]); lo, hi = Rm.min(0), Rm.max(0)
    for m in ("PD", "CSR", "NUC"):
        G = sim[(sim.archetype == A) & (sim.model == m)][gcols].to_numpy(float)
        gs = G.mean(0); cov = np.mean((gs >= lo) & (gs <= hi))  # su r = 0.5..30
        gsI = gs[I]; covI = np.mean((gsI >= lo[I]) & (gsI <= hi[I]))
        per_roi = [np.mean((G[i][I] >= lo[I]) & (G[i][I] <= hi[I])) for i in range(G.shape[0])]
        b2rows.append(dict(archetype=A, model=m, coverage=covI, coverage_per_roi_median=np.median(per_roi), coverage_per_roi_max=np.max(per_roi)))
        if m == "PD":
            add(A, "B2", f"{covI:.2f}", "PASS" if covI >= 0.8 else ("WARN" if covI >= 0.5 else "FAIL"), f"per ROI: mediana {np.median(per_roi):.2f}, max {np.max(per_roi):.2f}")
res = pd.DataFrame(rows)
S = pd.read_csv(f"{RES}/S1.2_checkB_summary.csv")
S["k"] = S.groupby(["archetype", "check"]).cumcount(); res["k"] = res.groupby(["archetype", "check"]).cumcount()
cmp = S.merge(res, on=["archetype", "check", "k"], how="outer")
cmp["value_match"] = cmp.value == cmp.value_rv; cmp["status_match"] = cmp.status == cmp.status_rv
cmp.drop(columns="k").to_csv(f"{RV}/rv_task2_summary_compare.csv", index=False)
pd.DataFrame(cp2rows).to_csv(f"{RV}/rv_task2_cp2.csv", index=False)
pd.DataFrame(cp1rows).to_csv(f"{RV}/rv_task2_cp1.csv", index=False)
pd.DataFrame(b2rows).to_csv(f"{RV}/rv_task2_b2.csv", index=False)
sim.drop(columns=gcols).to_csv(f"{RV}/rv_task2_sim_D.csv", index=False)
pd.DataFrame([dict(check=c[0], values=";".join(map(str, c[1:]))) for c in checks]).to_csv(f"{RV}/rv_task2_consistency.csv", index=False)
print(pd.DataFrame([dict(check=c[0], values=c[1:]) for c in checks]).to_string())
print(cmp[["archetype", "check", "value", "value_rv", "status", "status_rv", "value_match", "status_match"]].to_string())
print(cmp[["archetype", "check", "extra"]].to_string())
