"""tools/S1.5c_g4_mc.py - S1.5c, indagine sul cancello G-4 (parte a, fallito come scritto: 75/79 bin entro 0.10). Misura l'errore Monte Carlo delle curve
di calibrazione: con i parametri di tutte le ROI, 5 lotti indipendenti di 10 repliche dei campi (stesse funzioni di tools/S1.5c_calibration.py, semi diversi),
curva del rapporto di connessione per lotto. Confronta: la curva di Python con 3 repliche (file di calibrazione), quella di R con 10 repliche (G-4) e la
dispersione fra lotti. Uso: .venv/bin/python tools/S1.5c_g4_mc.py   Uscita: results/S1.5c/S1.5c_G4_mc.csv"""
import sys, importlib.util, numpy as np, pandas as pd
argv = sys.argv; sys.argv = [argv[0]]
s = importlib.util.spec_from_file_location("c15", "tools/S1.5c_calibration.py"); C = importlib.util.module_from_spec(s); s.loader.exec_module(C)
sys.argv = argv
G4 = pd.read_csv("results/S1.5c/S1.5c_G4_calcheck.csv"); G4 = G4[G4.part == "a"]
rows = []
for A in C.cal.ARCH:
    lab, _ = C.cal.load(A); lab = lab[~lab.type.isin(list(C.cal.UNTYPED))]
    types = sorted(lab.type.unique()); main_idx = [types.index(t) for t in C.MAIN[A]]
    f = lab.type.value_counts(normalize=True).reindex(types).values
    P0 = pd.read_csv(f"results/S1.5c/S1.5c_compartments_{A}.csv"); P0 = P0[P0.fold == "pooled"]
    ell = [C.ELL[0]] * len(types); sv = [0.0] * len(types)
    for _, r in P0.iterrows(): k = types.index(r.cell_type); ell[k] = int(r.ell_um); sv[k] = float(r.s)
    rng = np.random.default_rng(C.SEED); rois = []
    for r, e in lab.groupby("roi_id"):
        P = e[["x_um", "y_um"]].values; L = np.array([types.index(t) for t in e.type]); rois.append(C.ROI(P, L, types, rng))
    num_r = sum(ro.real[0] for ro in rois); den_r = sum(ro.real[1] for ro in rois); Rr = C.ratio(num_r, den_r, f)
    curves = []
    for b in range(5):
        ns = 0; ds = 0
        for k, ro in enumerate(rois):
            ro.Z = {}
            for rep in range(10):
                for t in main_idx:
                    if sv[t] > 0: ro.Z[(rep, t, ell[t])] = C.rff(ro.P, ell[t], np.random.default_rng([C.SEED, 1000 + 10 * b + rep, t, ell[t]]))
                Lb = C.labels(ro, f, ell, sv, ro.Z, rep, np.random.default_rng(C.SEED + 7919 * (1000 + 10 * b + rep) + 104729 * k))
                n_, d_, _ = ro.counts(Lb); ns = ns + n_; ds = ds + d_
        curves.append(C.ratio(ns, ds, f))
    cv = np.array(curves)
    for t in main_idx:
        for kb in range(C.NB):
            bn = f"{int(C.BINS[kb])}-{int(C.BINS[kb+1])}"; g = G4[(G4.archetype == A) & (G4.cell_type == types[t]) & (G4.bin == bn)]
            lb = np.log(cv[:, t, kb]); rows.append(dict(archetype=A, cell_type=types[t], bin=bn, real_pairs=num_r[t, kb], R_real=Rr[t, kb],
                 py10_mean=float(np.exp(lb.mean())), sd_log_10rep=float(lb.std(ddof=1)), sd_log_3rep=float(lb.std(ddof=1) * np.sqrt(10 / 3)),
                 R_fit_py3=float(g.R_fit_py.iloc[0]) if len(g) else np.nan, R_fit_R10=float(g.R_fit_R.iloc[0]) if len(g) else np.nan))
    print(A, "fatto", flush=True)
D = pd.DataFrame(rows)
D["z_R10_vs_py10"] = (np.log(D.R_fit_R10) - np.log(D.py10_mean)) / (D.sd_log_10rep * np.sqrt(2))
D["z_py3_vs_py10"] = (np.log(D.R_fit_py3) - np.log(D.py10_mean)) / np.sqrt(D.sd_log_3rep ** 2 + D.sd_log_10rep ** 2 / 5)
D.to_csv("results/S1.5c/S1.5c_G4_mc.csv", index=False); print(D.round(3).to_string())
