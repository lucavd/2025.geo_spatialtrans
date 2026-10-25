"""tools/S1.5d_size_resid_table.py - S1.5d: tabella dei residui veri del modello di taglia per gruppo (candidato M3, docs/steps/S1.5d.md).
Residuo della cella interna vera = log A_N - [2 log rho_g + (1 - beta_g) log(1e6 / density_model) + beta_g log A_T] (media del modello del
gruppo di S1.5b nel pacchetto, conversione «preserve_area»: non dipende dalla densita' del tipo). Tipi del gruppo = tipi del preset
(results/S1.5d/S1.5d_preset_types.csv, da tools/S1.5d_export_presets.R); nuclei veri con territorio (insieme di S1.5b/S1.5c), tutte le 5 ROI;
linfonodo e cuore (un gruppo ciascuno): celle su cui R3b ha stimato il modello del preset (tools/S1.5d_export_r3b_cells.R, DEC-006).
Decimi di log A_T assoluto (quantili 0.1...0.9 delle celle del gruppo); per decimo i 101 quantili (p = 0, 0.01, ..., 1) dei residui; decimo 0 =
tutti i residui (controprova CP-3). Controllo: media dei residui del gruppo |m| < 1e-2 (nota D-2: arrotondamento dei parametri del preset) e deviazione standard entro 1e-3 da sigma_g (stessa
stima di S1.5b), altrimenti errore. Uscita: inst/extdata/size_resid.csv (key = <preset>:<gruppo>) e results/S1.5d/S1.5d_size_resid_check.csv."""
import os, importlib.util, numpy as np, pandas as pd
s = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(s); s.loader.exec_module(cal)
PT = pd.read_csv("results/S1.5d/S1.5d_preset_types.csv"); R3B = ("A4", "A6"); CELLS = "/mnt/micron/geo_spatialtrans/S1.5d/size_resid_cells"; P = np.round(np.linspace(0, 1, 101), 2)
rows, chk = [], []
for A in cal.ARCH:
    pt = PT[PT.archetype == A]
    if A in R3B:                                                       # linfonodo e cuore: modello stimato sull'insieme di R3b (DEC-006), nota D-1
        e = pd.read_parquet(f"/mnt/micron/geo_spatialtrans/S1.5d/r3b_sdarea_{A}.parquet")
    else:
        _, d = cal.load(A); e = d[d.interior.astype(bool)]
    for g, gp in pt.groupby("size_group"):
        rho, beta, sig, dm = gp.rho_g.iloc[0], gp.beta_g.iloc[0], gp.sigma_g.iloc[0], gp.density_model.iloc[0]
        assert gp[["rho_g", "beta_g", "sigma_g"]].nunique().max() == 1, (A, g)
        c = e if A in R3B else e[e.type.isin(gp.cell_type)]; x = np.log(c.a_t.values)
        r = np.log(c.a_n.values) - (2 * np.log(rho) + (1 - beta) * np.log(1e6 / dm) + beta * x)
        ok = abs(r.mean()) < 1e-2 and abs(r.std(ddof=2) - sig) < 1e-3   # media: parametri del preset arrotondati a 3 cifre (cuore -0.004, nota D-2)
        chk.append(dict(archetype=A, preset=gp.preset.iloc[0], size_group=int(g), n=len(r), types="|".join(sorted(gp.cell_type)), resid_mean=r.mean(), resid_sd=r.std(ddof=2), sigma_g=sig, ok=ok))
        if not ok: raise SystemExit(f"{A} gruppo {g}: residui non coerenti con il modello del gruppo (media {r.mean():.4g}, sd {r.std(ddof=2):.4g}, sigma_g {sig:.4g})")
        key = f"{gp.preset.iloc[0]}:{int(g)}"; os.makedirs(CELLS, exist_ok=True)
        pd.DataFrame(dict(a_t=c.a_t.values, a_n=c.a_n.values)).to_parquet(f"{CELLS}/{key.replace(':', '_')}.parquet")   # celle usate (test e G-4)
        edges = np.quantile(x, np.linspace(0.1, 0.9, 9)); lo = np.r_[-np.inf, edges]; hi = np.r_[edges, np.inf]
        di = np.searchsorted(edges, x, side="left") + 1                # (lo, hi]: stesso lato di findInterval(left.open = TRUE) in R
        for k in range(1, 11):
            q = np.quantile(r[di == k], P)
            rows += [dict(key=key, decile=k, x_lo=lo[k - 1], x_hi=hi[k - 1], n=int((di == k).sum()), p=p, r=v) for p, v in zip(P, q)]
        q = np.quantile(r, P); rows += [dict(key=key, decile=0, x_lo=-np.inf, x_hi=np.inf, n=len(r), p=p, r=v) for p, v in zip(P, q)]
os.makedirs("inst/extdata", exist_ok=True)
pd.DataFrame(rows).to_csv("inst/extdata/size_resid.csv", index=False, float_format="%.10g")
pd.DataFrame(chk).to_csv("results/S1.5d/S1.5d_size_resid_check.csv", index=False); print(pd.DataFrame(chk).to_string())
