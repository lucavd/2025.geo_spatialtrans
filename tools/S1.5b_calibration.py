# tools/S1.5b_calibration.py - S1.5b: riferimenti per tipo dei preset (composizione, gruppi di taglia, modello di taglia per gruppo).
# Definizioni: docs/steps/S1.5b.md (regola dei gruppi di taglia, insieme di valutazione, densita' del tipo puro).
# Ingressi (sola lettura): etichette R6 PH-1 (DEC-013), territori reali esportati da R6 (insieme di calibrazione di S1.4/R3b),
#   tabella R6 dei riferimenti per tipo, scarto fra ROI dell'area nucleare mediana (R6, P-2), preset di S1.5.
# Uso (root del repo): .venv/bin/python tools/S1.5b_calibration.py [--a5 r6|corrected]
import sys, json, numpy as np, pandas as pd
M = "/mnt/micron/geo_spatialtrans/R6"; RES = "results/S1.5b"
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
# scarto fra ROI dell'area nucleare mediana (max/min dei 5 ROI), R6 P-2, da results/S1.5/S1.5_real_refs_roi.csv
MARGIN = {"A1": 1.095, "A2": 1.161, "A3": 1.216, "A4": 1.062, "A5": 2.096, "A6": 1.808}
# preset di S1.5 (R/04b0_tissue_presets.R): densita' manuale, densita' del modello di taglia, rho, beta, sigma
PRESET = {"A1": (12419, 10591.0114668242, 0.380, 0.114, 0.475), "A2": (8924, 7723.28633246204, 0.510, 0.320, 0.547),
          "A3": (3100, 4278.87001473923, 0.374, 0.147, 0.553), "A4": (28096, 27897.5, 0.679, 0.455, 0.210),
          "A5": (1185, 1039.7, 0.282, 0.621, 0.813), "A6": (988, 1062.4, 0.154, 0.122, 0.490)}
UNTYPED = {"non risolto", "non classificabile"}
EPI = {"epiteliale normale": "epiteliale", "epiteliale tumorale": "epiteliale"}   # DEC-014

def fit(an, at, aref):
    """Stesso stimatore di fit_nucleus_size_model(model = "full") in R/04b4_derive_cell_geometry.R."""
    y = np.log(an); x = np.log(at) - np.log(aref)
    b = np.cov(x, y, ddof=1)[0, 1] / np.var(x, ddof=1); a = y.mean() - b * x.mean(); r = y - a - b * x
    return float(np.sqrt(np.exp(a) / aref)), float(b), float(np.sqrt((r ** 2).sum() / (len(y) - 2))), int(len(y))

def load(A):
    lab = pd.read_parquet(f"{M}/labels/{A}_exp_umi20_seed20261022_nuclei_ph1.parquet")
    lab["type"] = lab.cell_type.replace(EPI)
    parts = []
    for r in sorted(lab.roi_id.unique()):
        T = pd.read_parquet(f"{M}/terr/{A}_{r}.parquet").assign(roi_id=r)
        T["xk"] = T.x.round(6); T["yk"] = T.y.round(6)
        L = lab[lab.roi_id == r].assign(xk=lambda d: d.x_um.round(6), yk=lambda d: d.y_um.round(6))
        m = T.merge(L[["xk", "yk", "type", "area_um2"]], on=["xk", "yk"], how="left", validate="one_to_one")
        if m.type.isna().any(): raise SystemExit(f"{A} {r}: territori senza etichetta nell'unione per coordinate")
        m.attrs["n_label_without_territory"] = len(L) - len(m)   # nuclei etichettati senza territorio esportato (maschera R3b)
        parts.append(m.assign(n_lab_no_terr=len(L) - len(m)))
    return lab, pd.concat(parts, ignore_index=True)

def size_groups(lab, A):
    """Regola dei gruppi di taglia (pre-registrata in docs/steps/S1.5b.md).
    1) Tipi principali = tipi con almeno il 5 % di tutti i nuclei del tessuto (soglia di R6, previsione P-2).
    2) Fra i principali, legame completo sul logaritmo dell'area nucleare mediana: due gruppi si uniscono (prima la coppia piu'
       vicina) se nel gruppo unito il rapporto fra la mediana piu' grande e la piu' piccola resta entro lo scarto fra ROI del tessuto.
    3) Ogni tipo minore entra nel gruppo principale con la mediana (dei nuclei del gruppo) piu' vicina in logaritmo."""
    med = lab[~lab.type.isin(UNTYPED)].groupby("type").area_um2.median()
    frac = lab.type.value_counts() / len(lab)
    main = sorted([t for t in med.index if frac.get(t, 0) >= 0.05]); minor = sorted([t for t in med.index if t not in main])
    groups = [[t] for t in main]; trace = []
    while len(groups) >= 2:
        best = None
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                v = med[groups[i] + groups[j]]; r = v.max() / v.min()
                if best is None or r < best[0]: best = (r, i, j)
        r, i, j = best; ok = r <= MARGIN[A]
        trace.append(dict(archetype=A, step=len(trace) + 1, kind="merge_main", a="|".join(groups[i]), b="|".join(groups[j]),
                          ratio=r, margin=MARGIN[A], accepted=ok))
        if not ok: break
        groups[i] = sorted(groups[i] + groups[j]); del groups[j]
    for t in minor:
        gm = [lab.area_um2[lab.type.isin(g)].median() for g in groups]
        k = int(np.argmin([abs(np.log(med[t] / x)) for x in gm]))
        trace.append(dict(archetype=A, step=len(trace) + 1, kind="assign_minor", a=t, b="|".join(groups[k]), ratio=med[t] / gm[k],
                          margin=MARGIN[A], accepted=True))
        groups[k] = sorted(groups[k] + [t])
    return groups, trace

def main():
    a5 = sys.argv[sys.argv.index("--a5") + 1] if "--a5" in sys.argv else "r6"
    rows, traces, checks = [], [], []
    for A in ARCH:
        D, Dm, rho0, b0, s0 = PRESET[A]
        lab, ev_all = load(A)
        ev = ev_all[ev_all.interior & ev_all.matched & (ev_all.a_n > 0)]
        aref = 1e6 / Dm
        pr = fit(ev.a_n.values, ev.a_t.values, aref)
        checks.append(dict(archetype=A, n_label_without_territory=int(ev_all.groupby("roi_id").n_lab_no_terr.first().sum()), rho_refit=pr[0], beta_refit=pr[1], sigma_refit=pr[2], n=pr[3], rho_preset=rho0, beta_preset=b0, sigma_preset=s0))
        groups, tr = size_groups(lab, A); traces += tr
        typed = lab[~lab.type.isin(UNTYPED)]
        n_t = typed.type.value_counts()
        comp = {"tissue": n_t / n_t.sum()}
        if A == "A4":
            for v, rr in (("follicle", ["f1", "f2"]), ("paracortex", ["p1", "p2", "p3"])):
                x = typed[typed.roi_id.isin(rr)].type.value_counts(); comp[v] = x / x.sum()
        # territorio mediano per tipo sull'insieme di valutazione (densita' del tipo puro; braccio ricentrato)
        m_t = ev[~ev.type.isin(UNTYPED)].groupby("type").a_t.median()
        for gi, g in enumerate(groups):
            if len(groups) == 1:
                rho_g, beta_g, sig_g, n_g, src = rho0, b0, s0, np.nan, "preset S1.5"
            else:
                e = ev[ev.type.isin(g)]; rho_g, beta_g, sig_g, n_g = fit(e.a_n.values, e.a_t.values, aref); src = "fit S1.5b"
            for t in g:
                rows.append(dict(archetype=A, cell_type=t, size_group=gi + 1, group_types="|".join(g), n_typed=int(n_t[t]),
                                 rho_group=rho_g, beta_group=beta_g, sigma_group=sig_g, n_fit=n_g, size_source=src,
                                 median_area_nuc=float(typed.area_um2[typed.type == t].median()), median_a_t_eval=float(m_t[t]),
                                 **{f"fraction_{k}": float(v.get(t, 0.0)) for k, v in comp.items()}))
    cal = pd.DataFrame(rows)
    # territorio mediano delle celle interne simulate in S1.5 (braccio preserve_area, 50 unita' per tessuto): ancora del braccio ricentrato
    P = pd.read_csv("results/S1.5/S1.5_panel_units.csv"); P = P[P.arm == "preserve_area"]
    m_sim = (np.pi * P.median_eq_r ** 2).groupby(P.archetype).mean()
    # densita' del tipo puro: rho_t = D * mbar / m_t, mbar = somma f_t m_t  ->  media armonica della miscela = D
    out = []
    for A, g in cal.groupby("archetype", sort=False):
        D, Dm = PRESET[A][0], PRESET[A][1]
        mbar = (g.fraction_tissue * g.median_a_t_eval).sum()
        g = g.assign(density_type=D * mbar / g.median_a_t_eval)
        # braccio (a) "group": modello del gruppo com'e' stimato, portato alla densita' del tipo puro (stessa funzione A_N(A_T))
        g = g.assign(rho_catalog_group=g.rho_group * (g.density_type / Dm) ** ((1 - g.beta_group) / 2))
        # braccio (b) "recentred": stessa pendenza, intercetta spostata perche' al territorio simulato tipico (m_sim) il tipo abbia
        # la taglia che ha al proprio territorio reale tipico (m_t): A_N'(A) = A_N(A) * (m_t / m_sim)^beta
        g = g.assign(m_sim=float(m_sim[A]), rho_catalog_recentred=g.rho_catalog_group * (g.median_a_t_eval / float(m_sim[A])) ** (g.beta_group / 2))
        g = g.assign(mix_density_check=1 / (g.fraction_tissue / g.density_type).sum())
        out.append(g)
    cal = pd.concat(out)
    import os; os.makedirs(RES, exist_ok=True)
    cal.to_csv(f"{RES}/S1.5b_type_calibration.csv", index=False)
    pd.DataFrame(traces).to_csv(f"{RES}/S1.5b_size_groups.csv", index=False)
    pd.DataFrame(checks).to_csv(f"{RES}/S1.5b_refit_check.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(traces).to_string(index=False))
    print(cal[["archetype", "cell_type", "size_group", "fraction_tissue", "rho_group", "beta_group", "sigma_group", "n_fit", "density_type", "rho_catalog_group", "rho_catalog_recentred", "mix_density_check"]].round(4).to_string(index=False))
    print(pd.DataFrame(checks).round(6).to_string(index=False))

if __name__ == "__main__":
    main()
