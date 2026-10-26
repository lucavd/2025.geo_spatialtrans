"""tools/S1.5e_g4.py - S1.5e, cancello G-4 (docs/steps/S1.5e.md).
(a) forma dei vuoti (allungamento mediano pesato, quota allungata, numero di componenti) ricalcolata con una seconda implementazione indipendente
    (skimage.measure.label + regionprops: area e autovalori del tensore d'inerzia) sui 30 ROI veri e su 30 unita' simulate (braccio XA, seme 1);
    stessi pixel vuoti (la soglia 3 x mediana della distanza dal vicino e' la definizione, non la parte ricalcolata): <= 1e-9;
(a2) controllo di coerenza: la forma dei vuoti veri di tools/S1.5e_metrics.py = valori dell'esplorazione E4 (results/S1.5e/S1.5e_explore_void_shape.csv);
(b) il banco legge per ogni ROI i parametri della sua piega (X, XA, XAD): parametri registrati nel layer (results/S1.5e/S1.5e_panel_units.csv,
    colonne het_*) = riga loo_<roi> di results/S1.5e/S1.5e_het_<famiglia>.csv, 30/30 per famiglia (tutti i semi).
Uscite: results/S1.5e/S1.5e_G4_check.csv. Uso: .venv/bin/python tools/S1.5e_g4.py"""
import importlib.util, numpy as np, pandas as pd
from skimage.measure import label, regionprops
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
m5e = mod("m5e", "tools/S1.5e_metrics.py"); cal = mod("cal", "tools/S1.5b_calibration.py")
RES = "results/S1.5e"; CELLS = "/mnt/micron/geo_spatialtrans/S1.5e/cells"
def shape_skimage(v, px):
    rs = [r for r in regionprops(label(v, connectivity=2)) if r.area * px * px >= m5e.MIN_COMP_UM2]
    if not rs: return dict(elong_wmed=1.0, frac_elong_gt3=np.nan, n_comp=0)
    a = np.array([r.area for r in rs], float) * px * px; ev = [np.sort(r.inertia_tensor_eigvals) for r in rs]
    e = np.array([np.sqrt(x[-1] / max(x[0], 1e-12)) for x in ev]); w = a / a.sum(); o = np.argsort(e); cw = np.cumsum(w[o])
    return dict(elong_wmed=float(e[o][np.searchsorted(cw, 0.5)]), frac_elong_gt3=float(w[e > 3].sum()), n_comp=len(rs))
def cmp(d1, d2):
    out = []
    for k in ["elong_wmed", "frac_elong_gt3", "n_comp"]:
        a, b = d1[k], d2[k]; out.append(0.0 if (np.isnan(a) and np.isnan(b)) else abs(a - b))
    return max(out)
rows = []; R = m5e.m5d.rois_table(); E4 = pd.read_csv(f"{RES}/S1.5e_explore_void_shape.csv"); E4 = E4[E4.source == "real"]
for A, roi in zip(R.archetype, R.roi_id):
    g = m5e.RoiGeomE(A, roi); _, d = cal.load(A); e = d[d.roi_id == roi]
    v = m5e.void_image(e.x.values, e.y.values, e.interior.values.astype(bool), g); a = m5e.summarise(m5e.components(v, g.px)); b = shape_skimage(v, g.px)
    rows.append(dict(check="a vero", unit=f"{A}_{roi}", value=cmp(a, b), pass_=cmp(a, b) <= 1e-9))
    x = E4[(E4.archetype == A) & (E4.roi_id == roi)].iloc[0]
    ref = dict(elong_wmed=x.elong_wmed if x.n_comp > 0 else 1.0, frac_elong_gt3=x.frac_area_elong_gt3, n_comp=int(x.n_comp))
    rows.append(dict(check="a2 vero = esplorazione E4", unit=f"{A}_{roi}", value=cmp(a, ref), pass_=cmp(a, ref) <= 1e-9))
    c = pd.read_parquet(f"{CELLS}/XA/{A}_{roi}_s01.parquet"); v = m5e.void_image(c.x.values, c.y.values, c.interior.values.astype(bool), g)
    a = m5e.summarise(m5e.components(v, g.px)); b = shape_skimage(v, g.px)
    rows.append(dict(check="a simulato (XA, seme 1)", unit=f"{A}_{roi}", value=cmp(a, b), pass_=cmp(a, b) <= 1e-9))
U = pd.read_csv(f"{RES}/S1.5e_panel_units.csv")
for fam in ["X", "XA", "XAD"]:
    H = pd.read_csv(f"{RES}/S1.5e_het_{fam}.csv"); u = U[U.arm == fam]
    for A, roi in zip(R.archetype, R.roi_id):
        h = H[(H.archetype == A) & (H.fold == f"loo_{roi}")].iloc[0]; z = u[(u.archetype == A) & (u.roi_id == roi)]
        ok = len(z) == 10 and np.allclose(z.het_void_frac, h.void_frac) and np.allclose(z.het_void_aspect, h.void_aspect) \
             and np.allclose(z.het_void_ell_um, h.void_ell_um) and np.allclose(z.het_dens_sigma, h.dens_sigma) and (z.het_void_mode == h.void_mode).all()
        rows.append(dict(check=f"b piega {fam}", unit=f"{A}_{roi}", value=float(ok), pass_=bool(ok)))
G = pd.DataFrame(rows).rename(columns={"pass_": "pass"}); G.to_csv(f"{RES}/S1.5e_G4_check.csv", index=False)
print(G.groupby("check")["pass"].agg(["sum", "count"]).to_string()); print("GATE G-4", "PASS" if G["pass"].all() else "FAIL")
