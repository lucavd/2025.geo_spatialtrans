"""tools/S1.6_lib.py - S1.6: valori per ROI (per seme in P-1) delle previsioni P-1...P-17 di S1.5e, per braccio e insieme di semi.
Stesse definizioni, stessi dati e stessi riferimenti di tools/S1.5e_analysis.py (verificato in tools/S1.6_explore.py: i verdetti di S1.5e si riproducono
riga per riga), salvo i riferimenti per tipo della taglia (P-4, P-8b con `rt_file`; P-14 con `q90_file`): S1.6 li usa sull'insieme del modello di taglia
(tools/S1.6_real_refs.py, BL-163, docs/steps/S1.6.md).
Uso: import da tools/S1.6_explore.py e tools/S1.6_analysis.py."""
import numpy as np, pandas as pd
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
CONTACT = {"A1": (0.366, 0.650), "A2": (0.217, 0.567), "A3": (0.083, 0.283), "A4": (0.350, 0.750), "A5": (0.000, 0.217), "A6": (0.000, 0.070)}  # B-127...B-132
UM_COLS = ["median_area_nuc", "median_nc", "frac_nc095", "frac_trunc", "cv_terr", "cv_terr_nofrag", "median_ecc_T"]
rd = lambda f: pd.read_csv(f, keep_default_na=False, na_values=[""])

def load_refs(q90_file="results/S1.5d/S1.5d_real_q90.csv", rt_file="results/S1.5b/S1.5b_real_refs_roi_type.csv"):
    RI = rd("results/S1.5c/S1.5c_real_intervals.csv"); RI[["cell_type", "bin"]] = RI[["cell_type", "bin"]].fillna("")
    DR = pd.read_csv("results/S1.5e/S1.5e_explore_density_ref.csv"); DR["est"] = DR.dens_seg_roi * DR.dens_manual_window / DR.dens_seg_window
    RT = pd.read_csv(rt_file)
    return dict(RI=RI, RD=pd.read_csv("results/S1.5d/S1.5d_real_intervals.csv"), Q90=pd.read_csv(q90_file),
                RR=pd.read_csv("results/S1.5/S1.5_real_refs_roi.csv"), RT=RT[RT.level == "type"], DR=DR,
                VS=pd.read_csv("results/S1.5e/S1.5e_real_void_shape.csv"), q90_file=q90_file)

def load_bench(res):
    """Tabelle per unita' di un banco (stessi nomi di file di S1.5e con il prefisso della cartella)."""
    p = res.rstrip("/").split("/")[-1]
    UM = rd(f"{res}/{p}_unit_metrics.csv"); UM[["cell_type", "bin"]] = UM[["cell_type", "bin"]].fillna("")
    return dict(U=pd.read_csv(f"{res}/{p}_panel_units.csv"), T=pd.read_csv(f"{res}/{p}_panel_types.csv"), UM=UM,
                W=pd.read_csv(f"{res}/{p}_windows.csv"))

def roi_values(B, R, arm, seeds=None):
    """Righe id, archetype, unit, level ('roi' o 'seed'), key, value, lo, hi, source per il braccio `arm` (semi `seeds`, interi 1-10; None = tutti)."""
    U, T, UM, W = B["U"], B["T"], B["UM"], B["W"]
    U = U[U.arm == arm]; T = T[T.arm == arm]; UM = UM[UM.arm == arm]; W = W[W.arm == arm]
    if seeds is not None:
        ss = set(int(s) for s in seeds)
        U = U[U.seed.isin(ss)]; T = T[T.seed.isin(ss)]; UM = UM[UM.seed.isin(ss)]
        W = W[W.seed.map(lambda s: int(str(s).lstrip("s"))).isin(ss)]
    out = []
    def put(pid, A, unit, v, lo, hi, src, level="roi"):
        for k, x in v.items():
            out.append(dict(id=pid, archetype=A, unit=unit, level=level, key=str(k), value=float(x), lo=float(lo), hi=float(hi), source=src))
    ps = W.groupby(["archetype", "seed"]).s1.mean().reset_index()
    RR, RT, RI, RD_, Q90, DR, VS = R["RR"], R["RT"], R["RI"], R["RD"], R["Q90"], R["DR"], R["VS"]
    u = U.groupby(["archetype", "roi_id"])[UM_COLS].mean().reset_index().merge(RR, on=["archetype", "roi_id"], suffixes=("", "_real"))
    um = UM.groupby(["archetype", "roi_id", "metric", "cell_type", "bin"]).value.mean().reset_index()
    tt = T.groupby(["archetype", "roi_id", "cell_type"])[["median_area_nuc", "median_nc"]].mean().reset_index()
    for A in sorted(MAIN):
        s = u[u.archetype == A].set_index("roi_id"); r = RR[RR.archetype == A]; ua = um[um.archetype == A]
        g = lambda mt, t="", b="": ua[(ua.metric == mt) & (ua.cell_type == t) & (ua.bin == b)].set_index("roi_id").value
        lo, hi = CONTACT[A]; put("P-1", A, "tessuto", ps[ps.archetype == A].set_index("seed").s1, lo, hi, "B-127...B-132", "seed")
        put("P-2", A, "tessuto", s.cv_terr / s.cv_terr_real - 1, -0.10, 0.10, "B-2 di S1.3")
        put("P-3", A, "tessuto", s.median_nc, r.median_nc.min(), r.median_nc.max(), "P-2 di S1.5")
        ri = RI[(RI.archetype == A) & (RI.metric == "cv_local4")].iloc[0]
        put("P-7", A, "tessuto", g("cv_local4"), ri["min"], ri["max"], "fra ROI reali (S1.5c)")
        put("P-8a", A, "tessuto", s.median_area_nuc, r.median_area_nuc.min(), r.median_area_nuc.max(), "P-1 di S1.5b")
        if A == "A2": put("P-8c", A, "tessuto", s.frac_nc095, 0, 0.0287, "P-5 di S1.5b")
        put("P-9", A, "tessuto", s.median_ecc_T - s.median_ecc_T_real, -0.05, 0.05, "B-4 di S1.3")
        put("P-10", A, "tessuto", s.frac_trunc, 0, 0.05, "S1.4")
        for pid, mt in [("P-11", "H_rel_3.0"), ("P-12", "xcv2_50"), ("P-13", "frac_gt3")]:
            ri = RD_[(RD_.archetype == A) & (RD_.metric == mt)].iloc[0]
            put(pid, A, "tessuto", g(mt), ri["min"], ri["max"], "fra ROI vere (S1.5d)")
        d = DR[DR.archetype == A]; x = g("roi_density"); x = x[x.index.isin(d.roi_id)]
        put("P-15", A, "tessuto", x, d.est.min(), d.est.max(), "fra ROI vere (E1 di S1.5e)")
        if A != "A4":
            v = VS[VS.archetype == A]
            put("P-16", A, "tessuto", g("frac_elong_gt3"), v.frac_elong_gt3.min(), v.frac_elong_gt3.max(), "fra ROI vere (E4 di S1.5e)")
            put("P-17", A, "tessuto", g("n_comp"), v.n_comp.min(), v.n_comp.max(), "fra ROI vere (E4 di S1.5e)")
        for t in MAIN[A]:
            rt = RT[(RT.archetype == A) & (RT.unit == t)]; x = tt[(tt.archetype == A) & (tt.cell_type == t)].set_index("roi_id")
            put("P-4", A, t, x.median_nc, rt.median_nt.min(), rt.median_nt.max(), "P-7 di S1.5b")
            put("P-8b", A, t, x.median_area_nuc, rt.median_area_nuc.min(), rt.median_area_nuc.max(), "P-2 di S1.5b")
            q = Q90[(Q90.archetype == A) & (Q90.cell_type == t)]
            put("P-14", A, t, g("q90_an", t), q.q90_min.iloc[0], q.q90_max.iloc[0], R["q90_file"].split("/")[-1])
            ri = RI[(RI.archetype == A) & (RI.metric == "nn_rel") & (RI.cell_type == t)]
            if len(ri): put("P-5", A, t, g("nn_rel", t), ri["min"].iloc[0], ri["max"].iloc[0], "fra ROI reali (S1.5c)")
            for b in ["10-20", "40-80"]:
                ri = RI[(RI.archetype == A) & (RI.metric == "markconn") & (RI.cell_type == t) & (RI.bin == b)]
                if len(ri): put("P-6", A, f"{t} {b}", g("markconn", t, b), ri["min"].iloc[0], ri["max"].iloc[0], "fra ROI reali (S1.5c)")
    return pd.DataFrame(out)

def est(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return (v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v)) if len(v) > 1 else (np.nan, np.nan, len(v))

def judge(m, se, lo, hi):
    """Esito di S1.5e (nota D-6: nessun taglio a [0, 1]): confermata se stima +- 2 ES dentro [lo, hi], smentita se tutta fuori, altrimenti indeterminata."""
    if not np.isfinite(m): return "non valutabile"
    a, b = m - 2 * se, m + 2 * se
    if a >= lo and b <= hi: return "confermata"
    if b < lo or a > hi: return "smentita"
    return "indeterminata"

def verdicts(rv, arm_label):
    rows = []
    for (pid, A, unit), g in rv.groupby(["id", "archetype", "unit"], sort=False):
        m, se, n = est(g.value)
        rows.append(dict(id=pid, archetype=A, unit=unit, arm=arm_label, est=m, se=se, n_units=n, lo=g.lo.iloc[0], hi=g.hi.iloc[0],
                         verdict=judge(m, se, g.lo.iloc[0], g.hi.iloc[0]), source=g.source.iloc[0]))
    return pd.DataFrame(rows)

def paired(ra, rb, frac=0.25):
    """Differenza appaiata (per ROI, per seme in P-1) fra due insiemi di valori per ROI; margine = +- frac x larghezza dell'intervallo vero."""
    k = ["id", "archetype", "unit", "level", "key"]
    m = ra.merge(rb[k + ["value"]], on=k, suffixes=("_a", "_b"))
    rows = []
    for (pid, A, unit), g in m.groupby(["id", "archetype", "unit"], sort=False):
        d, se, n = est(g.value_a - g.value_b); w = g.hi.iloc[0] - g.lo.iloc[0]
        rows.append(dict(id=pid, archetype=A, unit=unit, n_pairs=n, diff=d, se=se, width=w, margin=frac * w,
                         verdict=judge(d, se, -frac * w, frac * w)))
    return pd.DataFrame(rows)
