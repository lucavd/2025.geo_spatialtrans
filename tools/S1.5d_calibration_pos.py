"""tools/S1.5d_calibration_pos.py - S1.5d: calibrazione dei campi dei vuoti e della densita' (definizione in docs/steps/S1.5d.md).
Stadi: grid (griglia pre-registrata, results/S1.5d/S1.5d_het_grid.csv) | metrics [worker] (metriche di calibrazione di ogni unita' di
tools/S1.5d_calibration_pos.R e dei nuclei veri, stessa funzione tools/S1.5d_metrics.calibration_metrics) | fit (perdita per piega, scelta per
famiglia RV / RD / RVD, bordi della griglia; uscite results/S1.5d/S1.5d_het_<famiglia>.csv, S1.5d_het_loss.csv e, se serve, S1.5d_het_grid_ext.csv).
Perdita della piega k: somma sulle ROI j != k (tutte per «pooled») e sugli 8 bersagli di (simulato_j - vero_j)^2 / varianza fra le 5 ROI vere
del bersaglio nel tessuto; simulato_j = media delle 2 repliche; combinazioni con cellule non piazzate in una ROI usata: escluse.
Uso: .venv/bin/python tools/S1.5d_calibration_pos.py grid|metrics [n]|fit"""
import sys, os, glob, importlib.util, itertools, numpy as np, pandas as pd
from multiprocessing import Pool
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
met = mod("met", "tools/S1.5d_metrics.py"); cal = mod("cal", "tools/S1.5b_calibration.py")
RES = "results/S1.5d"; OUT = "/mnt/micron/geo_spatialtrans/S1.5d/calpos"; MET = "/mnt/micron/geo_spatialtrans/S1.5d/calpos_met"
PHI = [0, 0.025, 0.05, 0.10, 0.15, 0.20]; LV = [8, 16, 32, 64]; SD = [0, 0.3, 0.6]; LD = [10, 40]
def combo_name(pv, lv, sd, ld): return f"v{pv:g}_lv{lv:g}_d{sd:g}_ld{ld:g}"
def grid():
    rows = []
    for pv in PHI:
        for lv in (LV if pv > 0 else [np.nan]):
            for sd in SD:
                for ld in (LD if sd > 0 else [np.nan]):
                    rows.append(dict(combo=combo_name(pv, lv if pv > 0 else 0, sd, ld if sd > 0 else 0), void_frac=pv, void_ell_um=lv, dens_sigma=sd, dens_ell_um=ld, archetypes=""))
    G = pd.DataFrame(rows); G.to_csv(f"{RES}/S1.5d_het_grid.csv", index=False); print(len(G), "combinazioni")
def one_unit(f):
    combo = f.split("/")[-2]; A, roi, r = os.path.basename(f)[:-8].split("_"); o = f"{MET}/{combo}/{A}_{roi}_{r}.csv"
    if os.path.exists(o): return o
    g = met.RoiGeom(A, roi); e = pd.read_parquet(f); m = met.calibration_metrics(e.x.values, e.y.values, g)
    d = pd.DataFrame([dict(combo=combo, archetype=A, roi_id=roi, rep=int(r[1:]), n=len(e), n_failed=int(e.n_failed.iloc[0]) if len(e) else -1,
                           n_exceed=int(e.n_exceed.iloc[0]) if len(e) else -1, **m)])
    os.makedirs(os.path.dirname(o), exist_ok=True); d.to_csv(o + ".tmp", index=False); os.replace(o + ".tmp", o); return o
def real_one(args):
    A, roi = args; g = met.RoiGeom(A, roi); _, d = cal.load(A); e = d[d.roi_id == roi]
    return dict(archetype=A, roi_id=roi, **met.calibration_metrics(e.x.values, e.y.values, g))
def metrics(n):
    R = met.rois_table()
    with Pool(n) as p:
        real = pd.DataFrame(p.map(real_one, list(zip(R.archetype, R.roi_id))))
        fs = sorted(glob.glob(f"{OUT}/*/*.parquet")); outs = p.map(one_unit, fs, chunksize=8)
    real.to_csv(f"{RES}/S1.5d_het_real.csv", index=False)
    pd.concat([pd.read_csv(o) for o in outs]).to_csv(f"{MET}/all_units.csv", index=False); print(len(outs), "unita'")
FAM = {"RV": lambda c: c.dens_sigma == 0, "RD": lambda c: c.void_frac == 0, "RVD": lambda c: np.ones(len(c), bool)}
def fit():
    real = pd.read_csv(f"{RES}/S1.5d_het_real.csv"); U = pd.read_csv(f"{MET}/all_units.csv")
    G = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{RES}/S1.5d_het_grid*.csv"))]).drop_duplicates("combo")
    T = met.CAL_TARGETS; S = U.groupby(["combo", "archetype", "roi_id"]).agg({**{t: "mean" for t in T}, "n_failed": "max", "rep": "count"}).reset_index()
    rows, best_rows, ext = [], {k: [] for k in FAM}, []
    for A, ra in real.groupby("archetype"):
        var = ra[T].var(ddof=1); rlist = sorted(ra.roi_id)
        for fold in [f"loo_{r}" for r in rlist] + ["pooled"]:
            use = rlist if fold == "pooled" else [r for r in rlist if r != fold[4:]]
            s = S[(S.archetype == A) & S.roi_id.isin(use)].merge(ra[["roi_id"] + T], on="roi_id", suffixes=("", "_real"))
            s["L"] = sum((s[t] - s[t + "_real"]) ** 2 / var[t] for t in T)
            L = s.groupby("combo").agg(L=("L", "sum"), n_roi=("roi_id", "nunique"), n_rep=("rep", "min"), failed=("n_failed", "max")).reset_index().merge(G, on="combo")
            L = L[(L.n_roi == len(use)) & (L.n_rep == 2) & (L.failed == 0)]
            L0 = float(L[(L.void_frac == 0) & (L.dens_sigma == 0)].L.iloc[0])
            for _, x in L.iterrows(): rows.append(dict(archetype=A, fold=fold, combo=x.combo, loss=x.L, loss_null=L0))
            for fam, sel in FAM.items():
                Lf = L[sel(L)].sort_values(["L", "void_frac", "dens_sigma", "combo"]); b = Lf.iloc[0]
                pv, lv, sd, ld = b.void_frac, b.void_ell_um, b.dens_sigma, b.dens_ell_um
                edge = []
                if pv > 0 and pv >= max(Lf.void_frac): edge.append("void_frac_hi")
                if pv > 0 and lv >= Lf.void_ell_um.max(): edge.append("void_ell_hi")
                if pv > 0 and lv <= Lf.void_ell_um.min(): edge.append("void_ell_lo")
                if sd > 0 and sd >= Lf.dens_sigma.max(): edge.append("dens_sigma_hi")
                if sd > 0 and ld >= Lf.dens_ell_um.max(): edge.append("dens_ell_hi")
                if sd > 0 and ld <= Lf.dens_ell_um.min(): edge.append("dens_ell_lo")
                best_rows[fam].append(dict(archetype=A, fold=fold, combo=b.combo, void_frac=pv, void_ell_um=lv if pv > 0 else np.nan, dens_sigma=sd,
                                           dens_ell_um=ld if sd > 0 else np.nan, loss=b.L, loss_null=L0, edge="|".join(edge)))
                if fold == "pooled" and edge: ext.append((A, fam, edge, b))
    pd.DataFrame(rows).to_csv(f"{RES}/S1.5d_het_loss.csv", index=False)
    for fam, b in best_rows.items(): pd.DataFrame(b).to_csv(f"{RES}/S1.5d_het_{fam}.csv", index=False)
    new = []
    if os.path.exists(f"{RES}/S1.5d_het_grid_ext.csv"):                    # la griglia si estende una sola volta (pre-registrazione)
        print("griglia gia' estesa: ottimi al bordo dopo l'estensione:", sorted({(A, fam, '|'.join(e)) for A, fam, e, b in ext})); ext = []
    for A, fam, edge, b in ext:                                        # estensione di un passo nella direzione del bordo, per tutte le pieghe del tessuto
        pv = sorted(set(G.void_frac[G.void_frac > 0])); lv = sorted(set(G.void_ell_um.dropna())); sd = sorted(set(G.dens_sigma[G.dens_sigma > 0])); ld = sorted(set(G.dens_ell_um.dropna()))
        if "void_frac_hi" in edge: pv = pv + [round(max(pv) + 0.05, 3)]
        if "void_ell_hi" in edge: lv = lv + [max(lv) * 2]
        if "void_ell_lo" in edge: lv = [min(lv) / 2] + lv
        if "dens_sigma_hi" in edge: sd = sd + [round(max(sd) + 0.3, 3)]
        if "dens_ell_hi" in edge: ld = ld + [max(ld) * 2]
        if "dens_ell_lo" in edge: ld = [min(ld) / 2] + ld
        vv = [(0, 0)] + [(a, l) for a in pv for l in lv]; dd = [(0, 0)] + [(a, l) for a in sd for l in ld]
        if fam == "RV": dd = [(0, 0)]
        if fam == "RD": vv = [(0, 0)]
        for (a, l), (s_, m_) in itertools.product(vv, dd):
            c = combo_name(a, l, s_, m_)
            if c not in set(G.combo): new.append(dict(combo=c, void_frac=a, void_ell_um=l if a > 0 else np.nan, dens_sigma=s_, dens_ell_um=m_ if s_ > 0 else np.nan, archetypes=A))
    if new:
        N = pd.DataFrame(new).groupby("combo").agg({"void_frac": "first", "void_ell_um": "first", "dens_sigma": "first", "dens_ell_um": "first", "archetypes": lambda v: "|".join(sorted(set(v)))}).reset_index()
        N.to_csv(f"{RES}/S1.5d_het_grid_ext.csv", index=False); print("estensione:", len(N), "combinazioni per", sorted({a for a, *_ in ext}))
    elif not os.path.exists(f"{RES}/S1.5d_het_grid_ext.csv"): print("nessuna estensione")
    for fam in FAM: print(fam); print(pd.read_csv(f"{RES}/S1.5d_het_{fam}.csv").query("fold == 'pooled'").to_string(index=False))
if __name__ == "__main__":
    st = sys.argv[1]
    if st == "grid": grid()
    elif st == "metrics": metrics(int(sys.argv[2]) if len(sys.argv) > 2 else 10)
    elif st == "fit": fit()
