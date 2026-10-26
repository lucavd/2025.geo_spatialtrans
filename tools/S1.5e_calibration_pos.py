"""tools/S1.5e_calibration_pos.py - S1.5e: calibrazione dei vuoti che tolgono cellule (X isotropi, XA allungati) e della densita' con il tetto (XAD)
(definizioni: docs/steps/S1.5e.md, «Cosa si implementa»). Stadi:
  grid                  griglie X e XA (results/S1.5e/S1.5e_het_grid.csv)
  metrics [worker]      metriche di calibrazione di ogni unita' di tools/S1.5e_calibration_pos.R e dei nuclei veri (tools/S1.5e_metrics.calibration_metrics_e)
  fit                   perdita per piega, ottimo per X e XA, bordi; estensione della griglia una volta (S1.5e_het_grid_ext.csv); griglia di XAD con i vuoti
                        dell'ottimo di XA della stessa piega (S1.5e_het_grid_xad.csv); se le unita' di XAD ci sono, ottimo di XAD
Perdita della piega k: somma sulle ROI j != k (tutte per «pooled») e sui 9 bersagli (8 di S1.5d + allungamento dei vuoti; nel linfonodo senza
allungamento) di (simulato_j - vero_j)^2 / varianza fra le 5 ROI vere del bersaglio; simulato_j = media delle 2 repliche; combinazioni con cellule
non piazzate in una ROI usata: escluse. Ordine dei pareggi: perdita, poi meno vuoti, meno allungamento, meno densita', nome.
Uscite: results/S1.5e/S1.5e_het_<X|XA|XAD>.csv, S1.5e_het_loss.csv, S1.5e_het_real.csv.
Uso: .venv/bin/python tools/S1.5e_calibration_pos.py grid|metrics [n]|fit"""
import sys, os, glob, importlib.util, itertools, numpy as np, pandas as pd
from multiprocessing import Pool
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
met = mod("met_e", "tools/S1.5e_metrics.py"); cal = mod("cal", "tools/S1.5b_calibration.py")
RES = "results/S1.5e"; OUT = "/mnt/micron/geo_spatialtrans/S1.5e/calpos"; MET = "/mnt/micron/geo_spatialtrans/S1.5e/calpos_met"
PHI = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]; LV = [8, 16, 32, 64]; KAP = [2, 4, 8]; SD = [0.3, 0.6, 0.9]; LD = [10, 40, 80]
def cname(pv, lv, k, sd=0, ld=0, tag=""): return f"{tag}rm_v{pv:g}_lv{lv:g}_k{k:g}_d{sd:g}_ld{ld:g}"
def row(pv, lv, k, sd=0.0, ld=np.nan, arch="", rois="", tag=""):
    return dict(combo=cname(pv, lv, k, sd, 0 if not np.isfinite(ld) else ld, tag), void_frac=pv, void_ell_um=lv, void_aspect=k, void_mode="remove",
                dens_sigma=sd, dens_ell_um=ld, dens_cap=bool(sd > 0), archetypes=arch, rois=rois)
def grid():
    G = pd.DataFrame([row(pv, lv, k) for pv in PHI for lv in LV for k in [1] + KAP]); os.makedirs(RES, exist_ok=True)
    G.to_csv(f"{RES}/S1.5e_het_grid.csv", index=False); print(len(G), "combinazioni (X:", int((G.void_aspect == 1).sum()), ", XA:", int((G.void_aspect > 1).sum()), ")")
def one_unit(f):
    combo = f.split("/")[-2]; A, roi, r = os.path.basename(f)[:-8].split("_"); o = f"{MET}/{combo}/{A}_{roi}_{r}.csv"
    if os.path.exists(o): return o
    g = met.RoiGeomE(A, roi); e = pd.read_parquet(f); m = met.calibration_metrics_e(e.x.values, e.y.values, g)
    d = pd.DataFrame([dict(combo=combo, archetype=A, roi_id=roi, rep=int(r[1:]), n=len(e), n_failed=int(e.n_failed.iloc[0]) if len(e) else -1,
                           n_exceed=int(e.n_exceed.iloc[0]) if len(e) else -1, **m)])
    os.makedirs(os.path.dirname(o), exist_ok=True); d.to_csv(o + ".tmp", index=False); os.replace(o + ".tmp", o); return o
def real_one(args):
    A, roi = args; g = met.RoiGeomE(A, roi); _, d = cal.load(A); e = d[d.roi_id == roi]
    return dict(archetype=A, roi_id=roi, **met.calibration_metrics_e(e.x.values, e.y.values, g))
def metrics(n):
    R = met.m5d.rois_table()
    with Pool(n) as p:
        real = pd.DataFrame(p.map(real_one, list(zip(R.archetype, R.roi_id))))
        fs = sorted(glob.glob(f"{OUT}/*/*.parquet")); outs = p.map(one_unit, fs, chunksize=8)
    real.to_csv(f"{RES}/S1.5e_het_real.csv", index=False)
    pd.concat([pd.read_csv(o) for o in outs]).to_csv(f"{MET}/all_units.csv", index=False); print(len(outs), "unita'")
def targets(A): return [t for t in met.CAL_TARGETS_E if not (A == "A4" and t == "elong_wmed")]
def grids():
    fs = [p for p in sorted(glob.glob(f"{RES}/S1.5e_het_grid*.csv"))]
    G = pd.concat([pd.read_csv(p, keep_default_na=False, na_values=[""]) for p in fs]).drop_duplicates("combo"); G["archetypes"] = G.archetypes.fillna(""); G["rois"] = G.rois.fillna("")
    return G
def losses(real, S, G, A, fold, fam_sel):
    ra = real[real.archetype == A]; T = targets(A); var = ra[T].var(ddof=1); rlist = sorted(ra.roi_id)
    use = rlist if fold == "pooled" else [r for r in rlist if r != fold[4:]]
    s = S[(S.archetype == A) & S.roi_id.isin(use)].merge(ra[["roi_id"] + T], on="roi_id", suffixes=("", "_real"))
    s["L"] = sum((s[t] - s[t + "_real"]) ** 2 / var[t] for t in T)
    L = s.groupby("combo").agg(L=("L", "sum"), n_roi=("roi_id", "nunique"), n_rep=("rep", "min"), failed=("n_failed", "max")).reset_index().merge(G, on="combo")
    L = L[(L.n_roi == len(use)) & (L.n_rep == 2) & (L.failed == 0) & fam_sel(L, A, fold)]
    return L.sort_values(["L", "void_frac", "void_aspect", "dens_sigma", "combo"])
def fam_X(L, A, fold): return (L.void_aspect == 1) & (L.dens_sigma == 0) & (L.archetypes.isin(["", A]) | L.archetypes.str.contains(A)) & (L.rois == "")
def fam_XA(L, A, fold): return (L.void_aspect > 1) & (L.dens_sigma == 0) & (L.archetypes.isin(["", A]) | L.archetypes.str.contains(A)) & (L.rois == "")
def edges(b, Lf, fam):
    e = []
    # XAD: i vuoti sono fissati a quelli di XA della piega (un solo valore), quindi i bordi dei vuoti non si valutano
    if fam != "XAD" and b.void_frac >= Lf.void_frac.max(): e.append("void_frac_hi")
    if fam != "XAD" and b.void_ell_um >= Lf.void_ell_um.max(): e.append("void_ell_hi")
    if fam != "XAD" and b.void_ell_um <= Lf.void_ell_um.min(): e.append("void_ell_lo")
    if fam == "XA" and b.void_aspect >= Lf.void_aspect.max(): e.append("void_aspect_hi")
    if fam == "XAD" and b.dens_sigma > 0:
        if b.dens_sigma >= Lf.dens_sigma.max(): e.append("dens_sigma_hi")
        if b.dens_ell_um >= Lf.dens_ell_um.max(): e.append("dens_ell_hi")
        if b.dens_ell_um <= Lf[Lf.dens_sigma > 0].dens_ell_um.min(): e.append("dens_ell_lo")
    return "|".join(e)
COLS = ["combo", "void_frac", "void_ell_um", "void_aspect", "void_mode", "dens_sigma", "dens_ell_um", "dens_cap"]
def fit():
    real = pd.read_csv(f"{RES}/S1.5e_het_real.csv"); U = pd.read_csv(f"{MET}/all_units.csv"); G = grids()
    T = met.CAL_TARGETS_E; S = U.groupby(["combo", "archetype", "roi_id"]).agg({**{t: "mean" for t in T}, "n_failed": "max", "rep": "count"}).reset_index()
    best = {"X": [], "XA": [], "XAD": [], "XADs": []}; lrows = []; ext = []; xad = []
    folds = lambda A: [f"loo_{r}" for r in sorted(real[real.archetype == A].roi_id)] + ["pooled"]
    for A in sorted(real.archetype.unique()):
        for fold in folds(A):
            for fam, sel in [("X", fam_X), ("XA", fam_XA)]:
                Lf = losses(real, S, G, A, fold, sel); b = Lf.iloc[0]; ed = edges(b, Lf, fam)
                for _, x in Lf.iterrows(): lrows.append(dict(family=fam, archetype=A, fold=fold, combo=x.combo, loss=x.L))
                best[fam].append(dict(archetype=A, fold=fold, **{c: b[c] for c in COLS}, loss=b.L, edge=ed))
                if fold == "pooled" and ed: ext.append((A, fam, ed))
                if fam == "XA":                                          # griglia di XAD: vuoti dell'ottimo di XA della stessa piega
                    use = sorted(real[real.archetype == A].roi_id) if fold == "pooled" else [r for r in sorted(real[real.archetype == A].roi_id) if r != fold[4:]]
                    for sd, ld in itertools.product(SD, LD):
                        xad.append(row(b.void_frac, b.void_ell_um, b.void_aspect, sd, ld, A, "|".join(use), tag=f"xad_{A}_{fold}_"))
    pd.DataFrame(lrows).to_csv(f"{RES}/S1.5e_het_loss.csv", index=False)
    for fam in ["X", "XA"]: pd.DataFrame(best[fam]).to_csv(f"{RES}/S1.5e_het_{fam}.csv", index=False)
    if os.path.exists(f"{RES}/S1.5e_het_grid_ext.csv"):
        print("griglia gia' estesa; ottimi al bordo dopo l'estensione:", sorted(set(ext))); ext = []
    new = []
    for A, fam, ed in ext:                                                # un passo nella direzione del bordo, tutte le pieghe del tessuto
        pv = sorted(set(G.void_frac)); lv = sorted(set(G.void_ell_um)); kk = sorted(set(G.void_aspect[G.void_aspect > 1]))
        if "void_frac_hi" in ed: pv = pv + [round(max(pv) + 0.05, 3)]
        if "void_ell_hi" in ed: lv = lv + [max(lv) * 2]
        if "void_ell_lo" in ed: lv = [min(lv) / 2] + lv
        if "void_aspect_hi" in ed: kk = kk + [max(kk) * 2]
        for a, l, k in itertools.product(pv, lv, [1] if fam == "X" else kk):
            r_ = row(a, l, k, arch=A)
            if r_["combo"] not in set(G.combo): new.append(r_)
    if new:
        N = pd.DataFrame(new).groupby("combo").agg({**{c: "first" for c in COLS[1:]}, "archetypes": lambda v: "|".join(sorted(set(v))), "rois": "first"}).reset_index()
        N.to_csv(f"{RES}/S1.5e_het_grid_ext.csv", index=False); print("estensione:", len(N), "combinazioni per", sorted({a for a, *_ in ext})); return
    X = pd.DataFrame(xad).drop_duplicates("combo"); X.to_csv(f"{RES}/S1.5e_het_grid_xad.csv", index=False); print("griglia di XAD:", len(X), "combinazioni")
    if set(X.combo) <= set(S.combo):                                      # unita' di XAD presenti: ottimo di XAD (sigma 0 = XA della piega)
        G = grids(); XA = pd.DataFrame(best["XA"])
        for A in sorted(real.archetype.unique()):
            for fold in folds(A):
                xa = XA[(XA.archetype == A) & (XA.fold == fold)].iloc[0]
                sel = lambda L, A_, f_, xa=xa: ((L.combo == xa.combo) | L.combo.str.startswith(f"xad_{A_}_{f_}_"))
                Lf = losses(real, S, G, A, fold, sel); b = Lf.iloc[0]; ed = edges(b, Lf, "XAD")
                for _, x in Lf.iterrows(): lrows.append(dict(family="XAD", archetype=A, fold=fold, combo=x.combo, loss=x.L))
                best["XAD"].append(dict(archetype=A, fold=fold, **{c: b[c] for c in COLS}, loss=b.L, edge=ed))
        pd.DataFrame(best["XAD"]).to_csv(f"{RES}/S1.5e_het_XAD.csv", index=False); pd.DataFrame(lrows).to_csv(f"{RES}/S1.5e_het_loss.csv", index=False)
        bs = []                                                          # XADs (nota D-5): come pre-registrato, solo sigma > 0 (XAD qui sopra ammette anche sigma = 0)
        for A in sorted(real.archetype.unique()):
            for fold in folds(A):
                Lf = losses(real, S, G, A, fold, lambda L, A_, f_: L.combo.str.startswith(f"xad_{A_}_{f_}_")); b = Lf.iloc[0]
                bs.append(dict(archetype=A, fold=fold, **{c: b[c] for c in COLS}, loss=b.L, edge=edges(b, Lf, "XAD")))
        pd.DataFrame(bs).to_csv(f"{RES}/S1.5e_het_XADs.csv", index=False); best["XADs"] = bs
        xb = pd.DataFrame(best["XAD"]); xe = xb[(xb.fold == "pooled") & (xb.edge != "")]
        if len(xe) and not os.path.exists(f"{RES}/S1.5e_het_grid_xadext.csv"):   # un passo nella direzione del bordo della densita', tutte le pieghe del tessuto
            new = []
            for _, b in xe.iterrows():
                A = b.archetype; sd = list(SD); ld = list(LD)
                if "dens_sigma_hi" in b.edge: sd = sd + [round(max(sd) + 0.3, 3)]
                if "dens_ell_hi" in b.edge: ld = ld + [max(ld) * 2]
                if "dens_ell_lo" in b.edge: ld = [min(ld) / 2] + ld
                for fold in folds(A):
                    xa = XA[(XA.archetype == A) & (XA.fold == fold)].iloc[0]
                    use = sorted(real[real.archetype == A].roi_id) if fold == "pooled" else [r for r in sorted(real[real.archetype == A].roi_id) if r != fold[4:]]
                    for a_, l_ in itertools.product(sd, ld):
                        r_ = row(xa.void_frac, xa.void_ell_um, xa.void_aspect, a_, l_, A, "|".join(use), tag=f"xad_{A}_{fold}_")
                        if r_["combo"] not in set(G.combo): new.append(r_)
            pd.DataFrame(new).to_csv(f"{RES}/S1.5e_het_grid_xadext.csv", index=False); print("estensione di XAD:", len(new), "combinazioni per", sorted(set(xe.archetype)))
        elif len(xe): print("griglia di XAD gia' estesa; ottimi al bordo:", list(zip(xe.archetype, xe.edge)))
    for fam in best:
        if best[fam]: print(fam); print(pd.DataFrame(best[fam]).query("fold == 'pooled'")[["archetype"] + COLS[1:] + ["loss", "edge"]].to_string(index=False))
if __name__ == "__main__":
    st = sys.argv[1]
    if st == "grid": grid()
    elif st == "metrics": metrics(int(sys.argv[2]) if len(sys.argv) > 2 else 10)
    elif st == "fit": fit()
