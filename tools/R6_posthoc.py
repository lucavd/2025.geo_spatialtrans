"""R6 post hoc (declared after the results): PH-2 chance-corrected agreement (kappa) for CP-2 and P-6; PH-3 observed vs permutation null
for CP-1; PH-4 total epithelium in colon; PH-5 lymph node fractions among resolved nuclei; PH-6 log-area dispersion per type."""
import sys, numpy as np, pandas as pd
from sklearn.metrics import cohen_kappa_score
sys.path.insert(0, "tools")
import importlib.util; spec = importlib.util.spec_from_file_location("an", "tools/R6_analysis.py")
M = "/mnt/micron/geo_spatialtrans/R6"; OUT = "results/R6"; NON = {"non classificabile", "non risolto"}
rng = np.random.default_rng(20261022); B = 2000
# re-use loaders without running the analysis: copy the minimal functions
exec(open("tools/R6_analysis.py").read().split("# ================= main")[0].replace('VARIANT = sys.argv[1] if len(sys.argv) > 1 else "prereg"', 'VARIANT = "ph1"'))
rows = []
for V, suf in [("prereg", ""), ("ph1", "_ph1")]:
    SUF = suf
    for A in ["A5", "A6"]:
        lab = add_sr(load_labels(A, "exp", 20, 20261022), A); sh = add_sr(load_labels(A, "shift", 20, 20261022), A)
        for name, d in [("primary", lab), ("shifted", sh)]:
            d = d[d.sr_cell_id.notna() & d.sr_type.notna() & ~d.sr_lowumi.fillna(True).astype(bool) & (d.sr_type != "altro") & ~d.cell_type.isin(NON)]
            p = d.cell_type.value_counts(normalize=True); q = d.sr_type.value_counts(normalize=True); chance = sum(p.get(t, 0) * q.get(t, 0) for t in set(p.index) | set(q.index))
            rows.append(dict(id="PH-2", variant=V, archetype=A, footprint=name, n=len(d), agreement=(d.cell_type == d.sr_type).mean(), chance_agreement=chance, kappa=cohen_kappa_score(d.cell_type, d.sr_type)))
ph2 = pd.DataFrame(rows); ph2.to_csv(f"{OUT}/R6_posthoc_kappa.csv", index=False); print(ph2.round(3).to_string())
# PH-3: observed vs permutation null (from predictions and counterproofs files)
out = []
for V, suf in [("prereg", ""), ("ph1", "_ph1")]:
    P = pd.read_csv(f"{OUT}/R6_predictions{suf}.csv"); C = pd.read_csv(f"{OUT}/R6_counterproofs{suf}.csv"); C = C[C.id == "CP-1"]
    for _, c in C.iterrows():
        if "max/min" in c.quantity: obs = P[(P.id == "P-2") & (P.archetype == c.archetype)].iloc[0]
        elif "neuron" in c.quantity: obs = P[(P.id == "P-1")].iloc[0]
        else: obs = P[(P.id == "P-3")].iloc[0]
        out.append(dict(id="PH-3", variant=V, archetype=c.archetype, quantity=c.quantity, observed=obs.estimate, observed_lo=obs.lo, observed_hi=obs.hi, null_median=c.median_null, null_q975=c.q975_null,
                        beyond_null=(obs.lo > c.q975_null) if "ecc" not in c.quantity else (abs(obs.estimate) > c.q975_null)))
ph3 = pd.DataFrame(out); ph3.to_csv(f"{OUT}/R6_posthoc_null.csv", index=False); print(ph3.round(3).to_string())
# PH-4 / PH-5 / PH-6 on PH-1 labels
SUF = "_ph1"; res = []
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    lab = load_labels(A, "exp", 20, 20261022)
    if A in ("A2", "A3"):
        e = boot_stat(lab, lambda d: d.cell_type.str.startswith("epiteliale").mean(), two_stage=False); res.append(dict(id="PH-4", archetype=A, quantity="fraction epithelial (tumorale + normale)", estimate=e[0], lo=e[1], hi=e[2]))
    if A == "A4":
        r = lab[~lab.cell_type.isin(NON)]
        for ty in ["linfociti B", "linfociti T"]:
            e = boot_stat(r, lambda d, ty=ty: (d.cell_type == ty).mean(), two_stage=False); res.append(dict(id="PH-5", archetype=A, quantity=f"fraction {ty} among resolved nuclei", estimate=e[0], lo=e[1], hi=e[2]))
            e = boot_stat(r, lambda d, ty=ty: (d.cell_type == ty).mean(), two_stage=False, rois=["f1", "f2"]); res.append(dict(id="PH-5", archetype=A, quantity=f"fraction {ty} among resolved, follicular ROIs", estimate=e[0], lo=e[1], hi=e[2]))
            e = boot_stat(r, lambda d, ty=ty: (d.cell_type == ty).mean(), two_stage=False, rois=["p1", "p2", "p3"]); res.append(dict(id="PH-5", archetype=A, quantity=f"fraction {ty} among resolved, paracortical ROIs", estimate=e[0], lo=e[1], hi=e[2]))
        for roi, g in lab.groupby("roi_id"): res.append(dict(id="PH-5", archetype=A, quantity=f"ROI {roi}: fraction unresolved / B / T (all nuclei)", estimate=(g.cell_type == "non risolto").mean(), lo=(g.cell_type == "linfociti B").mean(), hi=(g.cell_type == "linfociti T").mean()))
    if A == "A5":
        g = lab[lab.cell_type.isin(["astrociti", "oligodendrociti", "microglia"])]; res.append(dict(id="PH-1 detail", archetype=A, quantity="median nuclear area glia pooled (um2)", estimate=g.area_um2.median(), lo=np.nan, hi=np.nan))
        for roi, gg in lab.groupby("roi_id"): res.append(dict(id="PH-5", archetype=A, quantity=f"ROI {roi}: fraction neurons / oligodendrocytes / astrocytes", estimate=(gg.cell_type == "neuroni").mean(), lo=(gg.cell_type == "oligodendrociti").mean(), hi=(gg.cell_type == "astrociti").mean()))
    for ty, g in lab[~lab.cell_type.isin(NON)].groupby("cell_type"):
        if len(g) >= 50:
            la = np.log(g.area_um2); res.append(dict(id="PH-6", archetype=A, quantity=f"{ty}: sigma of log nuclear area; q01, q99 of area (um2)", estimate=la.std(), lo=g.area_um2.quantile(.01), hi=g.area_um2.quantile(.99)))
ph = pd.DataFrame(res); ph.to_csv(f"{OUT}/R6_posthoc.csv", index=False); print(ph.round(3).to_string())
