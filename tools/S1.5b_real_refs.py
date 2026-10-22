# tools/S1.5b_real_refs.py - S1.5b: valori reali per ROI, per tipo e per gruppo di taglia (calcolati prima di ogni simulazione).
# Insieme di valutazione = celle interne appaiate con area del nucleo > 0 (insieme di calibrazione di S1.4/R3b), unione per coordinate.
# Soglia del cervello: media geometrica delle aree mediane di neuroni e glia (gruppo di taglia), arrotondata a 0.1 µm²; banda = soglia /÷ sqrt(2).
# Uscita: results/S1.5b/S1.5b_real_refs_roi_type.csv, results/S1.5b/S1.5b_real_refs_roi.csv
import numpy as np, pandas as pd, sys
import importlib.util
spec = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(spec); spec.loader.exec_module(cal)
RES = "results/S1.5b"
C = pd.read_csv(f"{RES}/S1.5b_type_calibration.csv")
rows, rroi = [], []
thr = None
for A in cal.ARCH:
    lab, d = cal.load(A)
    ev = d[d.interior & d.matched & (d.a_n > 0)].copy()
    g = C[C.archetype == A].set_index("cell_type").size_group
    ev["size_group"] = ev.type.map(g)          # NaN per non risolto / non classificabile
    ev["nt"] = ev.a_n / ev.a_t
    if A == "A5":
        mn = ev.a_n[ev.type == "neuroni"].median(); mg = ev.a_n[ev.size_group == g["astrociti"]].median()
        thr = round(float(np.sqrt(mn * mg)), 1)
    for r, e in ev.groupby("roi_id"):
        rec = dict(archetype=A, roi_id=r, n_eval=len(e), median_area_nuc=e.a_n.median(), median_nt=e.nt.median(),
                   frac_nt095=float((e.nt >= 0.95).mean()))
        if A == "A5":
            rec.update(small_threshold_um2=thr, frac_small=float((e.a_n < thr).mean()),
                       frac_band=float(((e.a_n > thr / np.sqrt(2)) & (e.a_n < thr * np.sqrt(2))).mean()),
                       ratio_neuron_glia=e.a_n[e.type == "neuroni"].median() / e.a_n[e.size_group == g["astrociti"]].median())
        rroi.append(rec)
        units = [("type", t, x) for t, x in e.groupby("type")] + [("size_group", int(s), x) for s, x in e.groupby("size_group")]
        for k, v, x in units:
            rows.append(dict(archetype=A, roi_id=r, level=k, unit=str(v), n=len(x), median_area_nuc=x.a_n.median(), median_nt=x.nt.median(),
                             median_a_t=x.a_t.median(), frac_nt095=float((x.nt >= 0.95).mean())))
pd.DataFrame(rows).to_csv(f"{RES}/S1.5b_real_refs_roi_type.csv", index=False)
pd.DataFrame(rroi).to_csv(f"{RES}/S1.5b_real_refs_roi.csv", index=False)
T = pd.DataFrame(rows); T = T[T.unit.isin(C.cell_type.astype(str)) | (T.level == "size_group")]
s = T.groupby(["archetype", "level", "unit"]).agg(n_roi=("roi_id", "size"), area_min=("median_area_nuc", "min"), area_max=("median_area_nuc", "max"),
                                                  nt_min=("median_nt", "min"), nt_max=("median_nt", "max")).reset_index()
pd.set_option("display.width", 250); print(s.round(4).to_string(index=False))
print(pd.DataFrame(rroi).round(4).to_string(index=False))
