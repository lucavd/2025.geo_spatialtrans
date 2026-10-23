"""tools/S1.5c_real_refs.py - S1.5c: valori reali per ROI delle metriche nuove (nn_rel, markconn, cv_local4), calcolati PRIMA di ogni simulazione.
Nuclei veri con tipo (etichette R6 PH-1, insieme di territori esportato da R6, come S1.5b). Uso: .venv/bin/python tools/S1.5c_real_refs.py
Uscita: results/S1.5c/S1.5c_real_refs.csv (una riga per ROI, metrica, tipo, bin) e S1.5c_real_intervals.csv ([min, max] dei 5 ROI)."""
import importlib.util, numpy as np, pandas as pd
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py"); met = mod("met", "tools/S1.5c_metrics.py")
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
rows = []
for A in cal.ARCH:
    _, d = cal.load(A)
    for r, e in d.groupby("roi_id"):
        o = met.roi_metrics(e.x.values, e.y.values, e.type.values, e.interior.values, MAIN[A])
        for k, v in o.items():
            if isinstance(k, tuple): rows.append(dict(archetype=A, roi_id=r, metric=k[0], cell_type=k[1], bin=k[2] if len(k) > 2 else "", value=v))
            else: rows.append(dict(archetype=A, roi_id=r, metric=k, cell_type="", bin="", value=v))
R = pd.DataFrame(rows); R.to_csv("results/S1.5c/S1.5c_real_refs.csv", index=False)
I = R.groupby(["archetype", "metric", "cell_type", "bin"], dropna=False).value.agg(["min", "max", "mean", "count"]).reset_index()
I.to_csv("results/S1.5c/S1.5c_real_intervals.csv", index=False); print(I.to_string())
