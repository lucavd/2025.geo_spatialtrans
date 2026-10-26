"""tools/S1.5d_posthoc_p8c.py - S1.5d, post hoc PH-4 (dopo il banco, non pre-registrato): frazione di cellule quasi senza citoplasma
(nucleo >= 0.95 del territorio, celle interne) per ROI vera del tumore del colon, stessa definizione di P-8c (massimo = 0.0287 di S1.5b P-5),
accanto alle stime dei bracci. Ingressi: nuclei veri con territorio (R6, via tools/S1.5b_calibration.py), results/S1.5d/S1.5d_verdicts.csv.
Uscita: results/S1.5d/S1.5d_posthoc_p8c_real.csv. Uso: .venv/bin/python -B tools/S1.5d_posthoc_p8c.py"""
import importlib.util, numpy as np, pandas as pd
s = importlib.util.spec_from_file_location("cal", "tools/S1.5b_calibration.py"); cal = importlib.util.module_from_spec(s); s.loader.exec_module(cal)
_, d = cal.load("A2"); e = d[d.interior.astype(bool) & d.a_t.notna() & d.a_n.notna()]
r = e.groupby("roi_id").apply(lambda g: pd.Series({"n_interior": len(g), "value": float(np.mean(g.a_n / g.a_t >= 0.95))})).reset_index()
r.insert(0, "source", "vero")
V = pd.read_csv("results/S1.5d/S1.5d_verdicts.csv"); v = V[V.id == "P-8c"][["arm", "est", "se", "verdict"]].rename(columns={"arm": "source", "est": "value"})
v["roi_id"] = "media 5 ROI"; out = pd.concat([r, v], ignore_index=True)
out["below_real_min"] = out.value < r.value.min(); out["above_real_max"] = out.value > r.value.max()
out.to_csv("results/S1.5d/S1.5d_posthoc_p8c_real.csv", index=False); print(out.round(4).to_string(index=False))
