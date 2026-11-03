"""tools/R5_ph2_synth_summary.py — riassunto della controprova sintetica di PH-2 (tools/R5_ph2_synth.py): sigma_in ritrovato
(interpolazione lineare dove previsto - osservato cambia segno) e r a quel sigma, contro i valori veri.
Uso: .venv/bin/python tools/R5_ph2_synth_summary.py <outdir> -> <outdir>/tables/R5_PH2_synth.csv"""
import sys, json, glob, numpy as np, pandas as pd
OUT = sys.argv[1]
D = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/ph2_synth/*.json"))]); D["d"] = D.pred_far - D.obs_far
rows = []
for (t, st), g in D.sort_values("sigma_in").groupby(["tag", "sigma_true"]):
    g = g[np.isfinite(g.d)]; s, d = g.sigma_in.values, g.d.values; star = np.nan
    for i in range(len(s) - 1):
        if d[i] * d[i + 1] < 0: star = s[i] + (s[i + 1] - s[i]) * (-d[i]) / (d[i + 1] - d[i]); break
    rows.append({"tag": t, "sigma_true": st, "r_true": g.r_true.iloc[0], "obs_star": g.obs_star.iloc[0], "obs_far": g.obs_far.iloc[0],
                 "sigma_in_recovered": star, "r_recovered": float(np.interp(star, s, g.r.values)) if np.isfinite(star) else np.nan})
S = pd.DataFrame(rows); S.to_csv(f"{OUT}/tables/R5_PH2_synth.csv", index=False); print(S.round(3).to_string(index=False))
