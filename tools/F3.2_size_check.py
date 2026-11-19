#!/usr/bin/env python3
"""tools/F3.2_size_check.py - F3.2: le misure dei metodi cambiano con il lato dell'area (1 mm contro 500 um)? Solo dati di taratura (F3.1, aree vere).
Per ogni area vera di R2: M-4 (singoletti), M-5 (concordanza con il nucleo), M-6 (coerenza spaziale) e UMI per blocco ricalcolati dai blocchi di RCTD
di F3.1 sull'area intera e sui 4 quarti di 500 um, con le definizioni di tools/F3.1_rctd.R. Scarto di dimensione = media dei 4 quarti meno area intera.
Regola decisa prima (la stessa che Luca ha scelto per la versione di Space Ranger, DEC-091): per tessuto e misura, |scarto medio| + 2 ES (fra aree)
entro un quarto della deviazione standard fra le 5 aree intere vere. Se regge, le aree di 500 um si confrontano con aree simulate di 500 um e con la
variabilita' fra quarti di R2.
Uso: .venv/bin/python tools/F3.2_size_check.py -> results/F3.2/F3.2_size_check_areas.csv, F3.2_size_check.csv
"""
import json, numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; M = Path("/mnt/micron/geo_spatialtrans/F3.1/main/real")
def meas(b):
    run = b.spot_class.notna(); ok = run & (b.spot_class != "reject")
    out = dict(n_blocks=len(b), M4=float((b.spot_class[ok] == "singlet").mean()) if ok.any() else np.nan, umi_block=float(b.umi.mean()))
    hn = ok & b.nuc_type.notna()
    out["M5"] = float((b["first"][hn] == b.nuc_type[hn]).mean()) if hn.any() else np.nan
    ft = {(r, c): f for r, c, f, o in zip(b.r8, b.c8, b["first"], ok) if o and isinstance(f, str)}
    same = tot = 0
    for (r, c), f in ft.items():
        for nb in ((r, c + 1), (r + 1, c)):
            if nb in ft: tot += 1; same += ft[nb] == f
    out["M6"] = same / tot if tot else np.nan
    return out
rows = []
for d in sorted(M.iterdir()):
    f = d / "rctd_blocks.csv"
    if not f.exists(): continue
    b = pd.read_csv(f); A = d.name.split("_")[0]
    rows.append(dict(area=d.name, A=A, part="whole", **meas(b)))
    for j in (0, 1):
        for i in (0, 1):
            q = b[(b.x >= i * 500) & (b.x < (i + 1) * 500) & (b.y >= j * 500) & (b.y < (j + 1) * 500)]
            rows.append(dict(area=d.name, A=A, part=f"q{j}{i}", **meas(q)))
t = pd.DataFrame(rows); t.to_csv(ROOT / "results/F3.2/F3.2_size_check_areas.csv", index=False)
print("x range", t.shape, flush=True)
out = []
for A, g in t.groupby("A"):
    for m in ["M4", "M5", "M6", "umi_block"]:
        w = g[g.part == "whole"].set_index("area")[m]; q = g[g.part != "whole"].groupby("area")[m].mean()
        diff = (q - w).dropna(); sd_w = w.std(ddof=1); est = diff.mean(); se = diff.std(ddof=1) / np.sqrt(len(diff))
        sd_q = g[g.part != "whole"][m].std(ddof=1)
        out.append(dict(A=A, measure=m, n_areas=len(diff), whole_mean=w.mean(), size_shift=est, se=se, sd_between_whole=sd_w, margin=sd_w / 4,
                        ok=bool(abs(est) + 2 * se <= sd_w / 4), sd_between_quarters=sd_q, sd_ratio_quarter_over_whole=sd_q / sd_w if sd_w > 0 else np.nan))
o = pd.DataFrame(out); o.to_csv(ROOT / "results/F3.2/F3.2_size_check.csv", index=False)
print(o.round(4).to_string(index=False)); print("regge:", int(o.ok.sum()), "su", len(o))
