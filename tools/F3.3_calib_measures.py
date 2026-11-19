"""tools/F3.3_calib_measures.py - F3.3: misure M-1...M-7 sulle uscite di F3.1 (30 aree vere di taratura, 60 simulate), con le
etichette di F3.1 (vero: R6; simulato: verita') e con le etichette trasferite (tools/F3_3_common.py); nel linfonodo anche i quattro
quarti di 500 um delle aree di follicolo (DEC-093, R-9). Solo dati di taratura e simulazioni gia' esistenti; nessun dato di verifica.
Uso: .venv-b2c/bin/python tools/F3.3_calib_measures.py <area_dir> <real|sim> <out_json>"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, "tools"); import F3_3_common as L
d, side, outp = sys.argv[1:4]; run = os.path.basename(d); tag = run.split("_k")[0]; A = tag[:2]; roi = tag.split("_")[1]
CACHE = "/mnt/micron/geo_spatialtrans/F3.3/cache"
bins, genes, X = L.load_area(d); exp = L.expand(d, bins, genes, X, f"{CACHE}/{side}_{os.path.basename(os.path.dirname(d))}_{run}_exp.npy")
cells = pd.read_parquet(f"{d}/cells.parquet")
orig = pd.Series(cells.cell_type.to_numpy(), index=cells.label.to_numpy())
tl, nref = L.transfer_labels(X, exp, cells.label.to_numpy(), A, genes, exclude_roi=roi if side == "real" else None)
tls = pd.Series(tl.cell_type.to_numpy(), index=tl.label.to_numpy())
os.makedirs(f"{CACHE}/labels", exist_ok=True); tl.to_parquet(f"{CACHE}/labels/{side}_{os.path.basename(os.path.dirname(d))}_{run}_tl.parquet", index=False)
parts = [("whole", None)]
if A == "A4" and roi.startswith("f"): parts += [(f"q{i}{j}", (i * 500.0, j * 500.0, (i + 1) * 500.0, (j + 1) * 500.0)) for j in (0, 1) for i in (0, 1)]
rows = []
for lab_name, cts in [("f31", orig), ("tl", tls)]:
    for part, box in parts:
        r = L.measures(d, A, cts, exp, bins, genes, X, box=box, with_aniso=(part == "whole"))
        rows.append({"side": side, "run": run, "tag": tag, "A": A, "labels": lab_name, "part": part, "n_ref_used": nref, **r})
# controllo G-1 d (con le etichette di F3.1 sull'area intera: valori di F3.1)
bj = json.load(open(f"{d}/b2c_label.json")); rj = json.load(open(f"{d}/rctd.json")); w = rows[0]
chk = {"M1": bj["share_assigned_exp"], "M2": bj["umi_per_cell_exp"], "M3": bj["foreign_exp"], "M4": rj["frac_singlet"], "M5": rj["agree_nucleus"], "M6": rj["neighbour_same"]}
dev = max(abs(w[k] - v) / max(abs(v), 1e-12) for k, v in chk.items())
for r in rows: r["g1d_max_rel_dev"] = dev
json.dump(rows, open(outp, "w"))
print(run, side, "dev", f"{dev:.2e}", {k: round(rows[len(parts)][k], 4) for k in ["M1", "M2", "M3", "M4", "M5", "M6"]}, flush=True)
