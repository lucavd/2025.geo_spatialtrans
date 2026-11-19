"""tools/F3.3_holdout_measures.py - F3.3: etichette trasferite e misure M-1...M-7 (tools/F3_3_common.py) su un'esecuzione della catena
delle aree nuove (/mnt/micron/geo_spatialtrans/F3.3/runs/<area>[_cg]). Aree di taratura rifatte (livello IN): profili fuori area e anche le
misure con le etichette di R6. Varianti di BL-264 (solo descrittive): M-2 e M-3 sui soli quadratini in tessuto; M-3 senza i marcatori
presenti nelle liste di due tipi. Uso: .venv-b2c/bin/python tools/F3.3_holdout_measures.py <run_dir> <out_json>"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, "tools"); import F3_3_common as L
d, outp = sys.argv[1:3]; run = os.path.basename(d); ex = json.load(open(f"{d}/export.json")); A = ex["archetype"]; level = ex["level"]
area = ex["area"]; roi = area.split("_")[1] if level == "IN" else None
CACHE = "/mnt/micron/geo_spatialtrans/F3.3/cache"
bins, genes, X = L.load_area(d); exp = L.expand(d, bins, genes, X, f"{CACHE}/runs_{run}_exp.npy")
cells = pd.read_parquet(f"{d}/cells.parquet")
tl, nref = L.transfer_labels(X, exp, cells.label.to_numpy(), A, genes, exclude_roi=roi)
os.makedirs(f"{CACHE}/labels", exist_ok=True); tl.to_parquet(f"{CACHE}/labels/runs_{run}_tl.parquet", index=False)
labs = [("tl", pd.Series(tl.cell_type.to_numpy(), index=tl.label.to_numpy()))]
if level == "IN": labs.append(("f31", pd.Series(cells.cell_type.to_numpy(), index=cells.label.to_numpy())))
rows = []
for nm, cts in labs:
    r = L.measures(d, A, cts, exp, bins, genes, X, with_aniso=True)
    rv = L.measures_variants(d, A, cts, exp, bins, genes, X)
    rows.append({"run": run, "area": area, "A": A, "level": level, "section": ex["section"], "genes": ex["genes"], "labels": nm, "part": "whole",
                 "n_ref_used": nref, "n_genes": len(genes), **r, **rv})
json.dump(rows, open(outp, "w")); print(run, {k: round(v, 4) for k, v in rows[0].items() if k.startswith("M")}, flush=True)
