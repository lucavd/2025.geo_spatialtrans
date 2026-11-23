"""tools/sonda/S3_measures.py - S3: etichette trasferite e misure M-1...M-6 di un'esecuzione della catena (uscite di tools/F3.1_b2c.py e
tools/F3.1_rctd.R in <run_dir>), con le funzioni invariate di tools/F3_3_common.py e il riferimento del polmone (sigla L, profili da tutte
le cellule del riferimento). Per B calcola anche le misure con i tipi veri (labels = true, descrittive).
Uso: .venv-b2c/bin/python -B tools/sonda/S3_measures.py <run_dir> <tumour|stroma> <out_json> arm=.. unit=.. run=.. set=.. k=..
La radice del riferimento si cambia con S3_REF."""
import sys, os, json, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools")); import F3_3_common as L
L.REF = os.environ.get("S3_REF", "/mnt/micron/geo_spatialtrans/S3_sonda/ref")
d, tissue, outp = sys.argv[1:4]; meta = dict(a.split("=", 1) for a in sys.argv[4:]); t0 = time.time()
bins, genes, X = L.load_area(d); exp = L.expand(d, bins, genes, X, f"{d}/exp_labels_expanded.npy")
cells = pd.read_parquet(f"{d}/cells.parquet")
tl, nref = L.transfer_labels(X, exp, cells.label.to_numpy(), "L", genes)
tl.to_parquet(f"{d}/labels_transferred.parquet", index=False)
labs = [("tl", pd.Series(tl.cell_type.to_numpy(), index=tl.label.to_numpy()))]
if cells.cell_type.notna().any(): labs.append(("true", pd.Series(cells.cell_type.to_numpy(), index=cells.label.to_numpy())))
bj = json.load(open(f"{d}/b2c_label.json")); rj = json.load(open(f"{d}/rctd.json")); rows = []
for nm, cts in labs:
    r = L.measures(d, "L", cts, exp, bins, genes, X)
    rows.append({**meta, "A": tissue, "labels": nm, "n_genes": int(len(genes)), "n_ref_used": nref, **r, "umi_reg": float(bj["umi_reg"]),
                 "rctd_n_blocks": rj["n_blocks"], "rctd_n_run": rj["n_run"], "rctd_n_genes": rj["n_genes_common"], "UMI_min_sigma": rj["UMI_min_sigma"],
                 "b2c_sec": bj.get("sec"), "rctd_sec": rj.get("sec"), "measures_sec": round(time.time() - t0, 1)})
json.dump(rows, open(outp, "w"), indent=1, default=float)
print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[0].items()}), flush=True)
