"""tools/riesame_panel_measures.py - riesame, Parte B: etichette trasferite e misure M-1...M-7 su un'esecuzione della catena del
pannello (uscite di tools/F3.1_b2c.py e tools/F3.1_rctd.R in <run_dir>), con le funzioni di tools/F3_3_common.py, come
tools/F3.3_holdout_measures.py e tools/F3.3_calib_measures.py: profili del riferimento di F3.1 sui soli geni dell'area (quindi gia'
ristretti all'insieme di geni), fuori area per le aree vere di taratura (roi), tutti altrove. Con etichette proprie dell'area
(vero di taratura: R6; simulato: verita') calcola anche le misure con quelle etichette e lo scarto relativo massimo dai JSON di
bin2cell e RCTD della stessa esecuzione (controllo interno, come G-1 d di F3.3).
Uso: .venv-b2c/bin/python tools/riesame_panel_measures.py <run_dir> <A> <roi|-> <out_json> arm=.. genes=.. unit=.. run=.. src=.. set=.. level=.. u2=..
Scrive un JSON con una riga per etichettatura (labels = tl, e f31 se disponibile)."""
import sys, os, json, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools")); import F3_3_common as L
d, A, roi, outp = sys.argv[1:5]; meta = dict(a.split("=", 1) for a in sys.argv[5:]); roi = None if roi == "-" else roi
t0 = time.time()
bins, genes, X = L.load_area(d); exp = L.expand(d, bins, genes, X, f"{d}/exp_labels_expanded.npy")
cells = pd.read_parquet(f"{d}/cells.parquet")
tl, nref = L.transfer_labels(X, exp, cells.label.to_numpy(), A, genes, exclude_roi=roi)
tl.to_parquet(f"{d}/labels_transferred.parquet", index=False)
labs = [("tl", pd.Series(tl.cell_type.to_numpy(), index=tl.label.to_numpy()))]
if "cell_type" in cells and cells.cell_type.notna().any(): labs.append(("f31", pd.Series(cells.cell_type.to_numpy(), index=cells.label.to_numpy())))
bj = json.load(open(f"{d}/b2c_label.json")); rj = json.load(open(f"{d}/rctd.json"))
rows = []
for nm, cts in labs:
    r = L.measures(d, A, cts, exp, bins, genes, X, with_aniso=True); rv = L.measures_variants(d, A, cts, exp, bins, genes, X)
    rows.append({**meta, "A": A, "labels": nm, "roi_excluded": roi, "n_genes": int(len(genes)), "n_ref_used": nref, **r, **rv,
                 "umi_reg": float(bj["umi_reg"]), "rctd_n_blocks": rj["n_blocks"], "rctd_n_run": rj["n_run"], "rctd_n_genes": rj["n_genes_common"],
                 "UMI_min_sigma": rj["UMI_min_sigma"], "b2c_sec": bj.get("sec"), "rctd_sec": rj.get("sec")})
if len(rows) > 1:
    w = rows[1]; chk = {"M1": bj["share_assigned_exp"], "M2": bj["umi_per_cell_exp"], "M3": bj["foreign_exp"], "M4": rj["frac_singlet"], "M5": rj["agree_nucleus"], "M6": rj["neighbour_same"]}
    dev = max(abs(w[k] - v) / max(abs(v), 1e-12) for k, v in chk.items())
else: dev = None
for r in rows: r["internal_max_rel_dev"] = dev; r["measure_sec"] = round(time.time() - t0, 1)
json.dump(rows, open(outp, "w"))
print(meta.get("run"), meta.get("genes"), "dev", dev, {k: round(rows[0][k], 4) for k in ["M1", "M2", "M3", "M4", "M5", "M6", "M7"]}, flush=True)
