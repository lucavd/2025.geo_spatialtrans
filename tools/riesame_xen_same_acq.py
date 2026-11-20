"""tools/riesame_xen_same_acq.py - Riesame Parte B: criterio XD-24 (docs/riesame/DECISIONI_B_XENIUM.md), lo zip GEO di P2 e' un ritaglio
della stessa acquisizione Xenium del file di 10x? (1) quota dei cell_id GEO presenti nel file di 10x dopo formato uniforme; (2) scarto
(x, y) costante per le cellule comuni (deviazione standard di ogni componente < 1 um); (3) conteggi per cellula uguali in almeno il 99 %.
Solo lettura. Uscita: results/riesame/xen_same_acquisition.json. Uso: .venv/bin/python -B tools/riesame_xen_same_acq.py"""
import json, hashlib, numpy as np, pandas as pd
XP = "/mnt/micron/geo_spatialtrans/data_xenium/P2"
def md5(p):
    h = hashlib.md5(); f = open(p, "rb")
    for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()
fi = [f"{XP}/tenx_alignment/Xenium_P2_cell_info.csv.gz", f"{XP}/outs/cells.parquet"]; before = [md5(p) for p in fi]
gh = pd.read_csv(fi[0]); ce = pd.read_parquet(fi[1])
def norm(s, drop_suffix):
    s = s.map(lambda v: v.decode() if isinstance(v, (bytes, bytearray)) else str(v)).str.strip().str.lower()
    return s.str.replace(r"-\d+$", "", regex=True) if drop_suffix else s
res = {"geo_columns": list(ce.columns), "github_columns": list(gh.columns), "geo_id_dtype": str(ce.cell_id.dtype), "github_id_dtype": str(gh.cell_id.dtype),
       "geo_id_examples": [str(v) for v in ce.cell_id.head(5)], "github_id_examples": [str(v) for v in gh.cell_id.head(5)],
       "n_geo": int(len(ce)), "n_github": int(len(gh)), "by_format": {}}
for ds in (False, True):
    g = norm(ce.cell_id, ds); h = norm(gh.cell_id, ds)
    res["by_format"]["drop_suffix" if ds else "as_is"] = {"share_geo_in_github": float(g.isin(set(h)).mean()), "n_common": int(g.isin(set(h)).sum()),
                                                          "geo_unique": bool(g.is_unique), "github_unique": bool(h.is_unique)}
best = max(res["by_format"], key=lambda k: res["by_format"][k]["share_geo_in_github"]); res["format_used"] = best
g = norm(ce.cell_id, best == "drop_suffix"); h = norm(gh.cell_id, best == "drop_suffix")
G = ce.assign(key=g).drop_duplicates("key").set_index("key"); H = gh.assign(key=h).drop_duplicates("key").set_index("key")
com = G.index.intersection(H.index); res["n_common"] = int(len(com)); res["share_geo_in_github"] = float(len(com) / len(ce))
c2 = {}
if len(com):
    d = H.loc[com, ["x_centroid", "y_centroid"]].to_numpy(float) - G.loc[com, ["x_centroid", "y_centroid"]].to_numpy(float)
    c2 = {"dx_mean": float(d[:, 0].mean()), "dy_mean": float(d[:, 1].mean()), "dx_sd": float(d[:, 0].std(ddof=1)) if len(com) > 1 else None,
          "dy_sd": float(d[:, 1].std(ddof=1)) if len(com) > 1 else None, "dist_q": {q: float(np.quantile(np.hypot(*d.T), q)) for q in (0.1, 0.5, 0.9)}}
    for col in ("transcript_counts", "total_counts"):
        if col in G.columns and col in H.columns: c2[f"share_equal_{col}"] = float((G.loc[com, col].to_numpy() == H.loc[com, col].to_numpy()).mean())
res["common_cells"] = c2
ok1 = res["share_geo_in_github"] >= 0.99
ok2 = bool(c2) and c2.get("dx_sd") is not None and c2["dx_sd"] < 1 and c2["dy_sd"] < 1
ok3 = bool(c2) and all(c2.get(f"share_equal_{c}", 0) >= 0.99 for c in ("transcript_counts", "total_counts"))
res["criterion"] = {"1_ids_99pct": ok1, "2_constant_vector": ok2, "3_counts_equal": ok3, "same_acquisition": bool(ok1 and ok2 and ok3)}
res["md5_unchanged"] = before == [md5(p) for p in fi]
json.dump(res, open("results/riesame/xen_same_acquisition.json", "w"), indent=1); print(json.dumps(res, indent=1))
