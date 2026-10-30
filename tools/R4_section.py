"""tools/R4_section.py — R4, stadio per sezione. Uso: .venv/bin/python tools/R4_section.py <KEY> <outdir>
G-3(a) medie per bin in tessuto a 2 e 8 um contro il riassunto 10x; G-3(b) bin grezzi a 8 um = somma dei 16 bin grezzi a 2 um;
D-4 medie della sezione; D-5 spettro dello sfondo (bin fuori dal tessuto a >= 50 um). Scrive <outdir>/section/<KEY>.json."""
import sys, os, json, time
import numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
from scipy import ndimage
sys.path.insert(0, "tools")
import R4_common as C

M = "/mnt/micron/geo_spatialtrans"; D = f"{M}/data_real"
KEY, OUT = sys.argv[1], sys.argv[2]
os.makedirs(f"{OUT}/section", exist_ok=True)
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
ds = DS[KEY]; sq = lambda s: f"{D}/{ds}/binned_outputs_x/binned_outputs/{s}"
fout = f"{OUT}/section/{KEY}.json"
if os.path.exists(fout):
    print(f"[{KEY}] skip"); sys.exit(0)
t0 = time.time()

def totals(path):
    """Totali UMI e geni rilevati per barcode (colonne CSC), tutte le feature."""
    with h5py.File(path, "r") as f:
        m = f["matrix"]; ip = m["indptr"][:]; dat = m["data"][:]
        bc = m["barcodes"][:].astype(str); ng = int(m["shape"][0]); ft = m["features"]["feature_type"][:].astype(str)
    cs = np.concatenate([[0], np.cumsum(dat, dtype=np.int64)])
    return bc, (cs[ip[1:]] - cs[ip[:-1]]).astype(np.int64), np.diff(ip).astype(np.int64), ng, sorted(set(ft))

met = pd.read_csv(f"{D}/{ds}/{ds}_metrics_summary.csv").iloc[0]
res = {"dataset": ds, "key": KEY}
# G-3(a)
ga = {}
for s, lab_u, lab_g, lab_n in (("square_002um", "Mean UMIs Under Tissue per Square 2 µm", "Mean Genes Under Tissue per Square 2 µm", "Number of Squares Under Tissue 2 µm"),
                               ("square_008um", "Mean UMIs Under Tissue per Bin 8 µm", "Mean Genes Under Tissue per Bin 8 µm", "Number of Bins Under Tissue 8 µm")):
    bc, umi, ngen, ng, ft = totals(f"{sq(s)}/filtered_feature_bc_matrix.h5")
    ours = {"n_bins": int(len(bc)), "mean_umi": float(umi.mean()), "mean_genes": float(ngen.mean()), "n_genes_matrix": ng, "feature_types": ft,
            "median_umi": float(np.median(umi)), "frac_zero": float((umi == 0).mean())}
    ref = {"n_bins": int(met[lab_n]), "mean_umi": float(met[lab_u]), "mean_genes": float(met[lab_g])}
    ga[s] = {"ours": ours, "ref": ref, "n_bins_equal": ours["n_bins"] == ref["n_bins"],
             "rel_diff_umi": abs(ours["mean_umi"] / ref["mean_umi"] - 1), "rel_diff_genes": abs(ours["mean_genes"] / ref["mean_genes"] - 1)}
    ga[s]["pass"] = bool(ga[s]["n_bins_equal"] and ga[s]["rel_diff_umi"] <= 1e-6 and ga[s]["rel_diff_genes"] <= 1e-6)
res["G3a"] = ga
print(f"[{KEY}] G3a {[(k, v['pass'], v['rel_diff_umi'], v['rel_diff_genes']) for k, v in ga.items()]} {time.time()-t0:.0f}s", flush=True)
# G-3(b): grezzi
bc2, u2, _, _, _ = totals(f"{sq('square_002um')}/raw_feature_bc_matrix.h5")
bc8, u8, _, _, _ = totals(f"{sq('square_008um')}/raw_feature_bc_matrix.h5")
tp2 = pq.read_table(f"{sq('square_002um')}/spatial/tissue_positions.parquet").to_pandas()
tp8 = pq.read_table(f"{sq('square_008um')}/spatial/tissue_positions.parquet").to_pandas()
p2 = tp2.merge(pd.DataFrame({"barcode": bc2, "umi": u2}), on="barcode", how="left"); p2["umi"] = p2.umi.fillna(0).astype(np.int64)
agg = C.aggregate_totals(p2.array_row.to_numpy(), p2.array_col.to_numpy(), p2.umi.to_numpy())
p8 = tp8.merge(pd.DataFrame({"barcode": bc8, "umi8": u8}), on="barcode", how="left"); p8["umi8"] = p8.umi8.fillna(0).astype(np.int64)
j = p8.merge(agg.rename(columns={"row": "array_row", "col": "array_col"}), on=["array_row", "array_col"], how="outer", indicator=True)
j["umi8"] = j.umi8.fillna(0).astype(np.int64); j["umi"] = j.umi.fillna(0).astype(np.int64)
d = j.umi8 - j.umi
res["G3b"] = {"n_bins8": int(len(p8)), "n_coarse_from_2um": int(len(agg)), "n_only_one_side": int((j._merge != "both").sum()),
              "n_only_one_side_nonzero": int(((j._merge != "both") & ((j.umi8 > 0) | (j.umi > 0))).sum()),
              "n_diff": int((d != 0).sum()), "max_abs_diff": int(d.abs().max()), "sum_umi8": int(j.umi8.sum()), "sum_umi2": int(j.umi.sum()),
              "pass": bool((d == 0).all())}
if not res["G3b"]["pass"]:
    dd = j[d != 0]
    res["G3b"]["diff_examples"] = dd.head(10)[["array_row", "array_col", "umi8", "umi"]].to_dict("records")
    res["G3b"]["frac_bins_diff"] = float((d != 0).mean()); res["G3b"]["rel_diff_total"] = float(abs(d.sum()) / j.umi.sum())
print(f"[{KEY}] G3b {res['G3b']['pass']} diff {res['G3b']['n_diff']} {time.time()-t0:.0f}s", flush=True)
# D-5: sfondo fuori dal tessuto a >= 50 um, tessere 64 x 64 interamente valide
R_, C_ = tp2.array_row.max() + 1, tp2.array_col.max() + 1
img = np.zeros((R_, C_)); T = np.zeros((R_, C_), bool)
img[p2.array_row, p2.array_col] = p2.umi.to_numpy(); T[p2.array_row, p2.array_col] = p2.in_tissue.to_numpy() == 1
dist = ndimage.distance_transform_edt(~T) * 2.0
valid = (~T) & (dist >= 50.0)
S, k = C.welch_spectrum(img, valid, tile=64, min_valid=1.0)
nanpk = {"ratio": float("nan"), "present": None, "fy": float("nan"), "fx": float("nan")}
spec = C.spectral_peak(S) if k >= 4 else dict(nanpk)
spec.update({f"offaxis_{kk}": v for kk, v in (C.spectral_peak(S, off_axis=True) if k >= 4 else nanpk).items()})
spec.update({"n_tiles": int(k), "mean_umi_bg": float(img[valid].mean()) if valid.any() else float("nan")})
St, kt = C.welch_spectrum(img, T, tile=64, min_valid=1.0)       # sensibilita': tessere interamente in tessuto
spec_t = C.spectral_peak(St) if kt >= 4 else dict(nanpk)
spec_t.update({f"offaxis_{kk}": v for kk, v in (C.spectral_peak(St, off_axis=True) if kt >= 4 else nanpk).items()})
spec_t["n_tiles"] = int(kt)
np.savez_compressed(f"{OUT}/section/{KEY}_spectrum.npz", bg=S, tissue=St)
res["D5_background"] = spec; res["D5_tissue_tiles"] = spec_t
print(f"[{KEY}] D5 bg {spec} tissue {spec_t} {time.time()-t0:.0f}s", flush=True)
res["elapsed_s"] = round(time.time() - t0, 1)
json.dump(res, open(fout, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(f"[{KEY}] done", flush=True)
