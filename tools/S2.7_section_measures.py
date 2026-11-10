"""tools/S2.7_section_measures.py - S2.7, P-2: tessuto e vetro della sezione intera a pezzi (S2.7) contro la sezione in un pezzo solo di S2.6
(motore di S2.5) e contro il vero (descrittivo). Definizioni di S2.6 (tools/S2.6_section_measures.py): blocchi di 8 um dell'area di cattura,
tessuto = maschera di Space Ranger, vetro per distanza dal tessuto (16-64, 64-256, 256-512, >= 512 um), esclusi i 50 um al bordo; rapporti
dell'UMI medio per blocco con intervallo percentile dal bootstrap appaiato delle tessere di 500 um (B = 2000, seme 20261221).
Uso (root del repo): .venv/bin/python tools/S2.7_section_measures.py <A> <s27_tag> [out.csv]"""
import sys, os, json, numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
from scipy import ndimage
A = sys.argv[1]; TAG = sys.argv[2]; OUT = sys.argv[3] if len(sys.argv) > 3 else f"results/S2.7/S2.7_section_{A}.csv"
SEC = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
       "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
M6 = "/mnt/micron/geo_spatialtrans/S2.6/q3"; M7 = "/mnt/micron/geo_spatialtrans/S2.7/sections"
D = f"/mnt/micron/geo_spatialtrans/data_real/{SEC[A]}/binned_outputs_x/binned_outputs/square_008um"
tp = pq.read_table(f"{D}/spatial/tissue_positions.parquet", columns=["in_tissue", "array_row", "array_col"]).to_pandas()
NR, NC = tp.array_row.max() + 1, tp.array_col.max() + 1
tis = np.zeros((NR, NC), bool); tis[tp.array_row[tp.in_tissue == 1], tp.array_col[tp.in_tissue == 1]] = True
dist = ndimage.distance_transform_edt(~tis) * 8.0
yc = (np.arange(NR) + 0.5) * 8; xc = (np.arange(NC) + 0.5) * 8
edge = (np.minimum(yc, NR * 8 - yc)[:, None] < 50) | (np.minimum(xc, NC * 8 - xc)[None, :] < 50)
band = np.full((NR, NC), "", object); band[tis] = "tissue"
for lo, hi, nm in [(16, 64, "glass_16_64"), (64, 256, "glass_64_256"), (256, 512, "glass_256_512"), (512, 1e9, "glass_512plus")]:
    band[(~tis) & (dist >= lo) & (dist < hi)] = nm
band[edge] = ""
tile = (np.floor(yc / 500)[:, None] * 100 + np.floor(xc / 500)[None, :]).astype(int)
def blocks(u2file, g):
    u2 = np.fromfile(u2file, dtype="<i4", count=g["n_row"] * g["n_col"]).reshape(g["n_row"], g["n_col"])
    n8r, n8c = g["n_row"] // 4, g["n_col"] // 4
    b8 = u2[: n8r * 4, : n8c * 4].reshape(n8r, 4, n8c, 4).sum(axis=(1, 3))
    o = int(round(-g["x0"] / 8)); oy = int(round(-g["y0"] / 8))
    assert abs(g["x0"] + 8 * o) < 1e-6 and abs(g["y0"] + 8 * oy) < 1e-6, "griglia non allineata ai blocchi di Space Ranger"
    return b8[oy: oy + NR, o: o + NC].astype(float)
S = {}
j7 = json.load(open(f"{M7}/{TAG}.json")); S["pieces"] = blocks(f"{M7}/{TAG}.bin", j7["grid"])
if os.path.exists(f"{M6}/one_{A}/run.json"):
    j6 = json.load(open(f"{M6}/one_{A}/run.json")); S["one"] = blocks(f"{M6}/one_{A}/u2.bin", j6["grid"])
f = h5py.File(f"{D}/raw_feature_bc_matrix.h5"); m = f["matrix"]
ip = m["indptr"][:]; dat = m["data"][:]; cs = np.diff(np.concatenate([[0], np.cumsum(dat, dtype=np.int64)])[ip])
bc = [b.decode() for b in m["barcodes"][:]]
rr = np.array([int(b.split("_")[2]) for b in bc]); cc = np.array([int(b.split("_")[3].split("-")[0]) for b in bc])
real = np.zeros((NR, NC)); real[rr, cc] = cs; S["real"] = real
rng = np.random.default_rng(20261221); B = 2000
rows = []
for bnd in ["tissue", "glass_16_64", "glass_64_256", "glass_256_512", "glass_512plus"]:
    w = band == bnd
    if w.sum() == 0: rows.append(dict(archetype=A, band=bnd, n_blocks=0)); continue
    t = tile[w]; ut, ti = np.unique(t, return_inverse=True); nt = len(ut)
    cnt = np.bincount(ti, minlength=nt).astype(float)
    T = {k: np.bincount(ti, weights=v[w], minlength=nt) for k, v in S.items()}
    r = dict(archetype=A, band=bnd, n_blocks=int(w.sum()), n_tiles=nt)
    for k in T: r[f"mean_{k}"] = T[k].sum() / cnt.sum()
    idx = rng.integers(0, nt, (B, nt))
    W_ = np.zeros((B, nt)); np.add.at(W_, (np.repeat(np.arange(B), nt), idx.ravel()), 1)
    for a_, b_ in [("pieces", "one"), ("pieces", "real")]:
        if a_ not in T or b_ not in T or T[b_].sum() == 0: continue
        den = W_ @ T[b_]; ok = den > 0
        bs = (W_ @ T[a_])[ok] / den[ok]
        r[f"ratio_{a_}_{b_}"] = T[a_].sum() / T[b_].sum(); r[f"lo_{a_}_{b_}"], r[f"hi_{a_}_{b_}"] = np.percentile(bs, [2.5, 97.5])
        r[f"n_boot_undefined_{a_}_{b_}"] = int((~ok).sum())
    rows.append(r)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
pd.DataFrame(rows).to_csv(OUT, index=False); print(pd.DataFrame(rows).to_string())
