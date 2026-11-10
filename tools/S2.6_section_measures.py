"""tools/S2.6_section_measures.py - S2.6, domanda 3: tessuto e vetro della sezione intera (a pezzi, in un pezzo solo, vero). Definizioni:
docs/steps/S2.6.md («Definizioni delle misure», P-6, P-7, misure descrittive). Blocchi di 8 um dell'area di cattura (838 x 838); tessuto =
maschera di Space Ranger; vetro per distanza fra centri dal blocco di tessuto piu' vicino: 16-64, 64-256, 256-512, >= 512 um; esclusi i blocchi
con il centro a meno di 50 um dal bordo dell'area di cattura. Rapporti dell'UMI medio per blocco con intervallo dal bootstrap appaiato delle
tessere di 500 um (B = 2000, seme 20261215). Lato vero: matrice raw di Space Ranger a 8 um (totali per blocco).
Uso (root del repo): .venv/bin/python tools/S2.6_section_measures.py <A> [out.csv]"""
import sys, os, json, numpy as np, pandas as pd, h5py, pyarrow.parquet as pq
from scipy import ndimage
A = sys.argv[1]; OUT = sys.argv[2] if len(sys.argv) > 2 else f"results/S2.6/S2.6_section_{A}.csv"
SEC = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
       "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
M = "/mnt/micron/geo_spatialtrans/S2.6/q3"; D = f"/mnt/micron/geo_spatialtrans/data_real/{SEC[A]}/binned_outputs_x/binned_outputs/square_008um"
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
def sim_blocks(tag):
    d = f"{M}/{tag}"
    if not os.path.exists(f"{d}/run.json"): return None, None
    j = json.load(open(f"{d}/run.json")); g = j["grid"]
    u2 = np.fromfile(f"{d}/u2.bin", dtype="<i4", count=g["n_row"] * g["n_col"]).reshape(g["n_row"], g["n_col"])
    n8r, n8c = g["n_row"] // 4, g["n_col"] // 4
    b8 = u2[: n8r * 4, : n8c * 4].reshape(n8r, 4, n8c, 4).sum(axis=(1, 3))
    o = int(round(-g["x0"] / 8)); oy = int(round(-g["y0"] / 8))
    assert abs(g["x0"] + 8 * o) < 1e-6 and abs(g["y0"] + 8 * oy) < 1e-6, "griglia non allineata ai blocchi di Space Ranger"
    cov = None
    if j["mode"] == "pieces":
        cv = np.zeros(g["n_row"] * g["n_col"], bool)
        for p in j["pieces"]:
            if p.get("empty"): continue
        # copertura: quadratini nell'area utile di pezzi con cellule (gli altri restano a zero per costruzione)
        U = [p["U"] for p in j["pieces"] if not p.get("empty")]
        yy = g["y0"] + (np.arange(g["n_row"]) + 0.5) * 2; xx = g["x0"] + (np.arange(g["n_col"]) + 0.5) * 2
        cvm = np.zeros((g["n_row"], g["n_col"]), bool)
        for u in U: cvm[np.ix_((yy >= u[1]) & (yy < u[3]), (xx >= u[0]) & (xx < u[2]))] = True
        cov = cvm[: n8r * 4, : n8c * 4].reshape(n8r, 4, n8c, 4).all(axis=(1, 3))[oy: oy + NR, o: o + NC]
    return b8[oy: oy + NR, o: o + NC].astype(float), cov
one, _ = sim_blocks(f"one_{A}"); pcs, cov = sim_blocks(f"pieces_{A}")
f = h5py.File(f"{D}/raw_feature_bc_matrix.h5"); m = f["matrix"]
ip = m["indptr"][:]; dat = m["data"][:]; cs = np.diff(np.concatenate([[0], np.cumsum(dat, dtype=np.int64)])[ip])   # somma per colonna (anche colonne vuote)
bc = [b.decode() for b in m["barcodes"][:]]
rr = np.array([int(b.split("_")[2]) for b in bc]); cc = np.array([int(b.split("_")[3].split("-")[0]) for b in bc])
real = np.zeros((NR, NC)); real[rr, cc] = cs
rng = np.random.default_rng(20261215); B = 2000
rows = []
for bnd in ["tissue", "glass_16_64", "glass_64_256", "glass_256_512", "glass_512plus"]:
    w = band == bnd
    if w.sum() == 0: rows.append(dict(archetype=A, band=bnd, n_blocks=0)); continue
    t = tile[w]; ut, ti = np.unique(t, return_inverse=True); nt = len(ut)
    def sums(x): return np.bincount(ti, weights=x[w], minlength=nt)
    cnt = np.bincount(ti, minlength=nt).astype(float)
    S = {k: (sums(v) if v is not None else None) for k, v in [("one", one), ("pieces", pcs), ("real", real)]}
    r = dict(archetype=A, band=bnd, n_blocks=int(w.sum()), n_tiles=nt, n_blocks_uncovered_pieces=int((~cov[w]).sum()) if cov is not None else None)
    for k in S: r[f"mean_{k}"] = S[k].sum() / cnt.sum() if S[k] is not None else np.nan
    idx = rng.integers(0, nt, (B, nt))
    W_ = np.zeros((B, nt)); np.add.at(W_, (np.repeat(np.arange(B), nt), idx.ravel()), 1)
    for a_, b_ in [("pieces", "one"), ("pieces", "real"), ("one", "real")]:
        if S[a_] is None or S[b_] is None or S[b_].sum() == 0: continue
        den = W_ @ S[b_]; ok = den > 0                       # ricampioni con denominatore nullo: rapporto non definito, esclusi (revisione G-5)
        bs = (W_ @ S[a_])[ok] / den[ok]
        r[f"ratio_{a_}_{b_}"] = S[a_].sum() / S[b_].sum(); r[f"lo_{a_}_{b_}"], r[f"hi_{a_}_{b_}"] = np.percentile(bs, [2.5, 97.5])
        r[f"n_boot_undefined_{a_}_{b_}"] = int((~ok).sum())
    rows.append(r)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
pd.DataFrame(rows).to_csv(OUT, index=False); print(pd.DataFrame(rows).to_string())
