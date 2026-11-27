"""Stage setup: chain, variant check, hold-out split, regions (own implementation),
region data caches (Visium UMIs and Xenium transcripts in 2 um array frame)."""
import json, sys, time, numpy as np, pandas as pd, pyarrow.parquet as pq, h5py
from g1_common import *
t0 = time.time(); res = {}
tp2 = pq.read_table(f"{BIN}/square_002um/spatial/tissue_positions.parquet").to_pandas()
tp8 = pq.read_table(f"{BIN}/square_008um/spatial/tissue_positions.parquet").to_pandas()
# ---- pixel -> array map (fit on a random 400k subsample, evaluated on all bins)
rng = np.random.default_rng(12345)
px = tp2[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
ar = tp2[["array_col", "array_row"]].to_numpy(float)
idx = rng.choice(len(tp2), 400000, replace=False)
B_aff = np.linalg.lstsq(to_h(px[idx]), ar[idx], rcond=None)[0]
e_aff = np.linalg.norm(to_h(px) @ B_aff - ar, axis=1)
H = fit_projective_lsq(px[idx], ar[idx])
e_pro = np.linalg.norm(proj_apply(H, px) - ar, axis=1)
res["px2arr_affine_max_resid_bins"] = float(e_aff.max())
res["px2arr_projective_max_resid_bins"] = float(e_pro.max())
res["px2arr_projective_p99_resid_bins"] = float(np.percentile(e_pro, 99))
M = np.loadtxt(f"{XEN}_he_imagealignment.csv", delimiter=",").reshape(3, 3)
res["M"] = M.tolist(); res["Minv"] = np.linalg.inv(M).tolist(); res["H_px2arr"] = H.tolist()
sf = json.load(open(f"{BIN}/square_002um/spatial/scalefactors_json.json"))
res["microns_per_pixel"] = sf["microns_per_pixel"]
# pixel size implied by the fitted map (um per H&E pixel along each axis), 2 um per bin
J = H[:2, :2]
res["um_per_px_from_map"] = (2.0 / np.linalg.norm(J, axis=1)).tolist()
res["um_per_px_from_M"] = float(XPIX_UM * np.sqrt(abs(np.linalg.det(M[:2, :2]))))
json.dump(res, open(f"{OUT}/chain.json", "w"), indent=1)
# ---- variant check (fraction of Xenium cell centroids on in-tissue 2 um bins / within grid)
cells = pq.read_table(f"{XEN}_cells.parquet", columns=["x_centroid", "y_centroid"]).to_pandas().to_numpy(float)
nr2 = int(tp2.array_row.max()) + 1; nc2 = int(tp2.array_col.max()) + 1
tis2 = np.zeros((nr2, nc2), bool); t = tp2[tp2.in_tissue == 1]; tis2[t.array_row, t.array_col] = True
var = {}
for dn, MM in [("direct", M), ("inverse", np.linalg.inv(M))]:
    for un, sc in [("px0.2125", 1 / XPIX_UM), ("um", 1.0)]:
        a = proj_apply(H, proj_apply(MM, cells * sc)); k = np.rint(a).astype(int)
        ins = (k[:, 0] >= 0) & (k[:, 0] < nc2) & (k[:, 1] >= 0) & (k[:, 1] < nr2)
        on = np.zeros(len(k), bool); on[ins] = tis2[k[ins, 1], k[ins, 0]]
        var[f"{dn}|{un}"] = dict(frac_in_grid=float(ins.mean()), frac_on_tissue_all=float(on.mean()),
                                 frac_on_tissue_given_grid=float(on[ins].mean()) if ins.any() else None)
res["variants"] = var
# ---- hold-out split and regions (8 um grid, x = (array_col+0.5)*8)
x8 = (tp8.array_col.to_numpy() + 0.5) * 8.0
m = float(np.median(x8[tp8.in_tissue.to_numpy() == 1])); res["split_m_um"] = m
nr8 = int(tp8.array_row.max()) + 1; nc8 = int(tp8.array_col.max()) + 1
a = xenium_to_array(cells, res)
k8 = np.rint(a).astype(int) // 4
ok = (k8[:, 0] >= 0) & (k8[:, 0] < nc8) & (k8[:, 1] >= 0) & (k8[:, 1] < nr8)
occ = np.zeros((nr8, nc8), bool); occ[k8[ok, 1], k8[ok, 0]] = True
occ3 = occ.copy()                       # 3x3 neighbourhood by explicit shifts
for dy in (-1, 0, 1):
    for dx in (-1, 0, 1):
        sh = np.zeros_like(occ)
        sh[max(dy, 0):nr8 + min(dy, 0), max(dx, 0):nc8 + min(dx, 0)] = occ[max(-dy, 0):nr8 + min(-dy, 0), max(-dx, 0):nc8 + min(-dx, 0)]
        occ3 |= sh
tis8 = np.zeros((nr8, nc8), bool); t = tp8[tp8.in_tissue == 1]; tis8[t.array_row, t.array_col] = True
cov = tis8 & occ3
def region_search(side, thr):
    nb = int(round(side / 8)); C = []
    for x0 in np.arange(0, nc8 * 8 - side + 1e-9, 100.0):
        if x0 + side > m - 50: continue
        for y0 in np.arange(0, nr8 * 8 - side + 1e-9, 100.0):
            c0, r0 = int(round(x0 / 8)), int(round(y0 / 8))
            f = cov[r0:r0 + nb, c0:c0 + nb].mean()
            if f >= thr: C.append((x0, y0, f))
    return np.array(C)
for side, thr in [(1000.0, 0.9), (1000.0, 0.8), (750.0, 0.8)]:
    C = region_search(side, thr)
    if len(C) >= 3: break
rr, cc = np.nonzero(cov); dev = (cc + 0.5) * 8 < m - 50
bary = np.array([((cc[dev] + 0.5) * 8).mean(), ((rr[dev] + 0.5) * 8).mean()])
ctr = C[:, :2] + side / 2
def pick(score, chosen):   # smallest score, ties: x then y
    order = np.lexsort((C[:, 1], C[:, 0], np.round(score, 9)))
    for i in order:
        if all(abs(C[i, 0] - C[j, 0]) >= side or abs(C[i, 1] - C[j, 1]) >= side for j in chosen):
            return int(i)
r1 = pick(np.linalg.norm(ctr - bary, axis=1), [])
r2 = pick(-np.linalg.norm(ctr - ctr[r1], axis=1), [r1])
r3 = pick(-np.minimum(np.linalg.norm(ctr - ctr[r1], axis=1), np.linalg.norm(ctr - ctr[r2], axis=1)), [r1, r2])
res["region_side_um"] = side; res["coverage_threshold"] = thr; res["n_candidates"] = int(len(C))
res["regions"] = [dict(name=f"R{i+1}", x0=float(C[j, 0]), y0=float(C[j, 1]), cover=float(C[j, 2])) for i, j in enumerate([r1, r2, r3])]
json.dump(res, open(f"{OUT}/chain.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("M", "Minv", "H_px2arr")}, indent=1), flush=True)
# ---- gene set
panel = json.load(open(f"{XEN}_gene_panel.json"))["payload"]["targets"]
sym2ens = {t["type"]["data"]["name"]: t["type"]["data"]["id"] for t in panel
           if t.get("type", {}).get("descriptor") == "gene" and str(t["type"]["data"].get("id", "")).startswith("ENSG")}
with h5py.File(f"{BIN}/square_002um/filtered_feature_bc_matrix.h5", "r") as f:
    g = f["matrix"]; fid = np.array([x.decode() for x in g["features/id"][:]])
    bcs = np.array([x.decode() for x in g["barcodes"][:]])
    indptr = g["indptr"][:].astype(np.int64); ind = g["indices"][:]; dat = g["data"][:]
genes = sorted(set(sym2ens.values()) & set(fid)); gidx = {e: i for i, e in enumerate(genes)}
print("panel genes with ENSG", len(sym2ens), "shared with Visium", len(genes), flush=True)
json.dump(dict(genes=genes, n_panel_ensg=len(sym2ens)), open(f"{OUT}/genes.json", "w"))
frow = np.array([gidx.get(e, -1) for e in fid])
# Visium entries -> (gene, col, row, n)
bc_pos = tp2.set_index("barcode").loc[bcs, ["array_col", "array_row"]].to_numpy(np.int64)
col_of_entry = np.repeat(np.arange(len(bcs)), np.diff(indptr))
gk = frow[ind]; keep = gk >= 0
Vg = gk[keep]; Vc = bc_pos[col_of_entry[keep], 0]; Vr = bc_pos[col_of_entry[keep], 1]; Vn = dat[keep]
print("visium shared-gene UMI total", int(Vn.sum()), flush=True)
# Xenium transcripts
tx = pq.read_table(f"{XEN}_transcripts.parquet", columns=["feature_name", "x_location", "y_location", "qv"]).to_pandas()
fn = tx.feature_name.astype(str)
ens = fn.map(sym2ens); keepx = (tx.qv.to_numpy() >= 20) & ens.isin(gidx.keys()).to_numpy()
res_tx = dict(n_all=int(len(tx)), n_qv20=int((tx.qv.to_numpy() >= 20).sum()), n_kept=int(keepx.sum()))
Xg = ens[keepx].map(gidx).to_numpy(np.int32)
Xa = xenium_to_array(tx.loc[keepx, ["x_location", "y_location"]].to_numpy(float), res)
print("xenium", res_tx, flush=True)
for R in res["regions"]:
    c0 = int(round(R["x0"] / 2)) - MARGIN_BINS; r0 = int(round(R["y0"] / 2)) - MARGIN_BINS
    n = int(round(side / 2)) + 2 * MARGIN_BINS
    kv = (Vc >= c0) & (Vc < c0 + n) & (Vr >= r0) & (Vr < r0 + n)
    kx = (Xa[:, 0] >= c0 - 0.5) & (Xa[:, 0] < c0 + n - 0.5) & (Xa[:, 1] >= r0 - 0.5) & (Xa[:, 1] < r0 + n - 0.5)
    np.savez_compressed(f"{OUT}/region_{R['name']}.npz", c0=c0, r0=r0, n=n,
                        vg=Vg[kv], vc=Vc[kv] - c0, vr=Vr[kv] - r0, vn=Vn[kv],
                        xg=Xg[kx], xc=Xa[kx, 0] - c0, xr=Xa[kx, 1] - r0)
    print(R["name"], "visium UMI", int(Vn[kv].sum()), "xenium tx", int(kx.sum()), flush=True)
print("setup done in", round(time.time() - t0), "s")
