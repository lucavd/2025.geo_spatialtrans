"""tools/sonda/S4_test.py - S4: test del simulatore dei bracci (tools/sonda/S4_sim_arms.py). PASS/FAIL per asserzione.
T1 griglia minuscola: dens V riproduce gli UMI per quadratino; T2 geo V riproduce i nuclei veri; T3 typ V = tipo della cellula stessa o della
tipizzata piu' vicina; T4 prof V usa solo i geni della cellula sorgente; T5 determinismo con lo stesso seme; T6 identita' con il banale:
tutti gli interruttori su B e il seme di B rigenerano bit per bit la simulazione B esistente (H2_tumour_1, seme 0).
Uso: .venv/bin/python -B tools/sonda/S4_test.py"""
import os, sys, json, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S4_sim_arms as S
res = []
def check(name, ok): res.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name, flush=True)
# griglia 30 x 30 quadratini da 2 um, angolo fuori tessuto, 4 nuclei veri 3 x 3
g = 30; yy, xx = np.divmod(np.arange(g * g), g); x = xx * 2.0 + 1; y = yy * 2.0 + 1
it = ~((xx < 4) & (yy < 4)); lab = np.zeros(g * g, np.int64)
for j, (a, b) in enumerate([(6, 6), (6, 20), (20, 6), (20, 20)]): lab[((xx >= a) & (xx < a + 3) & (yy >= b) & (yy < b + 3))] = j + 1
r0 = np.random.default_rng(1); umi = r0.poisson(3, g * g); ng = 7
Xexp = sp.csr_matrix(np.array([[30, 0, 5, 0, 0, 0, 0], [0, 40, 0, 0, 2, 0, 0], [0, 0, 0, 25, 0, 0, 0], [5, 0, 0, 0, 0, 10, 0]], float))
real_type = np.array([0, 1, -1, 0]); types = ["t0", "t1"]; P = np.full((2, ng), 1.0 / ng)
par = {"lambda_um2": 4 / (g * g * 4.0), "r_um": 3.0, "composition": {"t0": 0.5, "t1": 0.5}, "rho_per_bin": 2.0}
out = {}
for arm, spec in S.ARMS.items():
    out[arm] = S.simulate(x, y, it, lab, umi, Xexp, real_type, par, P, types, spec, np.random.default_rng([7, S.ARM_ID[arm]]))
for arm in ("U", "A", "AE"):
    label, rows, gg, ct, src, ids, d = out[arm]; per = np.bincount(rows, minlength=g * g)
    check(f"T1 {arm}: UMI per quadratino = V in tissue, 0 fuori", np.array_equal(per[it], umi[it]) and per[~it].sum() == 0)
for arm in ("G", "A", "AE"):
    label = out[arm][0]; check(f"T2 {arm}: nuclei = nuclei veri", np.array_equal(label > 0, lab > 0) and np.array_equal(label, lab))
for arm in ("A", "AE"):
    ct = out[arm][3]; check(f"T3 {arm}: tipi = tipo vero; non tipizzata -> tipizzata piu' vicina", list(ct) == [0, 1, int(real_type[out[arm][4][2]]), 0] and out[arm][4][2] != 2)
ct, src = out["Y"][3], out["Y"][4]; check("T3 Y: tipo = tipo della sorgente tipizzata", np.all(real_type[src] >= 0) and np.array_equal(ct, real_type[src]))
label, rows, gg, ct, src, ids, d = out["AE"]; cellbin = np.zeros(g * g, np.int64)
ok = True
from scipy.spatial import cKDTree
_, cxy, _ = S.real_cells(x, y, lab); _, own = cKDTree(cxy).query(np.c_[x, y]); cellbin = own; cellbin[lab > 0] = lab[lab > 0] - 1
for b_, gene in zip(rows, gg):
    s = src[cellbin[b_]]
    if Xexp.getrow(s).sum() >= S.MIN_UMI_PROF and Xexp[s, gene] == 0: ok = False; break
check("T4 AE: geni solo nel supporto della cellula sorgente", ok)
for arm, spec in S.ARMS.items():
    a2 = S.simulate(x, y, it, lab, umi, Xexp, real_type, par, P, types, spec, np.random.default_rng([7, S.ARM_ID[arm]]))
    check(f"T5 {arm}: determinismo", all(np.array_equal(u, v) for u, v in zip(out[arm][:6], a2[:6])))
# T6 identita' con B esistente
src_dir = f"{S.S3}/runs/V/H2_tumour_1"; bdir = f"{S.S3}/runs/B/H2_tumour_1_k0"
par = json.load(open(f"{S.PAR}/params_tumour.json")); z = np.load(f"{S.PAR}/params_tumour.npz"); P = z["profiles"]; types = list(z["types"])
bins = pd.read_parquet(f"{src_dir}/bins.parquet"); xb, ybb = bins.x_um.to_numpy(), bins.y_um.to_numpy(); itb = bins.in_tissue.to_numpy().astype(bool)
labv = bins.label.to_numpy().astype(np.int64); ids, _, _ = S.real_cells(xb, ybb, labv)
rt = np.zeros(len(ids), np.int64)                                  # irrilevante con tutti gli interruttori su B (serve solo una tipizzata)
spec = dict(geo="B", typ="B", dens="B", prof="B"); Xe = sp.csr_matrix((len(ids), len(types)))
label, rows, gg, ct, src, _, d = S.simulate(xb, ybb, itb, labv, np.zeros(len(bins), np.int64), Xe, rt, par, P, types, spec, np.random.default_rng([20261125, 0, 0]))
Xs = sp.csr_matrix((np.ones(len(rows), np.int64), (rows, gg)), shape=(len(bins), P.shape[1])); Xs.sum_duplicates()
TB = pd.read_parquet(f"{bdir}/X.parquet"); Xb = sp.csr_matrix((TB.n.to_numpy(), (TB.bin.to_numpy(), TB.gene.to_numpy())), shape=Xs.shape)
lb = pd.read_parquet(f"{bdir}/bins.parquet").label.to_numpy(); cb = pd.read_parquet(f"{bdir}/cells.parquet")
check("T6 tutti su B + seme di B = simulazione B esistente (conteggi, nuclei, tipi)", (Xs != Xb).nnz == 0 and np.array_equal(label, lb)
      and list(np.array(types, dtype=object)[ct]) == list(cb.cell_type))
print(f"{sum(res)}/{len(res)} PASS"); sys.exit(0 if all(res) else 1)
