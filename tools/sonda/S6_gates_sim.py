"""tools/sonda/S6_gates_sim.py - S6 (docs/sonda/REGOLA_S6.md): cancelli della simulazione.
G-MAP : la mappa H2 -> H1 di tools/sonda/S6_sim_types.py porta il centro di ogni area di V nel centro riportato di results/sonda/S3b/areas_H1_hom.csv
        (scarto <= 1e-6 um) e h1_to_h2(h2_to_h1(p)) = p (<= 1e-6 um) su una griglia di 10 000 punti per area; la trasformazione fra hires e' affine.
G-SELF: Yself (sorgente = area misurata) rigenera bit per bit Y di S4 (conteggi per quadratino e gene, nuclei, tipi): 10/10.
G-ID  : Yalt, YH1, Yshuf hanno nuclei identici a Y, UMI totali identici a Y, numero di cellule identico, tipo = tipo della sorgente: 30/30.
Uscita: results/sonda/S6/gates_sim.json. Uso: .venv/bin/python -B tools/sonda/S6_gates_sim.py"""
import os, sys, json, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S6_sim_types as M
REPO = M.REPO; Y4 = "/mnt/micron/geo_spatialtrans/S4_sonda/runs/Y"; OUT = f"{REPO}/results/sonda/S6"; os.makedirs(OUT, exist_ok=True)
def X_of(d, shape):
    T = pd.read_parquet(f"{d}/X.parquet"); return sp.csr_matrix((T.n.to_numpy(), (T.bin.to_numpy(), T.gene.to_numpy())), shape=shape)
g = {}
A = pd.read_csv(f"{REPO}/results/sonda/S3b/areas_H1_hom.csv"); dev_c, dev_r = [], []
for a in A.itertuples():
    c = M.h2_to_h1(np.array([[a.V_x0 + a.side / 2, a.V_y0 + a.side / 2]]))[0]; dev_c.append(float(np.hypot(c[0] - a.cx_H1, c[1] - a.cy_H1)))
    gx, gy = np.meshgrid(np.linspace(a.V_x0, a.V_x0 + a.side, 100), np.linspace(a.V_y0, a.V_y0 + a.side, 100)); p = np.c_[gx.ravel(), gy.ravel()]
    dev_r.append(float(np.abs(M.h1_to_h2(M.h2_to_h1(p)) - p).max()))
g["G-MAP"] = dict(affine_last_row=M.Wh[2].tolist(), max_centre_dev_um=max(dev_c), max_roundtrip_dev_um=max(dev_r),
                  pass_=bool(np.allclose(M.Wh[2], [0, 0, 1]) and max(dev_c) <= 1e-6 and max(dev_r) <= 1e-6))
self_ok, id_rows = [], []
for i, (t, k) in enumerate(M.S.UNITS):
    u = f"H2_{t}_{k}_k0"; yd = f"{Y4}/{u}"; by = pd.read_parquet(f"{yd}/bins.parquet"); genes = open(f"{yd}/genes.tsv").read().split()
    shape = (len(by), len(genes)); Xy = X_of(yd, shape); cy = pd.read_parquet(f"{yd}/cells.parquet"); umi_y = int(Xy.sum())
    for arm in M.ARMS6:
        d = f"{M.ROOT}/runs/{arm}/{u}"
        if not os.path.exists(f"{d}/S6_source.json"): id_rows.append(dict(arm=arm, unit=u, present=False)); continue
        b = pd.read_parquet(f"{d}/bins.parquet"); X = X_of(d, shape); c = pd.read_parquet(f"{d}/cells.parquet"); chk = json.load(open(f"{d}/S4_check.json"))
        r = dict(arm=arm, unit=u, present=True, nuclei_equal=bool(np.array_equal(b.label.to_numpy(), by.label.to_numpy())), umi_total_equal=int(X.sum()) == umi_y,
                 n_cells_equal=len(c) == len(cy), type_equals_source=chk.get("typ_identity_share") == 1.0, counts_equal=bool((X != Xy).nnz == 0),
                 types_equal=bool(len(c) == len(cy) and np.array_equal(c.cell_type.to_numpy(), cy.cell_type.to_numpy())),
                 share_types_equal_Y=float(np.mean(c.cell_type.to_numpy() == cy.cell_type.to_numpy())) if len(c) == len(cy) else None)
        id_rows.append(r)
R = pd.DataFrame(id_rows); R.to_csv(f"{OUT}/gates_sim_units.csv", index=False)
s = R[R.arm == "Yself"]; g["G-SELF"] = dict(n=int(len(s)), n_pass=int((s.present & s.nuclei_equal & s.counts_equal & s.types_equal).sum()))
g["G-SELF"]["pass_"] = g["G-SELF"]["n_pass"] == 10
o = R[R.arm != "Yself"]; okk = o.present & o.nuclei_equal & o.umi_total_equal & o.n_cells_equal & o.type_equals_source
g["G-ID"] = dict(n=int(len(o)), n_pass=int(okk.sum()), pass_=bool(okk.sum() == 30),
                 share_types_equal_Y_by_arm={a: float(o[o.arm == a].share_types_equal_Y.mean()) for a in ("Yalt", "YH1", "Yshuf") if (o.arm == a).any()})
json.dump(g, open(f"{OUT}/gates_sim.json", "w"), indent=1, default=float)
for kk, v in g.items(): print(kk, "PASS" if v["pass_"] else "FAIL", {a: b for a, b in v.items() if a != "pass_"})
sys.exit(0 if all(v["pass_"] for v in g.values()) else 1)
