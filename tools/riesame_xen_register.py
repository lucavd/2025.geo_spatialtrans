"""tools/riesame_xen_register.py - Riesame Parte B, braccio Xenium: registrazione Xenium di P2 -> coordinate del H&E a piena risoluzione di V
(Visium_HD_Human_Colon_Cancer). Scelte XD-3, XD-4, XD-5 di docs/riesame/DECISIONI_B_XENIUM.md.
1) Stesso sistema di coordinate fra il file di 10x (GitHub, Xenium_P2_cell_info.csv.gz) e lo zip GEO (cells.parquet): abbinamento per
   cell_id e, in ogni caso, distanza dal centroide GEO piu' vicino.
2) Affine (x, y Xenium in um) -> (colonna, riga in pixel di V) con minimi quadrati robusti (Huber, soglia 1.345 * 1.4826 * mediana della
   norma del residuo, 20 iterazioni); residui in um (pixel * um_per_px di V); se il 95 percentile supera 2 um: affine locale per area.
Uscite: /mnt/micron/geo_spatialtrans/riesame/B/registration.json, results/riesame/xen_registration.csv.
Uso: .venv/bin/python -B tools/riesame_xen_register.py"""
import json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
M = "/mnt/micron/geo_spatialtrans"; XP = f"{M}/data_xenium/P2"; OUT = f"{M}/riesame/B"
AREAS_CSV = "results/F3.3/F3.3_areas.csv"
gh = pd.read_csv(f"{XP}/tenx_alignment/Xenium_P2_cell_info.csv.gz")
ce = pd.read_parquet(f"{XP}/outs/cells.parquet")
A = pd.read_csv(AREAS_CSV); A = A[(A.level == "IN") & A.A.isin(["A2", "A3"])].reset_index(drop=True)
UPP = float(A.upp.iloc[0]); assert np.allclose(A.upp, UPP)
res = {"n_github": int(len(gh)), "n_geo_cells": int(len(ce)), "upp_um_per_px": UPP, "github_columns": list(gh.columns)}
# --- 1) stesso sistema di coordinate?
ids = np.intersect1d(gh.cell_id.astype(str), ce.cell_id.astype(str)); res["n_id_common"] = int(len(ids))
if len(ids):
    g1 = gh.set_index("cell_id").loc[ids, ["x_centroid", "y_centroid"]].to_numpy(); c1 = ce.set_index("cell_id").loc[ids, ["x_centroid", "y_centroid"]].to_numpy()
    dd = np.hypot(*(g1 - c1).T); res["id_match_centroid_dist_um_q"] = {q: float(np.quantile(dd, q)) for q in (0.5, 0.9, 0.99, 1.0)}
tr = cKDTree(ce[["x_centroid", "y_centroid"]].to_numpy()); dn, _ = tr.query(gh[["x_centroid", "y_centroid"]].to_numpy())
res["nn_github_to_geo_um_q"] = {q: float(np.quantile(dn, q)) for q in (0.1, 0.5, 0.9, 0.99)}
res["share_github_within_1um_of_geo"] = float((dn <= 1).mean())
tr2 = cKDTree(gh[["x_centroid", "y_centroid"]].to_numpy()); dn2, _ = tr2.query(ce[["x_centroid", "y_centroid"]].to_numpy())
res["share_geo_within_1um_of_github"] = float((dn2 <= 1).mean())
res["github_bbox_um"] = [float(gh.x_centroid.min()), float(gh.y_centroid.min()), float(gh.x_centroid.max()), float(gh.y_centroid.max())]
res["geo_bbox_um"] = [float(ce.x_centroid.min()), float(ce.y_centroid.min()), float(ce.x_centroid.max()), float(ce.y_centroid.max())]
# --- 2) affine robusta
def fit(X, U, it=20):
    D = np.c_[X, np.ones(len(X))]; w = np.ones(len(X))
    for _ in range(it):
        W = np.sqrt(w)[:, None]; B, *_ = np.linalg.lstsq(D * W, U * W, rcond=None)
        r = np.hypot(*(U - D @ B).T); s = 1.4826 * np.median(r); c = 1.345 * max(s, 1e-12)
        w = np.minimum(1.0, c / np.maximum(r, 1e-12))
    return B, r
X = gh[["x_centroid", "y_centroid"]].to_numpy(float); U = gh[["x_centroid_visium_scale", "y_centroid_visium_scale"]].to_numpy(float)
B, r = fit(X, U); r_um = r * UPP
def qstats(r_um):
    return {"median_um": float(np.median(r_um)), "p95_um": float(np.quantile(r_um, 0.95)), "max_um": float(r_um.max()),
            "share_within_2um": float((r_um <= 2).mean()), "share_within_5um": float((r_um <= 5).mean()), "share_within_10um": float((r_um <= 10).mean())}
res["global"] = {"B": B.tolist(), "n": int(len(X)), **qstats(r_um)}
# scala e rotazione: colonne di B = [[a, d], [b, e], [tx, ty]]: u = a x + b y + tx, v = d x + e y + ty
L = B[:2].T; sv = np.linalg.svd(L, compute_uv=False); res["global"]["px_per_um_singular"] = sv.tolist(); res["global"]["um_per_um_singular"] = (sv * UPP).tolist()
res["global"]["det_sign"] = float(np.sign(np.linalg.det(L)))
mode = "global" if res["global"]["p95_um"] <= 2.0 else "local"; res["mode"] = mode
rows = []; res["areas"] = {}
for _, a in A.iterrows():
    # riquadro dell'area in pixel di V; cellule GitHub nel riquadro (allargato di 250 um per la stima locale)
    m = 250.0 / UPP
    inb = (U[:, 0] >= a.c0) & (U[:, 0] < a.c1) & (U[:, 1] >= a.r0) & (U[:, 1] < a.r1)
    inm = (U[:, 0] >= a.c0 - m) & (U[:, 0] < a.c1 + m) & (U[:, 1] >= a.r0 - m) & (U[:, 1] < a.r1 + m)
    Bl, rl = fit(X[inm], U[inm]); rg = r_um[inb]
    rl_in = (np.hypot(*(U[inb] - np.c_[X[inb], np.ones(inb.sum())] @ Bl).T)) * UPP
    Buse = B if mode == "global" else Bl
    res["areas"][a.area] = {"B": Buse.tolist(), "B_local": Bl.tolist(), "n_cells_box": int(inb.sum()), "n_cells_fit_local": int(inm.sum())}
    rows.append({"area": a.area, "n_cells_box": int(inb.sum()), "cells_per_mm2_github": inb.sum() / ((a.side_px * UPP) / 1000) ** 2,
                 **{f"global_{k}": v for k, v in qstats(rg).items()}, **{f"local_{k}": v for k, v in qstats(rl_in).items()},
                 "local_vs_global_shift_um_at_center": float(np.hypot(*(np.r_[X[inb].mean(0), 1.0] @ (Bl - B))) * UPP)})
json.dump(res, open(f"{OUT}/registration.json", "w"), indent=1)
df = pd.DataFrame(rows); df.to_csv("results/riesame/xen_registration.csv", index=False)
print(json.dumps({k: v for k, v in res.items() if k not in ("areas",)}, indent=1)); print(df.round(4).to_string())
