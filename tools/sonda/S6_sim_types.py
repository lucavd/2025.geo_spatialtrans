"""tools/sonda/S6_sim_types.py - S6 (docs/sonda/REGOLA_S6.md): il braccio Y di S4 (geometria, densita' ed espressione del banale B; tipo della
cellula vera tipizzata piu' vicina) con la sorgente dei tipi che NON e' l'area misurata. Stesso codice di S4 (tools/sonda/S4_sim_arms.py,
argomento type_source), stesso seme di Y in S4: default_rng([20261129, 0, i, 0]), quindi geometria e UMI per cellula identici a Y.
  Yself : sorgente = l'area misurata (controllo: deve rigenerare bit per bit Y di S4, cancello G-SELF)
  Yalt  : tipi dall'area successiva dello stesso tessuto in V (k -> k mod 5 + 1), stesso riferimento locale (quadrati di 750 um)
  YH1   : tipi dall'area omologa di H1 (S3b_sonda/runs/T, tipizzata dalla sua catena); baricentri di H1 portati in H2 con la catena di
          coordinate di tools/sonda/S3b_areas_T.py (array H1 -> hires H1 -> hires H2 -> array H2), poi nel riferimento locale dell'area di V
  Yshuf : tipi delle cellule tipizzate dell'area misurata permutati a caso (composizione dell'area, disposizione distrutta);
          permutazione default_rng([20261124, i])
Uscite in $S6_ROOT/runs/<braccio>/H2_<tessuto>_<k>_k0/ come S4, piu' S6_source.json (sorgente, distanze, composizioni) e cells_src.parquet
con l'etichetta della cellula sorgente nella sua sezione. Uso: .venv/bin/python -B tools/sonda/S6_sim_types.py <W> <bracci,..> [i,..|all]"""
import sys, os, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from multiprocessing import Pool
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S4_sim_arms as S
REPO = os.path.expanduser("~/2025.geo_spatialtrans"); SEED_PERM = 20261124
ROOT = os.environ.get("S6_ROOT", "/mnt/micron/geo_spatialtrans/S6_sonda"); T_DIR = "/mnt/micron/geo_spatialtrans/S3b_sonda/runs/T"
ARMS6 = ("Yself", "Yalt", "YH1", "Yshuf")
S.ROOT = ROOT
for _a in ARMS6: S.ARMS[_a] = dict(S.ARMS["Y"]); S.ARM_ID[_a] = S.ARM_ID["Y"]      # stesso interruttore e stesso seme di Y
_I = json.load(open(f"{REPO}/results/sonda/S3b/overlap.json")); Wh = np.array(_I["H2hires_to_H1hires_colrow"])
co = {s: np.array(_I[s]["array_to_fullres"]).T for s in ("H1", "H2")}; sf = {s: _I[s]["hires_scalef"] for s in ("H1", "H2")}


def um_to_hires(sec, xy):  # copia di tools/sonda/S3b_areas_T.py
    r, c = xy[:, 1] / 8 - 0.5, xy[:, 0] / 8 - 0.5; P = np.c_[r, c, np.ones(len(r))] @ co[sec]; return np.c_[P[:, 1], P[:, 0]] * sf[sec]


def hires_to_um(sec, xy):  # copia di tools/sonda/S3b_areas_T.py
    A = co[sec][:2].T; b = co[sec][2]; rc = np.linalg.solve(A, (xy[:, ::-1] / sf[sec] - b).T).T; return np.c_[(rc[:, 1] + 0.5) * 8, (rc[:, 0] + 0.5) * 8]


def h2_to_h1(xy): return hires_to_um("H1", (Wh @ np.c_[um_to_hires("H2", xy), np.ones(len(xy))].T).T[:, :2])


def h1_to_h2(xy): return hires_to_um("H2", (np.linalg.inv(Wh) @ np.c_[um_to_hires("H1", xy), np.ones(len(xy))].T).T[:, :2])


def area_cells(d, types):
    """cellule vere di un'area esportata: ids, baricentri locali (um), tipo (indice in types o -1) come in S4_sim_arms.run, export.json."""
    b = pd.read_parquet(f"{d}/bins.parquet"); ids, cxy, _ = S.real_cells(b.x_um.to_numpy(), b.y_um.to_numpy(), b.label.to_numpy().astype(np.int64))
    tl = pd.read_parquet(f"{d}/labels_transferred.parquet").set_index("label").cell_type.reindex(ids)
    tix = {t: j for j, t in enumerate(types)}; rt = np.array([tix.get(t, -1) if isinstance(t, str) else -1 for t in tl.to_numpy()], np.int64)
    return ids, cxy, rt, json.load(open(f"{d}/export.json"))


def source(arm, i, types):
    """sorgente dei tipi per l'unita' i: (baricentri nel riferimento locale dell'area misurata, tipi, ids sorgente, descrizione)."""
    tissue, k = S.UNITS[i]; vdir = f"{S.S3}/runs/V/H2_{tissue}_{k}"
    if arm in ("Yself", "Yshuf"):
        ids, cxy, rt, _ = area_cells(vdir, types); d = vdir
        if arm == "Yshuf":
            rt = rt.copy(); t = rt >= 0; rt[t] = np.random.default_rng([SEED_PERM, i]).permutation(rt[t])
    elif arm == "Yalt":
        d = f"{S.S3}/runs/V/H2_{tissue}_{k % 5 + 1}"; ids, cxy, rt, _ = area_cells(d, types)
    elif arm == "YH1":
        d = f"{T_DIR}/H1_{tissue}_{k}"; ids, c1, rt, eT = area_cells(d, types); eV = json.load(open(f"{vdir}/export.json"))
        assert eT["section"] == "H1" and eT["side_um"] == eV["side_um"]
        cxy = h1_to_h2(c1 + np.array([eT["x0"], eT["y0"]])) - np.array([eV["x0"], eV["y0"]])
    else:
        raise ValueError(arm)
    return cxy, rt, ids, d


def run6(arm, i, k=0):
    tissue, kk = S.UNITS[i]; vdir = f"{S.S3}/runs/V/H2_{tissue}_{kk}"; dst = f"{ROOT}/runs/{arm}/H2_{tissue}_{kk}_k{k}"
    if os.path.exists(f"{dst}/S6_source.json"): return f"present {dst}"
    types = list(np.load(f"{S.PAR}/params_{tissue}.npz")["types"]); cxy, rt, sids, sdir = source(arm, i, types)
    keep = {}
    def sim(*a, **kw):
        out = S._simulate_orig(*a, **kw, type_source=(cxy, rt)); keep["src"] = out[4]; keep["ctype"] = out[3]
        return out[:5] + (sids,) + out[6:]                                # ids della sorgente: S4 scrive ids[src] in cells_src.parquet
    S.simulate = sim
    try:
        msg = S.run(arm, i, k)
    finally:
        S.simulate = S._simulate_orig
    src = keep["src"]; n = len(src)
    pd.DataFrame({"label": np.arange(1, n + 1), "src_real_label": sids[src], "src_dir": sdir}).to_parquet(f"{dst}/cells_src.parquet", index=False)
    b = pd.read_parquet(f"{vdir}/bins.parquet"); it = b.in_tissue.to_numpy().astype(bool); xy = np.c_[b.x_um.to_numpy(), b.y_um.to_numpy()][it]
    typed = rt >= 0; dist, _ = cKDTree(cxy[typed]).query(xy); side = json.load(open(f"{vdir}/export.json"))["side_um"]
    inside = (cxy[:, 0] >= 0) & (cxy[:, 0] <= side) & (cxy[:, 1] >= 0) & (cxy[:, 1] <= side)
    _, _, rtv, _ = area_cells(vdir, types)
    comp = lambda r: {types[t]: float(np.mean(r[r >= 0] == t)) for t in range(len(types))}
    info = {"arm": arm, "unit": os.path.basename(dst), "source_dir": sdir, "n_source_cells": int(len(rt)), "n_source_typed": int(typed.sum()),
            "share_source_inside_square": float(inside.mean()), "bin_to_source_dist_um": {"p50": float(np.median(dist)), "p95": float(np.percentile(dist, 95)),
            "share_gt_20um": float(np.mean(dist > 20))}, "composition_source_typed": comp(rt), "composition_measured_area_typed": comp(rtv),
            "composition_simulated": {types[t]: float(np.mean(keep["ctype"] == t)) for t in range(len(types))},
            "seed_sim": [S.SEED_S4, S.ARM_ID[arm], i, k], "seed_perm": [SEED_PERM, i] if arm == "Yshuf" else None}
    json.dump(info, open(f"{dst}/S6_source.json", "w"), indent=1)
    return msg


S._simulate_orig = S.simulate


def _one(a): return run6(*a)


if __name__ == "__main__":
    W = int(sys.argv[1]); arms = sys.argv[2].split(","); assert all(a in ARMS6 for a in arms)
    ii = list(range(10)) if len(sys.argv) < 4 or sys.argv[3] == "all" else [int(v) for v in sys.argv[3].split(",")]
    with Pool(W) as p:
        for r in p.imap_unordered(_one, [(a, i, 0) for a in arms for i in ii]): print(r, flush=True)
