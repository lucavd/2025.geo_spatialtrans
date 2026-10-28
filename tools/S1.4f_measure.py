#!/usr/bin/env python
"""tools/S1.4f_measure.py - S1.4f: misure del banco del nucleo orientato in parte (definizioni: docs/steps/S1.4f.md). Riusa le funzioni di
tools/S1.4e_measure.py (poly_shape, nn_dir, psi_deg, focal_mask, pair_counts, contacts, raster_of) e l'operatore unico di tools/S1.4c_raster.py.
Modi (root del repo):
  shape cal|B [ncores]  allineamento per classe di distanza dal vicino (< 5, 5-6, 6-8 um) ed eccentricita' mediana, nuclei interni
  pairs B|A|A_dd [nc]   coppie di nuclei a contatto (Chebyshev <= 2 px) per EP nel banco B, per K2, ET, EP nei banchi A e A_dd
  windows               P-1 e P-2 di EP (banco B), come S1.4e
  cp1                   contatto delle coppie giudicate da Luca nei banchi A e A_dd, per braccio
  g4                    cancello G-4 (asse del territorio, peso e uniformita' di F_rho(2 delta) ricalcolati in Python)
Uscite: results/S1.4f/*.csv; file per unita' in /mnt/micron/geo_spatialtrans/S1.4f."""
import importlib.util, json, os, sys, glob
import numpy as np, pandas as pd, shapely
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; SRC = Path(os.environ.get("S14F_PYSRC", ROOT))
_s = importlib.util.spec_from_file_location("m14e", ROOT / "tools/S1.4e_measure.py"); ME = importlib.util.module_from_spec(_s); _s.loader.exec_module(ME)
M = Path("/mnt/micron/geo_spatialtrans/S1.4f"); E14 = Path("/mnt/micron/geo_spatialtrans/S1.4e"); RES = ROOT / "results/S1.4f"
R2B = ROOT / "results/R2b"; R2C = ROOT / "results/R2c"
SEEDS = [f"s{k:02d}" for k in range(1, 11)]; LAWS = ["V0", "V1", "V2"]; ARMS_B = ["K2", "ET"] + LAWS
CLS = {"nn<5": (0, 5), "nn5-6": (5, 6), "nn6-8": (6, 8)}
LAWNAME = {"V0": "constant", "V1": "elongation", "V2": "compression"}

def profile(S, interior):
    """Quota P-3 per classe di distanza dal vicino (nuclei interni, allungamento >= 1.2, psi > 45 gradi) ed eccentricita' mediana."""
    d, j, u = ME.nn_dir(S[:, :2]); psi = ME.psi_deg(S[:, 2], u[:, 0], u[:, 1]); base = interior & (S[:, 3] >= 1.2); out = {}
    for c, (lo, hi) in CLS.items():
        m = base & (d >= lo) & (d < hi); out[c] = (int(m.sum()), float(np.mean(psi[m] > 45)) if m.any() else np.nan)
    return out, float(np.median(S[interior, 4]))

def shape_one(f, cells_f, tag):
    nuc = pd.read_parquet(f); cells = pd.read_parquet(cells_f, columns=["interior"])
    S = np.array([ME.poly_shape(g) for g in shapely.from_wkb(nuc.wkb.to_numpy())]); pr, p4 = profile(S, cells.interior.to_numpy())
    A, roi, sd = os.path.basename(f)[:-8].split("_")
    return [dict(tag=tag, archetype=A, roi_id=roi, seed=sd, cls=c, n=n, p3=p, p4=p4) for c, (n, p) in pr.items()]

def stage_shape(bench, nc):
    from multiprocessing import Pool
    jobs = []
    if bench == "cal":
        for law in LAWS:
            for kd in sorted(glob.glob(str(M / "cal" / "nuc" / law / "k*"))):
                for f in sorted(glob.glob(kd + "/*.parquet")):
                    jobs.append((f, str(M / "cal" / "cells" / os.path.basename(f)), f"{law}_{os.path.basename(kd)}"))
    else:
        for arm in ARMS_B:
            for f in sorted(glob.glob(str(M / "B" / "nuc" / arm / "*.parquet"))): jobs.append((f, str(M / "B" / "cells" / os.path.basename(f)), arm))
    with Pool(nc) as p: r = p.starmap(shape_one, jobs)
    D = pd.DataFrame([x for y in r for x in y]); D.to_csv(RES / f"S1.4f_shape_{bench}.csv", index=False); print(f"[shape {bench}] {len(jobs)} file, {len(D)} righe")

# ---------------------------------------------------------------------------------------------------- pairs
def stage_pairs(bench, nc):
    from multiprocessing import Pool
    rois = ME.rois_table(); jobs = []
    arms = ["EP"] if bench == "B" else ["K2", "ET", "EP"]
    wins = pd.read_csv(R2B / "R2b_windows.csv"); wr = set(wins.archetype + "_" + wins.roi_id)
    for arm in arms:
        od = M / bench / "pairs" / arm; od.mkdir(parents=True, exist_ok=True)
        for f in sorted(glob.glob(str(M / bench / "nuc" / arm / "*.parquet"))):
            b = os.path.basename(f)
            if bench == "B" and "_".join(b.split("_")[:2]) not in wr: continue
            jobs.append((f, str(od / b), rois))
    with Pool(nc) as p: r = p.starmap(ME.pairs_one, jobs)
    print(f"[pairs {bench}] {len(jobs)} unita', {sum(x == 'ok' for x in r)} nuove")

# ---------------------------------------------------------------------------------------------------- windows: P-1, P-2 (EP)
def stage_windows():
    wins = pd.read_csv(R2B / "R2b_windows.csv"); dm = pd.read_csv(R2C / "R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]; p1, p2 = [], []
    for wr in wins.itertuples():
        A = wr.archetype; dmax = float(dm.d_max_um[A]); ex = ME.K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
        for sd in SEEDS:
            key = f"{A}_{wr.roi_id}_{sd}.parquet"
            t = pd.read_parquet(M / "B" / "raster" / "EP" / key); P = t[["cx", "cy"]].to_numpy(); foc = ME.focal_mask(P, wr, dmax, u); c = t.n_c2.values > 0
            p1.append(dict(arm="EP", win_id=wr.win_id, archetype=A, roi=wr.roi_id, seed=sd, n=int(foc.sum()), n_contact=int(c[foc].sum()),
                           n_conflict=int(t.n_conflict.iloc[0]) if len(t) else 0))
            pr = pd.read_parquet(M / "B" / "pairs" / "EP" / key); M_ = int(max(t.label.max(), pr.label_j.max() if len(pr) else 0)) + 1
            cset = np.unique(pr.label_i.to_numpy().astype(np.int64) * M_ + pr.label_j.to_numpy()); ids = t.label.to_numpy().astype(np.int64)
            for c_, (n, k) in ME.pair_counts(P, ids, foc, cset, dmax, M_).items():
                p2.append(dict(arm="EP", win_id=wr.win_id, archetype=A, seed=sd, noise="none", draw=0, cls=c_, n=n, n_contact=k))
    a = pd.DataFrame(p1); a["s1"] = a.n_contact / a.n; a.to_csv(RES / "S1.4f_windows.csv", index=False)
    pd.DataFrame(p2).to_csv(RES / "S1.4f_pairs_B.csv", index=False); print(a.groupby("archetype").s1.mean().round(3).to_string())

# ---------------------------------------------------------------------------------------------------- cp1 (banchi A e A_dd)
def stage_cp1():
    it = pd.read_csv(R2C / "R2c_items_Luca.csv"); a = pd.DataFrame(json.load(open(R2C / "annotations/R2c_answers_Luca.json"))["answers"])
    J = it.merge(a, on="item", validate="1:1"); J = J[J.kind == "pair"]; wins = pd.read_csv(R2B / "R2b_windows.csv").set_index("win_id"); rows = []
    for bench in ["A", "A_dd"]:
        for wid, g in J.groupby("win_id"):
            w = wins.loc[wid]; mp = pd.read_parquet(M / bench / "map" / f"{wid}.parquet"); mm = ME.luca_to_cell(mp)
            for arm in ["K2", "ET", "EP"]:
                for sd in SEEDS:
                    pr = pd.read_parquet(M / bench / "pairs" / arm / f"{w.archetype}_{w.roi_id}_{wid}_{sd}.parquet")
                    Mx = int(max(mp.cell_id.max(), pr.label_j.max() if len(pr) else 0)) + 1; cset = set((pr.label_i.astype(np.int64) * Mx + pr.label_j).tolist())
                    for r in g.itertuples():
                        ci, cj = int(mm.loc[r.i]), int(mm.loc[r.j]); k = min(ci, cj) * Mx + max(ci, cj)
                        rows.append(dict(bench=bench, win_id=wid, archetype=w.archetype, arm=arm, seed=sd, i=r.i, j=r.j, d_um=r.d_um, answer=r.answer, contact=int(k in cset)))
    d = pd.DataFrame(rows); d.to_csv(RES / "S1.4f_pairs_A.csv", index=False)
    print(d[(d.d_um < 6) & (d.answer != "unsure")].groupby(["bench", "archetype", "arm"]).contact.mean().unstack().round(3).to_string())

# ---------------------------------------------------------------------------------------------------- g4
def g4_one(f):
    from scipy.stats import wrapcauchy
    t = pd.read_parquet(f); T = shapely.from_wkb(t.terr_wkb.to_numpy()); S = np.array([ME.poly_shape(g) for g in T])
    th_t, asp_t = S[:, 2], S[:, 3]
    dist = shapely.distance(shapely.points(t.x.to_numpy(), t.y.to_numpy()), shapely.boundary(T))
    law = t.law.iloc[0]; lam = float(t["lambda"].iloc[0])
    w = np.ones(len(t)) if law == "constant" else (asp_t - 1 if law == "elongation" else np.sqrt(t.nucleus_area_target.to_numpy() / np.pi) / dist)
    rho = 1 - np.exp(-lam * w) if np.isfinite(lam) else np.ones(len(t))
    dax = np.abs(((t.terr_axis_pkg.to_numpy() - th_t) + np.pi / 2) % np.pi - np.pi / 2)
    phi = 2 * (((t.nucleus_orient.to_numpy() - th_t) + np.pi / 2) % np.pi - np.pi / 2)
    x = np.mod(phi, 2 * np.pi); ok = rho < 1 - 1e-12; pit = np.where(rho > 1e-12, wrapcauchy.cdf(x, np.clip(rho, 1e-12, 1 - 1e-12)), x / (2 * np.pi))
    A, roi, sd = os.path.basename(f)[:-8].split("_")
    return dict(law=os.path.basename(os.path.dirname(f)), archetype=A, roi_id=roi, n=len(t), max_axis_diff=float(dax.max()),
                max_aspect_diff=float(np.abs(t.terr_aspect_pkg.to_numpy() - asp_t).max()), max_w_diff=float(np.abs(t.w_pkg.to_numpy() - w).max()),
                max_dist_diff=float(np.abs(t.center_dist.to_numpy() - dist).max()), pit=pit[ok])

def stage_g4():
    from multiprocessing import Pool
    from scipy.stats import kstest
    fs = sorted(glob.glob(str(M / "B" / "g4" / "*" / "*.parquet")))
    with Pool(15) as p: r = p.map(g4_one, fs)
    U = pd.DataFrame([{k: v for k, v in x.items() if k != "pit"} for x in r]); U.to_csv(RES / "S1.4f_G4_units.csv", index=False)
    ks = []
    for (law, A), idx in U.groupby(["law", "archetype"]).groups.items():
        pit = np.concatenate([r[i]["pit"] for i in idx]); ks.append(dict(law=law, archetype=A, n=len(pit), ks_p=float(kstest(pit, "uniform").pvalue)))
    K = pd.DataFrame(ks); K.to_csv(RES / "S1.4f_G4_ks.csv", index=False)
    out = pd.DataFrame([dict(check="territory_axis_max_diff_rad", value=float(U.max_axis_diff.max()), limit=1e-9),
                        dict(check="territory_aspect_max_diff", value=float(U.max_aspect_diff.max()), limit=1e-9),
                        dict(check="weight_max_diff", value=float(U.max_w_diff.max()), limit=1e-9),
                        dict(check="ks_uniform_min_p", value=float(K.ks_p.min()), limit=0.01)])
    out["pass_"] = [out.value[0] <= 1e-9, out.value[1] <= 1e-9, out.value[2] <= 1e-9, out.value[3] > 0.01]
    out.to_csv(RES / "S1.4f_G4_check.csv", index=False); print(out.to_string()); print(K.round(4).to_string())
    print("GATE G-4", "PASS" if out.pass_.all() else "FAIL")

if __name__ == "__main__":
    a = sys.argv[1:]; RES.mkdir(parents=True, exist_ok=True)
    if a[0] == "shape": stage_shape(a[1], int(a[2]) if len(a) > 2 else 20)
    elif a[0] == "pairs":
        ME.M = M; stage_pairs(a[1], int(a[2]) if len(a) > 2 else 20)
    elif a[0] == "windows": stage_windows()
    elif a[0] == "cp1": stage_cp1()
    elif a[0] == "g4": stage_g4()
    else: raise SystemExit("modo sconosciuto")
