#!/usr/bin/env python
"""tools/S1.4e_measure.py - S1.4e: misure del banco del disegno del nucleo (definizioni: docs/steps/S1.4e.md).
Operatore unico: raster di tools/S1.4c_raster.py (stessa rasterize(), contatto = Chebyshev <= 2 px). Modi (root del repo):
  pairs B|A [ncores]  coppie di nuclei a contatto (contact_pairs di S1.4c) per ogni unita' dei ROI con finestre (B) o di ogni finestra (A)
  windows             P-1 (frazione di focali a contatto, come tools/S1.5e_windows.py), P-2 e CP-2 (coppie con i centri a meno di 6 um)
  shape [ncores]      P-3 (allineamento) e P-4 (eccentricita' mediana) dai momenti dei poligoni, tutte le unita' del banco B
  cp1                 CP-1: contatto delle coppie giudicate da Luca nel banco A, per braccio
  g4                  cancello G-4 (seconda implementazione: skimage sul raster; contatto di coppia contro n_c2)
Uscite: results/S1.4e/*.csv; file per unita' in /mnt/micron/geo_spatialtrans/S1.4e/{A,B}/pairs/<braccio>/."""
import importlib.util, json, os, sys, glob
import numpy as np, pandas as pd, shapely
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("s14c_raster", ROOT / "tools/S1.4c_raster.py"); RAS = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(RAS)
_spec = importlib.util.spec_from_file_location("r2c_cand", ROOT / "tools/R2c_candidates.py"); K = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(K)
M = Path("/mnt/micron/geo_spatialtrans/S1.4e"); RES = ROOT / "results/S1.4e"; R2B = ROOT / "results/R2b"; R2C = ROOT / "results/R2c"
ARMS = ["K2", "ER", "ET"]; SEEDS = [f"s{k:02d}" for k in range(1, 11)]
CHEB = 2                                   # contatto: Chebyshev <= 2 px (operatore di S1.4c)
D_P2 = 6.0                                 # P-2: coppie con i centri a meno di 6 um
CLASSES = [(0, 6, "lt6"), (6, 10, "6to10")]
SIGMA_CP2 = {"A1": 0.672, "A2": 1.247, "A3": 1.552, "A4": 0.504}   # E6: deviazione standard per asse, punti di Luca contro maschera
P3_TISSUES = ["A1", "A2", "A3", "A4", "A6"]

# ---------------------------------------------------------------------------------------------------- geometria
def ring_moments(xy):
    """Momenti con segno di un anello (shoelace): area, Sx, Sy, Ixx, Iyy, Ixy (integrali di 1, x, y, x^2, y^2, xy)."""
    x, y = xy[:-1, 0], xy[:-1, 1]; x2, y2 = xy[1:, 0], xy[1:, 1]; cr = x * y2 - x2 * y
    return np.array([cr.sum() / 2, ((x + x2) * cr).sum() / 6, ((y + y2) * cr).sum() / 6, ((x * x + x * x2 + x2 * x2) * cr).sum() / 12,
                     ((y * y + y * y2 + y2 * y2) * cr).sum() / 12, ((x * y2 + 2 * x * y + 2 * x2 * y2 + x2 * y) * cr).sum() / 24])

def poly_shape(g):
    """Centroide, angolo dell'asse maggiore (rad, da +x verso +y), allungamento (sqrt degli autovalori) ed eccentricita' di un poligono."""
    s = np.zeros(6); parts = shapely.get_parts(g); x0, y0 = parts[0].exterior.coords[0]      # D-3: una sola origine per tutte le parti
    for p in parts:
        p = shapely.geometry.polygon.orient(p, 1.0)
        for r in [p.exterior] + list(p.interiors):
            xy = np.asarray(r.coords) - (x0, y0); s += ring_moments(xy)
    a = s[0]; cx, cy = s[1] / a, s[2] / a
    mxx = s[3] / a - cx * cx; myy = s[4] / a - cy * cy; mxy = s[5] / a - cx * cy
    th = 0.5 * np.arctan2(2 * mxy, mxx - myy)
    t = np.sqrt(((mxx - myy) / 2) ** 2 + mxy ** 2); l1 = (mxx + myy) / 2 + t; l2 = max((mxx + myy) / 2 - t, 1e-15)
    return cx + x0, cy + y0, th, float(np.sqrt(l1 / l2)), float(np.sqrt(max(0.0, 1 - l2 / l1)))

def psi_deg(theta, ux, uy):
    """Angolo (gradi, 0-90) fra l'asse maggiore (direzione cos t, sin t) e la direzione (ux, uy) del vicino: 90 = di fianco."""
    return np.degrees(np.arccos(np.clip(np.abs(np.cos(theta) * ux + np.sin(theta) * uy), 0, 1)))

def nn_dir(P):
    d, j = cKDTree(P).query(P, k=2); d, j = d[:, 1], j[:, 1]; u = (P[j] - P) / d[:, None]; return d, j, u

# ---------------------------------------------------------------------------------------------------- pairs
def rois_table():
    return pd.read_csv(ROOT / "results/R2/R2_rois_checked.csv")

def raster_of(f, rois):
    A, roi = os.path.basename(f).split("_")[:2]; rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]
    upp = float(rw.um_per_px); side = int(rw.side_px); t = pd.read_parquet(f)
    lab, ncf = RAS.rasterize(shapely.from_wkb(t["wkb"].to_numpy()), t["idx"].to_numpy(), (side, side), upp)
    return lab, ncf, upp

def contacts(lab):
    """Coppie di etichette a contatto: Chebyshev <= CHEB px (contact_pairs di S1.4c)."""
    return RAS.contact_pairs(lab, CHEB)

def pairs_one(f, out, rois):
    if os.path.exists(out): return "skip"
    lab, ncf, upp = raster_of(f, rois); pr = contacts(lab); pr["n_conflict"] = ncf
    pr.to_parquet(out + ".tmp", index=False); os.replace(out + ".tmp", out); return "ok"

def stage_pairs(bench, nc):
    from multiprocessing import Pool
    rois = rois_table(); wins = pd.read_csv(R2B / "R2b_windows.csv"); wr = set(wins.archetype + "_" + wins.roi_id)
    jobs = []
    for arm in ARMS:
        od = M / bench / "pairs" / arm; od.mkdir(parents=True, exist_ok=True)
        for f in sorted(glob.glob(str(M / bench / "nuc" / arm / "*.parquet"))):
            b = os.path.basename(f)
            if bench == "B" and "_".join(b.split("_")[:2]) not in wr: continue
            jobs.append((f, str(od / b), rois))
    with Pool(nc) as p: r = p.starmap(pairs_one, jobs)
    print(f"[pairs {bench}] {len(jobs)} unita', {sum(x == 'ok' for x in r)} nuove")

# ---------------------------------------------------------------------------------------------------- windows: P-1, P-2, CP-2
def focal_mask(P, wr, dmax, u):
    x0, y0 = float(wr.x0_um), float(wr.y0_um); s = int(wr.side_px) * float(wr.um_per_px)
    de = np.minimum.reduce([P[:, 0] - x0, x0 + s - P[:, 0], P[:, 1] - y0, y0 + s - P[:, 1]]); keep = de >= dmax
    if u is not None: keep &= shapely.distance(shapely.points(P), u) >= dmax
    return keep

def pair_counts(P, lab_ids, foc, cset, dmax, M_):
    """Coppie con almeno un estremo focale e distanza <= dmax: per classe, numero e numero a contatto (contatto da cset, chiavi i*M+j)."""
    t = cKDTree(P); pr = t.query_pairs(dmax, output_type="ndarray")
    if not len(pr): return {c: (0, 0) for _, _, c in CLASSES}
    pr = pr[foc[pr[:, 0]] | foc[pr[:, 1]]]; d = np.hypot(*(P[pr[:, 0]] - P[pr[:, 1]]).T)
    li, lj = lab_ids[pr[:, 0]], lab_ids[pr[:, 1]]; key = np.minimum(li, lj) * M_ + np.maximum(li, lj)
    con = np.isin(key, cset); out = {}
    for lo, hi, c in CLASSES:
        m = (d >= lo) & (d < hi); out[c] = (int(m.sum()), int(con[m].sum()))
    return out

def stage_windows():
    wins = pd.read_csv(R2B / "R2b_windows.csv"); dm = pd.read_csv(R2C / "R2c_dmax.csv").set_index("archetype")
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]; rng = np.random.default_rng(20261028)
    p1, p2 = [], []
    for arm in ARMS:
        for wr in wins.itertuples():
            A = wr.archetype; dmax = float(dm.d_max_um[A]); ex = K.exclusions_um(wjson, wr); u = shapely.union_all(ex) if ex else None
            for sd in SEEDS:
                key = f"{A}_{wr.roi_id}_{sd}.parquet"
                t = pd.read_parquet(M / "B" / "raster" / arm / key); P = t[["cx", "cy"]].to_numpy(); foc = focal_mask(P, wr, dmax, u)
                c = t.n_c2.values > 0
                p1.append(dict(arm=arm, win_id=wr.win_id, archetype=A, roi=wr.roi_id, seed=sd, n=int(foc.sum()), n_contact=int(c[foc].sum()),
                               n_conflict=int(t.n_conflict.iloc[0]) if len(t) else 0))
                pr = pd.read_parquet(M / "B" / "pairs" / arm / key); M_ = int(max(t.label.max(), pr.label_j.max() if len(pr) else 0)) + 1
                cset = np.unique(pr.label_i.to_numpy().astype(np.int64) * M_ + pr.label_j.to_numpy()); ids = t.label.to_numpy().astype(np.int64)
                for c_, (n, k) in pair_counts(P, ids, foc, cset, dmax, M_).items():
                    p2.append(dict(arm=arm, win_id=wr.win_id, archetype=A, seed=sd, noise="none", draw=0, cls=c_, n=n, n_contact=k))
                if A in SIGMA_CP2:                                    # CP-2: rumore dei punti manuali solo sulla distanza
                    for dr in range(10):
                        Pn = P + rng.normal(0, SIGMA_CP2[A], P.shape)
                        tq = cKDTree(Pn); pq = tq.query_pairs(dmax, output_type="ndarray")
                        pq = pq[foc[pq[:, 0]] | foc[pq[:, 1]]]; d = np.hypot(*(Pn[pq[:, 0]] - Pn[pq[:, 1]]).T)
                        li, lj = ids[pq[:, 0]], ids[pq[:, 1]]; con = np.isin(np.minimum(li, lj) * M_ + np.maximum(li, lj), cset)
                        for lo, hi, c_ in CLASSES:
                            m = (d >= lo) & (d < hi)
                            p2.append(dict(arm=arm, win_id=wr.win_id, archetype=A, seed=sd, noise="click", draw=dr + 1, cls=c_, n=int(m.sum()), n_contact=int(con[m].sum())))
    a = pd.DataFrame(p1); a["s1"] = a.n_contact / a.n; a.to_csv(RES / "S1.4e_windows.csv", index=False)
    pd.DataFrame(p2).to_csv(RES / "S1.4e_pairs_B.csv", index=False); print(a.groupby(["arm", "archetype"]).s1.mean().unstack().round(3).to_string())

# ---------------------------------------------------------------------------------------------------- shape: P-3, P-4
def shape_one(arm, key):
    nuc = pd.read_parquet(M / "B" / "nuc" / arm / f"{key}.parquet"); cells = pd.read_parquet(M / "B" / "cells" / f"{key}.parquet")
    G = shapely.from_wkb(nuc.wkb.to_numpy()); S = np.array([poly_shape(g) for g in G])
    P = S[:, :2]; d, j, u = nn_dir(P); it = cells.interior.to_numpy()
    psi = psi_deg(S[:, 2], u[:, 0], u[:, 1]); el = S[:, 3] >= 1.2; m = it & el & (d < 5)
    A, roi, sd = key.split("_")
    return dict(arm=arm, archetype=A, roi_id=roi, seed=sd, n_interior=int(it.sum()), n_p3=int(m.sum()),
                p3=float(np.mean(psi[m] > 45)) if m.any() else np.nan, p4=float(np.median(S[it, 4])),
                frac_elong=float(np.mean(el[it])), aspect_med=float(np.median(S[it, 3])))

def stage_shape(nc):
    from multiprocessing import Pool
    keys = sorted(os.path.basename(f)[:-8] for f in glob.glob(str(M / "B" / "cells" / "*.parquet")))
    with Pool(nc) as p: r = p.starmap(shape_one, [(a, k) for a in ARMS for k in keys])
    d = pd.DataFrame(r); d.to_csv(RES / "S1.4e_shape_B.csv", index=False)
    print(d.groupby(["arm", "archetype"])[["p3", "p4"]].mean().unstack(0).round(3).to_string())

# ---------------------------------------------------------------------------------------------------- cp1 (banco A)
def luca_to_cell(mp):
    """Mappa del banco A: indice del punto di Luca (R2b) -> cell_id del nucleo nato da quel punto."""
    mm = mp[mp.src == "manual"].set_index("luca_idx").cell_id; assert mm.index.is_unique; return mm

def stage_cp1():
    it = pd.read_csv(R2C / "R2c_items_Luca.csv"); a = pd.DataFrame(json.load(open(R2C / "annotations/R2c_answers_Luca.json"))["answers"])
    J = it.merge(a, on="item", validate="1:1"); J = J[J.kind == "pair"]; wins = pd.read_csv(R2B / "R2b_windows.csv").set_index("win_id")
    # D-2: le coppie «non giudicabile» restano nel file (colonna answer); CP-1 pre-registrata le esclude (come R-1 e P-2), in tools/S1.4e_analysis.py
    rows = []
    for wid, g in J.groupby("win_id"):
        w = wins.loc[wid]; mp = pd.read_parquet(M / "A" / "map" / f"{wid}.parquet"); mm = luca_to_cell(mp)
        for arm in ARMS:
            for sd in SEEDS:
                pr = pd.read_parquet(M / "A" / "pairs" / arm / f"{w.archetype}_{w.roi_id}_{wid}_{sd}.parquet")
                Mx = int(max(mp.cell_id.max(), pr.label_j.max() if len(pr) else 0)) + 1; cset = set((pr.label_i.astype(np.int64) * Mx + pr.label_j).tolist())
                for r in g.itertuples():
                    ci, cj = int(mm.loc[r.i]), int(mm.loc[r.j]); k = min(ci, cj) * Mx + max(ci, cj)
                    rows.append(dict(win_id=wid, archetype=w.archetype, arm=arm, seed=sd, i=r.i, j=r.j, d_um=r.d_um, answer=r.answer, contact=int(k in cset)))
    d = pd.DataFrame(rows); d.to_csv(RES / "S1.4e_pairs_A.csv", index=False)
    print(d[d.d_um < D_P2].groupby(["archetype", "arm"]).contact.mean().unstack().round(3).to_string())

# ---------------------------------------------------------------------------------------------------- g4
def g4_one(key, rois):
    from skimage.measure import regionprops
    arm = "ET"; f = M / "B" / "nuc" / arm / f"{key}.parquet"; lab, ncf, upp = raster_of(str(f), rois)
    G = shapely.from_wkb(pd.read_parquet(f).wkb.to_numpy()); S = np.array([poly_shape(g) for g in G]); cells = pd.read_parquet(M / "B" / "cells" / f"{key}.parquet")
    it = cells.interior.to_numpy(); P = S[:, :2]; d, j, u = nn_dir(P)
    th_r = np.full(len(G), np.nan); ecc_r = np.full(len(G), np.nan)
    for p in regionprops(lab):
        th_r[p.label - 1] = np.arctan2(np.cos(p.orientation), np.sin(p.orientation)); ecc_r[p.label - 1] = p.eccentricity   # asse (sin t, cos t) in (x, y)
    ok = it & np.isfinite(th_r) & (S[:, 3] >= 1.2)
    dpsi = np.abs(psi_deg(S[ok, 2], u[ok, 0], u[ok, 1]) - psi_deg(th_r[ok], u[ok, 0], u[ok, 1]))
    return dict(key=key, n=int(ok.sum()), dpsi_median=float(np.median(dpsi)), ecc_poly=float(np.median(S[it, 4])), ecc_raster=float(np.nanmedian(ecc_r[it])))

def stage_g4():
    from multiprocessing import Pool
    rois = rois_table(); keys = [f"{r.archetype}_{r.roi_id}_s01" for r in rois.itertuples()]
    with Pool(15) as p: a = pd.DataFrame(p.starmap(g4_one, [(k, rois) for k in keys]))
    a["decc"] = (a.ecc_poly - a.ecc_raster).abs()
    rows = []                                                      # contatto di coppia contro n_c2 del raster di S1.4c, tutte le unita' con pairs
    for arm in ARMS:
        for f in sorted(glob.glob(str(M / "B" / "pairs" / arm / "*.parquet"))):
            pr = pd.read_parquet(f); t = pd.read_parquet(M / "B" / "raster" / arm / os.path.basename(f), columns=["label", "n_c2"])
            any_c = np.isin(t.label.to_numpy(), np.r_[pr.label_i.to_numpy(), pr.label_j.to_numpy()])
            rows.append(dict(arm=arm, key=os.path.basename(f), n=len(t), agree=int((any_c == (t.n_c2.to_numpy() > 0)).sum())))
    b = pd.DataFrame(rows)
    out = pd.DataFrame([dict(check="psi_median_abs_diff_deg_max_over_units", value=float(a.dpsi_median.max()), limit=2.0, pass_=bool(a.dpsi_median.max() <= 2)),
                        dict(check="ecc_median_abs_diff_max_over_units", value=float(a.decc.max()), limit=0.02, pass_=bool(a.decc.max() <= 0.02)),
                        dict(check="pair_contact_vs_nc2_agreement", value=float(b.agree.sum() / b.n.sum()), limit=1.0, pass_=bool(b.agree.sum() == b.n.sum())),
                        dict(check="pair_contact_units", value=float(len(b)), limit=np.nan, pass_=True)])
    a.to_csv(RES / "S1.4e_G4_units.csv", index=False); out.to_csv(RES / "S1.4e_G4_check.csv", index=False); print(out.to_string())
    print("GATE G-4", "PASS" if out.pass_.all() else "FAIL")

if __name__ == "__main__":
    a = sys.argv[1:]; RES.mkdir(parents=True, exist_ok=True)
    if a[0] == "pairs": stage_pairs(a[1], int(a[2]) if len(a) > 2 else 20)
    elif a[0] == "windows": stage_windows()
    elif a[0] == "shape": stage_shape(int(a[1]) if len(a) > 1 else 20)
    elif a[0] == "cp1": stage_cp1()
    elif a[0] == "g4": stage_g4()
    else: raise SystemExit("modo sconosciuto")
