#!/usr/bin/env python
"""tools/R2c_candidates.py — R2c: d_max, coppie candidate, nuclei focali, check C-1/C-3 e liste degli item da annotare.

Pre-registrazione: results/R2c/R2c_preregistration.md (approvata da Luca 2026-10-16).
- Punti: punti manuali di Luca R2b (results/R2b/R2b_points.csv), coordinate in µm del ROI (x = colonna * upp).
- Contatto dei segmentatori: operatore di S1.4c (tools/S1.4c_raster.py, contact_pairs, Chebyshev <= 2 px) sulle maschere
  native del ROI filtrate come nello stadio cp2 (centroide del nucleo nel tessuto valido).
- d_max(archetipo) = q99.5 della distanza fra i due punti manuali appaiati alle etichette di ogni coppia a contatto
  (cellpose_rgb U stardist_he, coppie con entrambi i punti appaiati a etichette distinte) + 1 µm.
- Candidati: lati di Delaunay dei punti di Luca della finestra con d <= d_max.
- Focali: 15 per finestra fra i punti a distanza >= d_max dal bordo della finestra e dalle esclusioni; rng([SEED, k]),
  k = indice della finestra in R2b_windows.csv.
- Item: coppie candidate con almeno un estremo focale; ripetute 10 % (inserite >= 30 item dopo l'originale); trappole 5 %
  (coppie di punti della stessa finestra con 20 < d <= lato del ritaglio - 6 µm). Federica: gli item delle finestre rater2, ordine con seme proprio.
Output: results/R2c/R2c_{dmax,pairs,focal,coverage,c3,items_Luca,items_Federica}.csv
Uso (root del repo): .venv/bin/python tools/R2c_candidates.py
"""
import importlib.util, json
import numpy as np, pandas as pd
from pathlib import Path
from scipy.spatial import Delaunay
import shapely
from shapely.geometry import Polygon, box

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/R2c"; RES.mkdir(parents=True, exist_ok=True)
R2B = ROOT / "results/R2b"
MASKS = "/mnt/micron/geo_spatialtrans/R2/masks"
SEED = 20261016
N_FOCAL = 15
Q_DMAX, DMAX_PAD = 0.995, 1.0
REPEAT_FRAC, CATCH_FRAC, CATCH_DMIN, REPEAT_GAP = 0.10, 0.05, 20.0, 30
CROP_K, CROP_MIN_UM, CATCH_FIT_UM = 2.5, 30.0, 6.0   # lato del ritaglio = max(2.5 d_max, 30 µm); trappole con d <= lato - 6 µm
SEGS = ["cellpose_rgb", "stardist_he"]
SR_ARCH = {"A4", "A6"}


def load_s14c():
    spec = importlib.util.spec_from_file_location("s14c_raster", ROOT / "tools/S1.4c_raster.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


S = load_s14c()


def cp2_lab(A, roi, seg, rois):
    """Etichette native filtrate come nello stadio cp2 di S1.4c (centroide nel tessuto valido)."""
    rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]
    upp = float(rw.um_per_px); side_um = float(rw.side_px) * upp
    lab = np.load(f"{MASKS}/{A}_{roi}_{seg}_native.npz")["labels"].astype(np.int64)
    v = S.read_valid(A, roi); pxv = side_um / v.shape[1]
    r, c = np.nonzero(lab); l = lab[r, c]; M = int(lab.max()) + 1
    n = np.bincount(l, minlength=M); ids = np.nonzero(n)[0]
    cx = np.bincount(l, c.astype(float), M)[ids] / n[ids] * upp
    cy = np.bincount(l, r.astype(float), M)[ids] / n[ids] * upp
    ci = np.clip((cx / pxv).astype(int), 0, v.shape[1] - 1); ri = np.clip((cy / pxv).astype(int), 0, v.shape[0] - 1)
    keep = np.zeros(M, bool); keep[ids[v[ri, ci]]] = True
    lab[~keep[lab]] = 0
    return lab, upp


def delaunay_edges(xy):
    if len(xy) < 3:
        return np.zeros((0, 2), int)
    e = set()
    for s in Delaunay(xy).simplices:
        for a, b in ((0, 1), (1, 2), (0, 2)):
            i, j = sorted((int(s[a]), int(s[b]))); e.add((i, j))
    return np.array(sorted(e), int)


def exclusions_um(wins_json, wrow):
    """Poligoni di esclusione di Luca (R2b) in µm del ROI: pts immagine - margine -> finestra, * upp + origine."""
    wj = [w for w in wins_json if w["win_id"] == wrow.win_id][0]
    m = float(wj["margin_px"]); upp = float(wrow.um_per_px)
    ex = wj.get("exclusions") or []
    if isinstance(ex, str):
        import ast; ex = ast.literal_eval(ex)
    out = []
    for e in ex:
        p = (np.asarray(e["pts"], float) - m) * upp + np.array([float(wrow.x0_um), float(wrow.y0_um)])
        if len(p) >= 3:
            out.append(Polygon(p).buffer(0))
    return out


def seg_tables(wins, rois):
    """Per ROI x segmentatore: insieme delle coppie a contatto (label_i < label_j) e S-1 per etichetta; C-3 contro cp2."""
    need = sorted(set((r.archetype, r.roi_id) for r in wins.itertuples()))
    pairs, c3 = {}, []
    for A, roi in need:
        segs = SEGS + (["spaceranger"] if A in SR_ARCH else [])
        for seg in segs:
            lab, upp = cp2_lab(A, roi, seg, rois)
            pr = S.contact_pairs(lab, 2)
            pairs[(A, roi, seg)] = set(zip(pr.label_i.astype(int), pr.label_j.astype(int)))
            if seg in SEGS:                                      # C-3: identita' con lo stadio cp2 di S1.4c
                cp = pd.read_parquet(f"/mnt/micron/geo_spatialtrans/S1.4c/cp2/{A}_{roi}_{seg}_cp2.parquet")
                in_pair = set(pr.label_i.astype(int)) | set(pr.label_j.astype(int))
                a = set(cp.label[cp.n_c2 > 0].astype(int)); b = in_pair
                c3.append(dict(archetype=A, roi=roi, seg=seg, n_labels=len(cp), n_contact_cp2=len(a), n_contact_r2c=len(b),
                               n_diff=len(a ^ b), labels_equal=bool(set(cp.label.astype(int)) == set(np.unique(lab[lab > 0])))))
    return pairs, pd.DataFrame(c3)


def main():
    wins = pd.read_csv(R2B / "R2b_windows.csv")
    rois = S.rois_table()
    pts = pd.read_csv(R2B / "R2b_points.csv"); pts = pts[pts.rater == "Luca"].copy()
    wjson = json.load(open(R2B / "annotations/R2b_annotations_Luca.json"))["windows"]
    seg_pairs, c3 = seg_tables(wins, rois)
    c3.to_csv(RES / "R2c_c3.csv", index=False)

    # ---- coppie di punti manuali a contatto in ciascun segmentatore; lati di Delaunay
    rows = []
    for wr in wins.itertuples():
        g = pts[pts.win_id == wr.win_id].sort_values("idx")
        xy = g[["x_um", "y_um"]].to_numpy(); idx = g["idx"].to_numpy()
        E = delaunay_edges(xy)
        dset = set((int(idx[i]), int(idx[j])) for i, j in E)
        segs = SEGS + (["spaceranger"] if wr.archetype in SR_ARCH else [])
        lab = {seg: dict() for seg in segs}
        for seg in segs:
            mm = g[g[f"match_{seg}"].astype(str) == "True"]
            lab[seg] = dict(zip(mm["idx"].astype(int), mm[f"label_{seg}"].astype(int)))
        P = {int(k): v for k, v in zip(idx, xy)}
        allp = set(dset)
        for seg in segs:                                          # coppie a contatto fra punti appaiati
            inv = {}
            for k, l in lab[seg].items():
                inv.setdefault(l, []).append(k)
            for li, lj in seg_pairs[(wr.archetype, wr.roi_id, seg)]:
                if li in inv and lj in inv:
                    for a in inv[li]:
                        for b in inv[lj]:
                            allp.add(tuple(sorted((a, b))))
        for a, b in sorted(allp):
            row = dict(win_id=wr.win_id, archetype=wr.archetype, roi_id=wr.roi_id, i=a, j=b,
                       d_um=float(np.hypot(*(P[a] - P[b]))), delaunay=(a, b) in dset)
            for seg in ["cellpose_rgb", "stardist_he", "spaceranger"]:
                if seg not in segs:
                    row[f"both_{seg}"] = np.nan; row[f"contact_{seg}"] = np.nan; continue
                la, lb = lab[seg].get(a), lab[seg].get(b)
                both = la is not None and lb is not None and la != lb
                row[f"both_{seg}"] = both
                row[f"contact_{seg}"] = (tuple(sorted((la, lb))) in seg_pairs[(wr.archetype, wr.roi_id, seg)]) if both else np.nan
            rows.append(row)
    pr = pd.DataFrame(rows)

    # ---- d_max per archetipo (regola pre-registrata)
    dm = []
    for A, g in pr.groupby("archetype"):
        cu = g[(g.contact_cellpose_rgb == True) | (g.contact_stardist_he == True)]
        q = float(np.quantile(cu.d_um, Q_DMAX)) if len(cu) else np.nan
        dm.append(dict(archetype=A, n_contact_pairs_union=len(cu), q995_um=q, d_max_um=q + DMAX_PAD,
                       n_cp=int((g.contact_cellpose_rgb == True).sum()), n_sd=int((g.contact_stardist_he == True).sum())))
    dm = pd.DataFrame(dm); dm["crop_um"] = np.maximum(CROP_K * dm.d_max_um, CROP_MIN_UM)
    dm.to_csv(RES / "R2c_dmax.csv", index=False)
    crop = dict(zip(dm.archetype, dm.crop_um))
    dmax = dict(zip(dm.archetype, dm.d_max_um))
    pr["candidate"] = pr.delaunay & (pr.d_um <= pr.archetype.map(dmax))

    # ---- C-1 copertura
    cov = []
    for (A, seg), _ in [((A, s), 0) for A in sorted(dmax) for s in ["cellpose_rgb", "stardist_he", "spaceranger"]]:
        g = pr[(pr.archetype == A) & (pr[f"contact_{seg}"] == True)]
        if seg == "spaceranger" and A not in SR_ARCH:
            continue
        f = float(g.candidate.mean()) if len(g) else np.nan
        verdict = "NA" if not len(g) else ("PASS" if f >= 0.99 else "WARN" if f >= 0.97 else "FAIL")
        cov.append(dict(archetype=A, seg=seg, n_contact_pairs=len(g), n_not_delaunay=int((~g.delaunay).sum()),
                        n_beyond_dmax=int((g.d_um > dmax[A]).sum()), coverage=f, verdict=verdict,
                        role="descrittivo" if seg == "spaceranger" else "pre-registrato"))
    pd.DataFrame(cov).to_csv(RES / "R2c_coverage.csv", index=False)

    # ---- nuclei focali
    foc = []
    for k, wr in enumerate(wins.itertuples()):
        g = pts[pts.win_id == wr.win_id].sort_values("idx")
        upp = float(wr.um_per_px); x0, y0 = float(wr.x0_um), float(wr.y0_um); s = int(wr.side_px) * upp
        d = dmax[wr.archetype]
        dist_edge = np.asarray(np.minimum.reduce([g.x_um - x0, x0 + s - g.x_um, g.y_um - y0, y0 + s - g.y_um]), float)
        ex = exclusions_um(wjson, wr)
        if ex:
            u = shapely.union_all(ex)
            dist_ex = shapely.distance(shapely.points(g[["x_um", "y_um"]].to_numpy()), u)
        else:
            dist_ex = np.full(len(g), np.inf)
        elig = (dist_edge >= d) & (dist_ex >= d)
        rng = np.random.default_rng([SEED, k])
        cand_idx = g["idx"].to_numpy()[elig]
        sel = set(rng.choice(cand_idx, size=min(N_FOCAL, len(cand_idx)), replace=False).tolist())
        for i_, de, dx, el in zip(g["idx"].to_numpy(), dist_edge, dist_ex, elig):
            foc.append(dict(win_id=wr.win_id, archetype=wr.archetype, idx=int(i_), dist_edge_um=float(de),
                            dist_excl_um=float(dx), eligible=bool(el), focal=int(i_) in sel, n_eligible=int(elig.sum())))
    foc = pd.DataFrame(foc); foc.to_csv(RES / "R2c_focal.csv", index=False)
    F = set(zip(foc.win_id[foc.focal], foc.idx[foc.focal]))
    pr["focal_i"] = [(w, i) in F for w, i in zip(pr.win_id, pr.i)]
    pr["focal_j"] = [(w, j) in F for w, j in zip(pr.win_id, pr.j)]
    pr["judged"] = pr.candidate & (pr.focal_i | pr.focal_j)
    pr["pair_id"] = pr.win_id + ":" + pr.i.astype(str) + "-" + pr.j.astype(str)
    pr.to_csv(RES / "R2c_pairs.csv", index=False)

    # ---- item: coppie giudicate, trappole, ripetute
    rng = np.random.default_rng([SEED, 999])
    J = pr[pr.judged].copy()
    base = [dict(pair_id=p, win_id=w, i=i, j=j, d_um=d, kind="pair") for p, w, i, j, d in
            zip(J.pair_id, J.win_id, J.i, J.j, J.d_um)]
    n_catch = int(round(CATCH_FRAC * len(J)))
    wcount = J.win_id.value_counts()
    wlist, wprob = wcount.index.to_numpy(), (wcount / wcount.sum()).to_numpy()
    used = set()
    while len(used) < n_catch:
        w = rng.choice(wlist, p=wprob)
        g = pts[pts.win_id == w]
        a, b = rng.choice(g["idx"].to_numpy(), 2, replace=False)
        a, b = int(min(a, b)), int(max(a, b))
        pa = g[g.idx == a][["x_um", "y_um"]].to_numpy()[0]; pb = g[g.idx == b][["x_um", "y_um"]].to_numpy()[0]
        dd = float(np.hypot(*(pa - pb)))
        if CATCH_DMIN < dd <= crop[w[:2]] - CATCH_FIT_UM and (w, a, b) not in used:
            used.add((w, a, b))
            base.append(dict(pair_id=f"{w}:{a}-{b}", win_id=w, i=a, j=b, d_um=dd, kind="catch"))
    rng.shuffle(base)
    rep = rng.choice(len(J), size=int(round(REPEAT_FRAC * len(J))), replace=False)
    reps = [dict(base_row) for base_row in [dict(pair_id=J.pair_id.iloc[r], win_id=J.win_id.iloc[r], i=int(J.i.iloc[r]),
            j=int(J.j.iloc[r]), d_um=float(J.d_um.iloc[r]), kind="repeat") for r in rep]]
    order = list(base)
    for rr in reps:
        o = next(n for n, it in enumerate(order) if it["pair_id"] == rr["pair_id"] and it["kind"] == "pair")
        lo = o + REPEAT_GAP + 1
        pos = int(rng.integers(lo, len(order) + 1)) if lo <= len(order) else len(order)
        order.insert(pos, rr)
    items = pd.DataFrame(order); items.insert(0, "item", np.arange(1, len(items) + 1))
    items.to_csv(RES / "R2c_items_Luca.csv", index=False)
    r2w = set(wins.win_id[wins.rater2])
    it2 = items[items.win_id.isin(r2w)].drop(columns="item").sample(frac=1, random_state=SEED + 2).reset_index(drop=True)
    # ripetute di Federica a >= REPEAT_GAP dall'originale: riposizionate con la stessa regola
    b2 = it2[it2.kind != "repeat"].to_dict("records"); r2 = it2[it2.kind == "repeat"].to_dict("records")
    rng2 = np.random.default_rng([SEED, 2])
    for rr in r2:
        o = next(n for n, it in enumerate(b2) if it["pair_id"] == rr["pair_id"] and it["kind"] == "pair")
        lo = o + REPEAT_GAP + 1
        pos = int(rng2.integers(lo, len(b2) + 1)) if lo <= len(b2) else len(b2)
        b2.insert(pos, rr)
    it2 = pd.DataFrame(b2); it2.insert(0, "item", np.arange(1, len(it2) + 1))
    it2.to_csv(RES / "R2c_items_Federica.csv", index=False)

    print(dm.to_string(index=False))
    print(pd.DataFrame(cov).to_string(index=False))
    print("C-3 n_diff tot", int(c3.n_diff.sum()), "labels_equal", bool(c3.labels_equal.all()))
    print(foc.groupby("archetype").agg(n=("idx", "size"), elig=("eligible", "sum"), focal=("focal", "sum")).to_string())
    print(J.groupby("archetype").size().to_string())
    print("items Luca", len(items), items.kind.value_counts().to_dict(), "| Federica", len(it2), it2.kind.value_counts().to_dict())


if __name__ == "__main__":
    main()
