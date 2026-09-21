"""tools/R2b_match.py — R2b: appaiamento verita' manuale <-> segmentazioni R2 e metriche per finestra/archetipo.

Input : results/R2b/annotations/*.json   (export dell'annotatore HTML; uno o piu' annotatori)
        results/R2b/R2b_windows.csv, results/R2b/R2b_targets.csv
        /mnt/micron/geo_spatialtrans/R2/masks/<A>_<roi>_<method>_native.npz   (label image R2, scala nativa)
        results/R2/nuclei/<A>_<roi>_<method>_native.parquet                    (flag keep di R2)
Regole (pre-registrate, C-R2b.3): un punto manuale e' appaiato a un oggetto se cade dentro la sua maschera; altrimenti al
centroide piu' vicino entro DMAX_UM = 3 um; assegnazione 1:1 (Hungarian sul costo = 0 dentro maschera, distanza altrimenti).
Contano solo i punti dentro la finestra e fuori dalle esclusioni manuali; contano solo gli oggetti con flag `keep` di R2,
centroide dentro la finestra e fuori dalle esclusioni. Densita' primaria su area netta (finestra - esclusioni).
Output: results/R2b/R2b_manual_windows.csv, R2b_matching.csv, R2b_points.csv, R2b_areas.csv, R2b_interrater.csv,
        R2b_archetype_summary.csv; con --selftest: results/R2b/R2b_selftest.csv (PASS/FAIL)
Uso: .venv/bin/python tools/R2b_match.py [--selftest] [--primary Luca]
"""
import sys, json, glob, numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import linear_sum_assignment
from scipy.stats import chi2
from skimage import measure, draw
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT, load_labels, NUC_DIR

RES = ROOT / "results/R2b"; ANN_DIR = RES / "annotations"
METHODS = ["cellpose_rgb", "stardist_he", "spaceranger"]
DMAX_UM = 3.0; PAD_UM = 10.0; BIG = 1e6
TARGET_TOL_UM = 15.0

# ----------------------------------------------------------------------------------------------------- geometria di base
def poly_mask(polys, side):
    """Maschera booleana (side x side) dell'unione dei poligoni (coordinate finestra, px)."""
    m = np.zeros((side, side), bool)
    for pts in polys:
        if len(pts) < 3: continue
        # skimage testa i punti a coordinate intere: spostando i vertici di -0.5 il pixel (i, j) e' dentro se il suo CENTRO (i+0.5, j+0.5)
        # e' dentro il poligono -> area della maschera = area geometrica (+- bordo), coerente con il test floor(x) sui punti.
        rr, cc = draw.polygon(np.array(pts)[:, 1] - 0.5, np.array(pts)[:, 0] - 0.5, shape=m.shape); m[rr, cc] = True
    return m

def shoelace(pts):
    x, y = np.array(pts).T
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

def objects_table(labels):
    """Centroide (riga, colonna) e area (px) di ogni label dell'immagine completa."""
    p = measure.regionprops_table(labels, properties=("label", "centroid", "area"))
    return pd.DataFrame(p).rename(columns={"centroid-0": "cy", "centroid-1": "cx", "area": "area_px"})

def match_points(points, objs, labels_crop, upp, dmax_um=DMAX_UM):
    """points: (N,2) x,y px finestra; objs: DataFrame con label, cx, cy (px finestra). labels_crop: label image della finestra.
    Ritorna array (N) con l'indice riga di objs appaiato o -1, e il costo (0 = dentro maschera, altrimenti distanza in um)."""
    N, M = len(points), len(objs)
    assign = np.full(N, -1); cost_out = np.full(N, np.nan)
    if N == 0 or M == 0: return assign, cost_out
    cost = np.full((N, M), BIG)
    lab_idx = {int(l): j for j, l in enumerate(objs.label.values)}
    side = labels_crop.shape[0]
    for i, (x, y) in enumerate(points):
        d = np.hypot(objs.cx.values - x, objs.cy.values - y) * upp
        cost[i, d <= dmax_um] = d[d <= dmax_um]
        xi, yi = int(np.floor(x)), int(np.floor(y))
        if 0 <= xi < side and 0 <= yi < side:
            l = int(labels_crop[yi, xi])
            if l and l in lab_idx: cost[i, lab_idx[l]] = 0.0
    ri, cj = linear_sum_assignment(cost)
    for i, j in zip(ri, cj):
        if cost[i, j] < BIG: assign[i] = j; cost_out[i] = cost[i, j]
    return assign, cost_out

def window_objects(labels_roi, objs_full, c0, r0, side, upp, keep_labels=None):
    """Oggetti con centroide nella finestra allargata di PAD_UM (candidati all'appaiamento); colonna in_win = centroide dentro."""
    pad = PAD_UM / upp
    o = objs_full[(objs_full.cx >= c0 - pad) & (objs_full.cx < c0 + side + pad) & (objs_full.cy >= r0 - pad) & (objs_full.cy < r0 + side + pad)].copy()
    o["cx"] -= c0; o["cy"] -= r0
    o["in_win"] = (o.cx >= 0) & (o.cx < side) & (o.cy >= 0) & (o.cy < side)
    o["keep"] = True if keep_labels is None else o.label.isin(keep_labels)
    o["area_um2"] = o.area_px * upp ** 2
    return o.reset_index(drop=True)

def poisson_ci(n, area, alpha=0.05):
    lo = chi2.ppf(alpha / 2, 2 * n) / 2 if n > 0 else 0.0; hi = chi2.ppf(1 - alpha / 2, 2 * n + 2) / 2
    return lo / area, hi / area

# ----------------------------------------------------------------------------------------------------- una finestra
def analyse_window(wj, wrow, rater, label_sets, keep_sets, targets=None):
    """wj: finestra del JSON; wrow: riga di R2b_windows.csv; label_sets: {method: (labels_roi, objs_full)}; keep_sets: {method: set(label)}.
    Ritorna (riga manuale, righe matching per metodo, righe punti, righe aree)."""
    m, side, upp = int(wj["margin_px"]), int(wj["side_px"]), float(wj["um_per_px"])
    c0, r0 = int(wrow.c0_roi), int(wrow.r0_roi)
    pts_all = np.array(wj["points"], float).reshape(-1, 2) - m                       # coordinate finestra
    inside = (pts_all[:, 0] >= 0) & (pts_all[:, 0] < side) & (pts_all[:, 1] >= 0) & (pts_all[:, 1] < side) if len(pts_all) else np.zeros(0, bool)
    excl = [{"type": e["type"], "pts": (np.array(e["pts"], float) - m).tolist()} for e in wj.get("exclusions", [])]
    emask = poly_mask([e["pts"] for e in excl], side)
    def in_excl(xy):
        xi = np.clip(np.floor(xy[:, 0]).astype(int), 0, side - 1); yi = np.clip(np.floor(xy[:, 1]).astype(int), 0, side - 1)
        return emask[yi, xi]
    pts = pts_all[inside]; pe = in_excl(pts) if len(pts) else np.zeros(0, bool)
    pts_use = pts[~pe]
    area_gross = (side * upp) ** 2 / 1e6; area_net = area_gross - emask.sum() * upp ** 2 / 1e6
    n_decl = wj.get("n_declared")                                                 # solo annotatori automatici (R2b_import_model)
    count_only = bool(len(pts_use) == 0 and n_decl)                               # modello con conteggio ma senza posizioni
    n_manual = int(n_decl) if count_only else len(pts_use)
    row = dict(win_id=wj["win_id"], archetype=wj["archetype"], roi_id=wj["roi_id"], rater=rater, side_um=side * upp, area_gross_mm2=area_gross,
               area_net_mm2=area_net, excl_frac=emask.mean(), n_excl_unreadable=sum(e["type"] == "unreadable" for e in excl),
               n_excl_off=sum(e["type"] != "unreadable" for e in excl), n_points_total=len(pts_all), n_points_in_window=int(inside.sum()),
               n_points_in_excl=int(pe.sum()), n_manual=n_manual, n_declared=n_decl, count_only=count_only,
               density_manual_net=n_manual / area_net if area_net > 0 else np.nan,
               density_manual_gross=(n_manual if count_only else len(pts)) / area_gross, n_polygons=len(wj.get("polygons", [])), done=bool(wj.get("done")), note=wj.get("note", ""))
    row["ci_lo"], row["ci_hi"] = poisson_ci(n_manual, area_net) if area_net > 0 else (np.nan, np.nan)
    prow = pd.DataFrame({"win_id": wj["win_id"], "rater": rater, "idx": np.arange(len(pts_use)), "x_win_px": pts_use[:, 0], "y_win_px": pts_use[:, 1],
                         "x_um": float(wrow.x0_um) + pts_use[:, 0] * upp, "y_um": float(wrow.y0_um) + pts_use[:, 1] * upp})
    mrows, arows = [], []
    polys = [(np.array(p, float) - m) for p in wj.get("polygons", [])]
    for k, pg in enumerate(polys):
        if len(pg) < 3: continue
        cen = pg.mean(axis=0); a = dict(win_id=wj["win_id"], archetype=wj["archetype"], rater=rater, poly_idx=k, n_vertices=len(pg),
                                       cx_win_px=cen[0], cy_win_px=cen[1], area_manual_um2=shoelace(pg) * upp ** 2)
        a["eq_diam_manual_um"] = 2 * np.sqrt(a["area_manual_um2"] / np.pi)
        if targets is not None and len(targets):
            d = np.hypot(targets.x_win_px.values - cen[0], targets.y_win_px.values - cen[1]) * upp
            a["target"] = int(targets.target.values[d.argmin()]) if d.min() <= TARGET_TOL_UM else -1; a["target_dist_um"] = float(d.min())
        arows.append(a)
    for meth, (labels_roi, objs_full) in label_sets.items():
        crop = labels_roi[r0:r0 + side, c0:c0 + side]
        o = window_objects(labels_roi, objs_full, c0, r0, side, upp, keep_sets.get(meth))
        oe = in_excl(o[["cx", "cy"]].values) if len(o) else np.zeros(0, bool); o["in_excl"] = oe
        cand = o[o.keep].reset_index(drop=True)                                   # candidati: keep, anche fuori finestra (bordo)
        assign, cost = match_points(pts_use, cand, crop, upp)
        matched_obj = set(cand.label.values[assign[assign >= 0]]) if len(cand) else set()
        counted = cand[cand.in_win & ~cand.in_excl]
        tp = int(counted.label.isin(matched_obj).sum())                           # appaiamenti con oggetti conteggiati
        tp_any = int((assign >= 0).sum())                                          # punti appaiati (anche a oggetti di bordo)
        n_obj = len(counted)
        mrows.append(dict(win_id=wj["win_id"], archetype=wj["archetype"], roi_id=wj["roi_id"], rater=rater, method=meth, n_manual=len(pts_use), n_obj=n_obj,
                          n_obj_unkept=int((~o.keep & o.in_win).sum()), n_obj_in_excl=int((o.keep & o.in_win & o.in_excl).sum()),
                          tp=tp, tp_points=tp_any, recall=tp_any / len(pts_use) if len(pts_use) else np.nan, precision=tp / n_obj if n_obj else np.nan,
                          ghost_frac=1 - tp / n_obj if n_obj else np.nan, density_obj_net=n_obj / area_net if area_net > 0 else np.nan,
                          median_cost_um=float(np.nanmedian(cost)) if tp_any else np.nan, frac_inside_mask=float(np.mean(cost[assign >= 0] == 0)) if tp_any else np.nan))
        prow[f"match_{meth}"] = assign >= 0
        prow[f"label_{meth}"] = np.where(assign >= 0, cand.label.values[np.maximum(assign, 0)] if len(cand) else 0, 0)
        for a in arows:                                                            # area dell'oggetto che contiene il centroide del poligono
            xi, yi = int(np.floor(a["cx_win_px"])), int(np.floor(a["cy_win_px"]))
            l = int(crop[yi, xi]) if 0 <= xi < side and 0 <= yi < side else 0
            a[f"area_{meth}_um2"] = float(o.loc[o.label == l, "area_um2"].iloc[0]) if l and (o.label == l).any() else np.nan
    return row, mrows, prow, arows

# ----------------------------------------------------------------------------------------------------- inter-annotatore
def interrater(points_df, wins):
    rows = []
    for win_id, g in points_df.groupby("win_id"):
        raters = sorted(g.rater.unique())
        if len(raters) < 2: continue
        upp = float(wins.loc[win_id, "um_per_px"])
        for ia in range(len(raters)):
            for ib in range(ia + 1, len(raters)):
                a = g[g.rater == raters[ia]][["x_win_px", "y_win_px"]].values; b = g[g.rater == raters[ib]][["x_win_px", "y_win_px"]].values
                r = dict(win_id=win_id, archetype=win_id.split("_")[0], rater_a=raters[ia], rater_b=raters[ib], n_a=len(a), n_b=len(b),
                         rel_diff=abs(len(a) - len(b)) / max(1, (len(a) + len(b)) / 2))
                for dmax in (3.0, 5.0):
                    if len(a) and len(b):
                        d = np.hypot(a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1]) * upp
                        cost = np.where(d <= dmax, d, BIG); ri, cj = linear_sum_assignment(cost); tp = int((cost[ri, cj] < BIG).sum())
                    else: tp = 0
                    r[f"tp_{dmax:g}um"] = tp; r[f"f1_{dmax:g}um"] = 2 * tp / (len(a) + len(b)) if (len(a) + len(b)) else np.nan
                rows.append(r)
    return pd.DataFrame(rows)

# ----------------------------------------------------------------------------------------------------- pipeline reale
def load_label_sets(arch, roi_id, cache):
    key = (arch, roi_id)
    if key not in cache:
        sets, keeps = {}, {}
        for meth in METHODS:
            try: labels = load_labels(arch, roi_id, meth, 1)
            except FileNotFoundError: continue
            sets[meth] = (labels, objects_table(labels))
            f = NUC_DIR / f"{arch}_{roi_id}_{meth}_native.parquet"
            if f.exists():
                d = pd.read_parquet(f, columns=["label", "keep"]); keeps[meth] = set(d.label[d.keep].astype(int))
        cache[key] = (sets, keeps)
    return cache[key]

def main(primary="Luca"):
    wins = pd.read_csv(RES / "R2b_windows.csv").set_index("win_id")
    targets = pd.read_csv(RES / "R2b_targets.csv")
    files = sorted(glob.glob(str(ANN_DIR / "*.json"))); assert files, f"nessun JSON in {ANN_DIR}"
    cache = {}; man, mat, pts, areas = [], [], [], []
    seen = set()
    for f in files:
        obj = json.load(open(f)); rater = obj.get("rater", Path(f).stem)
        for wj in obj["windows"]:
            if (wj["win_id"], rater) in seen: continue                             # file piu' recente vince (ordine alfabetico = data)
            seen.add((wj["win_id"], rater))
            if wj["n_points_in_window"] == 0 and not wj.get("polygons") and not wj.get("n_declared"): continue  # non annotata
            wrow = wins.loc[wj["win_id"]]
            sets, keeps = load_label_sets(wj["archetype"], wj["roi_id"], cache)
            row, mrows, prow, arows = analyse_window(wj, wrow, rater, sets, keeps, targets[targets.win_id == wj["win_id"]])
            row["source_file"] = Path(f).name; man.append(row); mat += mrows; pts.append(prow); areas += arows
            print(f"{wj['win_id']:7s} {rater:10s} n_manual={row['n_manual']:4d} net={row['area_net_mm2']*1e6:8.0f} um2 " +
                  " ".join(f"{r['method'][:2]}: n={r['n_obj']:3d} R={r['recall']:.2f} P={r['precision']:.2f}" for r in mrows))
    man = pd.DataFrame(man); mat = pd.DataFrame(mat); pts = pd.concat(pts, ignore_index=True) if pts else pd.DataFrame(); areas = pd.DataFrame(areas)
    man.to_csv(RES / "R2b_manual_windows.csv", index=False); mat.to_csv(RES / "R2b_matching.csv", index=False)
    pts.to_csv(RES / "R2b_points.csv", index=False); areas.to_csv(RES / "R2b_areas.csv", index=False)
    ir = interrater(pts, wins) if len(pts) else pd.DataFrame(); ir.to_csv(RES / "R2b_interrater.csv", index=False)
    # confronto dei conteggi di ogni annotatore (umano o modello) con il primario, per finestra
    prim = man[man.rater == primary].set_index("win_id").n_manual
    cmp_rows = []
    for _, r in man[man.rater != primary].iterrows():
        if r.win_id in prim.index:
            cmp_rows.append(dict(win_id=r.win_id, archetype=r.archetype, rater=r.rater, n_rater=r.n_manual, n_primary=int(prim[r.win_id]),
                                 rel_err=(r.n_manual - prim[r.win_id]) / prim[r.win_id] if prim[r.win_id] else np.nan, count_only=r.count_only))
    pd.DataFrame(cmp_rows).to_csv(RES / "R2b_rater_vs_primary.csv", index=False)
    # sommario per archetipo (annotatore primario): densita' pooled su area netta, IC Poisson, media/SD fra finestre, confronto R2
    cons = pd.read_csv(ROOT / "results/R2/R2_consensus_density.csv")
    rows = []
    for arch, g in man[man.rater == primary].groupby("archetype"):
        n, area = g.n_manual.sum(), g.area_net_mm2.sum(); lo, hi = poisson_ci(n, area)
        r = dict(archetype=arch, rater=primary, n_windows=len(g), n_manual=n, area_net_mm2=area, density_manual=n / area, ci_lo=lo, ci_hi=hi,
                 density_mean_windows=g.density_manual_net.mean(), density_sd_windows=g.density_manual_net.std(ddof=1), density_manual_gross=g.n_points_in_window.sum() / g.area_gross_mm2.sum(),
                 consensus_R2_5roi=cons[cons.archetype == arch].density_consensus.mean(), consensus_R2_sampled_roi=cons[(cons.archetype == arch) & cons.roi_id.isin(g.roi_id)].density_consensus.mean(),
                 cellpose_R2_5roi=cons[cons.archetype == arch].density_cellpose.mean(), stardist_R2_5roi=cons[cons.archetype == arch].density_stardist.mean())
        mm = mat[(mat.archetype == arch) & (mat.rater == primary)]
        for meth, gm in mm.groupby("method"):
            r[f"n_obj_{meth}"] = gm.n_obj.sum(); r[f"density_{meth}_windows"] = gm.n_obj.sum() / area
            r[f"recall_{meth}"] = gm.tp_points.sum() / gm.n_manual.sum(); r[f"precision_{meth}"] = gm.tp.sum() / gm.n_obj.sum() if gm.n_obj.sum() else np.nan
            r[f"ghost_{meth}"] = 1 - r[f"precision_{meth}"]
        pp = pts[(pts.rater == primary) & pts.win_id.str.startswith(arch + "_")]
        cols = [c for c in ("match_cellpose_rgb", "match_stardist_he") if c in pp]
        r["common_fn_frac_cp_sd"] = float((~pp[cols].any(axis=1)).mean()) if len(pp) and cols else np.nan
        cols3 = [c for c in pp.columns if c.startswith("match_")]
        r["common_fn_frac_all"] = float((~pp[cols3].any(axis=1)).mean()) if len(pp) and cols3 else np.nan
        aa = areas[(areas.archetype == arch) & (areas.rater == primary)]
        r["n_polygons"] = len(aa); r["area_manual_median_um2"] = aa.area_manual_um2.median() if len(aa) else np.nan
        for meth in METHODS:
            c = f"area_{meth}_um2"
            if c in aa and aa[c].notna().any():
                q = aa[aa[c].notna()]; r[f"area_{meth}_median_paired"] = q[c].median(); r[f"area_manual_median_paired_{meth}"] = q.area_manual_um2.median()
                r[f"area_ratio_manual_over_{meth}"] = (q.area_manual_um2 / q[c]).median(); r[f"n_paired_{meth}"] = len(q)
        rows.append(r)
    summ = pd.DataFrame(rows); summ.to_csv(RES / "R2b_archetype_summary.csv", index=False)
    if summ.empty: print(f"nessuna finestra dell'annotatore primario '{primary}': sommario per archetipo vuoto"); return
    print(summ[["archetype", "n_manual", "density_manual", "ci_lo", "ci_hi", "consensus_R2_5roi"] + [c for c in summ if c.startswith("recall_") or c.startswith("precision_")]].round(3).to_string())

# ----------------------------------------------------------------------------------------------------- autotest sintetico
def selftest():
    rng = np.random.default_rng(7); upp = 0.274; S = 1200
    def make_labels(n, rmin, rmax):
        lab = np.zeros((S, S), np.int32); cents = []
        while len(cents) < n:
            x, y, r = rng.uniform(30, S - 30), rng.uniform(30, S - 30), rng.uniform(rmin, rmax)
            if all(np.hypot(x - cx, y - cy) > r + cr + 8 for cx, cy, cr in cents):
                cents.append((x, y, r)); rr, cc = draw.disk((y, x), r, shape=lab.shape); lab[rr, cc] = len(cents)
        return lab, np.array(cents)
    results = []
    def rec(name, ok, detail=""): results.append(dict(check=name, result="PASS" if ok else "FAIL", detail=detail)); print(("PASS" if ok else "FAIL"), name, detail)
    # Scenario A: punti dentro la maschera (jitter 1.5 px), 80 % degli oggetti + 30 falsi positivi lontani; esclusione in un angolo
    lab, cents = make_labels(300, 6, 12); objs = objects_table(lab)
    side = 1000; c0 = r0 = 100
    det = rng.random(300) < 0.8
    pts = cents[det, :2] + rng.normal(0, 1.5, (det.sum(), 2))
    fp = []
    while len(fp) < 30:
        p = rng.uniform(c0 + 5, c0 + side - 5, 2)
        if np.hypot(cents[:, 0] - p[0], cents[:, 1] - p[1]).min() * upp > 6: fp.append(p)
    pts = np.vstack([pts, np.array(fp)])
    m = 37; excl = [{"type": "unreadable", "pts": [[m, m], [300 + m, m], [300 + m, 300 + m], [m, 300 + m]]}]   # coordinate immagine = finestra + margine
    wj = dict(win_id="T_w1", archetype="T", roi_id="r1", margin_px=m, side_px=side, um_per_px=upp, points=(pts - [c0, r0] + m).tolist(), polygons=[], exclusions=excl, done=True, note="")
    wrow = pd.Series(dict(c0_roi=c0, r0_roi=r0, x0_um=c0 * upp, y0_um=r0 * upp))
    row, mrows, prow, _ = analyse_window(wj, wrow, "test", {"synthetic": (lab, objs)}, {})
    in_win = (cents[:, 0] >= c0) & (cents[:, 0] < c0 + side) & (cents[:, 1] >= r0) & (cents[:, 1] < r0 + side)
    in_ex = (cents[:, 0] < c0 + 300) & (cents[:, 1] < r0 + 300)
    exp_obj = int((in_win & ~in_ex).sum())
    ptsw = pts - [c0, r0]; pin = (ptsw[:, 0] >= 0) & (ptsw[:, 0] < side) & (ptsw[:, 1] >= 0) & (ptsw[:, 1] < side); pex = (ptsw[:, 0] < 300) & (ptsw[:, 1] < 300)
    exp_manual = int((pin & ~pex).sum()); exp_tp = int((det & in_win & ~in_ex).sum())   # punti dei rilevati dentro finestra, fuori esclusione
    mr = mrows[0]
    rec("C-R2b.3a n_manual = punti dentro finestra e fuori esclusione", row["n_manual"] == exp_manual, f"{row['n_manual']} vs {exp_manual}")
    rec("C-R2b.3b n_obj = oggetti keep con centroide dentro finestra e fuori esclusione", mr["n_obj"] == exp_obj, f"{mr['n_obj']} vs {exp_obj}")
    rec("C-R2b.3c recall recuperato entro 1 %", abs(mr["recall"] - exp_tp / exp_manual) <= 0.01, f"{mr['recall']:.4f} vs {exp_tp/exp_manual:.4f}")
    rec("C-R2b.3d precisione recuperata entro 1 %", abs(mr["precision"] - exp_tp / exp_obj) <= 0.01, f"{mr['precision']:.4f} vs {exp_tp/exp_obj:.4f}")
    rec("C-R2b.3e area netta = finestra - esclusione", abs(row["area_net_mm2"] - ((side ** 2 - 300 ** 2) * upp ** 2 / 1e6)) < 1e-9, f"{row['area_net_mm2']:.6f}")
    rec("C-R2b.3f appaiamenti tutti dentro maschera (costo 0)", mr["frac_inside_mask"] > 0.97, f"{mr['frac_inside_mask']:.3f}")
    # Scenario B: oggetti piccoli (r 2-3 px), punti a 1.5-2.5 um dal centroide FUORI maschera -> ramo distanza; un oggetto con due punti
    lab2, cents2 = make_labels(200, 2, 3); objs2 = objects_table(lab2)
    ang = rng.uniform(0, 2 * np.pi, 200); dist = rng.uniform(1.5, 2.5, 200) / upp
    pts2 = cents2[:, :2] + np.c_[np.cos(ang), np.sin(ang)] * dist[:, None]
    pts2 = np.vstack([pts2, pts2[:1] + [2, 2]])
    pts3_ang = ang                                     # doppione: solo uno dei due si appaia
    wj2 = dict(win_id="T_w2", archetype="T", roi_id="r2", margin_px=0, side_px=S, um_per_px=upp, points=pts2.tolist(), polygons=[], exclusions=[], done=True, note="")
    row2, mrows2, _, _ = analyse_window(wj2, pd.Series(dict(c0_roi=0, r0_roi=0, x0_um=0, y0_um=0)), "test", {"synthetic": (lab2, objs2)}, {})
    mr2 = mrows2[0]
    rec("C-R2b.3g ramo distanza (<= 3 um): recall = 200/201", abs(mr2["recall"] - 200 / 201) < 1e-9, f"{mr2['recall']:.4f}")
    rec("C-R2b.3h assegnazione 1:1: precisione = 1 con un punto doppio", abs(mr2["precision"] - 1.0) < 1e-9 and mr2["tp"] == 200, f"P={mr2['precision']:.4f} tp={mr2['tp']}")
    # Scenario C: punti oltre 3 um da tutto -> nessun appaiamento
    pts3 = cents2[:, :2] + np.c_[np.cos(ang), np.sin(ang)] * (3.5 / upp)
    far = np.array([np.hypot(cents2[:, 0] - p[0], cents2[:, 1] - p[1]).min() * upp > 3.0 for p in pts3])
    wj3 = dict(wj2, win_id="T_w3", points=pts3[far].tolist())
    _, mrows3, _, _ = analyse_window(wj3, pd.Series(dict(c0_roi=0, r0_roi=0, x0_um=0, y0_um=0)), "test", {"synthetic": (lab2, objs2)}, {})
    # alcuni punti a 3.5 um possono cadere dentro la maschera di un oggetto vicino? no: raggio <= 3 px = 0.8 um e distanza minima fra dischi > 8 px
    rec("C-R2b.3i punti a > 3 um da ogni centroide e fuori maschera: recall = 0", mrows3[0]["tp_points"] == 0, f"tp={mrows3[0]['tp_points']} su {int(far.sum())}")
    # Scenario D: poligono manuale -> area shoelace vs area disco
    pg = [[c0 + m + 500 + 20 * np.cos(t), r0 + m + 500 + 20 * np.sin(t)] for t in np.linspace(0, 2 * np.pi, 200, endpoint=False)]
    wj4 = dict(wj, win_id="T_w4", points=[], polygons=[pg], exclusions=[])
    _, _, _, arows = analyse_window(wj4, wrow, "test", {"synthetic": (lab, objs)}, {})
    rec("C-R2b.3j area del poligono (cerchio r=20 px) entro 0.1 % di pi*r^2", abs(arows[0]["area_manual_um2"] - np.pi * 400 * upp ** 2) / (np.pi * 400 * upp ** 2) < 1e-3, f"{arows[0]['area_manual_um2']:.3f} um2")
    df = pd.DataFrame(results); df.to_csv(RES / "R2b_selftest.csv", index=False)
    print("ALL PASS" if (df.result == "PASS").all() else f"{(df.result=='FAIL').sum()} FAIL")
    return (df.result == "PASS").all()

if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    prim = sys.argv[sys.argv.index("--primary") + 1] if "--primary" in sys.argv else "Luca"
    main(prim)
