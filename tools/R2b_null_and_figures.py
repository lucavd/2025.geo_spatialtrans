"""tools/R2b_null_and_figures.py — R2b: (1) CP-R2b.4 modello nullo dei punti casuali; (2) figure di sovrapposizione per finestra.

(1) Per ogni finestra e annotatore con posizioni: N punti uniformi nel tessuto valido e fuori dalle esclusioni (N = punti dell'annotatore),
    20 repliche (seed fisso), appaiati agli oggetti keep di ogni metodo con le stesse regole di R2b_match (punto dentro maschera,
    altrimenti centroide <= 3 um, 1:1). Output: results/R2b/R2b_null.csv (per finestra, metodo, annotatore: F1 osservato, F1 nullo medio/sd,
    precisione osservata vs nulla) e R2b_null_summary.csv per archetipo.
(2) results/R2b/figures/overlay_<win_id>.png: H&E + punti dell'annotatore primario (rossi: appaiati a Cellpose, gialli: persi da Cellpose e
    StarDist) + contorni Cellpose (verde) e StarDist (blu); pannello di sommario per archetipo.
Uso: .venv/bin/python tools/R2b_null_and_figures.py [--primary Luca] [--no-figures]
"""
import sys, json, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from skimage import segmentation
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT, load_roi, masks_for
from R2b_match import RES, ANN_DIR, METHODS, load_label_sets, window_objects, match_points, poly_mask

FIG = RES / "figures"; FIG.mkdir(exist_ok=True)
N_REP = 20; SEED = 20260921

def f1(tp, n_a, n_b): return 2 * tp / (n_a + n_b) if (n_a + n_b) else np.nan

def main(primary="Luca", figures=True):
    wins = pd.read_csv(RES / "R2b_windows.csv").set_index("win_id")
    anns = {}
    for f in sorted(glob.glob(str(ANN_DIR / "*.json"))):
        obj = json.load(open(f)); rater = obj.get("rater", Path(f).stem)
        for wj in obj["windows"]:
            if wj["n_points_in_window"] > 0: anns[(wj["win_id"], rater)] = wj
    cache = {}; rows = []; rng = np.random.default_rng(SEED)
    for (win_id, rater), wj in sorted(anns.items()):
        w = wins.loc[win_id]; m, side, upp = int(wj["margin_px"]), int(wj["side_px"]), float(wj["um_per_px"]); c0, r0 = int(w.c0_roi), int(w.r0_roi)
        rgb = load_roi(w.archetype, w.roi_id); t, b = masks_for(rgb, upp, 1, w.archetype, w.roi_id)
        valid = t[r0:r0 + side, c0:c0 + side] & ~b[r0:r0 + side, c0:c0 + side]
        emask = poly_mask([(np.array(e["pts"], float) - m).tolist() for e in wj.get("exclusions", [])], side); valid &= ~emask
        pts = np.array(wj["points"], float).reshape(-1, 2) - m
        inside = (pts[:, 0] >= 0) & (pts[:, 0] < side) & (pts[:, 1] >= 0) & (pts[:, 1] < side); pts = pts[inside]
        xi, yi = np.clip(np.floor(pts[:, 0]).astype(int), 0, side - 1), np.clip(np.floor(pts[:, 1]).astype(int), 0, side - 1)
        pts = pts[~emask[yi, xi]]; N = len(pts)
        vy, vx = np.nonzero(valid)
        sets, keeps = load_label_sets(w.archetype, w.roi_id, cache)
        for meth, (labels_roi, objs_full) in sets.items():
            crop = labels_roi[r0:r0 + side, c0:c0 + side]
            cand = window_objects(labels_roi, objs_full, c0, r0, side, upp, keeps.get(meth)); cand = cand[cand.keep].reset_index(drop=True)
            cx = np.clip(np.floor(cand.cx).astype(int), 0, side - 1); cy = np.clip(np.floor(cand.cy).astype(int), 0, side - 1)
            counted = cand[cand.in_win & ~emask[cy, cx]]; n_obj = len(counted)
            a_obs, _ = match_points(pts, cand, crop, upp)
            ml = cand.label.values[a_obs[a_obs >= 0]] if (a_obs >= 0).any() else np.array([]); tp_obs = int(counted.label.isin(ml).sum())
            nulls = []
            for k in range(N_REP):
                idx = rng.integers(0, len(vy), N); rp = np.c_[vx[idx] + rng.random(N), vy[idx] + rng.random(N)]
                a, _ = match_points(rp, cand, crop, upp); tp = int(counted.label.isin(cand.label.values[a[a >= 0]]).sum()) if (a >= 0).any() else 0
                nulls.append(tp)
            nulls = np.array(nulls)
            rows.append(dict(win_id=win_id, archetype=w.archetype, rater=rater, method=meth, n_points=N, n_obj=n_obj, tp_obs=tp_obs, f1_obs=f1(tp_obs, N, n_obj),
                             precision_obs=tp_obs / n_obj if n_obj else np.nan, tp_null_mean=nulls.mean(), f1_null_mean=np.mean([f1(x, N, n_obj) for x in nulls]),
                             f1_null_sd=np.std([f1(x, N, n_obj) for x in nulls]), precision_null_mean=nulls.mean() / n_obj if n_obj else np.nan,
                             f1_gain=f1(tp_obs, N, n_obj) - np.mean([f1(x, N, n_obj) for x in nulls])))
        if figures and rater == primary and "cellpose_rgb" in sets:
            pdf = pd.read_csv(RES / "R2b_points.csv"); pp = pdf[(pdf.win_id == win_id) & (pdf.rater == primary)]
            img = rgb[r0:r0 + side, c0:c0 + side].copy()
            fig, ax = plt.subplots(1, 1, figsize=(9, 9), facecolor="white")
            ax.imshow(img)
            for meth, col in (("cellpose_rgb", "#00ff00"), ("stardist_he", "#3399ff")):
                if meth in sets:
                    cr = sets[meth][0][r0:r0 + side, c0:c0 + side]
                    bd = segmentation.find_boundaries(cr, mode="inner"); yy, xx = np.nonzero(bd); ax.scatter(xx, yy, s=0.15, c=col, marker=".", linewidths=0)
            if len(pp):
                both = ~(pp.get("match_cellpose_rgb", False) | pp.get("match_stardist_he", False))
                ax.scatter(pp.x_win_px[~both], pp.y_win_px[~both], s=14, c="red", edgecolors="white", linewidths=0.4, label=f"{primary}: appaiato ({(~both).sum()})")
                ax.scatter(pp.x_win_px[both], pp.y_win_px[both], s=22, c="yellow", edgecolors="black", linewidths=0.5, marker="^", label=f"perso da CP e SD ({both.sum()})")
            for e in wj.get("exclusions", []):
                q = np.array(e["pts"], float) - m; ax.fill(q[:, 0], q[:, 1], color="orange", alpha=0.25, lw=1, edgecolor="orange")
            ax.set_xlim(0, side); ax.set_ylim(side, 0); ax.set_axis_off()
            ax.set_title(f"{win_id} ({w.roi_id}, {w.side_um:.0f} µm) — verde Cellpose, blu StarDist; n {primary} = {N}", fontsize=11)
            ax.legend(loc="lower left", fontsize=8, framealpha=0.85); fig.tight_layout(); fig.savefig(FIG / f"overlay_{win_id}.png", dpi=110); plt.close(fig)
        print(win_id, rater, N, " ".join(f"{r['method'][:2]}:F1 {r['f1_obs']:.2f} vs null {r['f1_null_mean']:.2f}" for r in rows[-len(sets):]))
    df = pd.DataFrame(rows); df.to_csv(RES / "R2b_null.csv", index=False)
    summ = df.groupby(["archetype", "rater", "method"]).apply(lambda g: pd.Series(dict(
        f1_obs=f1(g.tp_obs.sum(), g.n_points.sum(), g.n_obj.sum()), f1_null=f1(g.tp_null_mean.sum(), g.n_points.sum(), g.n_obj.sum()),
        precision_obs=g.tp_obs.sum() / g.n_obj.sum(), precision_null=g.tp_null_mean.sum() / g.n_obj.sum(), n_points=g.n_points.sum(), n_obj=g.n_obj.sum())), include_groups=False).reset_index()
    summ["f1_gain"] = summ.f1_obs - summ.f1_null; summ["CP-R2b.4"] = np.where(summ.f1_gain > 0.20, "PASS", "FAIL")
    summ.to_csv(RES / "R2b_null_summary.csv", index=False)
    print(summ[summ.method == "cellpose_rgb"].round(2).to_string(index=False))

if __name__ == "__main__":
    prim = sys.argv[sys.argv.index("--primary") + 1] if "--primary" in sys.argv else "Luca"
    main(prim, figures="--no-figures" not in sys.argv)
