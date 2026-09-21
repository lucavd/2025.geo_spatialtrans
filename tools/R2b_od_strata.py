"""tools/R2b_od_strata.py — R2b: CP-R2b.5, evidenza di ematossilina dei punti manuali appaiati vs persi (stesso arbitro di R2).
Per ogni punto dell'annotatore primario: OD ematossilina (Ruifrok, canale H di rgb2hed) media nel disco r <= 1.5 um meno media nell'anello
3-5 um; classe = 'matched' se appaiato ad almeno uno fra cellpose_rgb / stardist_he (colonne match_* di R2b_points.csv), 'missed_both'
altrimenti; 'matched_any3' include anche spaceranger. Controllo: OD degli stessi dischi centrati su punti casuali nel tessuto (fondo).
Output: results/R2b/R2b_od_points.csv (per punto), results/R2b/R2b_od_strata.csv (per archetipo), results/R2b/figures/od_strata.png
Uso: .venv/bin/python tools/R2b_od_strata.py [--primary Luca]
"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from skimage import color
sys.path.insert(0, str(Path(__file__).resolve().parent))
from R2_common import ROOT, load_roi, masks_for
RES = ROOT / "results/R2b"; R_IN, R_OUT1, R_OUT2 = 1.5, 3.0, 5.0

def od_feature(h, x, y, upp):
    """h: canale ematossilina (float, OD); disco e anello attorno a (x, y) in px."""
    r2 = int(np.ceil(R_OUT2 / upp)); H, W = h.shape
    x0, x1, y0, y1 = max(0, int(x) - r2), min(W, int(x) + r2 + 1), max(0, int(y) - r2), min(H, int(y) + r2 + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]; d = np.hypot(xx - x, yy - y) * upp; sub = h[y0:y1, x0:x1]
    inner = sub[d <= R_IN]; ring = sub[(d >= R_OUT1) & (d <= R_OUT2)]
    return float(inner.mean()) if inner.size else np.nan, float(ring.mean()) if ring.size else np.nan

def main(primary="Luca"):
    wins = pd.read_csv(RES / "R2b_windows.csv").set_index("win_id")
    pts = pd.read_csv(RES / "R2b_points.csv"); pts = pts[pts.rater == primary].copy()
    rng = np.random.default_rng(20260921); rows = []
    for win_id, g in pts.groupby("win_id"):
        w = wins.loc[win_id]; upp = float(w.um_per_px); side = int(w.side_px); c0, r0 = int(w.c0_roi), int(w.r0_roi)
        rgb = load_roi(w.archetype, w.roi_id); img = rgb[r0:r0 + side, c0:c0 + side]
        h = color.rgb2hed(img)[..., 0]
        t, b = masks_for(rgb, upp, 1, w.archetype, w.roi_id); valid = t[r0:r0 + side, c0:c0 + side] & ~b[r0:r0 + side, c0:c0 + side]
        for _, p in g.iterrows():
            oi, orr = od_feature(h, p.x_win_px, p.y_win_px, upp)
            m_cp, m_sd = bool(p.get("match_cellpose_rgb", False)), bool(p.get("match_stardist_he", False)); m_sr = bool(p.get("match_spaceranger", False)) if "match_spaceranger" in p and pd.notna(p.get("match_spaceranger")) else False
            rows.append(dict(win_id=win_id, archetype=w.archetype, kind="manual", cls="matched" if (m_cp or m_sd) else "missed_both",
                             matched_any3=(m_cp or m_sd or m_sr), od_in=oi, od_ring=orr, od_contrast=oi - orr, x=p.x_win_px, y=p.y_win_px))
        vy, vx = np.nonzero(valid); idx = rng.integers(0, len(vy), len(g))
        for xr, yr in zip(vx[idx] + rng.random(len(g)), vy[idx] + rng.random(len(g))):
            oi, orr = od_feature(h, xr, yr, upp)
            rows.append(dict(win_id=win_id, archetype=w.archetype, kind="random", cls="random", matched_any3=False, od_in=oi, od_ring=orr, od_contrast=oi - orr, x=xr, y=yr))
        print(win_id, len(g))
    df = pd.DataFrame(rows); df.to_csv(RES / "R2b_od_points.csv", index=False)
    out = []
    for arch, g in df.groupby("archetype"):
        m = g[g.cls == "matched"].od_contrast; x = g[g.cls == "missed_both"].od_contrast; r = g[g.cls == "random"].od_contrast
        out.append(dict(archetype=arch, n_matched=len(m), n_missed=len(x), n_random=len(r), od_matched_median=m.median(), od_missed_median=x.median(), od_random_median=r.median(),
                        od_matched_q25=m.quantile(.25), od_missed_q75=x.quantile(.75), ratio_missed_over_matched=x.median() / m.median() if m.median() else np.nan,
                        frac_missed_above_matched_q25=float((x > m.quantile(.25)).mean()) if len(x) else np.nan,
                        frac_missed_above_random_q75=float((x > r.quantile(.75)).mean()) if len(x) else np.nan))
    s = pd.DataFrame(out); s["missed_paler_25pct"] = s.ratio_missed_over_matched <= 0.75; s.to_csv(RES / "R2b_od_strata.csv", index=False)
    print(s.round(3).to_string(index=False))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), facecolor="white", sharey=False)
    for ax, (arch, g) in zip(axes.ravel(), df.groupby("archetype")):
        data = [g[g.cls == c].od_contrast.dropna() for c in ("matched", "missed_both", "random")]
        ax.boxplot(data, tick_labels=[f"appaiati\n(n={len(data[0])})", f"persi CP&SD\n(n={len(data[1])})", f"casuali\n(n={len(data[2])})"], showfliers=False, widths=0.6)
        ax.set_title(arch); ax.set_ylabel("OD ematossilina: disco 1.5 µm − anello 3–5 µm"); ax.axhline(0, color="gray", lw=0.5)
    fig.suptitle(f"CP-R2b.5 — evidenza di ematossilina dei punti di {primary}: appaiati da ≥1 segmentatore vs persi da entrambi vs punti casuali nel tessuto", fontsize=11)
    fig.tight_layout(); (RES / "figures").mkdir(exist_ok=True); fig.savefig(RES / "figures/od_strata.png", dpi=110)

if __name__ == "__main__":
    main(sys.argv[sys.argv.index("--primary") + 1] if "--primary" in sys.argv else "Luca")
