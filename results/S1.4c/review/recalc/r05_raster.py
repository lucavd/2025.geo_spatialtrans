#!/usr/bin/env python
"""results/S1.4c/review/recalc/r05_raster.py — revisione avversariale S1.4c (recalc).
Rasterizzazione indipendente dei poligoni simulati (gen_in) con shapely.contains_xy sui centri dei pixel
(x = c*upp, y = r*upp), invece di skimage.draw.polygon; confronto di n_px, centroide, conflitti per nucleo con sim_out;
poi n_bnd/n_c2 con il ciclo per etichetta di r04 sul raster proprio e S-1 sull'insieme ev.
Uso: .venv/bin/python r05_raster.py <file gen_in> ...
"""
import sys, os, numpy as np, pandas as pd, shapely
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r04_masks import per_label
B = "/mnt/micron/geo_spatialtrans/S1.4c/"; D = "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/dump/"
rois = pd.read_csv("results/R2/R2_rois_checked.csv")
out = []
for name in sys.argv[1:]:
    A, roi = name.split("_")[:2]; rw = rois[(rois.archetype == A) & (rois.roi_id == roi)].iloc[0]
    upp = float(rw.um_per_px); side = int(rw.side_px)
    t = pd.read_parquet(B + "gen_in/" + name); G = shapely.from_wkb(t.wkb.to_numpy()); ids = t.idx.to_numpy()
    lab = np.zeros((side, side), np.int64); cnt = np.zeros((side, side), np.uint8)
    for g, i in zip(G, ids):
        if g is None or shapely.is_empty(g): continue
        x0, y0, x1, y1 = g.bounds
        c0 = max(int(np.ceil(x0 / upp)) - 1, 0); c1 = min(int(np.floor(x1 / upp)) + 2, side)
        r0 = max(int(np.ceil(y0 / upp)) - 1, 0); r1 = min(int(np.floor(y1 / upp)) + 2, side)
        if c1 <= c0 or r1 <= r0: continue
        cc, rr = np.meshgrid(np.arange(c0, c1), np.arange(r0, r1))
        m = shapely.contains_xy(g, cc * upp, rr * upp)
        sub = lab[r0:r1, c0:c1]; sub[m] = i; cnt[r0:r1, c0:c1] += m
    df, _ = per_label(lab, 2); df = df.set_index("label")
    df["cx"] = df.sc / df.n_px * upp; df["cy"] = df.sr / df.n_px * upp
    so = pd.read_parquet(B + "sim_out/" + name).set_index("label")
    allid = np.union1d(df.index, so.index)
    a = df.n_px.reindex(allid, fill_value=0); b = so.n_px.reindex(allid, fill_value=0)
    dif = (a - b).to_numpy()
    cc_ = pd.read_parquet(D + f"cells_{A}_{roi}.parquet"); ev = cc_.idx.to_numpy()[cc_.ev_prep.to_numpy().astype(bool)]
    ev_m = np.intersect1d(ev, df.index); ev_s = np.intersect1d(ev, so.index)
    com = np.intersect1d(df.index, so.index)
    out.append(dict(file=name, n_poly=len(G), n_lab_mine=len(df), n_lab_simout=len(so), n_px_tot_mine=int(a.sum()), n_px_tot_simout=int(b.sum()),
                    n_nuc_npx_diff=int((dif != 0).sum()), max_abs_npx_diff=int(np.abs(dif).max()), conflicts_mine=int((cnt > 1).sum()),
                    conflicts_simout=int(so.n_conflict.iloc[0]), max_dcentroid_um=float(np.max(np.hypot(df.cx[com] - so.cx[com], df.cy[com] - so.cy[com]))),
                    n_bnd_diff=int((df.n_bnd[com] != so.n_bnd[com]).sum()), n_c2_diff=int((df.n_c[com] != so.n_c2[com]).sum()),
                    s1_ev_mine=float((df.n_c[ev_m] > 0).mean()), s1_ev_simout=float((so.n_c2[ev_s] > 0).mean())))
    print(out[-1], flush=True)
pd.DataFrame(out).to_csv("/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/r05_raster.csv", index=False)
