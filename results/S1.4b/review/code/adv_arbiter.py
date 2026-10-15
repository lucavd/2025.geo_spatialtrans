#!/usr/bin/env python
"""review/code/adv_arbiter.py — arbitri indipendenti per i casi avversari di power_diagram() (revisore "code").
(1) BF: cella_i = intersezione di TUTTI i semipiani pow_i <= pow_j (j != i), Sutherland-Hodgman numpy proprio, poi
    shapely intersection con la regione; confronto per cella con le celle di R (area, differenza simmetrica, vuote).
(2) RASTER: centri di pixel di passo h; argmin della potenza per forza bruta (nessun KD-tree, nessun sollevamento);
    pixel assegnato male se la cella di R che lo contiene != argmin e la distanza dal bordo di potenza > TOL.
Uso: python adv_arbiter.py <dir> <out.csv> [h]"""
import sys, glob, os
import numpy as np, pandas as pd, shapely

def clip(P, a, b):
    if len(P) == 0: return P
    s = P @ a - b; ins = s <= 0
    if ins.all(): return P
    if not ins.any(): return P[:0]
    out = []; m = len(P)
    for i in range(m):
        j = (i + 1) % m
        if ins[i]: out.append(P[i])
        if ins[i] != ins[j]:
            t = s[i] / (s[i] - s[j]); out.append(P[i] + t * (P[j] - P[i]))
    return np.array(out)

def bf_cells(x, y, w, box):
    mx, my = x.mean(), y.mean(); xc, yc = x - mx, y - my; q = xc**2 + yc**2 - w
    B = np.array([[box[0], box[1]], [box[2], box[1]], [box[2], box[3]], [box[0], box[3]]]) - [mx, my]
    out = []
    for i in range(len(x)):
        P = B.copy()
        # ordine: prima i vicini euclidei (solo efficienza; si usano comunque TUTTI i j)
        for j in np.argsort((xc - xc[i])**2 + (yc - yc[i])**2):
            if j == i: continue
            P = clip(P, 2 * np.array([xc[j] - xc[i], yc[j] - yc[i]]), q[j] - q[i])
            if len(P) == 0: break
        out.append(shapely.Polygon(P + [mx, my]) if len(P) >= 3 else shapely.Polygon())
    return out

def run(f, h):
    d = pd.read_parquet(f); reg = shapely.from_wkb(open(f.replace(".parquet", "_region.wkb"), "rb").read())
    cells = shapely.from_wkb(d["wkb"].to_numpy()); x, y, w = (d[k].to_numpy() for k in ("x", "y", "w"))
    n = len(x); bb = reg.bounds; pad = 50 + 0.5 * max(bb[2] - bb[0], bb[3] - bb[1])
    box = (min(bb[0], x.min()) - pad, min(bb[1], y.min()) - pad, max(bb[2], x.max()) + pad, max(bb[3], y.max()) + pad)
    bf = shapely.intersection(np.array(bf_cells(x, y, w, box), dtype=object), reg)
    aR = shapely.area(cells); aB = shapely.area(bf)
    sym = shapely.area(shapely.symmetric_difference(cells, bf))
    emR = aR <= 0; emB = aB <= 0
    # coppie esattamente coincidenti con pesi uguali: partizione ambigua -> confronto sull'unione
    key = pd.Series(list(zip(x, y, w))); dup = key.duplicated(keep=False).to_numpy()
    nd = ~dup
    # raster
    xs = np.arange(bb[0] + h / 2, bb[2], h); ys = np.arange(bb[1] + h / 2, bb[3], h)
    X, Y = np.meshgrid(xs, ys); X = X.ravel(); Y = Y.ravel(); inr = shapely.contains_xy(reg, X, Y); X, Y = X[inr], Y[inr]
    mx, my = x.mean(), y.mean()
    best = np.empty(len(X), int); gapd = np.empty(len(X))
    for s0 in range(0, len(X), 20000):
        Xc = X[s0:s0 + 20000, None] - mx; Yc = Y[s0:s0 + 20000, None] - my
        Pw = (Xc - (x - mx))**2 + (Yc - (y - my))**2 - w
        o = np.argsort(Pw, axis=1)[:, :2]; b1, b2 = o[:, 0], o[:, 1]; r = np.arange(len(b1))
        dc = np.hypot(x[b1] - x[b2], y[b1] - y[b2])
        gapd[s0:s0 + 20000] = np.where(dc > 0, (Pw[r, b2] - Pw[r, b1]) / (2 * np.maximum(dc, 1e-300)), np.inf * (Pw[r, b2] > Pw[r, b1]))
        best[s0:s0 + 20000] = b1
    lab = np.full(len(X), -1); cnt = np.zeros(len(X), int)
    for i in range(n):
        if emR[i]: continue
        x0, y0, x1, y1 = cells[i].bounds; s = (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)
        if not s.any(): continue
        ii = np.where(s)[0]; hit = shapely.contains_xy(cells[i], X[ii], Y[ii]); lab[ii[hit]] = i; cnt[ii[hit]] += 1
    TOL = 1e-4
    far = gapd > TOL
    okdup = dup[best] & dup[np.maximum(lab, 0)] & (lab >= 0)   # pixel delle coppie coincidenti: qualunque dei due
    mis = far & (lab >= 0) & (lab != best) & ~okdup
    return dict(case=os.path.basename(f).replace(".parquet", ""), n=n, n_dup=int(dup.sum()),
                bf_emptyR_vs_emptyB=f"{int(emR[nd].sum())}/{int(emB[nd].sum())}", empty_mismatch=int((emR[nd] != emB[nd]).sum()),
                max_abs_dA=float(np.max(np.abs(aR - aB)[nd])), max_rel_dA=float(np.max((np.abs(aR - aB) / np.maximum(aB, 1e-300))[nd & ~emB])) if (nd & ~emB).any() else 0.0,
                max_symdiff=float(np.max(sym[nd])), n_symdiff_gt_1e6=int((sym[nd] > 1e-6).sum()),
                dup_union_dA=float(abs(aR[dup].sum() - aB[dup].sum()) / 2) if dup.any() else 0.0,
                px=len(X), h=h, px_mis=int(mis.sum()), px_uncov=int(((lab < 0) & far).sum()), px_double=int(((cnt >= 2) & far).sum()),
                px_mis_wrongcells=int(len(np.unique(best[mis]))))

if __name__ == "__main__":
    a = sys.argv[1:]; h = float(a[2]) if len(a) > 2 else 0.1
    rows = []
    for f in sorted(glob.glob(os.path.join(a[0], "*.parquet"))):
        try: rows.append(run(f, h))
        except Exception as e: rows.append(dict(case=os.path.basename(f), error=str(e)[:200]))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(a[1], index=False)
