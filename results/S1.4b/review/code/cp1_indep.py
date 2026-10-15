#!/usr/bin/env python
"""review/code/cp1_indep.py — revisore "code": CP-1 ricalcolato senza il sollevamento 3D.
Pixel "fuori" se esiste j con pow_j(p) < pow_own(p) (pow = |p-c|^2 - w). Candidati: K vicini euclidei; esatto quando
d_K^2 >= pow_own + W (nessun altro j puo' battere il proprio); altrimenti forza bruta sul pixel. Stesse permutazioni del
banco (default_rng(20261016 + indice di riga del ROI), 20 permutation() in sequenza). Confronto con S1.4b_cp1_raw.csv."""
import os, sys, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from multiprocessing import Pool
R3 = "/mnt/micron/geo_spatialtrans/R3"; MASKS = "/mnt/micron/geo_spatialtrans/R2/masks"; S14B = "/mnt/micron/geo_spatialtrans/S1.4b"
SUF = {"cellpose_rgb": "cellpose_rgb_native", "spaceranger": "spaceranger_native"}
PRIM = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "spaceranger", "A5": "cellpose_rgb", "A6": "spaceranger"}
K = 24
def out_mask(px, py, own, gx, gy, w, tree):
    W = w.max(); po = (px - gx[own])**2 + (py - gy[own])**2 - w[own]
    dk, ik = tree.query(np.column_stack([px, py]), k=K, workers=4)
    pw = dk**2 - w[ik]; better = (pw < po[:, None]) & (ik != own[:, None])
    out = better.any(1); cert = dk[:, -1]**2 >= po + W
    nb = np.where(~cert & ~out)[0]
    for t in nb:                                   # forza bruta sui pixel non certificati
        pp = (px[t] - gx)**2 + (py[t] - gy)**2 - w; pp[own[t]] = np.inf; out[t] = pp.min() < po[t]
    return out, len(nb)
def roi(args):
    i, A, r = args; m = PRIM[A]
    gen = pd.read_parquet(f"{R3}/real/{A}_{r}_{m}_gen.parquet"); it = pd.read_parquet(f"{S14B}/cp1_in/{A}_{r}_{m}.parquet")
    lab = np.load(f"{MASKS}/{A}_{r}_{SUF[m]}.npz")["labels"]; upp = float(gen.um_per_px.iloc[0])
    # mappa etichetta -> generatore costruita in modo diverso (dizionario) per controllare lab2gen
    d = {int(l): k for k, l in enumerate(gen.label.to_numpy())}; assert len(d) == len(gen), "etichette duplicate"
    rr, cc = np.nonzero(lab); L = lab[rr, cc]; own = np.array([d.get(int(l), -1) for l in np.unique(L)])
    u, inv = np.unique(L, return_inverse=True); own = own[inv]; k = own >= 0; rr, cc, own = rr[k], cc[k], own[k]
    px = cc * upp; py = rr * upp; gx = gen.x.to_numpy(); gy = gen.y.to_numpy(); a = gen.area_nuc.to_numpy()
    inter = it.interior.to_numpy().astype(bool); tree = cKDTree(np.column_stack([gx, gy]))
    rng = np.random.default_rng(20261016 + i)
    conds = [("vor", 0.0, None), ("lag", 1.0, None)] + [(f"perm{p + 1:02d}", 1.0, rng.permutation(len(a))) for p in range(20)]
    G = len(gx); tot = np.bincount(own, minlength=G); rows = []; nbf = 0
    for nm, lam, perm in conds:
        w = lam * (a[perm] if perm is not None else a) / np.pi
        o, nb = out_mask(px, py, own, gx, gy, w, tree); nbf += nb
        bad = np.bincount(own[o], minlength=G); sel = inter & (tot > 0); f = np.where(tot > 0, bad / np.maximum(tot, 1), np.nan)
        rows.append(dict(archetype=A, roi=r, cond=nm, m1=bad[sel].sum() / tot[sel].sum(), m2=float(np.mean(f[sel] > 0.05)), n_nuc=int(sel.sum())))
    # controlli sul dato: pixel con etichetta non tra i generatori; generatori senza pixel
    rows.append(dict(archetype=A, roi=r, cond="_diag", m1=float(np.mean(~k)), m2=float(np.mean(tot == 0)), n_nuc=nbf))
    return rows
if __name__ == "__main__":
    rois = pd.read_csv("results/R2/R2_rois_checked.csv"); jobs = [(i, rw.archetype, rw.roi_id) for i, rw in rois.iterrows()]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 15) as p: R = [x for rr in p.map(roi, jobs) for x in rr]
    me = pd.DataFrame(R); me.to_csv("results/S1.4b/review/code/cp1_indep_raw.csv", index=False)
    ref = pd.read_csv("results/S1.4b/S1.4b_cp1_raw.csv")
    j = me.merge(ref[["archetype", "roi", "cond", "m1", "m2", "n_nuc"]], on=["archetype", "roi", "cond"], suffixes=("", "_ref"))
    print("conds compared", len(j), "max|dm1|", np.abs(j.m1 - j.m1_ref).max(), "max|dm2|", np.abs(j.m2 - j.m2_ref).max(), "n_nuc equal", (j.n_nuc == j.n_nuc_ref).all())
    print(me[me.cond == "_diag"].describe().to_string())
    b = []
    for (A, r), z in me[me.cond != "_diag"].groupby(["archetype", "roi"]):
        pm = z[z.cond.str.startswith("perm")]; g = z.set_index("cond")
        b.append(dict(A=A, roi=r, beats=bool(g.loc["lag", "m1"] < pm.m1.min() and g.loc["lag", "m2"] < pm.m2.min()), lag_beats_vor=bool(g.loc["lag", "m1"] < g.loc["vor", "m1"] and g.loc["lag", "m2"] < g.loc["vor", "m2"]),
                      vor_beats_permmin=bool(g.loc["vor", "m1"] < pm.m1.min() and g.loc["vor", "m2"] < pm.m2.min()), m2_red=1 - g.loc["lag", "m2"] / g.loc["vor", "m2"]))
    b = pd.DataFrame(b); print(b.groupby("A").agg(beats=("beats", "sum"), lag_beats_vor=("lag_beats_vor", "sum"), vor_beats_perm=("vor_beats_permmin", "sum")).to_string())
    b.to_csv("results/S1.4b/review/code/cp1_indep_roi.csv", index=False)
