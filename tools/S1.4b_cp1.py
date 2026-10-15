#!/usr/bin/env python
"""tools/S1.4b_cp1.py — S1.4b, CP-1 (controprova senza modello): i nuclei reali stanno meglio nelle celle di Laguerre
costruite con le loro aree che nel Voronoi dei centroidi, e meglio che con le stesse aree permutate?

Per ogni ROI (segmentatore primario di S1.4): pixel delle maschere native (R2), generatori = centroidi (R3 *_gen),
pesi w_i = lambda * A_own,i / pi (A_own = area del nucleo-generatore). Ogni pixel nucleare va al generatore di potenza
minima (KD-tree 3D sui generatori sollevati, esatto). frac_out = frazione dei pixel del nucleo fuori dalla propria cella.
Metriche sui nuclei interni (S1.4 `interior`): m1 = pixel fuori / pixel totali; m2 = frazione di nuclei con frac_out > 0.05.
Condizioni: vor (lambda = 0), lag (lambda = 1), perm01..perm20 (pesi permutati nel ROI), lam0.5, lam1.5, lam2.
Controllo: con lambda = 0 frac_out deve coincidere con R3 (`*_nuc.parquet`).
Uso: .venv/bin/python tools/S1.4b_cp1.py <A> <roi> <method> <outdir> | all <outdir>
"""
import sys, os, glob
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

R3 = "/mnt/micron/geo_spatialtrans/R3"; MASKS = "/mnt/micron/geo_spatialtrans/R2/masks"
S14B = os.environ.get("S14B_OUT", "/mnt/micron/geo_spatialtrans/S1.4b")
SUFFIX = {"cellpose_rgb": "cellpose_rgb_native", "spaceranger": "spaceranger_native", "stardist_he": "stardist_he_native"}
NPERM = 20; SEED = 20261016


def frac_out(r, c, own, upp, gx, gy, w):
    W = w.max(); tree = cKDTree(np.column_stack([gx, gy, np.sqrt(W - w)]))
    pts = np.column_stack([c * upp, r * upp, np.zeros(len(r))])
    _, near = tree.query(pts, k=1, workers=8)
    G = len(gx); tot = np.bincount(own, minlength=G); bad = np.bincount(own[near != own], minlength=G)
    with np.errstate(invalid="ignore", divide="ignore"):
        return bad / tot, tot, bad


def run_roi(A, roi, method, outdir, roi_index):
    gen = pd.read_parquet(f"{R3}/real/{A}_{roi}_{method}_gen.parquet")
    nuc = pd.read_parquet(f"{R3}/real/{A}_{roi}_{method}_nuc.parquet")
    it = pd.read_parquet(f"{S14B}/cp1_in/{A}_{roi}_{method}.parquet")
    assert (it["idx"].to_numpy() == np.arange(1, len(gen) + 1)).all() and len(it) == len(gen)
    labels = np.load(f"{MASKS}/{A}_{roi}_{SUFFIX[method]}.npz")["labels"]
    upp = float(gen["um_per_px"].iloc[0])
    lab = gen["label"].to_numpy().astype(np.int64)
    lab2gen = np.full(int(labels.max()) + 1, -1, np.int64); lab2gen[lab] = np.arange(len(lab))
    r, c = np.nonzero(labels); own = lab2gen[labels[r, c]]; k = own >= 0; r, c, own = r[k], c[k], own[k]
    gx = gen["x"].to_numpy(); gy = gen["y"].to_numpy(); a_own = gen["area_nuc"].to_numpy()
    inter = it["interior"].to_numpy().astype(bool)
    rng = np.random.default_rng(SEED + roi_index)
    conds = [("vor", 0.0, None), ("lag", 1.0, None)] + [(f"perm{p + 1:02d}", 1.0, rng.permutation(len(a_own))) for p in range(NPERM)] \
            + [("lam0.5", 0.5, None), ("lam1.5", 1.5, None), ("lam2", 2.0, None)]
    rows = []; f_vor = None
    for name, lam, perm in conds:
        w = lam * (a_own[perm] if perm is not None else a_own) / np.pi
        f, tot, bad = frac_out(r, c, own, upp, gx, gy, w)
        if name == "vor": f_vor = f
        sel = inter & (tot > 0)
        rows.append(dict(archetype=A, roi=roi, method=method, cond=name, lam=lam, n_nuc=int(sel.sum()),
                         m1=bad[sel].sum() / tot[sel].sum(), m2=float(np.mean(f[sel] > 0.05)),
                         m_any=float(np.mean(f[sel] > 0))))
    ref = nuc["frac_out"].to_numpy()
    both = np.isfinite(ref) & np.isfinite(f_vor)
    rows.append(dict(archetype=A, roi=roi, method=method, cond="check_R3", lam=0.0, n_nuc=int(both.sum()),
                     m1=float(np.max(np.abs(ref[both] - f_vor[both]))), m2=float(np.mean(ref[both] == f_vor[both])), m_any=np.nan))
    out = pd.DataFrame(rows); out.to_csv(os.path.join(outdir, f"{A}_{roi}_{method}.csv"), index=False)
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "all":
        outdir = a[1]; os.makedirs(outdir, exist_ok=True)
        rois = pd.read_csv("results/R2/R2_rois_checked.csv")
        prim = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "spaceranger", "A5": "cellpose_rgb", "A6": "spaceranger"}
        for i, rw in rois.iterrows():
            o = run_roi(rw["archetype"], rw["roi_id"], prim[rw["archetype"]], outdir, i)
            print(rw["archetype"], rw["roi_id"], "done", flush=True)
    else:
        os.makedirs(a[3], exist_ok=True); print(run_roi(a[0], a[1], a[2], a[3], int(a[4]) if len(a) > 4 else 0).to_string())
