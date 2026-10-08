
# results/R3/review/claims/rv_cp3_masks.py — CP-3: perche' Space Ranger 53% e StarDist 1% in A4?
# Contatto fra maschere (nuclei che toccano un'altra etichetta) ed erosione dei bordi (k pixel) prima del contenimento.
import numpy as np, pandas as pd, sys
from scipy.spatial import cKDTree
OUT = "/mnt/micron/geo_spatialtrans/R3"; MASKS = "/mnt/micron/geo_spatialtrans/R2/masks"; RV = "results/R3/review/claims"
SUF = dict(cellpose_rgb="cellpose_rgb_native", spaceranger="spaceranger_native", stardist_he="stardist_he_native")
cases = [("A4", r, m) for r in ["f1", "f2", "p1", "p2", "p3"] for m in ["spaceranger", "stardist_he"]] + \
        [("A6", r, m) for r in ["r1", "r3"] for m in ["spaceranger", "stardist_he"]] + [("A1", "r1", "cellpose_rgb"), ("A2", "r1", "cellpose_rgb")]
cells = pd.read_parquet(f"{OUT}/cells_all.parquet", columns=["archetype", "roi_id", "method", "idx", "interior", "frac_out", "area", "area_nuc"])
def erode_labels(L):
    E = L.copy()
    for sh in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        S = np.roll(L, sh, axis=(0, 1)); E[S != L] = 0
    return E
def touching(L):
    t = np.zeros(int(L.max()) + 1, bool)
    for sh in [(1, 0), (0, 1)]:
        S = np.roll(L, sh, axis=(0, 1)); m = (L > 0) & (S > 0) & (S != L)
        t[L[m]] = True; t[S[m]] = True
    return t
rows = []
for A, roi, me in cases:
    gen = pd.read_parquet(f"{OUT}/real/{A}_{roi}_{me}_gen.parquet"); upp = float(gen.um_per_px.iloc[0])
    L = np.load(f"{MASKS}/{A}_{roi}_{SUF[me]}.npz")["labels"].astype(np.int64)
    lab = gen.label.to_numpy().astype(np.int64); lab2gen = np.full(int(L.max()) + 1, -1, np.int64); lab2gen[lab] = np.arange(len(lab))
    tree = cKDTree(gen[["x", "y"]].to_numpy())
    cc = cells[(cells.archetype == A) & (cells.roi_id == roi) & (cells.method == me)].set_index("idx").loc[gen.idx]
    inter = cc.interior.to_numpy()
    tch = touching(L)[lab]
    E = L; res = {}
    for k in range(0, 5):
        if k > 0: E = erode_labels(E)
        r_, c_ = np.nonzero(E); Lk = E[r_, c_]; own = lab2gen[Lk]; ok = own >= 0
        _, near = tree.query(np.column_stack([c_[ok] * upp, r_[ok] * upp]), k=1, workers=4)
        tot = np.bincount(own[ok], minlength=len(lab)); bad = np.bincount(own[ok][near != own[ok]], minlength=len(lab))
        with np.errstate(invalid="ignore", divide="ignore"): f = bad / tot
        valid = inter & (tot > 0)
        res[k] = ((f[valid] > .05).mean(), valid.sum())
    f0 = cc.frac_out.to_numpy()
    rows.append(dict(archetype=A, roi_id=roi, method=me, upp=upp, n_gen=len(lab), n_interior=int(inter.sum()),
                     coverage_nuclear_px=(L > 0).mean(), frac_touching=tch.mean(), frac_touching_interior=tch[inter].mean(),
                     cut5_k0_check=(f0[inter] > .05).mean(),
                     cut5_touching=(f0[inter & tch] > .05).mean() if (inter & tch).any() else np.nan,
                     cut5_not_touching=(f0[inter & ~tch] > .05).mean() if (inter & ~tch).any() else np.nan,
                     **{f"cut5_erode{k}px": res[k][0] for k in range(1, 5)}, n_valid_erode4=res[4][1],
                     nuc_area_median=np.median(cc.area_nuc[inter]), terr_area_median=np.median(cc.area[inter])))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f"{RV}/rv_cp3_masks.csv", index=False)

