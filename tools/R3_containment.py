#!/usr/bin/env python
"""tools/R3_containment.py — R3: contenimento dei nuclei reali nel proprio territorio di Voronoi (CP-3),
momenti nucleari dai pixel (eccentricita' e orientazione nelle coordinate del ROI) e oggetto hed che contiene
il generatore (D-1). Pre-registrazione: results/R3/R3_preregistration.md (commit 3fac4e7).

Uso (root del repo):
  .venv/bin/python tools/R3_containment.py selftest [--mutant M3]     -> results/R3/R3_selftest_py[_M3].csv
  .venv/bin/python tools/R3_containment.py roi <A> <roi> <method>      -> <R3>/real/<A>_<roi>_<method>_nuc.parquet
  .venv/bin/python tools/R3_containment.py all                          -> tutti i *_gen.parquet presenti

Convenzioni: coordinate del pixel (r, c) = (c*upp, r*upp) µm, come i centroidi di R2 (R2_common.py: centroid * upp).
Generatore piu' vicino = cella di Voronoi (esatto). Orientazione: rad da +x verso +y (y in basso); regionprops
`orientation` si converte con theta = pi/2 - orientation (assiale), verificato in selftest (C-R3.3c).
"""
import sys, os, glob
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

R3 = "/mnt/micron/geo_spatialtrans/R3"
MASKS = "/mnt/micron/geo_spatialtrans/R2/masks"
SUFFIX = {"cellpose_rgb": "cellpose_rgb_native", "spaceranger": "spaceranger_native",
          "stardist_he": "stardist_he_native", "cellpose_hed": "cellpose_hed_native"}


def pixel_moments(labels, upp, labs):
    """Area, centroide, eccentricita', orientazione per le etichette `labs` (momenti dei centri dei pixel)."""
    r, c = np.nonzero(labels); L = labels[r, c]
    mx = int(labels.max()) + 1
    x = c * upp; y = r * upp
    n = np.bincount(L, minlength=mx).astype(float)
    sx = np.bincount(L, x, mx); sy = np.bincount(L, y, mx)
    sxx = np.bincount(L, x * x, mx); syy = np.bincount(L, y * y, mx); sxy = np.bincount(L, x * y, mx)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx = sx / n; cy = sy / n
        mxx = sxx / n - cx ** 2; myy = syy / n - cy ** 2; mxy = sxy / n - cx * cy
        h = np.sqrt(((mxx - myy) / 2) ** 2 + mxy ** 2); t2 = (mxx + myy) / 2
        l1 = t2 + h; l2 = np.maximum(t2 - h, 0)
        ecc = np.sqrt(np.maximum(1 - l2 / l1, 0)); th = 0.5 * np.arctan2(2 * mxy, mxx - myy)
    labs = np.asarray(labs)
    return pd.DataFrame({"label": labs, "n_px": n[labs], "area_px_um2": n[labs] * upp ** 2, "cx_px": cx[labs],
                         "cy_px": cy[labs], "ecc_N": ecc[labs], "theta_N": th[labs]})


def containment(labels, upp, gen_xy, lab2gen, swap=False):
    """Frazione dei pixel di ogni nucleo-generatore il cui generatore piu' vicino non e' il proprio."""
    r, c = np.nonzero(labels); L = labels[r, c]
    own = lab2gen[L]; k = own >= 0; r, c, own = r[k], c[k], own[k]
    pts = np.column_stack([r * upp, c * upp]) if swap else np.column_stack([c * upp, r * upp])
    _, near = cKDTree(gen_xy).query(pts, k=1, workers=8)
    G = len(gen_xy)
    tot = np.bincount(own, minlength=G); bad = np.bincount(own[near != own], minlength=G)
    with np.errstate(invalid="ignore", divide="ignore"):
        return bad / tot, tot


def selftest(mutant=None):
    rows = []
    # C-R3.3c: ellisse rasterizzata, momenti dei pixel vs regionprops
    from skimage import measure
    upp = 0.1; H = 400; yy, xx = np.mgrid[0:H, 0:H]; X = xx * upp; Y = yy * upp
    for phi_deg in (0, 30, 60, -30, 75):
        phi = np.deg2rad(phi_deg); a, b = 6.0, 3.0
        u = (X - 20) * np.cos(phi) + (Y - 20) * np.sin(phi); v = -(X - 20) * np.sin(phi) + (Y - 20) * np.cos(phi)
        Lm = ((u / a) ** 2 + (v / b) ** 2 <= 1).astype(np.uint32)
        pm = pixel_moments(Lm, upp, [1]).iloc[0]; rp = measure.regionprops(Lm)[0]
        th_rp = np.pi / 2 - rp.orientation
        d_th = abs((pm.theta_N - th_rp + np.pi / 2) % np.pi - np.pi / 2)
        d_an = abs((pm.theta_N - phi + np.pi / 2) % np.pi - np.pi / 2)
        rows.append(dict(check="C-R3.3c", case=f"ellipse_phi{phi_deg}", value=max(abs(pm.ecc_N - rp.eccentricity), 0),
                         value2=np.rad2deg(d_th), value3=np.rad2deg(d_an),
                         passed=bool(abs(pm.ecc_N - rp.eccentricity) < 0.02 and np.rad2deg(d_th) < 2 and np.rad2deg(d_an) < 2)))
    # C-R3.5: dischi tangenti di raggio diverso; dischi separati
    upp = 0.05; H = 800; yy, xx = np.mgrid[0:H, 0:H]; X = xx * upp; Y = yy * upp
    for case, c2 in (("tangent", 15.0), ("separate", 20.0)):
        Lm = np.zeros((H, H), np.uint32)
        Lm[(X - 10) ** 2 + (Y - 10) ** 2 <= 9] = 1; Lm[(X - c2) ** 2 + (Y - 10) ** 2 <= 4] = 2
        pm = pixel_moments(Lm, upp, [1, 2]); gen = pm[["cx_px", "cy_px"]].values
        lab2gen = np.full(3, -1); lab2gen[1] = 0; lab2gen[2] = 1
        f, _ = containment(Lm, upp, gen, lab2gen, swap=(mutant == "M3"))
        h = (gen[1, 0] - gen[0, 0]) / 2; R = 3.0
        expect1 = (R ** 2 * np.arccos(h / R) - h * np.sqrt(R ** 2 - h ** 2)) / (np.pi * R ** 2) if h < R else 0.0
        ok = abs(f[0] - expect1) < 0.01 and (f[1] == 0 if case == "tangent" else (f[0] == 0 and f[1] == 0))
        rows.append(dict(check="C-R3.5", case=case, value=f[0], value2=f[1], value3=expect1, passed=bool(ok)))
    out = pd.DataFrame(rows)
    os.makedirs("results/R3", exist_ok=True)
    out.to_csv(f"results/R3/R3_selftest_py{'_' + mutant if mutant else ''}.csv", index=False)
    print(out.to_string()); return out


def run_roi(A, roi, method, upp=None):
    gen = pd.read_parquet(f"{R3}/real/{A}_{roi}_{method}_gen.parquet")
    labels = np.load(f"{MASKS}/{A}_{roi}_{SUFFIX[method]}.npz")["labels"]
    upp = float(gen["um_per_px"].iloc[0])
    lab = gen["label"].to_numpy().astype(np.int64)
    lab2gen = np.full(int(labels.max()) + 1, -1, np.int64); lab2gen[lab] = np.arange(len(lab))
    pm = pixel_moments(labels, upp, lab)
    frac, tot = containment(labels, upp, gen[["x", "y"]].to_numpy(), lab2gen)
    hed = np.load(f"{MASKS}/{A}_{roi}_{SUFFIX['cellpose_hed']}.npz")["labels"]
    ri = np.clip(np.round(gen["y"].to_numpy() / upp).astype(int), 0, hed.shape[0] - 1)
    ci = np.clip(np.round(gen["x"].to_numpy() / upp).astype(int), 0, hed.shape[1] - 1)
    hl = hed[ri, ci]; ha = np.bincount(hed.ravel(), minlength=int(hed.max()) + 1) * upp ** 2
    out = pd.DataFrame({"idx": gen["idx"].to_numpy(), "label": lab}).merge(pm, on="label", how="left")
    out["frac_out"] = frac; out["n_px_own"] = tot
    out["hed_label"] = hl; out["hed_area_um2"] = np.where(hl > 0, ha[hl], np.nan)
    out.to_parquet(f"{R3}/real/{A}_{roi}_{method}_nuc.parquet", index=False)
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "selftest":
        selftest("M3" if "--mutant" in a and a[a.index("--mutant") + 1] == "M3" else None)
    elif a[0] == "roi":
        run_roi(a[1], a[2], a[3])
    elif a[0] == "all":
        for f in sorted(glob.glob(f"{R3}/real/*_gen.parquet")):
            A, roi, *m = os.path.basename(f).replace("_gen.parquet", "").split("_")
            run_roi(A, roi, "_".join(m)); print("done", os.path.basename(f), flush=True)
