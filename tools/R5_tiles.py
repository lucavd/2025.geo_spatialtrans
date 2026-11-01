"""tools/R5_tiles.py — R5, tessere di bordo: selezione (regola pre-registrata), maschere H&E sulla sottogriglia di 0.5 µm, fondo del vetro lontano.
Uso: .venv/bin/python tools/R5_tiles.py <KEY> <outdir>. Legge <outdir>/section/<KEY>.npz; scrive <outdir>/tiles/<KEY>/t<r>_<c>.npz,
<outdir>/tiles/<KEY>_candidates.csv, <outdir>/tiles/<KEY>_background.json.
Regola (docs/steps/R5.md, «Finestre di bordo»): tessere di 150 x 150 bin allineate alla griglia; in_tissue di Space Ranger fra 0.2 e 0.8; tutta dentro
l'immagine; a >= 100 µm dal bordo dell'area di cattura; tessuto H&E fra 0.2 e 0.8; bolle < 1 %; >= 2000 bin di vetro a 10-100 µm dal tessuto."""
import sys, os, json, time
import numpy as np, pandas as pd
from scipy import ndimage
sys.path.insert(0, "tools")
import R5_common as C
from R2_common import tissue_mask, bubble_mask, downsample
from R1_verify import read_region
from skimage import color, filters, morphology

D = "/mnt/micron/geo_spatialtrans/data_real"
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
TB, MB = 150, 75                      # tessera e margine (bin)
BORDER_BINS = 50                      # 100 µm dal bordo dell'area di cattura
FAR_UM, FAR_MIN_BINS = 300.0, 10000   # fondo del vetro lontano (D-2)

OD_THR = 0.15     # D-4: vetro all'H&E OD 0.075-0.125 nelle cinque sezioni (results/R5/R5_D4_od.txt), tessuto pallido sopra

def od_sum(rgb):
    return -np.log10((rgb.astype(float) + 1.0) / 256.0).sum(axis=2)

def tissue_mask_r5(rgb, upp, smooth_um=2.0, close_um=5.0):
    """D-4 (dopo E2, prima di ogni adattamento): maschera di R2 unita alla densita' ottica totale > OD_THR. La maschera di R2
    (saturazione > 0.06 e valore < 0.97) perde il tessuto rosa pallido (neuropilo del cervello, miocardio pallido: S1.3c)."""
    hsv = color.rgb2hsv(rgb)
    s = max(0.5, smooth_um / upp)
    sat = filters.gaussian(hsv[..., 1], sigma=s); val = filters.gaussian(hsv[..., 2], sigma=s)
    od = filters.gaussian(od_sum(rgb), sigma=s)
    m = ((sat > 0.06) & (val < 0.97)) | (od > OD_THR)
    if close_um > 0:
        m = morphology.closing(m, morphology.disk(max(1, int(round(close_um / upp)))))
    m = morphology.remove_small_objects(m, max_size=int(500 / upp ** 2))
    return morphology.remove_small_holes(m, max_size=int(100 / upp ** 2))

def tissue_mask_alt(rgb, upp):
    """S-1: come tissue_mask di R2 ma lisciatura 0.5 µm e senza chiusura."""
    hsv = color.rgb2hsv(rgb)
    sat = filters.gaussian(hsv[..., 1], sigma=max(0.5, 0.5 / upp)); val = filters.gaussian(hsv[..., 2], sigma=max(0.5, 0.5 / upp))
    m = (sat > 0.06) & (val < 0.97)
    m = morphology.remove_small_objects(m, max_size=int(500 / upp ** 2))
    return morphology.remove_small_holes(m, max_size=int(100 / upp ** 2))

def img_path(ds):
    import glob
    c = [p for p in glob.glob(f"{D}/{ds}/{ds}_tissue_image.*") if p.endswith((".btf", ".tif", ".tiff"))]
    assert len(c) == 1, c; return c[0]

def sub_points(R0, C0, nbins):
    """Coordinate di array (indici continui) dei centri dei sottopixel di una regione che parte dal bin (R0, C0)."""
    j = (np.arange(nbins * C.SUB) + 0.5) / C.SUB - 0.5
    return R0 + j, C0 + j

def region_masks(p, H, W, upp, M, t, Rr, Cr, nb):
    """Maschere H&E (R2 primaria, S-1 alternativa, bolle) e copertura dell'immagine sulla sottogriglia di una regione di nb x nb bin."""
    rr, cc = sub_points(Rr, Cr, nb)
    RR, CC = np.meshgrid(rr, cc, indexing="ij")
    PR = M[0, 0] * RR + M[0, 1] * CC + t[0]; PC = M[1, 0] * RR + M[1, 1] * CC + t[1]
    r0 = int(max(0, np.floor(PR.min()) - 4)); r1 = int(min(H, np.ceil(PR.max()) + 4))
    c0 = int(max(0, np.floor(PC.min()) - 4)); c1 = int(min(W, np.ceil(PC.max()) + 4))
    rgb = read_region(p, r0, r1, c0, c1)
    rgb2 = downsample(rgb, 2); u2 = upp * 2
    if os.environ.get("R5_MASK_R2") == "1":
        tm, ta = tissue_mask(rgb2, u2), tissue_mask_alt(rgb2, u2)
    else:
        tm, ta = tissue_mask_r5(rgb2, u2), tissue_mask_r5(rgb2, u2, smooth_um=0.5, close_um=0.0)
    bm = bubble_mask(rgb2, u2)
    if os.environ.get("R5_MASK_R2") != "1":
        tm, ta = tm & ~bm, ta & ~bm          # gli anelli delle bolle hanno densita' ottica da tessuto: non sono sorgenti
    ir = np.floor((PR - r0) / 2).astype(int); ic = np.floor((PC - c0) / 2).astype(int)
    ok = (PR >= 0) & (PR < H) & (PC >= 0) & (PC < W) & (ir >= 0) & (ir < tm.shape[0]) & (ic >= 0) & (ic < tm.shape[1])
    irc, icc = np.clip(ir, 0, tm.shape[0] - 1), np.clip(ic, 0, tm.shape[1] - 1)
    f = np.where(ok, tm[irc, icc], False); fa = np.where(ok, ta[irc, icc], False); bub = np.where(ok, bm[irc, icc], False)
    return f, fa, bub, ok

def pixel_grids(ds, nr, nc):
    """Coordinate nell'immagine (pixel) del centro di ogni bin a 2 µm, dalla registrazione di Space Ranger."""
    import pyarrow.parquet as pq
    tp = pq.read_table(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/tissue_positions.parquet",
                       columns=["array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]).to_pandas()
    PR = np.full((nr, nc), np.nan); PC = np.full((nr, nc), np.nan)
    PR[tp.array_row.to_numpy(), tp.array_col.to_numpy()] = tp.pxl_row_in_fullres.to_numpy()
    PC[tp.array_row.to_numpy(), tp.array_col.to_numpy()] = tp.pxl_col_in_fullres.to_numpy()
    return PR, PC

def local_affine(PXR, PXC, Rr, Cr, nb):
    """Affine locale griglia -> immagine stimata sui bin della regione (la registrazione di Space Ranger non e' esattamente affine
    sull'intera area: residuo massimo globale fino a 0.94 pixel)."""
    nr, nc = PXR.shape
    a0, a1, b0, b1 = max(Rr, 0), min(Rr + nb, nr), max(Cr, 0), min(Cr + nb, nc)
    rr, cc = np.mgrid[a0:a1, b0:b1]
    X = np.c_[rr.ravel(), cc.ravel(), np.ones(rr.size)]
    pr, pc = PXR[a0:a1, b0:b1].ravel(), PXC[a0:a1, b0:b1].ravel()
    ok = np.isfinite(pr) & np.isfinite(pc)
    cr, *_ = np.linalg.lstsq(X[ok], pr[ok], rcond=None); cc_, *_ = np.linalg.lstsq(X[ok], pc[ok], rcond=None)
    res = max(np.max(np.abs(X[ok] @ cr - pr[ok])), np.max(np.abs(X[ok] @ cc_ - pc[ok])))
    return np.array([[cr[0], cr[1]], [cc_[0], cc_[1]]]), np.array([cr[2], cc_[2]]), float(res)

def image_info(ds):
    import tifffile
    p = img_path(ds)
    with tifffile.TiffFile(p) as tf:
        sh = tf.series[0].shape; H, W = (sh[1], sh[2]) if tf.series[0].axes[0] in "SC" else (sh[0], sh[1])
    upp = json.load(open(f"{D}/{ds}/binned_outputs_x/binned_outputs/square_002um/spatial/scalefactors_json.json"))["microns_per_pixel"]
    return p, H, W, upp

if __name__ == "__main__":
    KEY, OUT = sys.argv[1], sys.argv[2]; ds = DS[KEY]
    os.makedirs(f"{OUT}/tiles/{KEY}", exist_ok=True)
    t0 = time.time()
    z = np.load(f"{OUT}/section/{KEY}.npz"); U, IT, M, t = z["U"], z["in_tissue"], z["M"], z["t"]
    nr, nc = U.shape
    PXR, PXC = pixel_grids(ds, nr, nc)
    p, H, W, upp = image_info(ds)
    to_px = lambda r, c: (M[0, 0] * r + M[0, 1] * c + t[0], M[1, 0] * r + M[1, 1] * c + t[1])
    # fondo del vetro lontano (D-2): bin fuori dal tessuto di Space Ranger a >= 300 µm da esso e a >= 100 µm dal bordo dell'area
    dist = ndimage.distance_transform_edt(~IT) * C.BIN_UM
    inner = np.zeros_like(IT); inner[BORDER_BINS:-BORDER_BINS, BORDER_BINS:-BORDER_BINS] = True
    far = (~IT) & inner & (dist >= FAR_UM)
    rule = "distance >= 300 um"
    if far.sum() < FAR_MIN_BINS:
        cand = (~IT) & inner
        thr = np.sort(dist[cand])[-FAR_MIN_BINS]; far = cand & (dist >= thr); rule = f"farthest {FAR_MIN_BINS} bins (>= {thr:.0f} um)"
    bg = {"dataset": ds, "b_far": float(U[far].mean()), "n_bins": int(far.sum()), "rule": rule, "median_distance_um": float(np.median(dist[far])),
          "b_far_by_parity": [float(U[far & (C.parity_class(*np.indices(U.shape)) == k)].mean()) for k in range(4)]}
    json.dump(bg, open(f"{OUT}/tiles/{KEY}_background.json", "w"), indent=1)
    rows = []
    for i in range(nr // TB):
        for j in range(nc // TB):
            R0, C0 = i * TB, j * TB
            rec = {"tile": f"t{R0}_{C0}", "R0": R0, "C0": C0}
            fr = IT[R0:R0 + TB, C0:C0 + TB].mean(); rec["sr_frac"] = float(fr)
            if not (0.2 <= fr <= 0.8):
                rec["status"] = "sr_frac"; rows.append(rec); continue
            if R0 < BORDER_BINS or C0 < BORDER_BINS or R0 + TB > nr - BORDER_BINS or C0 + TB > nc - BORDER_BINS:
                rec["status"] = "array_border"; rows.append(rec); continue
            cr = [to_px(r, c) for r in (R0 - 0.5, R0 + TB - 0.5) for c in (C0 - 0.5, C0 + TB - 0.5)]
            if not all(0 <= a < H and 0 <= b < W for a, b in cr):
                rec["status"] = "outside_image"; rows.append(rec); continue
            # regione = tessera + margine; maschere H&E sulla sottogriglia
            Rr, Cr, nb = R0 - MB, C0 - MB, TB + 2 * MB
            Ml, tl_, rec["affine_local_maxres_px"] = local_affine(PXR, PXC, Rr, Cr, nb)
            f, fa, bub, ok = region_masks(p, H, W, upp, Ml, tl_, Rr, Cr, nb)
            # regione in bin, con i bin fuori dall'array marcati
            Ureg = np.zeros((nb, nb), np.int32); inarr = np.zeros((nb, nb), bool)
            a0, a1, b0, b1 = max(Rr, 0), min(Rr + nb, nr), max(Cr, 0), min(Cr + nb, nc)
            Ureg[a0 - Rr:a1 - Rr, b0 - Cr:b1 - Cr] = U[a0:a1, b0:b1]; inarr[a0 - Rr:a1 - Rr, b0 - Cr:b1 - Cr] = True
            occ = C.bin_occupancy(f.astype(float)); dist_b = C.edge_distance_bins(f.astype(float))
            bub_bin = C.bin_occupancy(bub.astype(float)) > 0
            inimg_bin = C.bin_occupancy(ok.astype(float)) >= 1.0
            tl = (slice(MB, MB + TB), slice(MB, MB + TB))
            rec["he_frac"] = float(f[MB * 4:(MB + TB) * 4, MB * 4:(MB + TB) * 4].mean())
            rec["bubble_frac"] = float(bub[MB * 4:(MB + TB) * 4, MB * 4:(MB + TB) * 4].mean())
            g = (occ[tl] == 0) & (dist_b[tl] >= 10) & (dist_b[tl] <= 100)
            rec["n_glass_10_100"] = int(g.sum())
            rec["he_frac_alt"] = float(fa[MB * 4:(MB + TB) * 4, MB * 4:(MB + TB) * 4].mean())
            if not (0.2 <= rec["he_frac"] <= 0.8):
                rec["status"] = "he_frac"
            elif rec["bubble_frac"] >= 0.01:
                rec["status"] = "bubbles"
            elif rec["n_glass_10_100"] < 2000:
                rec["status"] = "glass_10_100"
            else:
                rec["status"] = "admitted"
                np.savez_compressed(f"{OUT}/tiles/{KEY}/{rec['tile']}.npz", U=Ureg, f_sub=f.astype(np.uint8), f_alt=fa.astype(np.uint8),
                                    bubble_bin=bub_bin, inimg_bin=inimg_bin, inarr_bin=inarr, R0=Rr, C0=Cr, tile_R0=R0, tile_C0=C0, M_local=Ml, t_local=tl_)
            rows.append(rec)
    T = pd.DataFrame(rows); T.to_csv(f"{OUT}/tiles/{KEY}_candidates.csv", index=False)
    print(f"[{KEY}] candidates {len(T)} status {T.status.value_counts().to_dict()} b_far {bg['b_far']:.3f} ({bg['rule']}) {time.time()-t0:.0f}s", flush=True)
