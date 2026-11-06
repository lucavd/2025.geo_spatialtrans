"""tools/S2.3_real.py — S2.3: statistiche delle 30 aree vere (lato vero del banco e del pannello; esplorazione E-2).
Definizioni: docs/steps/S2.3.md («Il banco della taratura», «La controprova sulle aree simulate»). Maschere lette allineate con lo
scarto della sezione (regola C2, parametri corretti di S2.0 run1_wfix). Nessuna simulazione.
Uso: .venv/bin/python tools/S2.3_real.py <outdir> [A_roi ...] -> <outdir>/real/<A>_<roi>.json (processi dal file <outdir>/WORKERS)"""
import sys, os, json, time, traceback
import numpy as np, pandas as pd, pyarrow.parquet as pq
import shapely
from multiprocessing import Pool
sys.path.insert(0, "tools")
R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
S20 = "/mnt/micron/geo_spatialtrans/S2.0/run1_wfix"
LAYERS = "/mnt/micron/geo_spatialtrans/S2.3/layers"
_saved = list(sys.argv); sys.argv = ["R5_nc.py", R5RUN, f"{R5RUN}/kernels.json"]
import importlib.util
_spec = importlib.util.spec_from_file_location("r5nc", "tools/R5_nc.py"); NC = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(NC)
sys.argv = _saved
import R5_common as C
import S2_0_common as S0
import S2_3_common as S
import R5_tiles as T5
from R2_common import bubble_mask, downsample
from R1_verify import read_region
from R6_common import footprint_labels

FRAME_UM = 25.0
LUMEN_MIN_UM2 = 2500.0
BANDS = [(2.0, 4.0), (4.0, 8.0)]

def roi_tissue(Rw, Dd):
    """Copia di tools/S2.0_disp.py roi_tissue (non importabile: lo script legge argv all'import). Tessuto H&E della ROI sulla
    griglia di 0.5 um (maschera di R5: R2 unita a densita' ottica > 0.15, bolle escluse)."""
    p, H, W, upp_img = T5.image_info(NC.DS[Rw.archetype])
    rgb = read_region(p, int(Rw.r0), int(Rw.r1), int(Rw.c0), int(Rw.c1))
    assert rgb.shape[:2] == Dd["mask"].shape, (rgb.shape, Dd["mask"].shape)
    rgb2 = downsample(rgb, 2); u2 = Dd["upp"] * 2
    tm = T5.tissue_mask_r5(rgb2, u2) & ~bubble_mask(rgb2, u2)
    H0, W0 = Dd["mask"].shape; step = S0.SUB_UM
    n0, n1 = int(np.floor(H0 * Dd["upp"] / step)), int(np.floor(W0 * Dd["upp"] / step))
    ir = np.minimum(((np.arange(n0) + 0.5) * step / Dd["upp"] / 2).astype(int), tm.shape[0] - 1)
    ic = np.minimum(((np.arange(n1) + 0.5) * step / Dd["upp"] / 2).astype(int), tm.shape[1] - 1)
    return tm[np.ix_(ir, ic)].astype(float)

def box_bins(Rw, in_tissue_only):
    sq = f"{NC.D}/{NC.DS[Rw.archetype]}/binned_outputs_x/binned_outputs/square_002um"
    tp = pq.read_table(f"{sq}/spatial/tissue_positions.parquet").to_pandas()
    m = (tp.pxl_row_in_fullres >= Rw.r0) & (tp.pxl_row_in_fullres < Rw.r1) & (tp.pxl_col_in_fullres >= Rw.c0) & (tp.pxl_col_in_fullres < Rw.c1)
    if in_tissue_only: m &= tp.in_tissue == 1
    b = tp[m]
    pos = np.c_[(b.pxl_row_in_fullres.to_numpy() - Rw.r0) * Rw.um_per_px, (b.pxl_col_in_fullres.to_numpy() - Rw.c0) * Rw.um_per_px]
    return b.array_row.to_numpy(), b.array_col.to_numpy(), pos

def stripes_c(A, rows, cols):
    Fz = np.load(f"{S20}/stripes/{NC.KEY[A]}_factors.npz")
    cf = np.exp(np.nan_to_num(Fz["alpha"])[rows] + np.nan_to_num(Fz["beta"])[cols])
    return cf / cf.mean()

def tissue_lookup(tm, pos_al, step=0.5):
    q = np.floor(pos_al / step).astype(np.int64)
    return S.lookup(tm.astype(np.int8), q[:, 0], q[:, 1], 0).astype(bool)

def region_geom(tag):
    t = pd.read_csv(f"{LAYERS}/{tag}_regions.tsv", sep="\t")
    g = shapely.from_wkb([bytes.fromhex(h) for h in t.wkb])
    return shapely.union_all(g)

def far_flags(Mp, rr, cc, lab):
    """Come R5 far_of: nel nucleo distanza dal bordo dentro >= 2 um, fuori distanza dal nucleo >= 2 um."""
    din = S.lookup(Mp["d_in"], rr, cc, 0.0); dout = S.lookup(Mp["d_out"], rr, cc, 0.0)
    ins = (rr >= 0) & (rr < Mp["lab"].shape[0]) & (cc >= 0) & (cc < Mp["lab"].shape[1])
    return ins & np.where(lab > 0, din >= 2.0, dout >= 2.0)

def lumen_stats(R, pos_al, u, inreg, side):
    """Buchi delle regioni interamente nel riquadro (non toccano il bordo) di area >= 2500 um^2. Coordinate dell'area: x = colonna, y = riga."""
    holes = []
    for poly in getattr(R, "geoms", [R]):
        for ring in poly.interiors:
            h = shapely.Polygon(ring)
            if h.area < LUMEN_MIN_UM2: continue
            x0, y0, x1, y1 = h.bounds
            if x0 <= 0 or y0 <= 0 or x1 >= side or y1 >= side: continue
            holes.append(h)
    x, y = pos_al[:, 1], pos_al[:, 0]
    out = {"n_lumens": len(holes), "ref_sum": 0.0, "ref_n": 0}
    for lo, hi in BANDS: out[f"band_{lo:g}_{hi:g}_sum"] = 0.0; out[f"band_{lo:g}_{hi:g}_n"] = 0
    for h in holes:
        x0, y0, x1, y1 = h.bounds
        near = np.flatnonzero((x >= x0 - 8) & (x <= x1 + 8) & (y >= y0 - 8) & (y <= y1 + 8))
        if len(near) == 0: continue
        pts = shapely.points(x[near], y[near])
        inside = shapely.contains_xy(h, x[near], y[near])
        d = shapely.distance(pts, h.exterior)
        ref = (~inside) & inreg[near] & (d < 8.0)
        out["ref_sum"] += float(u[near][ref].sum()); out["ref_n"] += int(ref.sum())
        for lo, hi in BANDS:
            m = inside & (d >= lo) & (d < hi)
            out[f"band_{lo:g}_{hi:g}_sum"] += float(u[near][m].sum()); out[f"band_{lo:g}_{hi:g}_n"] += int(m.sum())
    return out

def run(Rw, OUT):
    tag = f"{Rw.archetype}_{Rw.roi_id}"; fo = f"{OUT}/real/{tag}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        A, upp = Rw.archetype, Rw.um_per_px
        PJ = json.load(open(f"{S20}/params_flat.json")); prm = PJ[NC.KEY[A]]["C2"]
        Dd = NC.roi_data(Rw)                                  # bin in tessuto (Space Ranger), maschera, impronte di R4
        di = NC.d_image(Dd, prm)
        rows, cols, pos = box_bins(Rw, True)
        assert len(rows) == len(Dd["u"]) and np.allclose(pos, Dd["pos"])
        c2 = stripes_c(A, rows, cols)
        Mp = S.nuclear_maps(Dd["mask"], upp)
        res = {"tag": tag, "archetype": A, "roi_id": Rw.roi_id, "d_img_um": di.tolist(), "n_bins_in_tissue": int(len(rows))}
        # G-4 (c): psi osservato del nucleo come S2.0 PH-2 (maschera non allineata, fattori della sezione)
        lab0 = S0.deep_labels(Dd["mask"], Dd["rr"], Dd["cc"], upp, 1.0)
        res["psi_nuc_unaligned_S20def"] = S0.psi_index(Dd["u"], lab0, c2, Dd["keep"])[0]
        # banco: gruppi allineati
        rr, cc = S.aligned_px(Dd["pos"], di, upp)
        tm = roi_tissue(Rw, Dd)
        pos_al = Dd["pos"] - di[None, :]
        tis = tissue_lookup(tm, pos_al)
        G = S.groups_from_maps(Mp, rr, cc, Dd["keep"], tis, rows, cols)
        for k, g in G.items():
            st = S.group_stats(Dd["u"], g, c2, rows, cols)
            for kk, v in st.items(): res[f"bench_{k}_{kk}"] = v
        # pannello: tutti i bin del riquadro, tessuto = regioni dello Step 1 al centro allineato, fuori 25 um dal riquadro
        rowsA, colsA, posA = box_bins(Rw, False)
        z = np.load(f"{R5RUN}/section/{NC.KEY[A]}.npz"); uA = z["U"][rowsA, colsA].astype(float)
        posA_al = posA - di[None, :]
        side = float(Rw.side_px * upp)
        R = region_geom(tag)
        inreg = shapely.contains_xy(R, posA_al[:, 1], posA_al[:, 0])
        fr = (posA_al[:, 0] >= FRAME_UM) & (posA_al[:, 0] <= side - FRAME_UM) & (posA_al[:, 1] >= FRAME_UM) & (posA_al[:, 1] <= side - FRAME_UM)
        P = inreg & fr
        rrA, ccA = S.aligned_px(posA, di, upp)
        L = {"nuc": footprint_labels(Dd["mask"], rrA, ccA, upp, 8.0, "nuc"), "exp": footprint_labels(Dd["mask"], rrA, ccA, upp, 8.0, "exp")}
        s = int(round(20.0 / upp)); ms = np.roll(np.roll(Dd["mask"], s, 0), s, 1)
        L["nuc_shift"] = footprint_labels(ms, rrA, ccA, upp, 8.0, "nuc"); L["exp_shift"] = footprint_labels(ms, rrA, ccA, upp, 8.0, "exp")
        Ms = S.nuclear_maps(ms, upp)
        far = {"far": far_flags(Mp, rrA, ccA, L["nuc"]), "far_shift": far_flags(Ms, rrA, ccA, L["nuc_shift"])}
        LP = {k: np.where(P, v, 0) for k, v in L.items()}
        FP = {k: v & P for k, v in far.items()}
        res["panel_star"] = C.star_fast(uA, C.stat_index(LP, Dd["keep"]))
        res["panel_star_far"] = C.star_fast(uA, C.stat_index(LP, Dd["keep"], FP))
        dn = S.lookup(Mp["d_out"], rrA, ccA, np.inf)
        res["panel_radial_R"] = S.profile_ratio(dn[P], uA[P])
        cA = stripes_c(A, rowsA, colsA)
        GA = S.groups_from_maps(Mp, rrA, ccA, Dd["keep"], inreg, rowsA, colsA)
        for k in ("nucleus", "corona"):
            g = np.where(P, GA[k], -1)
            st = S.group_stats(uA, g, cA, rowsA, colsA)
            for kk, v in st.items(): res[f"panel_{k}_{kk}"] = v
        res["panel_n_bins"] = int(P.sum()); res["panel_mean_umi"] = float(uA[P].mean())
        res.update({f"lumen_{k}": v for k, v in lumen_stats(R, posA_al, uA, inreg, side).items()})
        res["sec"] = round(time.time() - t0, 1)
    except Exception:
        res = {"tag": tag, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

if __name__ == "__main__":
    OUT = sys.argv[1]; ONLY = sys.argv[2:]
    os.makedirs(f"{OUT}/real", exist_ok=True)
    rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    jobs = [Rw for _, Rw in rois.iterrows() if (not ONLY or f"{Rw.archetype}_{Rw.roi_id}" in ONLY)]
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 4
    print(f"real jobs {len(jobs)} workers {nw}", flush=True); t0 = time.time()
    with Pool(min(nw, len(jobs))) as p:
        for i, f in enumerate(p.imap_unordered(run, jobs, chunksize=1) if False else p.starmap(run, [(j, OUT) for j in jobs]), 1):
            pass
    print(f"done {time.time() - t0:.0f}s", flush=True)
