"""tools/S2.3_panel_stats.py — S2.3: statistiche del pannello sui quadratini simulati (stesse funzioni del lato vero, tools/S2.3_real.py).
Uso: .venv/bin/python tools/S2.3_panel_stats.py <panel_dir> <variante> -> <panel_dir>/<variante>/stats/<A>_<roi>_k<k>.json (processi da WORKERS)
Nuclei simulati rasterizzati a 0.25 um (etichetta = posizione nel livello + 1); posizioni vere (nessuno scarto); fattori di cattura c = 1."""
import sys, os, json, glob, time, traceback
import numpy as np, pandas as pd, shapely
from multiprocessing import Pool
from skimage.draw import polygon as skpoly
sys.path.insert(0, "tools")
import importlib.util
_spec = importlib.util.spec_from_file_location("s23real", "tools/S2.3_real.py"); RL = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(RL)
S, C = RL.S, RL.C
from R6_common import footprint_labels
UPP = 0.25

def nuclei_mask(tag, side):
    t = pd.read_csv(f"{RL.LAYERS}/{tag}_nuclei.tsv", sep="\t")
    g = shapely.from_wkb([bytes.fromhex(h) for h in t.wkb])
    n = int(np.ceil(side / UPP)); M = np.zeros((n, n), np.int32)
    for k, p in enumerate(g):
        if p.is_empty: continue
        for q in getattr(p, "geoms", [p]):
            xy = np.asarray(q.exterior.coords)
            # centro del pixel (riga r, colonna c) = ((r + 0.5) UPP, (c + 0.5) UPP): coordinate in pixel meno 0.5
            rr, cc = skpoly(xy[:, 1] / UPP - 0.5, xy[:, 0] / UPP - 0.5, shape=M.shape)
            M[rr, cc] = k + 1
    return M, len(g)

def run(args):
    f, Rw, OUT = args
    tag_k = os.path.basename(f)[:-5]; fo = f"{OUT}/stats/{tag_k}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        meta = json.load(open(f)); tag = meta["tag"]; g = meta["grid"]
        side = float(Rw.side_px * Rw.um_per_px)
        tot = np.fromfile(f[:-5] + ".bin", dtype=np.int32).reshape(g["n_row"], g["n_col"]).astype(float)
        i, j = np.divmod(np.arange(g["n_row"] * g["n_col"]), g["n_col"])
        x = g["x0"] + (j + 0.5) * g["bin_um"]; y = g["y0"] + (i + 0.5) * g["bin_um"]
        box = (x >= 0) & (x < side) & (y >= 0) & (y < side)
        x, y, i, j = x[box], y[box], i[box], j[box]; u = tot.ravel()[box]
        pos = np.c_[y, x]
        mask, nn = nuclei_mask(tag, side)
        Mp = S.nuclear_maps(mask, UPP)
        rr, cc = np.floor(y / UPP).astype(np.int64), np.floor(x / UPP).astype(np.int64)
        R = RL.region_geom(tag)
        inreg = shapely.contains_xy(R, x, y)
        fr = (y >= RL.FRAME_UM) & (y <= side - RL.FRAME_UM) & (x >= RL.FRAME_UM) & (x <= side - RL.FRAME_UM)
        P = inreg & fr
        keep = np.arange(1, nn + 1)
        L = {"nuc": footprint_labels(mask, rr, cc, UPP, 8.0, "nuc"), "exp": footprint_labels(mask, rr, cc, UPP, 8.0, "exp")}
        s = int(round(20.0 / UPP)); ms = np.roll(np.roll(mask, s, 0), s, 1)
        L["nuc_shift"] = footprint_labels(ms, rr, cc, UPP, 8.0, "nuc"); L["exp_shift"] = footprint_labels(ms, rr, cc, UPP, 8.0, "exp")
        Ms = S.nuclear_maps(ms, UPP)
        far = {"far": RL.far_flags(Mp, rr, cc, L["nuc"]), "far_shift": RL.far_flags(Ms, rr, cc, L["nuc_shift"])}
        LP = {k: np.where(P, v, 0) for k, v in L.items()}; FP = {k: v & P for k, v in far.items()}
        res = {"tag": tag, "archetype": meta["archetype"], "k": meta["k"], "variant": meta["variant"]}
        res["panel_star"] = C.star_fast(u, C.stat_index(LP, keep))
        res["panel_star_far"] = C.star_fast(u, C.stat_index(LP, keep, FP))
        dn = S.lookup(Mp["d_out"], rr, cc, np.inf)
        res["panel_radial_R"] = S.profile_ratio(dn[P], u[P])
        GA = S.groups_from_maps(Mp, rr, cc, keep, inreg, i, j)
        c1 = np.ones(len(u))
        for k in ("nucleus", "corona"):
            st = S.group_stats(u, np.where(P, GA[k], -1), c1, i, j)
            for kk, v in st.items(): res[f"panel_{k}_{kk}"] = v
        res["panel_n_bins"] = int(P.sum()); res["panel_mean_umi"] = float(u[P].mean())
        res.update({f"lumen_{k}": v for k, v in RL.lumen_stats(R, pos, u, inreg, side).items()})
        res["sec"] = round(time.time() - t0, 1)
    except Exception:
        res = {"file": f, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

if __name__ == "__main__":
    PD, VAR = sys.argv[1], sys.argv[2]; OUT = f"{PD}/{VAR}"
    os.makedirs(f"{OUT}/stats", exist_ok=True)
    rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois["tag"] = rois.archetype + "_" + rois.roi_id
    files = sorted(glob.glob(f"{OUT}/*_k*.json"))
    jobs = [(f, rois.set_index("tag").loc[json.load(open(f))["tag"]], OUT) for f in files]
    nw = int(open(f"{PD}/WORKERS").read().strip()) if os.path.exists(f"{PD}/WORKERS") else 4
    print(f"stats {len(jobs)} workers {nw}", flush=True); t0 = time.time()
    with Pool(min(nw, max(1, len(jobs)))) as p:
        for k_, fo in enumerate(p.imap_unordered(run, jobs, chunksize=1), 1): print(f"{k_}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
