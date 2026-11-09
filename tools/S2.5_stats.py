"""tools/S2.5_stats.py — S2.5: misure del pannello simulato con le definizioni del lato vero (tools/S2.5_real.py, definizione reg).
Quadratini con il centro nelle regioni dello Step 1 e a >= 25 um dal bordo dell'area; blocchi di 8 um con almeno 8 dei 16 quadratini;
UMI medio, frazione di vuoti, geni rilevati a 2 e 8 um; anisotropia dei blocchi di 8 um (tools/S2_5_common.py); rapporto nucleo/corona con la funzione del pannello di S2.3 (panel_star).
Uso: .venv/bin/python tools/S2.5_stats.py <panel_dir> <variante> -> <panel_dir>/stats_<variante>/<A>_<roi>_k<k>.json (processi da WORKERS)"""
import sys, os, json, glob, time, traceback
import numpy as np, pandas as pd, shapely
from multiprocessing import Pool
sys.path.insert(0, "tools")
import S2_5_common as CM
import importlib.util
def _load(name, path):
    sp = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m
FRAME_UM = 25.0; MIN_CHILD8 = 8

stats = CM.stats

def run(args):
    f, Rw, OUT = args
    tag_k = os.path.basename(f)[:-5]; fo = f"{OUT}/{tag_k}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    try:
        RL = _load("s23real", "tools/S2.3_real.py")
        meta = json.load(open(f)); tag = meta["tag"]; g = meta["grid"]; nr, nc, bu, bl = g["n_row"], g["n_col"], g["bin_um"], g["block"]
        side = float(Rw.side_px * Rw.um_per_px)
        base = f[:-5]
        u2 = np.fromfile(base + ".bin", dtype=np.int32); g2 = np.fromfile(base + ".g2.bin", dtype=np.int32)
        u8 = np.fromfile(base + ".u8.bin", dtype=np.int32); g8 = np.fromfile(base + ".g8.bin", dtype=np.int32)
        assert len(u2) == nr * nc and len(u8) == (nr // bl) * (nc // bl)
        i, j = np.divmod(np.arange(nr * nc), nc)
        x = g["x0"] + (j + 0.5) * bu; y = g["y0"] + (i + 0.5) * bu
        R = RL.region_geom(tag)
        fr = (y >= FRAME_UM) & (y <= side - FRAME_UM) & (x >= FRAME_UM) & (x <= side - FRAME_UM)
        P = np.zeros(nr * nc, bool); P[fr] = shapely.contains_xy(R, x[fr], y[fr])
        res = {"tag": tag, "archetype": meta["archetype"], "k": meta["k"], "variant": meta["variant"]}
        for k, v in stats(u2[P], g2[P]).items(): res[f"reg_{k}2"] = v
        b8 = (i // bl) * (nc // bl) + j // bl
        nP = np.bincount(b8, weights=P.astype(float), minlength=len(u8))
        sel = nP >= MIN_CHILD8
        for k, v in stats(u8[sel], g8[sel]).items(): res[f"reg_{k}8"] = v
        res["reg_aniso8"], res["reg_angle8"] = CM.aniso(u8.reshape(nr // bl, nc // bl), sel.reshape(nr // bl, nc // bl))
        res["n_tissue_no_territory"] = meta.get("n_tissue_no_territory"); res["n_in_tissue"] = meta.get("n_in_tissue")
        # rapporto nucleo/corona: funzione del pannello di S2.3 sullo stesso file di totali (scrive in OUT/s23/)
        PS = _load("s23ps", "tools/S2.3_panel_stats.py")
        os.makedirs(f"{OUT}/s23/stats", exist_ok=True)
        tmpd = f"{OUT}/s23"
        fj = PS.run((f, Rw, tmpd)); j23 = json.load(open(fj))
        if "error" in j23: raise RuntimeError(j23["error"])
        res["reg_star"] = float(j23["panel_star"]); res["reg_star_n_bins"] = int(j23["panel_n_bins"])
        res["sec"] = round(time.time() - t0, 1)
    except Exception:
        res = {"file": f, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

if __name__ == "__main__":
    PD, VAR = sys.argv[1], sys.argv[2]; OUT = f"{PD}/stats_{VAR}"
    os.makedirs(OUT, exist_ok=True)
    rois = pd.read_csv("results/R2/R2_rois_checked.csv"); rois["tag"] = rois.archetype + "_" + rois.roi_id
    files = sorted(f for f in glob.glob(f"{PD}/{VAR}/*_k*.json"))
    jobs = [(f, rois.set_index("tag").loc[json.load(open(f))["tag"]], OUT) for f in files]
    nw = int(open(f"{PD}/WORKERS").read().strip()) if os.path.exists(f"{PD}/WORKERS") else 4
    print(f"stats {len(jobs)} workers {nw}", flush=True); t0 = time.time()
    with Pool(min(nw, max(1, len(jobs)))) as p:
        for k_, fo in enumerate(p.imap_unordered(run, jobs, chunksize=1), 1): print(f"{k_}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    S = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/*_k*.json"))])
    if "error" in S.columns and S["error"].notna().any(): print(S.loc[S["error"].notna(), "error"].iloc[0]); sys.exit(1)
    S.to_csv(f"{PD}/S2.5_sim_{VAR}.csv", index=False); print(S.round(4).to_string(index=False))
