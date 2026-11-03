"""tools/R5_ph5_narrowcore.py — R5, post hoc PH-5 (dopo PH-2): i bordi del tessuto con la parte gaussiana di K2 fissata alla
sfocatura interna (sigma = SIG µm, default 0.7; coda w, ell e scarto liberi). Confronto di verosimiglianza con K2 libero, stessa tessera,
stessi bin, stesso stimatore (b libero, rho da tutti i bin di tessuto). Non modifica i json delle tessere.
Uso: R5_PH5_SIG=0.7 .venv/bin/python tools/R5_ph5_narrowcore.py <outdir> <workers> -> <outdir>/ph5_s<SIG>/<KEY>_<tile>.json"""
import sys, os, json, glob
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R5_common as C
import importlib.util
OUT, NW = sys.argv[1], int(sys.argv[2])
SIG = float(os.environ.get("R5_PH5_SIG", "0.7"))
ODIR = f"{OUT}/ph5_s{SIG:g}"
spec = importlib.util.spec_from_file_location("r5fit", "tools/R5_fit.py")
_orig_bounds = C._bounds
def _bounds_fixed(model):
    b = _orig_bounds(model)
    if model == "K2": b = [(np.log(SIG), np.log(SIG))] + b[1:]
    return b

def one(path):
    key = path.split("/")[-2]; tile = os.path.basename(path)[:-4]; fo = f"{ODIR}/{key}_{tile}.json"
    if os.path.exists(fo): return fo
    fj = path.replace(".npz", ".json"); d = json.load(open(fj))
    if "error" in d: return fo
    sys.argv = ["R5_fit.py", OUT]; F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
    C._bounds = _bounds_fixed
    z = np.load(path); y = z["U"].astype(float); f = z["f_sub"].astype(float)
    R0, C0 = int(z["R0"]), int(z["C0"]); rows, cols = np.indices(y.shape); rows = rows + R0; cols = cols + C0
    valid = z["inarr_bin"] & z["inimg_bin"] & ~z["bubble_bin"]
    P0 = F.region_P(y, f, rows, cols, valid); sel = C.likelihood_mask(P0["occ"], P0["dist"], F.TILE) & valid
    k2 = d["K2"]
    starts = [dict(k2, sigma=SIG), dict(k2, sigma=SIG, w=0.2, ell=20.0), dict(k2, sigma=SIG, w=0.05, ell=60.0)]
    r = C.fit_model("K2", P0["Tft"], P0["kr"], P0["kc"], P0["pshape"], P0["region_bins"], y, P0["c"], sel, starts)
    C._bounds = _orig_bounds
    res = {"key": key, "tile": tile, "sigma_fixed": SIG, "fit": r, "nll_K2_free": k2["nll"], "sigma_K2_free": k2["sigma"], "w_K2_free": k2["w"],
           "delta_nll": r["nll"] - k2["nll"], "n_bins": r["n_bins"]}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo

if __name__ == "__main__":
    os.makedirs(ODIR, exist_ok=True)
    jobs = sorted(p for key in ("A1", "A23", "A4", "A5", "A6") for p in glob.glob(f"{OUT}/tiles/{key}/t*.npz") if os.path.exists(p.replace(".npz", ".json")))
    with Pool(NW) as pool:
        for i, f in enumerate(pool.imap_unordered(one, jobs, chunksize=1), 1):
            if i % 20 == 0: print(i, len(jobs), flush=True)
    print("PH5 done", len(jobs), flush=True)
