"""tools/R5_fit_s0.py — R5: ricalcola la sola sensibilita' S-0 (fondo fissato dal vetro lontano) nelle tessere in cui il primo calcolo
l'aveva eseguita per errore con il fondo libero (sostituzione di testo in tools/R5_fit.py, scoperta durante il calcolo del 2026-10-31).
Stesso avvio (la stima K2 primaria della tessera) e stessi dati di tools/R5_fit.py: risultato identico a quello del codice corretto.
Uso: .venv/bin/python tools/R5_fit_s0.py <outdir>. Processi da <outdir>/WORKERS."""
import sys, os, json, glob
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R5_common as C
import importlib.util
spec = importlib.util.spec_from_file_location("r5fit", "tools/R5_fit.py")
OUT = sys.argv[1]

def one(args):
    path, bfar = args
    fj = path.replace(".npz", ".json"); d = json.load(open(fj))
    if "error" in d or d.get("S0_fixed_ok"): return fj
    sys.argv = ["R5_fit.py", OUT]; F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
    z = np.load(path); y = z["U"].astype(float); f = z["f_sub"].astype(float)
    R0, C0 = int(z["R0"]), int(z["C0"]); rows, cols = np.indices(y.shape); rows = rows + R0; cols = cols + C0
    valid = z["inarr_bin"] & z["inimg_bin"] & ~z["bubble_bin"]
    P = F.region_P(y, f, rows, cols, valid); sel = C.likelihood_mask(P["occ"], P["dist"], F.TILE) & valid
    d["S0_bfixed"] = C.fit_model("K2", P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, [d["K2"]], bfar)
    d["S0_fixed_ok"] = True
    tmp = fj + ".tmp"; json.dump(d, open(tmp, "w")); os.replace(tmp, fj)
    return fj

if __name__ == "__main__":
    jobs = []
    for key in ("A1", "A23", "A4", "A5", "A6"):
        bfar = json.load(open(f"{OUT}/tiles/{key}_background.json"))["b_far"]
        jobs += [(p, bfar) for p in sorted(glob.glob(f"{OUT}/tiles/{key}/t*.npz")) if os.path.exists(p.replace(".npz", ".json"))]
    nw = int(open(f"{OUT}/WORKERS").read().strip())
    with Pool(nw) as pool: list(pool.imap_unordered(one, jobs, chunksize=1))
    print("S0 patch done", len(jobs), flush=True)
