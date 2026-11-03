"""tools/R5_fit_s2.py — R5: ricalcola la sola sensibilita' S-2 (parita' stimata dal vetro) con la densita' del tessuto come nella stima primaria
(tutti i bin interamente di tessuto). La prima versione di tools/R5_fit.py cambiava insieme parita' e densita' (rilievo F13 della revisione del codice).
Stesso avvio (stima K2 primaria della tessera). Uso: .venv/bin/python tools/R5_fit_s2.py <outdir> <workers>."""
import sys, os, json, glob
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R5_common as C
import importlib.util
OUT, NW = sys.argv[1], int(sys.argv[2])
spec = importlib.util.spec_from_file_location("r5fit", "tools/R5_fit.py")

def one(path):
    fj = path.replace(".npz", ".json"); d = json.load(open(fj))
    if "error" in d or d.get("S2_fixed_ok"): return fj
    sys.argv = ["R5_fit.py", OUT]; F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
    z = np.load(path); y = z["U"].astype(float); f = z["f_sub"].astype(float)
    R0, C0 = int(z["R0"]), int(z["C0"]); rows, cols = np.indices(y.shape); rows = rows + R0; cols = cols + C0
    valid = z["inarr_bin"] & z["inimg_bin"] & ~z["bubble_bin"]
    P0 = F.region_P(y, f, rows, cols, valid); sel = C.likelihood_mask(P0["occ"], P0["dist"], F.TILE) & valid
    Pg = F.region_P(y, f, rows, cols, valid, pf_from_glass=True)
    d["S2_parity_glass"] = C.fit_model("K2", Pg["Tft"], Pg["kr"], Pg["kc"], Pg["pshape"], Pg["region_bins"], y, Pg["c"], sel, [d["K2"]], None)
    d["pf_glass"] = Pg["pf"].tolist(); d["S2_fixed_ok"] = True
    tmp = fj + ".tmp"; json.dump(d, open(tmp, "w")); os.replace(tmp, fj)
    return fj

if __name__ == "__main__":
    jobs = sorted(p for key in ("A1", "A23", "A4", "A5", "A6") for p in glob.glob(f"{OUT}/tiles/{key}/t*.npz") if os.path.exists(p.replace(".npz", ".json")))
    with Pool(NW) as pool: list(pool.imap_unordered(one, jobs, chunksize=1))
    print("S2 patch done", len(jobs), flush=True)
