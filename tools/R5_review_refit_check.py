"""tools/R5_review_refit_check.py — R5, controllo di accettazione della correzione dell'ottimizzatore (revisione numerica F-1, F-3):
per le 37 tessere adattate da zero dal revisore, K2 primario con la nuova procedura (C.fit_k2_screened, w in [0, BOUNDS]) sui file di run1;
NLL nuova - NLL salvata (prima della correzione) confrontata con dnll del revisore (sua NLL - salvata, nella sua implementazione).
PASS se nuova - salvata <= dnll_revisore + 0.01 in ogni tessera.
Uso: [R5_WMAX=0.5] .venv/bin/python tools/R5_review_refit_check.py <run1> <workers> -> results/R5/review/R5_review_refit_check[_w05].csv"""
import sys, os, json
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R5_common as C
import importlib.util
OUT, NW = sys.argv[1], int(sys.argv[2])
spec = importlib.util.spec_from_file_location("r5fit", "tools/R5_fit.py")

def one(a):
    key, tile = a
    path = f"{OUT}/tiles/{key}/{tile}.npz"; d = json.load(open(path.replace(".npz", ".json")))
    sys.argv = ["R5_fit.py", OUT]; F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
    z = np.load(path); y = z["U"].astype(float); f = z["f_sub"].astype(float)
    R0, C0 = int(z["R0"]), int(z["C0"]); rows, cols = np.indices(y.shape); rows = rows + R0; cols = cols + C0
    valid = z["inarr_bin"] & z["inimg_bin"] & ~z["bubble_bin"]
    P = F.region_P(y, f, rows, cols, valid); sel = C.likelihood_mask(P["occ"], P["dist"], F.TILE) & valid
    old = d["K2"]
    # NLL della vecchia soluzione con la profilazione corretta di (a, b) (F-3)
    g = C.expected_bins(P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], old["sigma_r"], old["sigma_c"], old["w"], old["ell"], old["dr"], old["dc"])[sel]
    nll_old_fixedab = C.fit_ab(y[sel].astype(float), g, P["c"][sel].astype(float))[2]
    k2, k1 = C.fit_k2_screened(P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, None)
    return {"key": key, "tile": tile, "nll_saved": old["nll"], "nll_saved_refit_ab": nll_old_fixedab, "nll_new": k2["nll"],
            "dnll_new_vs_saved": k2["nll"] - old["nll"], "sigma_saved": old["sigma"], "sigma_new": k2["sigma"], "w_saved": old["w"], "w_new": k2["w"],
            "ell_saved": old["ell"], "ell_new": k2["ell"], "dr_saved": old["dr"], "dr_new": k2["dr"], "dc_saved": old["dc"], "dc_new": k2["dc"],
            "b_new": k2["b"], "converged_new": k2["converged"], "wmax": C.BOUNDS["w"][1]}

if __name__ == "__main__":
    rv = pd.read_csv("results/R5/review/R5_review_numeric_refit_summary.csv")
    with Pool(NW) as p: rows = p.map(one, list(zip(rv.key, rv.tile)), chunksize=1)
    T = pd.DataFrame(rows).merge(rv[["key", "tile", "dnll", "sigma_r", "w_r", "ell_r", "dr_r", "dc_r"]].rename(columns={"dnll": "dnll_reviewer"}), on=["key", "tile"])
    T["pass"] = T.dnll_new_vs_saved <= T.dnll_reviewer + 0.01
    sfx = "" if C.BOUNDS["w"][1] == 1.0 else f"_w{str(C.BOUNDS['w'][1]).replace('.', '')}"
    os.makedirs("results/R5/review", exist_ok=True); T.to_csv(f"results/R5/review/R5_review_refit_check{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(T[["key", "tile", "dnll_reviewer", "dnll_new_vs_saved", "sigma_saved", "sigma_new", "sigma_r", "w_saved", "w_new", "w_r", "pass"]].round(3).to_string(index=False))
    print(f"PASS {int(T['pass'].sum())}/{len(T)} (wmax={C.BOUNDS['w'][1]})")
