"""tools/R5_nc.py — R5, rapporto nucleo/corona prima della diffusione (P-5, BL-196) e CP-1, per ROI.
Uso: .venv/bin/python tools/R5_nc.py <outdir> <kernels.json> [A_roi ...]   (processi dal file <outdir>/WORKERS; un json per ROI e nucleo)
kernels.json: per sezione {"K2": {...mediane...}, "K1": {...}, "draws": [{...}, ...]} dal bootstrap delle tessere.
Per ogni ROI e ogni nucleo (primario K2, S-5 K1, J estrazioni K2): r_vero per bisezione; statistica senza i bin entro 2 µm dal bordo
del nucleo (CP-1) simulata con r_vero. G-3 (c): totali per cellula e strati uguali alle tabelle di R4. Scrive <outdir>/nc/<A>_<roi>.json."""
import sys, os, json, time, glob, traceback
import numpy as np, pandas as pd, pyarrow.parquet as pq
from multiprocessing import Pool
from scipy import ndimage
sys.path.insert(0, "tools")
import R5_common as C
from R6_common import footprint_labels

M5 = "/mnt/micron/geo_spatialtrans"; D = f"{M5}/data_real"
OUT, KJ = sys.argv[1], sys.argv[2]; ONLY = sys.argv[3:]
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A3": "Visium_HD_Human_Colon_Cancer",
      "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
KEY = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
SRC = {"A1": "cellpose_rgb", "A2": "cellpose_rgb", "A3": "cellpose_rgb", "A4": "consensus", "A5": "cellpose_rgb", "A6": "consensus"}
D_EXP = 8.0
NREP = 5

def roi_data(Rw):
    A, roi, upp = Rw.archetype, Rw.roi_id, Rw.um_per_px; tag = f"{A}_{roi}"
    mask = np.load(f"{M5}/{'R3b' if SRC[A]=='consensus' else 'R2'}/masks/{tag}_{SRC[A]}_native.npz")["labels"]
    nuc = pd.read_parquet(f"{M5}/R3b/nuclei/{tag}_consensus_native.parquet" if SRC[A] == "consensus" else f"results/R2/nuclei/{tag}_{SRC[A]}_native.parquet")
    keep = nuc.loc[nuc.keep, "label"].to_numpy()
    sq = f"{D}/{DS[A]}/binned_outputs_x/binned_outputs/square_002um"
    tp = pq.read_table(f"{sq}/spatial/tissue_positions.parquet").to_pandas()
    box = tp[(tp.pxl_row_in_fullres >= Rw.r0) & (tp.pxl_row_in_fullres < Rw.r1) & (tp.pxl_col_in_fullres >= Rw.c0) & (tp.pxl_col_in_fullres < Rw.c1) & (tp.in_tissue == 1)]
    z = np.load(f"{OUT}/section/{KEY[A]}.npz")
    u = z["U"][box.array_row.to_numpy(), box.array_col.to_numpy()].astype(float)
    rr = np.floor(box.pxl_row_in_fullres.to_numpy() - Rw.r0).astype(int); cc = np.floor(box.pxl_col_in_fullres.to_numpy() - Rw.c0).astype(int)
    pos = np.c_[(box.pxl_row_in_fullres.to_numpy() - Rw.r0) * upp, (box.pxl_col_in_fullres.to_numpy() - Rw.c0) * upp]
    s = int(round(20.0 / upp)); ms = np.roll(np.roll(mask, s, 0), s, 1)
    L = {"nuc": footprint_labels(mask, rr, cc, upp, D_EXP, "nuc"), "exp": footprint_labels(mask, rr, cc, upp, D_EXP, "exp"),
         "nuc_shift": footprint_labels(ms, rr, cc, upp, D_EXP, "nuc"), "exp_shift": footprint_labels(ms, rr, cc, upp, D_EXP, "exp")}
    def far_of(mm, ln):      # come R4 S-2: nel nucleo distanza dal bordo dentro, fuori distanza dal nucleo
        d_in = ndimage.distance_transform_edt(mm > 0) * upp; d_out = ndimage.distance_transform_edt(mm == 0) * upp
        ins = (rr >= 0) & (rr < mm.shape[0]) & (cc >= 0) & (cc < mm.shape[1])
        f = np.zeros(len(rr), bool); f[ins] = np.where(ln[ins] > 0, d_in[rr[ins], cc[ins]] >= 2.0, d_out[rr[ins], cc[ins]] >= 2.0)
        return f
    far = {"far": far_of(mask, L["nuc"]), "far_shift": far_of(ms, L["nuc_shift"])}
    rows, cols = box.array_row.to_numpy(), box.array_col.to_numpy()
    cls = C.parity_class(rows, cols); pf = C.parity_factors(u, cls, np.ones(len(u), bool)); c = pf[cls]
    return {"A": A, "roi": roi, "tag": tag, "upp": upp, "mask": mask, "keep": keep, "u": u, "pos": pos, "L": L, "far": far, "c": c, "pf": pf,
            "M": z["M"], "rr": rr, "cc": cc}

def g3c(Dd):
    """Totali per cellula (nucleo, corona; maschera vera e spostata) contro R4 (run1/posthoc)."""
    t = pd.read_parquet(f"{M5}/R4/run1/posthoc/{Dd['tag']}_cells.parquet")
    ok = np.array_equal(t.label.to_numpy(), Dd["keep"])
    for mk in ("", "_shift"):
        pc = C.per_cell_density(Dd["L"]["nuc" + mk], Dd["L"]["exp" + mk], Dd["keep"], Dd["u"])
        ok &= np.array_equal(pc["u_nuc"], t["umi_nuc" + mk].to_numpy()) and np.array_equal(pc["n_nuc"], t["n_nuc" + mk].to_numpy())
        ok &= np.array_equal(pc["u_ring"], t["umi_ring" + mk].to_numpy()) and np.array_equal(pc["n_ring"], t["n_ring" + mk].to_numpy())
    return bool(ok)

def d_image(Dd, kern):
    """Scarto della sezione (µm, righe e colonne della griglia) nel riferimento dell'immagine della ROI (µm)."""
    d_arr = np.array([kern["dr"], kern["dc"]]) / C.BIN_UM
    return (Dd["M"] @ d_arr) * Dd["upp"]

_CACHE = {}
def prep(tag, Rw):
    """Dati della ROI e quantita' che non dipendono dal nucleo di diffusione (in cache nel processo)."""
    if tag in _CACHE: return _CACHE[tag]
    Dd = roi_data(Rw)
    IX, IXf = C.stat_index(Dd["L"], Dd["keep"]), C.stat_index(Dd["L"], Dd["keep"], Dd["far"])
    nuc_g = C.labels_to_grid(Dd["mask"], Dd["upp"]); exp_g = C.expand_labels_grid(nuc_g, D_EXP)
    Dobs = C.observed_cell_density(Dd["L"], Dd["keep"], Dd["u"])
    am = Dd["L"]["exp"] == 0
    acell = float(Dd["u"][am].mean()) if am.any() else 0.0
    S = C.nc_setup(nuc_g, exp_g, Dobs, float(np.median(list(Dobs.values()))), acell)
    U = np.random.default_rng([C.SEED, int(Rw.name)]).uniform(size=(NREP, len(Dd["u"])))
    base = {"tag": tag, "archetype": Dd["A"], "roi_id": Dd["roi"], "n_bins": int(len(Dd["u"])), "n_cells": int(len(Dd["keep"])),
            "obs_star": C.star_fast(Dd["u"], IX), "obs_star_far": C.star_fast(Dd["u"], IXf), "G3c_identical_to_R4": g3c(Dd), "pf": Dd["pf"].tolist()}
    _CACHE.clear(); _CACHE[tag] = (Dd, IX, IXf, Dobs, am, acell, S, U, base)
    return _CACHE[tag]

def run_job(args):
    tag, Rw, name, kern = args
    fout = f"{OUT}/nc/{tag}__{name}.json"
    if os.path.exists(fout):
        return fout
    t0 = time.time()
    try:
        Dd, IX, IXf, Dobs, am, acell, S, U, base = prep(tag, Rw)
        kt = (kern["sigma_r"], kern["sigma_c"], kern["w"], kern["ell"]); di = d_image(Dd, kern)
        r, Sr = C.invert_r_refined(base["obs_star"], S, Dd["pos"], kt, di, Dd["c"], U, IX, Dd["L"], Dd["keep"], Dobs, acell, am)
        pred_far = float("nan")
        if np.isfinite(r):
            Y = C.nc_simulate(Sr, r, Dd["pos"], kt, di, Dd["c"], U)
            pred_far = float(np.mean([C.star_fast(y, IXf) for y in Y]))
        res = {**base, "kernel": name, "r_true": r, "pred_star_far": pred_far, "d_img": di.tolist(),
               **{k: kern[k] for k in ("sigma_r", "sigma_c", "w", "ell", "dr", "dc")}, "sec": round(time.time() - t0, 1)}
    except Exception:
        res = {"tag": tag, "kernel": name, "error": traceback.format_exc()}
    tmp = fout + ".tmp"; json.dump(res, open(tmp, "w")); os.replace(tmp, fout)
    return fout

if __name__ == "__main__":
    os.makedirs(f"{OUT}/nc", exist_ok=True)
    K = json.load(open(KJ))
    rois = pd.read_csv("results/R2/R2_rois_checked.csv")
    jobs = []
    for idx, Rw in rois.iterrows():
        tag = f"{Rw.archetype}_{Rw.roi_id}"
        if ONLY and tag not in ONLY: continue
        k = K[KEY[Rw.archetype]]
        for name, kern in [("K2", k["K2"]), ("K1", k["K1"])] + [(f"draw{j:02d}", d) for j, d in enumerate(k["draws"])]:
            jobs.append((tag, Rw, name, kern))          # in ordine di ROI: ogni processo riusa la ROI in cache
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 20
    todo = [j for j in jobs if not os.path.exists(f"{OUT}/nc/{j[0]}__{j[2]}.json")]
    print(f"jobs {len(jobs)} todo {len(todo)} workers {nw}", flush=True)
    t0 = time.time()
    with Pool(min(nw, max(1, len(todo)))) as pool:
        for i, f in enumerate(pool.imap(run_job, todo, chunksize=11), 1):
            if i % 22 == 0 or i == len(todo): print(f"{i}/{len(todo)} {time.time()-t0:.0f}s", flush=True)
    print("done", flush=True)
