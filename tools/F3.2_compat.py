#!/usr/bin/env python3
"""tools/F3.2_compat.py - F3.2, cancello G-3: compatibilita' di una sezione con il flusso del progetto.
Controlli (file dello step): (a) matrici a 2 e 8 um leggibili e coerenti con le posizioni (stesso numero di quadratini, barcode uguali);
(b) scala da scalefactors_json in [0.2735, 0.2745] um/px; (c) passo della griglia a 2 um stimato dalle posizioni (minimi quadrati
pixel ~ array, metodo di R1) in [1.995, 2.005] um; (d) H&E a piena risoluzione leggibile, >= 99 % dei quadratini in tessuto dentro l'immagine;
(e) StarDist su una tessera di prova di 250 um (centrata sulla mediana dei quadratini in tessuto) senza errori e con almeno un nucleo.
Guasti costruiti (devono dare FAIL): scale10 (scala * 1.10), rot90 (posizioni ruotate di 90 gradi), trunc (immagine troncata a meta' altezza).
Uso: .venv/bin/python tools/F3.2_compat.py <label> [fault]  -> results/F3.2/compat/<label>[_<fault>].json
La tessera di StarDist si segmenta con tools/py_stardist.sh in un passo a parte (stardist <label>), conteggio dei nuclei soltanto.
"""
import sys, json, math, subprocess
import numpy as np, pandas as pd, h5py, tifffile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util
_s = importlib.util.spec_from_file_location('f32areas', str(Path(__file__).resolve().parent / 'F3.2_areas.py'))
A = importlib.util.module_from_spec(_s); _s.loader.exec_module(A)

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/F3.2/compat"; OUT.mkdir(parents=True, exist_ok=True)

def n_h5(sq):
    f = h5py.File(sq / "filtered_feature_bc_matrix.h5", "r")["matrix"]
    return int(f["shape"][1]), [b.decode() for b in f["barcodes"][:]]

def img_shape(p):
    with tifffile.TiffFile(p) as tf:
        s = tf.series[0]; sh = s.shape; ax = s.axes
        H = sh[ax.index("Y")]; W = sh[ax.index("X")]
        return int(H), int(W)

def pitch_um(pos, upp):
    """minimi quadrati pixel = a + B [row, col] sugli array; passo = media delle norme delle colonne di B, in um"""
    s = pos.sample(min(len(pos), 200000), random_state=1)
    Xd = np.c_[np.ones(len(s)), s.array_row.values, s.array_col.values]
    br, *_ = np.linalg.lstsq(Xd, s.pxl_row_in_fullres.values, rcond=None)
    bc, *_ = np.linalg.lstsq(Xd, s.pxl_col_in_fullres.values, rcond=None)
    v_row = np.array([br[1], bc[1]]); v_col = np.array([br[2], bc[2]])
    return float((np.linalg.norm(v_row) + np.linalg.norm(v_col)) / 2 * upp), float(np.linalg.norm(v_row) / np.linalg.norm(v_col))

def run(label, fault=None):
    P = A.paths(label); res = dict(label=label, fault=fault or "", dataset=P["ds"]); ok = {}
    try:
        n2, bc2 = n_h5(P["sq2"]); n8, bc8 = n_h5(P["sq8"])
        p2 = pd.read_parquet(P["sq2"] / "spatial/tissue_positions.parquet"); p8 = pd.read_parquet(P["sq8"] / "spatial/tissue_positions.parquet")
        it2 = set(p2.barcode[p2.in_tissue == 1]); it8 = set(p8.barcode[p8.in_tissue == 1])
        ok["a"] = (n2 == len(it2)) and (set(bc2) == it2) and (n8 == len(it8)) and (set(bc8) == it8)
        res.update(n_bins2_matrix=n2, n_bins2_in_tissue=len(it2), n_bins8_matrix=n8, n_bins8_in_tissue=len(it8))
    except Exception as e:
        ok["a"] = False; res["a_error"] = repr(e); p2 = None
    try:
        upp = A.upp_of(P) * (1.10 if fault == "scale10" else 1.0)
        res["um_per_px"] = upp; ok["b"] = 0.2735 <= upp <= 0.2745
    except Exception as e:
        ok["b"] = False; res["b_error"] = repr(e); upp = np.nan
    try:
        q = p2.copy()
        if fault == "rot90":
            q["pxl_row_in_fullres"], q["pxl_col_in_fullres"] = q.pxl_col_in_fullres.values.copy(), -q.pxl_row_in_fullres.values.copy()
        pt, an = pitch_um(q, upp); res.update(pitch_um=pt, anisotropy=an); ok["c"] = 1.995 <= pt <= 2.005
    except Exception as e:
        ok["c"] = False; res["c_error"] = repr(e); q = p2
    try:
        H, W = img_shape(P["img"])
        if fault == "trunc": H = H // 2
        it = q[q.in_tissue == 1]
        inside = ((it.pxl_col_in_fullres >= 0) & (it.pxl_col_in_fullres < W) & (it.pxl_row_in_fullres >= 0) & (it.pxl_row_in_fullres < H)).mean()
        res.update(img_H=H, img_W=W, frac_in_tissue_inside=float(inside))
        r0 = int(it.pxl_row_in_fullres.median()); c0 = int(it.pxl_col_in_fullres.median())
        reg = A.read_region(P["img"], r0, r0 + 64, c0, c0 + 64)
        readable = reg.size > 0 and (fault != "trunc" or r0 + 64 <= H)
        ok["d"] = bool(inside >= 0.99 and readable)
    except Exception as e:
        ok["d"] = False; res["d_error"] = repr(e)
    sd = ROOT / "results/F3.2/compat" / f"{label}_stardist.json"
    if fault is None and sd.exists():
        s = json.load(open(sd)); res.update(stardist_n_nuclei=s.get("n_nuclei"), stardist_error=s.get("error", ""))
        ok["e"] = (s.get("n_nuclei", 0) >= 1) and not s.get("error")
    else:
        ok["e"] = None if fault is None else True                     # i guasti costruiti non toccano StarDist
    res["checks"] = ok; res["PASS"] = all(v for v in ok.values() if v is not None) and ok["e"] is not None
    json.dump(res, open(OUT / (f"{label}_{fault}.json" if fault else f"{label}.json"), "w"), indent=1, default=str)
    print(label, fault or "", "PASS" if res["PASS"] else "FAIL", ok, flush=True)

def tile(label):
    """ritaglio della tessera di prova (250 um, centrata sulla mediana dei quadratini in tessuto) per StarDist"""
    P = A.paths(label); upp = A.upp_of(P)
    p2 = pd.read_parquet(P["sq2"] / "spatial/tissue_positions.parquet"); it = p2[p2.in_tissue == 1]
    h = int(round(125 / upp)); r0 = int(it.pxl_row_in_fullres.median()) - h; c0 = int(it.pxl_col_in_fullres.median()) - h
    reg = A.read_region(P["img"], r0, r0 + 2 * h, c0, c0 + 2 * h)
    if reg.ndim == 3 and reg.shape[0] in (3, 4) and reg.shape[-1] not in (3, 4): reg = np.moveaxis(reg, 0, -1)
    out = Path("/mnt/micron/geo_spatialtrans/F3.2_sealed/stardist_tiles"); out.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out / f"{label}.tif", np.ascontiguousarray(reg[..., :3]).astype(np.uint8)); print(out / f"{label}.tif")

if __name__ == "__main__":
    if sys.argv[1] == "tile": tile(sys.argv[2])
    else: run(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
