#!/usr/bin/env python
"""tools/R3b_nuc.py — R3b: momenti dei pixel dei nuclei reali (eccentricita', orientazione) e contenimento per i generatori di
R3b, copia di tools/R3_containment.py::run_roi con i percorsi di R3b (maschera CO in R3b/masks, StarDist in R2/masks).
Uso (root del repo): .venv/bin/python tools/R3b_nuc.py   -> /mnt/micron/geo_spatialtrans/R3b/R3/real/*_nuc.parquet
"""
import importlib.util, glob, os
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r3c", ROOT / "tools/R3_containment.py"); R3C = importlib.util.module_from_spec(spec); spec.loader.exec_module(R3C)
OUT = Path("/mnt/micron/geo_spatialtrans/R3b"); MASKS = Path("/mnt/micron/geo_spatialtrans/R2/masks")


def run_roi(A, roi, method):
    o = OUT / "R3/real" / f"{A}_{roi}_{method}_nuc.parquet"
    if o.exists():
        return "skip"
    gen = pd.read_parquet(OUT / "R3/real" / f"{A}_{roi}_{method}_gen.parquet")
    mf = OUT / "masks" / f"{A}_{roi}_consensus_native.npz" if method == "consensus" else MASKS / f"{A}_{roi}_{method}_native.npz"
    labels = np.load(mf)["labels"]
    upp = float(gen["um_per_px"].iloc[0]); lab = gen["label"].to_numpy().astype(np.int64)
    lab2gen = np.full(int(labels.max()) + 1, -1, np.int64); lab2gen[lab] = np.arange(len(lab))
    pm = R3C.pixel_moments(labels, upp, lab)
    frac, tot = R3C.containment(labels, upp, gen[["x", "y"]].to_numpy(), lab2gen)
    hed = np.load(MASKS / f"{A}_{roi}_cellpose_hed_native.npz")["labels"]
    ri = np.clip(np.round(gen["y"].to_numpy() / upp).astype(int), 0, hed.shape[0] - 1)
    ci = np.clip(np.round(gen["x"].to_numpy() / upp).astype(int), 0, hed.shape[1] - 1)
    hl = hed[ri, ci]; ha = np.bincount(hed.ravel(), minlength=int(hed.max()) + 1) * upp ** 2
    out = pd.DataFrame({"idx": gen["idx"].to_numpy(), "label": lab}).merge(pm, on="label", how="left")
    out["frac_out"] = frac; out["n_px_own"] = tot
    out["hed_label"] = hl; out["hed_area_um2"] = np.where(hl > 0, ha[hl], np.nan)
    out.to_parquet(str(o) + ".tmp", index=False); os.replace(str(o) + ".tmp", o)
    return "ok"


if __name__ == "__main__":
    fs = sorted(glob.glob(str(OUT / "R3/real/*_gen.parquet")))
    for f in fs:
        A, roi, *m = os.path.basename(f).replace("_gen.parquet", "").split("_")
        print("[nuc]", A, roi, "_".join(m), run_roi(A, roi, "_".join(m)), flush=True)
    print(f"[nuc] {len(fs)} file")
