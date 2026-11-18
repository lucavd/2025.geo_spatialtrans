"""tools/F3.2_stardist.py - F3.2, cancello G-3 (e): StarDist 2D_versatile_he sulla tessera di prova; scrive solo il numero di nuclei.
Uso: bash tools/py_stardist.sh tools/F3.2_stardist.py <label> [<label> ...]"""
import sys, json, traceback
from pathlib import Path
import numpy as np, tifffile
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/F3.2/compat"; OUT.mkdir(parents=True, exist_ok=True)
from stardist.models import StarDist2D
from csbdeep.utils import normalize
m = StarDist2D.from_pretrained("2D_versatile_he")
for lab in sys.argv[1:]:
    rec = dict(label=lab)
    try:
        img = tifffile.imread(f"/mnt/micron/geo_spatialtrans/F3.2_sealed/stardist_tiles/{lab}.tif")
        lbl, _ = m.predict_instances(normalize(img, 1, 99.8, axis=(0, 1)))
        rec["n_nuclei"] = int(lbl.max()); rec["tile_shape"] = list(img.shape)
    except Exception as e:
        rec["error"] = repr(e); rec["trace"] = traceback.format_exc()[-800:]
    json.dump(rec, open(OUT / f"{lab}_stardist.json", "w")); print(rec.get("n_nuclei"), rec.get("error", ""), lab, flush=True)
