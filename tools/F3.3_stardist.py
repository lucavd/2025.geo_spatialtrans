"""tools/F3.3_stardist.py - F3.3: StarDist 2D_versatile_he sulle aree 'consensus' (linfonodo, cuore) con i parametri di
tools/R2_segment_stardist.py (normalizzazione 1-99.8, predict_instances_big a blocchi di 2048, seme 0). Legge il ritaglio scritto da
tools/F3.3_segment.py. Uso: tools/py_stardist.sh tools/F3.3_stardist.py <area,area,...|IN|holdout>"""
import sys, json, time, numpy as np, pandas as pd, tifffile
from pathlib import Path
sys.path.insert(0, "tools"); import R2_common as C
import tensorflow as tf
from stardist.models import StarDist2D
from csbdeep.utils import normalize
import os; OUT = Path(os.environ.get("F33_ROOT", "/mnt/micron/geo_spatialtrans/F3.3")) / "seg"
A = pd.read_csv("results/F3.3/F3.3_areas.csv"); A = A[A.seg == "consensus"]; sel = sys.argv[1]
A = A[A.level == "IN"] if sel == "IN" else A[A.level != "IN"] if sel == "holdout" else A[A.area.isin(sel.split(","))]
if os.environ.get("F33_SD_DET") == "1": tf.config.experimental.enable_op_determinism()   # G-3: senza, le aree di 500 um (via senza tasselli) non tornano identiche
tf.random.set_seed(0); np.random.seed(0); model = StarDist2D.from_pretrained("2D_versatile_he")
for _, r in A.iterrows():
    o = OUT / r.area
    if (o / "stardist_he_native.npz").exists(): print(r.area, "presente", flush=True); continue
    t0 = time.time(); img = tifffile.imread(o / "roi.tif"); x = normalize(img, 1, 99.8, axis=(0, 1))
    if max(img.shape[:2]) > 2048: lab, _ = model.predict_instances_big(x, axes="YXC", block_size=2048, min_overlap=128, context=128, n_tiles=(2, 2, 1))
    else: lab, _ = model.predict_instances(x)
    np.savez_compressed(o / "stardist_he_native.npz", labels=lab.astype(np.uint32))
    print(r.area, int(lab.max()), round(time.time() - t0, 1), flush=True)
