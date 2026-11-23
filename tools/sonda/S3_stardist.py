"""tools/sonda/S3_stardist.py - S3: StarDist 2D_versatile_he sui ritagli di tools/sonda/S3_export.py (crop), parametri di
tools/F3.3_stardist.py (normalizzazione 1 e 99,8; predict_instances_big a blocchi di 2048, min_overlap 128, context 128, n_tiles (2,2,1);
seme 0; determinismo delle operazioni attivo). Uso: tools/py_stardist.sh tools/sonda/S3_stardist.py <unit,unit,...|all>"""
import sys, os, glob, json, time, numpy as np, tifffile
import tensorflow as tf
from stardist.models import StarDist2D
from csbdeep.utils import normalize
ROOT = os.environ.get("S3_ROOT", "/mnt/micron/geo_spatialtrans/S3_sonda")
tf.config.experimental.enable_op_determinism(); tf.random.set_seed(0); np.random.seed(0)
model = StarDist2D.from_pretrained("2D_versatile_he")
sel = sys.argv[1]; dirs = sorted(glob.glob(f"{ROOT}/seg/*")) if sel == "all" else [f"{ROOT}/seg/{u}" for u in sel.split(",")]
for o in dirs:
    if os.path.exists(f"{o}/stardist.npz") or not os.path.exists(f"{o}/roi.tif"): continue
    t0 = time.time(); img = tifffile.imread(f"{o}/roi.tif"); x = normalize(img, 1, 99.8, axis=(0, 1))
    if max(img.shape[:2]) > 2048: lab, _ = model.predict_instances_big(x, axes="YXC", block_size=2048, min_overlap=128, context=128, n_tiles=(2, 2, 1))
    else: lab, _ = model.predict_instances(x)
    np.savez_compressed(f"{o}/stardist.npz", labels=lab.astype(np.uint32))
    json.dump({"n_nuclei": int(lab.max()), "sec": round(time.time() - t0, 1), "shape": list(lab.shape)}, open(f"{o}/stardist.json", "w"))
    print(os.path.basename(o), int(lab.max()), round(time.time() - t0, 1), flush=True)
