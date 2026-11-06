"""tools/S2.3_test_proto.py — S2.3, G-1 (e): lambda del prototipo di S2.0 (make_blur, regola C2) e trasformata di distanza sulla
maschera passata da R/testing/test_S2.3.R. Uso: .venv/bin/python tools/S2.3_test_proto.py <dir> (mass.bin, mask.bin, dist.bin, pos.bin,
meta.json) -> <dir>/lam.bin; stampa DIST_IDENTICAL True|False."""
import sys, json, os
os.environ["S20_F64"] = "1"
import numpy as np
from scipy import ndimage
sys.path.insert(0, "tools")
import S2_0_common as S0
d = sys.argv[1]; m = json.load(open(f"{d}/meta.json")); nr, nc = m["nr"], m["nc"]
mass = np.fromfile(f"{d}/mass.bin").reshape(nr, nc); mask = np.fromfile(f"{d}/mask.bin").reshape(nr, nc)
distR = np.fromfile(f"{d}/dist.bin").reshape(nr, nc); pos = np.fromfile(f"{d}/pos.bin").reshape(-1, 2)
dP = ndimage.distance_transform_edt(mask > 0.5) * 0.5
okd = np.allclose(np.where(mask > 0.5, dP, 0), np.where(mask > 0.5, distR, 0), rtol=0, atol=1e-12)
print("DIST_IDENTICAL", bool(okd), float(np.max(np.abs(np.where(mask > 0.5, dP - distR, 0)))))
# make_blur restituisce la densita' per sottopixel integrata sul quadratino divisa per l'area del sottopixel: si passa la massa per sottopixel
# e si moltiplica per (2 / 0.5)^2 = 16 sottopixel per quadratino; la riflessione ai bordi della griglia riguarda solo vetro (massa 0)
prm = {"sigma_e": m["sigma_e"], "D": m["D"], "w": m["w"], "ell": m["ell"]}
p = 80
massP = np.pad(mass, p); maskP = np.pad(mask, p)
blur = S0.make_blur("C2", prm, pos + p * 0.5, np.zeros(2), maskP, pad_um=20.0)
lam = blur(massP) * 16.0
lam.astype(np.float64).tofile(f"{d}/lam.bin")
print("LAM_SUM", float(lam.sum()), float(mass.sum()))
