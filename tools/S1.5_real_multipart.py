# tools/S1.5_real_multipart.py - S1.5: quanti nuclei reali sono in piu' pezzi? (dato per la scelta su BL-139, prima di ogni simulazione)
# Per ogni ROI e maschera (primaria: cellpose_rgb nativa in A1, A2, A3, A5; maschera unica CO di R3b in A4, A6; secondaria: stardist_he nativa)
# conta le etichette i cui pixel formano piu' di una componente connessa (8-vicinato e 4-vicinato).
# Uscita: results/S1.5/S1.5_real_multipart.csv.  Uso: .venv/bin/python tools/S1.5_real_multipart.py
import csv, glob, os, numpy as np
from scipy import ndimage as ndi
M = "/mnt/micron/geo_spatialtrans"
rois = list(csv.DictReader(open("results/R2/R2_rois_checked.csv")))
def load(p):
    z = np.load(p); k = [k for k in z.files][0]; return z[k], z.files
rows = []
for r in rois:
    A, roi = r["archetype"], r["roi_id"]
    prim = f"{M}/R3b/masks/{A}_{roi}_consensus_native.npz" if A in ("A4", "A6") else f"{M}/R2/masks/{A}_{roi}_cellpose_rgb_native.npz"
    for meth, p in (("primary", prim), ("stardist_he", f"{M}/R2/masks/{A}_{roi}_stardist_he_native.npz")):
        lab, keys = load(p)
        sl = ndi.find_objects(lab)
        n = n8 = n4 = 0; frac_area_minor = []
        for i, s in enumerate(sl, start=1):
            if s is None: continue
            n += 1; b = lab[s] == i
            l8, k8 = ndi.label(b, structure=np.ones((3, 3), int)); l4, k4 = ndi.label(b)
            if k8 > 1:
                n8 += 1; sz = np.bincount(l8.ravel())[1:]; frac_area_minor.append(1 - sz.max() / sz.sum())
            if k4 > 1: n4 += 1
        rows.append(dict(archetype=A, roi_id=roi, mask=meth, file=os.path.basename(p), npz_key=keys[0], n_labels=n, n_multipart_8=n8, n_multipart_4=n4,
                         max_minor_area_frac_8=(max(frac_area_minor) if frac_area_minor else 0.0)))
os.makedirs("results/S1.5", exist_ok=True)
with open("results/S1.5/S1.5_real_multipart.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
import collections
t = collections.defaultdict(lambda: [0, 0, 0])
for x in rows:
    k = (x["archetype"], x["mask"]); t[k][0] += x["n_labels"]; t[k][1] += x["n_multipart_8"]; t[k][2] += x["n_multipart_4"]
for k, v in sorted(t.items()): print(k, "labels", v[0], "multipart8", v[1], "multipart4", v[2])
print("keys", sorted(set(x["npz_key"] for x in rows)))
