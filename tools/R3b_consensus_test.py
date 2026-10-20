#!/usr/bin/env python
"""tools/R3b_consensus_test.py — R3b, cancello G-1 (parte Python): test del costruttore della maschera CO su etichette sintetiche.
Stampa una riga PASS/FAIL per asserzione e l'ultima riga «G-1 python: k/n PASS». Mutanti: R3B_MUTANT=M1..M5 (cancello G-2).
Uso (root del repo): .venv/bin/python tools/R3b_consensus_test.py
"""
import importlib.util, os, sys
import numpy as np
from pathlib import Path
from skimage import measure
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("r3b", ROOT / "tools/R3b_consensus.py"); B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
spec = importlib.util.spec_from_file_location("r3c", ROOT / "tools/R3_containment.py"); R3C = importlib.util.module_from_spec(spec); spec.loader.exec_module(R3C)
RES = []


def ck(name, ok, val=""):
    RES.append(bool(ok)); print(f"{'PASS' if ok else 'FAIL'} {name} {val}")


H = 120; yy, xx = np.mgrid[0:H, 0:H]
def disk(cx, cy, r): return (xx - cx) ** 2 + (yy - cy) ** 2 <= r * r
def rect(x0, y0, x1, y1): return (xx >= x0) & (xx < x1) & (yy >= y0) & (yy < y1)


def paint(*items):
    lab = np.zeros((H, H), np.int64)
    for l, m in items: lab[m] = l
    return lab


# T1 uno a uno: forma = StarDist
sr = paint((5, disk(30, 30, 6))); sd = paint((9, disk(31, 30, 5)))
co, tab = B.build_consensus(sr, sd)
ck("T1 un nucleo", len(tab) == 1 and tab.src.iloc[0] == "sd")
ck("T1 pixel = StarDist", np.array_equal(co > 0, sd > 0))
ck("T1 matched_pixels_identical = 0", B.matched_pixels_identical(co, sd, tab) == 0)

# T2 due StarDist dentro uno Space Ranger: appaiato il piu' sovrapposto, l'altro (solo StarDist) escluso
sr = paint((1, rect(10, 10, 50, 30))); sd = paint((3, rect(12, 12, 35, 28)), (4, rect(40, 14, 48, 26)))
co, tab = B.build_consensus(sr, sd)
ck("T2 un solo nucleo CO", len(tab) == 1 and int(tab.sd_label.iloc[0]) == 3, f"n={len(tab)}")
ck("T2 nessun pixel dello StarDist non appaiato", not np.any(co[sd == 4] > 0))
ck("T2 nessuna origine sd_only", "sd_only" not in set(tab.src))

# T3 uno StarDist su due Space Ranger: appaiato al piu' sovrapposto; l'altro resta con i pixel non coperti
sr = paint((1, rect(10, 10, 30, 30)), (2, rect(30, 10, 50, 30))); sd = paint((7, rect(14, 12, 36, 28)))
co, tab = B.build_consensus(sr, sd)
r1 = tab[tab.sr_label == 1]; r2 = tab[tab.sr_label == 2]
ck("T3 SR 1 appaiato a SD 7", len(r1) == 1 and r1.src.iloc[0] == "sd" and int(r1.sd_label.iloc[0]) == 7)
ck("T3 SR 2 solo Space Ranger", len(r2) == 1 and r2.src.iloc[0] == "sr")
l2 = int(r2.label.iloc[0]) if len(r2) else -1
ck("T3 SR 2 = pixel SR non coperti da SD", np.array_equal(co == l2, (sr == 2) & (sd != 7)))
ck("T3 vince StarDist sulla sovrapposizione", np.all(co[(sr == 2) & (sd == 7)] == int(r1.label.iloc[0]) if len(r1) else False))

# T4 Space Ranger coperto del tutto dallo StarDist appaiato a un altro: tolto
sr = paint((1, rect(10, 10, 40, 40)), (2, rect(42, 20, 46, 24))); sd = paint((5, rect(12, 12, 50, 38)))
co, tab = B.build_consensus(sr, sd)
ck("T4 SR coperto tolto", len(tab) == 1 and int(tab.sr_label.iloc[0]) == 1, f"n={len(tab)}")
ck("T4 conteggio etichette = presenti", len(np.unique(co[co > 0])) == len(tab))
ck("T4 attributo dei tolti", tab.attrs.get("n_sr_only_dropped") == 1)

# T5 Space Ranger senza StarDist: tenuto con i suoi pixel
sr = paint((3, disk(20, 20, 5)), (8, disk(80, 80, 6))); sd = paint((1, disk(20, 21, 4)))
co, tab = B.build_consensus(sr, sd)
r8 = tab[tab.sr_label == 8]
ck("T5 SR senza StarDist tenuto", len(r8) == 1 and r8.src.iloc[0] == "sr")
ck("T5 pixel SR", len(r8) == 1 and np.array_equal(co == int(r8.label.iloc[0]), sr == 8))

# T6 identificativi 1..n in ordine di etichetta SR
sr = paint((9, disk(20, 20, 5)), (2, disk(60, 20, 5)), (5, disk(20, 60, 5)), (7, disk(60, 60, 5)))
sd = paint((1, disk(60, 21, 4)), (2, disk(20, 61, 4)))
co, tab = B.build_consensus(sr, sd)
ck("T6 etichette 1..n", np.array_equal(np.unique(co[co > 0]), np.arange(1, 5)) and np.array_equal(tab.label, np.arange(1, 5)))
ck("T6 ordine SR", list(tab.sr_label) == [2, 5, 7, 9])

# T7 conservazione dei pixel: unione (SR solo + SD appaiati)
sr = paint((1, rect(10, 10, 30, 30)), (2, rect(30, 10, 50, 30)), (3, disk(90, 90, 8))); sd = paint((7, rect(14, 12, 36, 28)), (8, disk(70, 30, 5)))
co, tab = B.build_consensus(sr, sd)
exp = (sd == 7) | ((sr == 2) & (sd != 7)) | (sr == 3)
ck("T7 pixel CO = SD appaiati U SR solo non coperti", np.array_equal(co > 0, exp))
ck("T7 n_px della tabella", all(int(t.n_px) == int(np.sum(co == t.label)) for t in tab.itertuples()))

# T8 pareggio: due StarDist con la stessa sovrapposizione -> identificativo piu' piccolo
sr = paint((1, rect(10, 10, 30, 20))); sd = paint((4, rect(10, 10, 20, 20)), (3, rect(20, 10, 30, 20)))
co, tab = B.build_consensus(sr, sd)
ck("T8 pareggio all'identificativo piu' piccolo", len(tab) == 1 and int(tab.sd_label.iloc[0]) == 3)

# T9 reciprocita': SR 1 preferisce SD 5, ma SD 5 preferisce SR 2 -> SR 1 resta solo Space Ranger
sr = paint((1, rect(10, 10, 20, 30)), (2, rect(20, 10, 40, 30))); sd = paint((5, rect(14, 12, 40, 28)))
co, tab = B.build_consensus(sr, sd)
ck("T9 appaiamento reciproco", list(tab.src) == ["sr", "sd"] and int(tab[tab.src == "sd"].sr_label.iloc[0]) == 2)

# T10 filtro del tessuto valido (cp2): etichetta con centroide fuori tolta
lab = paint((1, disk(20, 20, 4)), (2, disk(100, 100, 4)))
valid = np.zeros((30, 30), bool); valid[:15, :15] = True                   # meta' superiore sinistra valida
f = B.cp2_filter(lab, 1.0, valid, float(H))
ck("T10 dentro tenuta, fuori tolta", np.any(f == 1) and not np.any(f == 2))

# T11 tabella dei nuclei: momenti dei pixel (R3) contro regionprops su un'ellisse CO
upp = 0.25; a, b_, phi = 18.0, 8.0, np.deg2rad(30)
u = (xx - 60) * np.cos(phi) + (yy - 60) * np.sin(phi); v = -(xx - 60) * np.sin(phi) + (yy - 60) * np.cos(phi)
sr = paint((4, (u / a) ** 2 + (v / b_) ** 2 <= 1)); sd = paint((6, (u / (a - 1)) ** 2 + (v / (b_ - 1)) ** 2 <= 1))
co, tab = B.build_consensus(sr, sd)
pm = R3C.pixel_moments(co, upp, [1]).iloc[0]; rp = measure.regionprops(co.astype(np.int32))[0]
ck("T11 area", abs(pm.area_px_um2 - rp.area * upp ** 2) < 1e-12)
ck("T11 centroide", abs(pm.cx_px - rp.centroid[1] * upp) < 1e-9 and abs(pm.cy_px - rp.centroid[0] * upp) < 1e-9)
ck("T11 eccentricita'", abs(pm.ecc_N - rp.eccentricity) < 0.02, f"{pm.ecc_N:.4f} vs {rp.eccentricity:.4f}")
dth = abs((pm.theta_N - (np.pi / 2 - rp.orientation) + np.pi / 2) % np.pi - np.pi / 2)
ck("T11 orientazione", np.rad2deg(dth) < 2, f"{np.rad2deg(dth):.3f} deg")

print(f"G-1 python [{os.environ.get('R3B_MUTANT', 'none')}]: {sum(RES)}/{len(RES)} PASS")
sys.exit(0 if all(RES) else 1)
