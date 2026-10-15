"""review/exact/f5_clipper_versions.py — punto (a): stessi interi (troncati come polyclip / arrotondati come dice la sua documentazione)
dati a Clipper 6.4.0 (sorgente di polyclip 1.10-7), 6.4.0 + sola guardia AddLocalMinPoly del 6.4.2, Clipper 6.4.2 (sorgente di pyclipper 1.4.0, C++ diretto),
pyclipper 1.4.0 (Python) e Clipper2 2.0.1 (pyclipr 0.1.8). Area in unità intere esatta (shoelace intero) -> µm² = area * 1e-12 (Fraction).
Confronto con l'area esatta (truth, verificata da rx in review_viol.py). Uscita: review_exact/f5/f5_versions.csv"""
import subprocess, sys, os
from fractions import Fraction as F
import pandas as pd, numpy as np
O = "/mnt/micron/geo_spatialtrans/S1.3d/review_exact/f5"; B = "/mnt/micron/geo_spatialtrans/S1.3d/review_exact/build"
G = pd.read_csv(f"{O}/f5_polyclip.csv"); I = pd.read_csv(f"{O}/f5_ints.csv", dtype={"xt": str, "yt": str, "xr": str, "yr": str})
I = I.sort_values(["case", "poly", "ring", "k"])
def paths(case, poly, v):
    d = I[(I.case == case) & (I.poly == poly)]
    return [[(int(a), int(b)) for a, b in zip(g["x" + v], g["y" + v])] for _, g in d.groupby("ring")]
P = {(cs, v): (paths(cs, "A", v), paths(cs, "B", v)) for cs in G.case for v in ("t", "r")}
txt = []
for (cs, v), (A, Bp) in P.items():
    txt.append(f"C {cs} {v} {len(A)} {len(Bp)}")
    for p in A + Bp: txt.append(f"R {len(p)}"); txt += [f"{x} {y}" for x, y in p]
inp = "\n".join(txt) + "\n"
res = {}
for drv in ("drv640", "drv640g", "drv642"):
    out = subprocess.run([f"{B}/{drv}"], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
    for ln in out:
        if ln.strip():
            cs, v, a2, npth, nv = ln.split(); res[(int(cs), v, drv)] = F(int(a2), 2) * F(1, 10**12)
import pyclipper, pyclipr
def area_int(paths_):
    s = F(0)
    for p in paths_:
        n = len(p)
        for k in range(n): s += F(int(p[k][0]) * int(p[(k + 1) % n][1]) - int(p[(k + 1) % n][0]) * int(p[k][1]), 2)
    return s * F(1, 10**12)
for (cs, v), (A, Bp) in P.items():
    pc = pyclipper.Pyclipper(); pc.AddPaths(A, pyclipper.PT_SUBJECT, True); pc.AddPaths(Bp, pyclipper.PT_CLIP, True)
    res[(cs, v, "pyclipper642")] = area_int(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
    c2 = pyclipr.Clipper(); c2.scaleFactor = 1
    for p in A: c2.addPath(np.array(p, dtype=np.float64), pyclipr.Subject, False)
    for p in Bp: c2.addPath(np.array(p, dtype=np.float64), pyclipr.Clip, False)
    sol = c2.execute(pyclipr.Intersection, pyclipr.FillRule.NonZero)
    res[(cs, v, "clipper2")] = area_int([[(round(float(x)), round(float(y))) for x, y in s_] for s_ in sol])
rows = []
for _, g in G.iterrows():
    tr = F(float(g.truth)); r = dict(case=g.case, key=g.key, i=g.i, ix=g.ix, iy=g.iy, s=g.s, truth=g.truth, cand_reported=g.cand_reported, pc_default=g.pc_default, nvert=g.nvert)
    for v in ("t", "r"):
        for drv in ("drv640", "drv640g", "drv642", "pyclipper642", "clipper2"):
            a = res[(g.case, v, drv)]; r[f"{drv}_{v}"] = float(a); r[f"{drv}_{v}_err"] = float(abs(a - tr))
    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(f"{O}/f5_versions.csv", index=False)
# riassunto: errore / tau con tau = g * perimetro del pezzo esatto (dal viol); classi: >10 tau = topologico
V = pd.read_csv("/mnt/micron/geo_spatialtrans/S1.3d/review_exact/f5_viol_tau.csv") if os.path.exists("/mnt/micron/geo_spatialtrans/S1.3d/review_exact/f5_viol_tau.csv") else None
print("drv640_t riproduce polyclip (|d| <= 1e-9 rel):", int((abs(R.drv640_t - R.pc_default) <= 1e-9 * R.pc_default.abs() + 1e-15).sum()), "/", len(R))
for col in ["drv640", "drv640g", "drv642", "pyclipper642", "clipper2"]:
    for v in ("t", "r"):
        e = R[f"{col}_{v}_err"]; print(f"{col:13s} {v}: errore max {e.max():.3e}  n(err > 1e-3 µm²) = {(e > 1e-3).sum():3d}  n(err > 1e-5) = {(e > 1e-5).sum():3d}")
