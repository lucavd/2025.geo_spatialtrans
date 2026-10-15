"""review/exact/f5_broad.py — punto (a), controprova estesa: TUTTI i pezzi esatti a 2 µm di un caso, interi troncati (come polyclip) o
arrotondati, dati a Clipper 6.4.0 (polyclip), 6.4.0+guardia, 6.4.2. Conta |errore| > tau (tau = g * perimetro esatto del pezzo) e > 1e-3 µm².
Limite dichiarato: si testano solo i pezzi della verità (i pezzi spuri, area esatta 0, non sono provati qui).
Uso: python f5_broad.py <key>   -> review_exact/f5/broad_<key>.csv (riassunto) e broad_<key>_bad.csv (pezzi oltre tau)"""
import sys, os, json, subprocess
from fractions import Fraction as F
import numpy as np, pandas as pd, pyarrow.parquet as pq
D = "/mnt/micron/geo_spatialtrans/S1.3d"; O = f"{D}/review_exact/f5"; B = f"{D}/review_exact/build"
key = sys.argv[1]; s = 2; g = 1e-6
gs = json.load(open(f"{D}/geom/{key}_meta.json"))["gs"]
T = pq.read_table(f"{D}/truth/{key}_o2_s2.parquet").to_pandas()
t = pq.read_table(f"{D}/geom/{key}_cells.parquet").to_pandas().sort_values(["cell", "part", "ring", "k"])
def sarea(x, y): return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
rings = {}
for (ci, pa, ri), d in t.groupby(["cell", "part", "ring"], sort=True):
    x, y = d.x.values, d.y.values
    if len(x) > 1 and x[0] == x[-1] and y[0] == y[-1]: x, y = x[:-1], y[:-1]
    a = sarea(x, y)
    if (ri == 1 and a < 0) or (ri > 1 and a > 0): x, y = x[::-1], y[::-1]
    rings.setdefault(int(ci), []).append((x, y))
def ints(v, mode):   # come interface.cpp: (cInt)((x - x0)/eps) con x0 = 0 -> troncamento verso zero del quoziente double
    q = v / 1e-6
    return np.trunc(q).astype(np.int64) if mode == "t" else np.round(q).astype(np.int64)
res = {}
for mode in ("t", "r"):
    lines = []
    for q, (ci, ix, iy) in enumerate(zip(T.i.values, T.ix.values, T.iy.values)):
        R = rings[int(ci)]; x0 = gs["X0"] + s * int(ix); y0 = gs["Y0"] + s * int(iy)
        bx = ints(np.array([x0, x0 + s, x0 + s, x0], float), mode); by = ints(np.array([y0, y0, y0 + s, y0 + s], float), mode)
        lines.append(f"C {q} {mode} {len(R)} 1")
        for x, y in R:
            X, Y = ints(x, mode), ints(y, mode); lines.append(f"R {len(X)}"); lines += [f"{a} {b}" for a, b in zip(X, Y)]
        lines.append("R 4"); lines += [f"{a} {b}" for a, b in zip(bx, by)]
    inp = "\n".join(lines) + "\n"
    for drv in ("drv640", "drv640g", "drv642"):
        out = subprocess.run([f"{B}/{drv}"], input=inp, capture_output=True, text=True, check=True).stdout.split()
        a2 = np.array(out[2::5], dtype=object)
        res[(drv, mode)] = a2
rows = []; bad = []
tau = g * T.perim.values; tr = T.area.values
for (drv, mode), a2 in res.items():
    err = np.array([float(abs(F(int(v), 2 * 10**12) - F(float(w)))) for v, w in zip(a2, tr)])
    v = err > tau
    rows.append(dict(key=key, drv=drv, mode=mode, n_pieces=len(T), n_over_tau=int(v.sum()), n_over_1e3=int((err > 1e-3).sum()), max_err=float(err.max()), max_ratio=float((err / tau).max())))
    for q in np.nonzero(v)[0]: bad.append(dict(key=key, drv=drv, mode=mode, i=int(T.i.values[q]), ix=int(T.ix.values[q]), iy=int(T.iy.values[q]), truth=float(tr[q]), clip=float(F(int(a2[q]), 2 * 10**12)), err=float(err[q]), tau=float(tau[q])))
pd.DataFrame(rows).to_csv(f"{O}/broad_{key}.csv", index=False); pd.DataFrame(bad).to_csv(f"{O}/broad_{key}_bad.csv", index=False)
print(pd.DataFrame(rows).to_string())
