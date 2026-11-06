"""tools/S2.3_repro_compare.py — S2.3, cancello G-3: confronto dei contenuti fra la corsa principale e quella rifatta (tempi esclusi;
NaN uguali a NaN: il confronto diretto dei dizionari contava come diverse 3 aree del linfonodo con correlazioni NaN).
Uso: .venv/bin/python tools/S2.3_repro_compare.py -> results/S2.3/S2.3_G3.txt"""
import json, glob, pickle, hashlib, numpy as np
M = "/mnt/micron/geo_spatialtrans/S2.3"; R = f"{M}/repro"
SKIP = {"sec", "elapsed_function_s", "elapsed_total_s", "params_file"}
def eq(x, y):
    if isinstance(x, dict): return x.keys() == y.keys() and all(eq(x[k], y[k]) for k in x if k not in SKIP)
    if isinstance(x, list): return len(x) == len(y) and all(eq(a, b) for a, b in zip(x, y))
    if isinstance(x, np.ndarray): return np.array_equal(x, y, equal_nan=True)
    if isinstance(x, float): return (np.isnan(x) and np.isnan(y)) or x == y
    return x == y
out = []
def cnt(files, other, load):
    n = sum(eq(load(f), load(other(f))) for f in files); return n, len(files)
n, t = cnt(sorted(glob.glob(f"{R}/bench/grid/*.json")), lambda f: f.replace(R + "/bench", M + "/bench"), lambda f: json.load(open(f)))
out.append(f"banco, verifica: {n}/{t} file identici (tempi esclusi)")
n, t = cnt(sorted(glob.glob(f"{R}/bench/prep/*.pkl")), lambda f: f.replace(R + "/bench", M + "/bench"), lambda f: pickle.load(open(f, "rb")))
out.append(f"banco, preparazione: {n}/{t} aree identiche in tutti i campi (tempi esclusi, NaN uguali)")
fs = sorted(glob.glob(f"{R}/panel/main/*.bin"))
n = sum(hashlib.md5(open(f, "rb").read()).hexdigest() == hashlib.md5(open(f.replace(R + "/panel", M + "/panel"), "rb").read()).hexdigest() for f in fs)
out.append(f"pannello: {n}/{len(fs)} file dei conteggi con md5 identico")
n, t = cnt(sorted(glob.glob(f"{R}/panel/main/stats/*.json")), lambda f: f.replace(R + "/panel", M + "/panel"), lambda f: json.load(open(f)))
out.append(f"pannello, statistiche: {n}/{t} file identici (tempi esclusi)")
open("results/S2.3/S2.3_G3.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
