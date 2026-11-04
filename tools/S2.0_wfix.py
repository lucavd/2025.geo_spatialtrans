"""tools/S2.0_wfix.py — S2.0, post hoc PH-6 (dalla revisione G-5 b, dubbio D1 del revisore dei bordi): le soluzioni salvate con la massa
della coda w = 0 possono essere minimi di bordo (tutte le partenze ereditano w dal K2 di R5 della finestra; con w = 0 la lunghezza
della coda non conta e l'ottimizzatore resta fermo). Per ogni finestra e regola con w < 1e-6: nuovo adattamento con la soluzione
salvata come partenza tenuta e 9 partenze a coda accesa (w in {0.01, 0.03, 0.06} x ell in {50, 100, 200} µm, gli altri parametri
salvati); stesso ottimizzatore (S.fit_rule). Uscite: <run>/edges_wfix/<sezione>_<tessera>.json; poi <run>_wfix/edges/ con le finestre
unite (la regola si sostituisce solo se la NLL scende di oltre 1e-6). Uso: .venv/bin/python tools/S2.0_wfix.py <run> <nproc>"""
import sys, os, json, glob, time, traceback, importlib.util
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "tools")
import S2_0_common as S
sp = importlib.util.spec_from_file_location("ed", "tools/S2.0_edges.py"); ED = importlib.util.module_from_spec(sp)
sys_argv = sys.argv; sys.argv = ["x", sys_argv[1], "none"]; sp.loader.exec_module(ED); sys.argv = sys_argv
RUN = sys.argv[1]; NP = int(sys.argv[2]) if len(sys.argv) > 2 else 8
R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
PK = {"C0": ("sigma",), "C1": ("q", "sigma_m"), "C2": ("sigma_e", "D"), "C3": ("sigma_g",)}
def job(f):
    d = json.load(open(f)); key, tile = d["key"], d["tile"]; fo = f"{RUN}/edges_wfix/{key}_{tile}.json"
    if os.path.exists(fo): return fo
    t0 = time.time(); out = {"key": key, "tile": tile}
    try:
        y, fsub, rows, cols, valid = ED.load_tile(f"{R5RUN}/tiles/{key}/{tile}.npz")
        P, X, sel = ED.setup(y, fsub, rows, cols, valid)
        for rl in ("C0", "C1", "C2", "C3"):
            s = d[rl]
            if s["w"] >= 1e-6: continue
            base = {k: float(s[k]) for k in PK[rl] + ("w", "ell", "dr", "dc")}
            starts = [dict(base, keep=True)] + [dict(base, w=w, ell=l) for w in (0.01, 0.03, 0.06) for l in (50.0, 100.0, 200.0)]
            fit = S.fit_rule(rl, P, X, y, sel, starts)
            out[rl] = {"saved_nll": float(s["nll"]), "new": fit, "improvement": float(s["nll"] - fit["nll"])}
    except Exception:
        out["error"] = traceback.format_exc()
    out["sec"] = round(time.time() - t0, 1)
    json.dump(out, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo); return fo
def merge():
    os.makedirs(f"{RUN}_wfix/edges", exist_ok=True); n = 0
    for f in sorted(glob.glob(f"{RUN}/edges/*.json")):
        d = json.load(open(f)); g = f"{RUN}/edges_wfix/{os.path.basename(f)}"
        if os.path.exists(g) and "error" not in d:
            w = json.load(open(g))
            for rl in ("C0", "C1", "C2", "C3"):
                if rl in w and w[rl]["improvement"] > 1e-6:
                    keep = {k: v for k, v in d[rl].items() if k in ("nll_recomputed",)}
                    d[rl] = {**w[rl]["new"], **keep, "wfix_improvement": w[rl]["improvement"]}; n += 1
        json.dump(d, open(f"{RUN}_wfix/edges/{os.path.basename(f)}", "w"))
    print("merged; rules replaced:", n)
if __name__ == "__main__":
    os.makedirs(f"{RUN}/edges_wfix", exist_ok=True)
    todo = []
    for f in sorted(glob.glob(f"{RUN}/edges/*.json")):
        d = json.load(open(f))
        if "error" not in d and any(d[r]["w"] < 1e-6 for r in ("C0", "C1", "C2", "C3")): todo.append(f)
    print("windows", len(todo), flush=True)
    with Pool(NP) as p:
        for i, _ in enumerate(p.imap_unordered(job, todo)): print(i + 1, flush=True)
    merge()
