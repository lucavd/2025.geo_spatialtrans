"""tools/S2.0_edges.py — S2.0, M2 (BL-197): regole d'uscita C1, C2, C3 sulle 314 finestre di bordo di R5 (run3). Pre-registrazione
docs/steps/S2.0.md (afef06b, a8690b3, f2a23d7). Stessa regione, stessi bin, stessa verosimiglianza di R5 (R5_fit.py: prepare_region,
likelihood_mask sulla tessera centrale, bin validi). Partenze: S.starts_for (K2 di R5 e PH-5 di R5 della finestra).
Uso: .venv/bin/python tools/S2.0_edges.py <outdir> real [N_MAX]      -> <outdir>/edges/<KEY>_<tile>.json   (processi dal file <outdir>/WORKERS)
     .venv/bin/python tools/S2.0_edges.py <outdir> e1                -> <outdir>/e1/<KEY>_<tile>__<gen>.json (E1: dati sintetici sulla geometria vera)"""
import sys, os, json, glob, time, traceback
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R5_common as C
import S2_0_common as S

R5RUN = "/mnt/micron/geo_spatialtrans/R5/run3"
TB, MB = 150, 75
TILE = (slice(MB, MB + TB), slice(MB, MB + TB))
KEYS = ("A1", "A23", "A4", "A5", "A6")

def load_tile(path):
    z = np.load(path); y = z["U"].astype(float); f = z["f_sub"].astype(float)
    R0, C0 = int(z["R0"]), int(z["C0"]); rows, cols = np.indices(y.shape); rows = rows + R0; cols = cols + C0
    valid = z["inarr_bin"] & z["inimg_bin"] & ~z["bubble_bin"]
    return y, f, rows, cols, valid

def setup(y, f, rows, cols, valid):
    P = C.prepare_region(y, f, rows, cols)
    sel = C.likelihood_mask(P["occ"], P["dist"], TILE) & valid
    return P, S.region_extra(P, f), sel

def r5_fits(key, tile):
    d = json.load(open(f"{R5RUN}/tiles/{key}/{tile}.json"))
    ph = f"{R5RUN}/ph5_s0.7/{key}_{tile}.json"
    p5 = json.load(open(ph))["fit"] if os.path.exists(ph) else None
    return d, p5

def one_real(path):
    key = path.split("/")[-2]; tile = os.path.basename(path)[:-4]; fo = f"{OUT}/edges/{key}_{tile}.json"
    if os.path.exists(fo): return fo
    if os.environ.get("S20_IGNORE_LOCKS") != "1":      # prenotazione: con piu' esecutori nessuna finestra si calcola due volte
        try:
            os.close(os.open(fo + ".lock", os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        except FileExistsError:
            return fo
    t0 = time.time()
    try:
        d, p5 = r5_fits(key, tile)
        if "error" in d: raise RuntimeError("R5 tile error")
        y, f, rows, cols, valid = load_tile(path)
        P, X, sel = setup(y, f, rows, cols, valid)
        k2 = d["K2"]
        a, b, nll0 = S.nll_at("C0", {"sigma": k2["sigma"], "w": k2["w"], "ell": k2["ell"], "dr": k2["dr"], "dc": k2["dc"]}, P, X, y, sel)
        res = {"key": key, "tile": tile, "n_bins": int(sel.sum()), "C0": {**{k: k2[k] for k in ("sigma", "w", "ell", "dr", "dc", "a", "b", "nll")},
               "aic": 2 * S.K_AIC["C0"] + 2 * k2["nll"], "nll_recomputed": nll0}, "PH5": p5}
        for rule in ("C1", "C2", "C3"):
            t1 = time.time()
            res[rule] = S.fit_rule(rule, P, X, y, sel, S.starts_for(rule, k2, p5)); res[rule]["sec"] = round(time.time() - t1, 1)
        res["sec"] = round(time.time() - t0, 1)
    except Exception:
        res = {"key": key, "tile": tile, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo)
    return fo

# ------------------------------------------------------------------ E1: si possono distinguere le regole? (dichiarata prima del calcolo vero)
E1_TILES = 2       # finestre per sezione (le prime in ordine di nome fra quelle con almeno 3000 bin nella verosimiglianza)
def e1_params(gen, k2):
    base = {"w": float(np.clip(k2["w"], 0.03, 0.2)), "ell": float(np.clip(k2["ell"], 10, 40)), "dr": float(k2["dr"]), "dc": float(k2["dc"])}
    s = float(np.clip(k2["sigma"], 2.5, 5.5))
    return {"C1": dict(base, q=0.6, sigma_m=s), "C2": dict(base, sigma_e=s, D=5.0), "C3": dict(base, sigma_g=s)}[gen]

def one_e1(args):
    path, gen = args
    key = path.split("/")[-2]; tile = os.path.basename(path)[:-4]; fo = f"{OUT}/e1/{key}_{tile}__{gen}.json"
    if os.path.exists(fo): return fo
    try:
        d, p5 = r5_fits(key, tile)
        y, f, rows, cols, valid = load_tile(path)
        P, X, sel = setup(y, f, rows, cols, valid)
        k2 = d["K2"]; prm = e1_params(gen, k2)
        lam = P["c"] * (k2["a"] * np.maximum(S.g_rule(gen, prm, P, X), 0) + k2["b"])
        seed = int(abs(hash_str(f"{key}_{tile}_{gen}")) % (2 ** 31))
        ys = np.random.default_rng([S.SEED, seed]).poisson(lam).astype(float)
        Ps, Xs, sels = setup(ys, f, rows, cols, valid)
        k2s, _ = C.fit_k2_screened(Ps["Tft"], Ps["kr"], Ps["kc"], Ps["pshape"], Ps["region_bins"], ys, Ps["c"], sels, None)
        res = {"key": key, "tile": tile, "gen": gen, "true": prm, "C0": dict(k2s, aic=2 * S.K_AIC["C0"] + 2 * k2s["nll"])}
        for rule in ("C1", "C2", "C3"):
            res[rule] = S.fit_rule(rule, Ps, Xs, ys, sels, S.starts_for(rule, k2s, None))
    except Exception:
        res = {"key": key, "tile": tile, "gen": gen, "error": traceback.format_exc()}
    json.dump(res, open(fo + ".tmp", "w")); os.replace(fo + ".tmp", fo)
    return fo

def hash_str(s):
    """Seme deterministico da una stringa (md5; non hash() di Python, BL-202)."""
    import hashlib
    return int(hashlib.md5(s.encode()).hexdigest()[:8], 16)

def tiles_all():
    return sorted(p for k in KEYS for p in glob.glob(f"{R5RUN}/tiles/{k}/t*.npz") if os.path.exists(p.replace(".npz", ".json")))

if __name__ == "__main__":
    OUT, MODE = sys.argv[1], sys.argv[2]
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 8
    if MODE == "real":
        os.makedirs(f"{OUT}/edges", exist_ok=True)
        jobs = tiles_all()
        if len(sys.argv) > 3:          # E2: prova dei tempi su un sottoinsieme (2 per sezione)
            n = int(sys.argv[3]); jobs = [p for k in KEYS for p in [q for q in jobs if f"/{k}/" in q][:n]]
        jobs = [j for j in jobs if not os.path.exists(f"{OUT}/edges/{j.split('/')[-2]}_{os.path.basename(j)[:-4]}.json")]
        if os.environ.get("S20_REVERSE") == "1": jobs = jobs[::-1]          # secondo esecutore dalla fine (stesso risultato per finestra)
        if os.environ.get("S20_MID") == "1": jobs = jobs[len(jobs) // 2:] + jobs[:len(jobs) // 2]   # terzo esecutore dalla meta'

        if os.environ.get("S20_NW"): nw = int(os.environ["S20_NW"])
        print(f"real jobs {len(jobs)} workers {nw}", flush=True); t0 = time.time()
        with Pool(nw) as pool:
            for i, f in enumerate(pool.imap_unordered(one_real, jobs, chunksize=1), 1):
                print(f"{i}/{len(jobs)} {time.time() - t0:.0f}s {os.path.basename(f)}", flush=True)
    elif MODE == "e1":
        os.makedirs(f"{OUT}/e1", exist_ok=True)
        jobs = []
        for k in KEYS:
            if k == "A6": continue
            cand = []
            for p in [q for q in tiles_all() if f"/{k}/" in q]:
                d = json.load(open(p.replace(".npz", ".json")))
                if "error" not in d and d["K2"]["n_bins"] >= 3000: cand.append(p)
                if len(cand) == E1_TILES: break
            jobs += [(p, g) for p in cand for g in ("C1", "C2", "C3")]
        print(f"e1 jobs {len(jobs)} workers {nw}", flush=True); t0 = time.time()
        with Pool(nw) as pool:
            for i, f in enumerate(pool.imap_unordered(one_e1, jobs, chunksize=1), 1):
                print(f"{i}/{len(jobs)} {time.time() - t0:.0f}s {os.path.basename(f)}", flush=True)
    print("done", flush=True)
