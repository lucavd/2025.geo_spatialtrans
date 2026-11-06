"""tools/S2.4_R3_boot.py — S2.4, incertezza dei parametri di R-3 (dalla revisione avversariale, rilievo del revisore del vetro).
Bootstrap delle tessere di 300 um del vetro (B = 200, seme 20261115; B ridotto rispetto alle 2000 della pre-registrazione per il costo,
dichiarato): per ogni replica, (b, w2) adattati esattamente (tools/S2.4_R3_glass.fit_bw) su 17 valori di ell2 log-spaziati fra 0.4 e 2.5
volte l'ottimo, e scelto il migliore; intervalli percentili 2.5-97.5 %, frazione di repliche al bordo della griglia.
Uso: .venv/bin/python tools/S2.4_R3_boot.py <outdir> <R3dir> -> results/S2.4/S2.4_R3_boot.csv"""
import sys, os, json
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
import importlib.util
_s = importlib.util.spec_from_file_location("r3", "tools/S2.4_R3_glass.py"); R3 = importlib.util.module_from_spec(_s); _s.loader.exec_module(R3)
from scipy import ndimage, fft as sfft
GL, C = R3.GL, R3.C; B = 200; SEED = 20261115; KEYS = R3.KEYS
def run(KEY, OUT, R3D):
    fo = f"{OUT}/boot_{KEY}.json"
    if os.path.exists(fo): return fo
    fit = json.load(open(f"{R3D}/R3_{KEY}.json"))
    z = np.load(f"{GL.R5RUN}/section/{KEY}.npz"); U = z["U"].astype(float); IT = z["in_tissue"].astype(bool)
    prm = json.load(open(f"{GL.S20}/params_flat.json"))[KEY]["C2"]; bfar = fit["b_far"]; P = R3.PADB
    Up = np.pad(U, P); ITp = np.pad(IT, P); shape = Up.shape
    inner = np.zeros(shape, bool); inner[P + GL.BORDER:P + U.shape[0] - GL.BORDER, P + GL.BORDER:P + U.shape[1] - GL.BORDER] = True
    din = ndimage.distance_transform_edt(ITp) * GL.BIN; pi = np.where(ITp, np.exp(-din / prm["D"]), 0.0); del din
    dg = ndimage.distance_transform_edt(~ITp) * GL.BIN
    K = GL.kft(shape); kr, kc, box = K
    T = np.where(ITp, np.maximum(Up - bfar, 0.0), 0.0); Pm, Tf = GL.model(T, pi, prm, shape, K); del pi
    fm = (~ITp) & inner & (dg >= R3.DMIN); y = Up[fm]; core = Pm[fm]
    tiles = (((np.arange(shape[0])[:, None] - P) // R3.TILE) * 100000 + ((np.arange(shape[1])[None, :] - P) // R3.TILE))[fm]
    ut, inv = np.unique(tiles, return_inverse=True); order = np.argsort(inv); starts = np.searchsorted(inv[order], np.arange(len(ut)))
    ends = np.r_[starts[1:], len(order)]
    ells = fit["ell2"] * np.exp(np.linspace(np.log(0.4), np.log(2.5), 17))
    E2 = [sfft.irfft2(Tf * C.exp_ft(kr, kc, l) * box, s=shape, workers=1)[fm] for l in ells]
    rng = np.random.default_rng([SEED, 77, KEYS.index(KEY)]); out = []
    for r in range(B):
        pick = rng.integers(0, len(ut), len(ut)); cnt = np.bincount(pick, minlength=len(ut))
        wt = np.zeros(len(y)); 
        for t in np.flatnonzero(cnt): wt[order[starts[t]:ends[t]]] = cnt[t]
        sel = wt > 0; ys, cs, ww = y[sel], core[sel], wt[sel]
        best = None
        for k, l in enumerate(ells):
            e2 = E2[k][sel]; dl = e2 - cs
            def f(th):
                b, w = th; lam = np.maximum(b + cs + w * dl, 1e-12); rr = ww * (1.0 - ys / lam)
                return float(np.sum(ww * (lam - ys * np.log(lam)))), np.array([np.sum(rr), np.sum(rr * dl)])
            from scipy import optimize
            res = optimize.minimize(f, np.array([fit["b"] + 1e-3, max(fit["w2"], 1e-3)]), jac=True, method="L-BFGS-B", bounds=[(0.0, 50.0), (0.0, 0.5)], options={"ftol": 1e-12, "maxiter": 300})
            if best is None or res.fun < best[0]: best = (res.fun, k, res.x)
        out.append({"rep": r, "ell2": float(ells[best[1]]), "b": float(best[2][0]), "w2": float(best[2][1]), "edge": best[1] in (0, len(ells) - 1)})
    D = pd.DataFrame(out); res = {"key": KEY, "B": B, "edge_frac": float(D.edge.mean())}
    for q in ("b", "w2", "ell2"): res[f"{q}_hat"] = fit[q]; res[f"{q}_lo"] = float(D[q].quantile(0.025)); res[f"{q}_hi"] = float(D[q].quantile(0.975))
    json.dump(res, open(fo, "w")); return fo
if __name__ == "__main__":
    OUT, R3D = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
    with Pool(5) as p: p.starmap(run, [(k, OUT, R3D) for k in KEYS], chunksize=1)
    D = pd.DataFrame([json.load(open(f"{OUT}/boot_{k}.json")) for k in KEYS]); D.to_csv("results/S2.4/S2.4_R3_boot.csv", index=False)
    print(D.round(4).to_string(index=False))
