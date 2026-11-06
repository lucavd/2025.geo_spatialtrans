"""tools/S2.4_R2_slow.py — S2.4, regola R-2 e controprova CP-3 (docs/steps/S2.4.md): parte lenta delle strisce stimata sul vetro.
Per sezione: lambda del vetro con il modello di R-3 (b, w2, ell2 da <R3dir>/R3_<KEY>.json; stessa costruzione di tools/S2.4_R3_glass.py),
quadratini di vetro a >= DMIN um dal tessuto, esclusi i 50 al bordo dell'area di cattura. y_ij ~ Poisson(exp(a_i + g_j) lambda_ij):
stima esatta alternata (20 iterazioni). Parita' tolta come S2.0 (periodic_split, indice della griglia di Space Ranger); parte lenta s =
media mobile su 21 righe (colonne) del residuo (media dei soli valori finiti, almeno 10 su 21), come il complemento di PH-5.
Meta': quadratini di vetro di ogni riga divisi alla mediana delle colonne (per le colonne: delle righe); a stimata nelle due meta' con g
fisso; riga valida se ogni meta' ha >= 200 quadratini e >= 50 UMI. Criteri (i) corr(s1, s2) >= 0.5, (ii) sqrt(cov(s1, s2)) >= 0.05,
(iii) >= 300 righe con entrambe, (iv) corr(s a >= 16 um, s a >= 100 um) >= 0.5. W del generatore da corr a 10 righe (1 - 10/W).
CP-3: conteggi del vetro sostituiti da Poisson(lambda exp(alpha_i + beta_j)), parte fine come il generatore della sezione (R-1,
capture), piu' (positivo) parte lenta = media mobile su 41 righe di normali, SD 0.08, o (negativo) nessuna; seme 20261114.
Uso: .venv/bin/python tools/S2.4_R2_slow.py <outdir> <R3dir>  -> <outdir>/R2_<KEY>.json; results/S2.4/S2.4_R2_slow.csv, S2.4_CP3.csv"""
import sys, os, json, time
import numpy as np, pandas as pd
from scipy import ndimage, fft as sfft
sys.path.insert(0, "tools")
import importlib.util
_sp = importlib.util.spec_from_file_location("s23glass", "tools/S2.3_glass.py"); GL = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(GL)
C = GL.C
import S2_0_common as S0
KEYS = ["A1", "A23", "A4", "A5", "A6"]; PADB = 1500; W = 21; MINB = 200; MINU = 50; NIT = 20; SEED = 20261114
GEN = pd.read_csv("results/S2.4/S2.4_R1_generator.csv")
SEC = pd.read_csv("results/S2.0/S2.0_PH5_sections.csv"); SEC = SEC[SEC.n_iter == 10]

def lam_model(KEY, R3):
    z = np.load(f"{GL.R5RUN}/section/{KEY}.npz"); U = z["U"].astype(float); IT = z["in_tissue"].astype(bool)
    prm = json.load(open(f"{GL.S20}/params_flat.json"))[KEY]["C2"]
    bfar = float(pd.read_csv("results/R5/R5_kernel_descriptive.csv").set_index("key").loc[KEY, "b_far"])
    Up = np.pad(U, PADB); ITp = np.pad(IT, PADB); shape = Up.shape
    din = ndimage.distance_transform_edt(ITp) * GL.BIN; pi = np.where(ITp, np.exp(-din / prm["D"]), 0.0); del din
    K = GL.kft(shape); kr, kc, box = K
    T = np.where(ITp, np.maximum(Up - bfar, 0.0), 0.0)
    Pm, Tf = GL.model(T, pi, prm, shape, K)
    E2 = sfft.irfft2(Tf * C.exp_ft(kr, kc, R3["ell2"]) * box, s=shape, workers=1)
    lam = R3["b"] + (1 - R3["w2"]) * Pm + R3["w2"] * E2
    sl = (slice(PADB, PADB + U.shape[0]), slice(PADB, PADB + U.shape[1]))
    dg = ndimage.distance_transform_edt(~IT) * GL.BIN
    inner = np.zeros(U.shape, bool); inner[GL.BORDER:-GL.BORDER, GL.BORDER:-GL.BORDER] = True
    return U, lam[sl], dg, (~IT) & inner

def fit_ag(Y, L, M):
    a = np.zeros(Y.shape[0]); g = np.zeros(Y.shape[1]); Ym = np.where(M, Y, 0.0); Lm = np.where(M, L, 0.0)
    for _ in range(NIT):
        num = Ym.sum(1); den = (Lm * np.exp(g)[None, :]).sum(1); a = np.where((num > 0) & (den > 0), np.log(np.maximum(num, 1e-300) / np.maximum(den, 1e-300)), np.nan)
        ea = np.nan_to_num(np.exp(a)); num = Ym.sum(0); den = (Lm * ea[:, None]).sum(0); g = np.where((num > 0) & (den > 0), np.log(np.maximum(num, 1e-300) / np.maximum(den, 1e-300)), np.nan)
        g = np.nan_to_num(g - np.nanmean(g))
    return a, g

def half_a(Y, L, M, g):
    """a per riga nelle due meta' (colonne del vetro divise alla mediana della riga), g fisso."""
    nr, nc = Y.shape; jj = np.arange(nc)[None, :]
    cnt = M.sum(1); csum = np.cumsum(M, 1); half = M & (csum <= (cnt[:, None] + 1) // 2)
    out = []
    for H in (half, M & ~half):
        num = np.where(H, Y, 0).sum(1); den = (np.where(H, L, 0) * np.exp(g)[None, :]).sum(1)
        ok = (H.sum(1) >= MINB) & (num >= MINU) & (den > 0)
        out.append(np.where(ok, np.log(np.maximum(num, 1e-300) / np.maximum(den, 1e-300)), np.nan))
    return out

def slow(f):
    f = np.asarray(f, float); ok = np.isfinite(f); k = np.ones(W)
    idx = np.arange(len(f)); mu, r = S0.periodic_split(f, idx)
    num = np.convolve(np.where(ok, r, 0.0), k, "same"); den = np.convolve(ok.astype(float), k, "same")
    with np.errstate(invalid="ignore", divide="ignore"): m = num / den
    return np.where(ok & (den >= W // 2), m, np.nan), np.where(ok & (den >= W // 2), r - m, np.nan), mu

def cc(x, y, lag=0):
    a, b = (x[:-lag], y[lag:]) if lag else (x, y); g = np.isfinite(a) & np.isfinite(b)
    if g.sum() < 10: return np.nan, np.nan, int(g.sum())
    cv = np.cov(a[g], b[g]); return float(cv[0, 1] / np.sqrt(cv[0, 0] * cv[1, 1])), float(cv[0, 1]), int(g.sum())

def analyse(Y, L, M, M100, axis):
    if axis == "cols": Y, L, M, M100 = Y.T, L.T, M.T, M100.T
    a, g = fit_ag(Y, L, M); h1, h2 = half_a(Y, L, M, g)
    s1, f1, _ = slow(h1); s2, f2, _ = slow(h2); s, fpart, mu = slow(a)
    r, cv, n = cc(s1, s2); sd = float(np.sqrt(cv)) if np.isfinite(cv) and cv > 0 else 0.0
    a100, _ = fit_ag(Y, L, M100); s100, _, _ = slow(a100); r100, _, n100 = cc(s, s100)
    c10 = 0.5 * (cc(s1, s2, 10)[1] + cc(s2, s1, 10)[1]); rho10 = c10 / cv if (np.isfinite(cv) and cv > 0) else np.nan
    Wg = 10.0 / (1.0 - rho10) if np.isfinite(rho10) and 0 < rho10 < 1 else np.nan
    rf, _, _ = cc(f1, f2)
    ident = bool((r >= 0.5) and (sd >= 0.05) and (n >= 300) and (r100 >= 0.5))
    return {"corr_halves_slow": r, "sd_repro_slow": sd, "n_both": n, "corr_slow_16_100": r100, "n_16_100": n100, "rho10": rho10, "W_generator": Wg,
            "corr_halves_fine": rf, "parity_glass": float(np.exp(mu[1] - mu[0])), "identified": ident}, fpart

def run(KEY, OUT, R3D):
    fo = f"{OUT}/R2_{KEY}.json"
    if os.path.exists(fo): return fo
    t0 = time.time()
    R3 = json.load(open(f"{R3D}/R3_{KEY}.json"))
    U, L, dg, G0 = lam_model(KEY, R3)
    M16 = G0 & (dg >= 16); M100 = G0 & (dg >= 100)
    res = {"key": KEY, "n_glass16": int(M16.sum()), "n_glass100": int(M100.sum()), "umi_glass16": float(U[M16].sum())}
    Fz = np.load(f"{GL.S20}/stripes/{KEY}_factors.npz")
    for ax in ("rows", "cols"):
        o, fpart = analyse(U, L, M16, M100, ax)
        tf = np.asarray(Fz["alpha" if ax == "rows" else "beta"], float); _, tfine, _ = slow(tf)
        o["corr_fine_glass_tissue"] = cc(fpart, tfine)[0]
        res[ax] = o
    # CP-3: conteggi sintetici sul vetro con fattori noti
    rng = np.random.default_rng([SEED, 30, KEYS.index(KEY)]); res["cp3"] = {}
    for lab, sd_slow in (("positive", 0.08), ("negative", 0.0)):
        fac = {}
        for ax, n in (("rows", U.shape[0]), ("cols", U.shape[1])):
            gr = GEN[(GEN.key == KEY) & (GEN.axis == ax)].iloc[0]; sdf = float(SEC[(SEC.key == KEY) & (SEC.axis == ax)].iloc[0].sd_repro_hp)
            rho = float(gr.rho_generator); th = (1 - np.sqrt(1 - 4 * rho ** 2)) / (2 * rho) if rho != 0 else 0.0
            e = rng.normal(size=n + 1); x = sdf * (e[1:] + th * e[:-1]) / np.sqrt(1 + th ** 2)
            if sd_slow > 0:
                e2 = rng.normal(size=n + 40); csum = np.concatenate([[0], np.cumsum(e2)]); x = x + sd_slow * np.sqrt(41) * (csum[41:] - csum[:-41]) / 41
            fac[ax] = x
        Ys = np.where(G0, rng.poisson(L * np.exp(fac["rows"][:, None] + fac["cols"][None, :])), 0).astype(float)
        res["cp3"][lab] = {ax: analyse(Ys, L, M16, M100, ax)[0] for ax in ("rows", "cols")}
    res["sec"] = round(time.time() - t0, 1)
    json.dump(res, open(fo + ".tmp", "w"), default=float); os.replace(fo + ".tmp", fo); return fo

def summary(OUT):
    rows, cp = [], []
    for k in KEYS:
        f = f"{OUT}/R2_{k}.json"
        if not os.path.exists(f): continue
        d = json.load(open(f))
        for ax in ("rows", "cols"): rows.append({"key": k, "axis": ax, **d[ax]})
        for lab, v in d["cp3"].items():
            for ax in ("rows", "cols"): cp.append({"key": k, "control": lab, "axis": ax, **v[ax]})
    A = pd.DataFrame(rows); B = pd.DataFrame(cp)
    A.to_csv("results/S2.4/S2.4_R2_slow.csv", index=False); B.to_csv("results/S2.4/S2.4_CP3.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["key", "axis", "corr_halves_slow", "sd_repro_slow", "n_both", "corr_slow_16_100", "W_generator", "corr_halves_fine", "corr_fine_glass_tissue", "parity_glass", "identified"]
    print(A[cols].round(3).to_string(index=False))
    print(B[["key", "control", "axis", "corr_halves_slow", "sd_repro_slow", "n_both", "corr_slow_16_100", "identified"]].round(3).to_string(index=False))
    for ax in ("rows", "cols"): print(ax, "identificata in", int(A[A.axis == ax].identified.sum()), "sezioni su", int((A.axis == ax).sum()))

if __name__ == "__main__":
    OUT, R3D = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
    if len(sys.argv) > 3 and sys.argv[3] == "summary": summary(OUT); sys.exit(0)
    from multiprocessing import Pool
    nw = int(open(f"{OUT}/WORKERS").read().strip()) if os.path.exists(f"{OUT}/WORKERS") else 3
    with Pool(min(nw, len(KEYS))) as p: p.starmap(run, [(k, OUT, R3D) for k in KEYS], chunksize=1)
    summary(OUT)
