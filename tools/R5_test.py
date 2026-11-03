"""tools/R5_test.py — R5, cancello G-1: batteria di test su dati sintetici con parametri noti, piu' 5 mutanti.
Uso (root del repo): .venv/bin/python tools/R5_test.py [--workers N] [--mutant M1..M5] [--quick]
Stampa PASS/FAIL per asserzione; senza --mutant scrive results/R5/R5_test_output.txt e results/R5/R5_G1_recovery.csv.
Pre-registrazione: docs/steps/R5.md (c5c9fce, ceca77b)."""
import sys, os, argparse, time, json
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
import R5_common as C
from scipy import ndimage, signal

ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=20); ap.add_argument("--mutant", default="")
ap.add_argument("--quick", action="store_true"); ap.add_argument("--nrep", type=int, default=30); A = ap.parse_args()

# ---------------------------------------------------------------- mutanti (applicati al modulo prima dei test)
GEOM_CODE = {"curve": 1, "blob": 2, "lumen": 3, "straight": 4}
def synth_seed(geom, s_r, s_c, w, dr, dc, rep):
    return [C.SEED, int(s_r * 100), int(s_c * 100), int(w * 100), int((dr + 5) * 10), int((dc + 5) * 10), rep, GEOM_CODE[geom]]

def t12_seed_reproducible():
    # T12: la stessa tessera sintetica generata in due processi con PYTHONHASHSEED diversi e' identica
    import subprocess
    code = ("import sys, numpy as np; sys.path.insert(0, 'tools'); import R5_common as C; "
            "rng = np.random.default_rng([C.SEED, 200, 200, 5, 65, 40, 3, 2]); "
            "S = C.synth_tile(rng, 2.0, 2.0, 0.05, 30.0, 1.5, -1.0, geometry='blob'); print(float(S['y'].sum()), float(S['f_sub'].sum()))")
    outs = [subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED=h)).stdout.strip() for h in ("1", "2")]
    ok = synth_seed("blob", 2.0, 2.0, 0.05, 1.5, -1.0, 3) == [C.SEED, 200, 200, 5, 65, 40, 3, 2]
    check("T12 synthetic tile identical across processes (PYTHONHASHSEED 1 vs 2)", ok and outs[0] == outs[1] and outs[0] != "", outs[0] + " | " + outs[1])

def apply_mutant(m):
    if m == "M1":     # nucleo non normalizzato: manca (1 - w)
        C.kernel_ft = lambda kr, kc, sr, sc, w=0.0, ell=30.0: C.gauss_ft(kr, kc, sr, sc) + (w * C.exp_ft(kr, kc, ell) if w > 0 else 0)
    elif m == "M2":   # segno dello scarto invertito
        C.shift_ft = lambda kr, kc, dr, dc: np.exp(+2j * np.pi * (kr * dr + kc * dc))
    elif m == "M3":   # parita' ignorata
        C.parity_factors = lambda U, cls, sel: np.ones(4)
    elif m == "M4":   # distanze in bin invece che in µm
        _fg = C.freq_grid; C.freq_grid = lambda shape, d=C.SUB_UM: _fg(shape, d / C.BIN_UM)
    elif m == "X1":   # convenzione di campionamento sbagliata sulle sole righe (+0.25 µm; rilievo C10 della revisione)
        _eb = C.expected_bins
        def eb(Tft, kr, kc, pshape, region_bins, sr, sc, w, ell, dr, dc): return _eb(Tft, kr, kc, pshape, region_bins, sr, sc, w, ell, dr + 0.5, dc)
        C.expected_bins = eb
    elif m == "M5":   # coda ignorata nel modello in avanti del rapporto
        _ba = C.blur_at; C.blur_at = lambda F, pos, sr, sc, w, ell, d, **k: _ba(F, pos, sr, sc, 0.0, ell, d, **k)
if A.mutant:
    apply_mutant(A.mutant)

RES = []
def check(name, ok, detail=""):
    RES.append((name, bool(ok), detail)); print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)

TILE = (slice(75, 225), slice(75, 225))

# ---------------------------------------------------------------- T1 normalizzazione del nucleo
def t1():
    shape = (2400, 2400); kr, kc = C.freq_grid(shape)
    for p in ((1.0, 1.0, 0.0, 30.0), (2.0, 3.0, 0.05, 30.0), (0.5, 0.5, 0.2, 80.0)):
        k = np.fft.irfft2(C.kernel_ft(kr, kc, *p), s=shape)
        check(f"T1 kernel sums to 1 {p}", abs(k.sum() - 1.0) < 1e-6 and abs(C.kernel_ft(kr, kc, *p)[0, 0] - 1.0) < 1e-12, f"sum={k.sum():.8f}")

# ---------------------------------------------------------------- T2 semipiano sfocato contro la formula
def t2():
    from scipy.stats import norm
    nb = (60, 60); n = (nb[0] * C.SUB, nb[1] * C.SUB); X0 = 60.0        # bordo verticale a 60 µm (spigolo di bin)
    x = (np.arange(n[1]) + 0.5) * C.SUB_UM
    T = np.tile((x < X0).astype(float), (n[0], 1))
    T = np.pad(T, ((0, 0), (0, 0)), mode="constant")
    for s in (0.5, 2.0):
        ps = C.pad_shape(n, 0)                                        # periodico: si confronta lontano dai lati
        kr, kc = C.freq_grid(ps)
        g = C.expected_bins(C.source_ft(T, ps), kr, kc, ps, nb, s, s, 0.0, 30.0, 0.0, 0.0)
        xb = (np.arange(nb[1]) + 0.5) * C.BIN_UM
        q = np.linspace(-1, 1, 2001)
        ref = np.array([norm.cdf((X0 - (xc + q)) / s).mean() for xc in xb])
        sel = (xb > 40) & (xb < 80)
        err = np.max(np.abs(g[30, sel] - ref[sel]))
        tol = 5e-3 if s < 1 else 1e-3      # campionamento del gradino a 0.5 µm: 0.5 % del salto a sigma = 0.5 equivale a < 0.01 µm di posizione
        check(f"T2 half-plane blur matches Phi, sigma={s}", err < tol, f"max err={err:.2e} tol={tol:.0e}")

# ---------------------------------------------------------------- T3 convenzione dello scarto
def t3():
    rng = np.random.default_rng(5); nb = (40, 40); n = (160, 160)
    T = ndimage.gaussian_filter(rng.random(n), 3)
    ps = C.pad_shape(n, 80); kr, kc = C.freq_grid(ps)
    g_shift = C.expected_bins(C.source_ft(T, ps), kr, kc, ps, nb, 1.5, 1.5, 0.0, 30.0, 1.0, -0.5)
    Tr = np.zeros_like(T); Tr[2:, :] = T[:-2, :]; Tr = np.roll(Tr, -1, axis=1); Tr[:, -1] = 0   # +1 µm righe (2 sottopixel), -0.5 µm colonne
    g_ref = C.expected_bins(C.source_ft(Tr, ps), kr, kc, ps, nb, 1.5, 1.5, 0.0, 30.0, 0.0, 0.0)
    err = np.max(np.abs(g_shift[5:-5, 5:-5] - g_ref[5:-5, 5:-5]))
    check("T3 shift delta = source moved by +delta", err < 1e-6, f"max err={err:.2e}")

# ---------------------------------------------------------------- T4 a, b profilati
def t4():
    rng = np.random.default_rng(6); g = rng.gamma(2.0, 2.0, 5000); c = rng.uniform(0.5, 1.5, 5000)
    lam = c * (1.7 * g + 0.25)
    a, b, f = C.fit_ab(lam, g, c)
    check("T4 fit_ab recovers a, b on exact data", abs(a - 1.7) < 1e-4 and abs(b - 0.25) < 1e-4, f"a={a:.6f} b={b:.6f}")

# ---------------------------------------------------------------- T9 fattori di parita'
def t9():
    rng = np.random.default_rng(9)
    S = C.synth_tile(rng, 2.0, 2.0, 0.05, 30.0, 0.0, 0.0, geometry="curve")
    P = C.prepare_region(S["y"], S["f_sub"], S["rows"], S["cols"])
    true = np.array([0.6, 1.0, 0.9, 1.5]); true = true / true.mean()
    err = np.max(np.abs(P["pf"] / true - 1))
    check("T9 parity factors recovered within 3 %", err < 0.03, f"est={np.round(P['pf'], 3).tolist()} max rel err={err:.3f}")

# ---------------------------------------------------------------- T5, T5b, T6 recupero dei parametri
def one_fit(args):
    geom, s_r, s_c, w, dr, dc, rep, model = args
    # seme deterministico: la prima versione usava hash(geom), che cambia a ogni processo (PYTHONHASHSEED): tessere non riproducibili
    rng = np.random.default_rng(synth_seed(geom, s_r, s_c, w, dr, dc, rep))
    S = C.synth_tile(rng, s_r, s_c, w, 30.0, dr, dc, geometry=geom)
    t0 = time.time()
    P = C.prepare_region(S["y"], S["f_sub"], S["rows"], S["cols"])
    # stimatore pre-registrato: fondo b libero (default dopo la revisione); R5_BFIXED=1 = variante scartata (b dal vetro lontano)
    bfar = C.far_glass_background(S["y"], S["f_sub"], P["c"]) if os.environ.get("R5_BFIXED") == "1" else None
    fit, P, sel = C.fit_region(S["y"], S["f_sub"], S["rows"], S["cols"], TILE, model="K2", P=P, b_fixed=bfar)
    if model == "K3":
        st = dict(fit); st["sigma_r"] = st["sigma_c"] = fit["sigma"]
        fit = C.fit_model("K3", P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], S["y"], P["c"], sel, [st], bfar)
    return {"geometry": geom, "model": model, "sigma_r_true": s_r, "sigma_c_true": s_c, "w_true": w, "dr_true": dr, "dc_true": dc, "rep": rep,
            **{k: fit[k] for k in ("sigma", "sigma_r", "sigma_c", "w", "ell", "dr", "dc", "a", "b", "nll", "converged")}, "sec": time.time() - t0}

def t5(workers, nrep):
    jobs = []
    for s in (0.5, 1.0, 2.0, 3.0, 5.0):
        for w in (0.0, 0.05):
            for (dr, dc) in ((0.0, 0.0), (1.5, -1.0)):
                jobs += [("curve", s, s, w, dr, dc, k, "K2") for k in range(nrep)]
    jobs += [("blob", 2.0, 2.0, 0.05, 1.5, -1.0, k, "K2") for k in range(nrep)]
    jobs += [("blob", 2.0, 2.0, 0.05, 0.0, 0.0, k, "K3") for k in range(nrep)] + [("blob", 3.0, 2.0, 0.05, 0.0, 0.0, k, "K3") for k in range(nrep)]
    with Pool(workers) as pool:
        out = pd.DataFrame(pool.map(one_fit, jobs, chunksize=1))
    return out

def judge_recovery(D):
    g = D.groupby(["geometry", "model", "sigma_r_true", "sigma_c_true", "w_true", "dr_true", "dc_true"])
    for key, d in g:
        geom, model, sr, sc, w, dr, dc = key
        lab = f"{geom} {model} sigma=({sr},{sc}) w={w} delta=({dr},{dc}) n={len(d)}"
        # stima come sulla sezione: mediana delle tessere (docs/steps/R5.md, «Stima per sezione»)
        if model == "K2":
            es = d.sigma.median() - sr; tol = 0.15 + 0.05 * sr
            check(f"T5 sigma recovered [{lab}]", abs(es) <= tol, f"bias={es:+.3f} tol={tol:.3f} sd_tile={d.sigma.std():.3f}")
            ew = d.w.median() - w
            check(f"T5 w recovered [{lab}]", abs(ew) <= 0.01, f"bias={ew:+.4f} sd_tile={d.w.std():.4f}")
            for comp, tv in (("dr", dr), ("dc", dc)):
                e = d[comp].median() - tv
                check(f"T5 delta {comp} recovered [{lab}]", abs(e) <= 0.25, f"bias={e:+.3f} sd_tile={d[comp].std():.3f}")
        else:
            for comp, tv in (("sigma_r", sr), ("sigma_c", sc)):
                e = d[comp].median() - tv; tol = 0.15 + 0.05 * tv
                check(f"T6 K3 {comp} recovered [{lab}]", abs(e) <= tol, f"bias={e:+.3f} tol={tol:.3f} sd_tile={d[comp].std():.3f}")
            rr = (d.sigma_r / d.sigma_c).median()
            check(f"T6 K3 ratio sigma_r/sigma_c [{lab}]", abs(rr - sr / sc) <= 0.1 * sr / sc, f"ratio={rr:.3f} true={sr/sc:.3f}")

# ---------------------------------------------------------------- T7 inversione del rapporto (generatore indipendente, nello spazio reale)
def ref_blur_bins(F, pos_um, sigma, w, ell, d_img, step=C.SUB_UM):
    """Implementazione di riferimento: nucleo campionato nello spazio reale, convoluzione diretta (FFT di scipy.signal), media sui 4 x 4
    sottopixel del bin, scarto multiplo del passo."""
    R = int(np.ceil(max(4 * sigma, 6 * ell if w > 0 else 0) / step))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1] * step
    rr2 = yy ** 2 + xx ** 2
    g = np.exp(-rr2 / (2 * sigma ** 2)); g /= g.sum()
    K = (1 - w) * g
    if w > 0:
        e = np.exp(-np.sqrt(rr2) / ell); e /= e.sum(); K = K + w * e
    p = R
    G = np.pad(F, p, mode="reflect")
    H = signal.fftconvolve(G, K, mode="same")[p:-p, p:-p]
    sh = np.round(np.asarray(d_img) / step).astype(int)
    H = np.roll(H, tuple(sh), axis=(0, 1))
    n0, n1 = H.shape[0] // C.SUB, H.shape[1] // C.SUB
    Hb = H[: n0 * C.SUB, : n1 * C.SUB].reshape(n0, C.SUB, n1, C.SUB).mean(axis=(1, 3))
    bi = np.floor(np.asarray(pos_um) / C.BIN_UM).astype(int)
    return Hb[bi[:, 0], bi[:, 1]]

def t7(quick):
    rng = np.random.default_rng(77)
    Sx = C.synth_roi(rng, side_um=300.0 if quick else 600.0)
    nuc_g = C.labels_to_grid(Sx["mask"], Sx["upp"]); exp_g = C.expand_labels_grid(nuc_g, 8.0)
    kern = (2.0, 2.0, 0.05, 30.0); d_img = np.array([0.5, -0.5])
    c = np.array([0.6, 1.0, 0.9, 1.5])[C.parity_class(Sx["rows"], Sx["cols"])]; c = c / c.mean()
    labs = np.unique(exp_g[exp_g > 0]); Dtrue = {int(l): float(v) for l, v in zip(labs, rng.lognormal(np.log(15), 0.5, len(labs)))}
    St = C.nc_setup(nuc_g, exp_g, Dtrue, 15.0, 3.0)
    for r_true in (0.7, 1.0, 2.0, 4.0):
        lam = c * ref_blur_bins(C.nc_field(St, r_true), Sx["pos"], 2.0, 0.05, 30.0, d_img)
        y = rng.poisson(lam).astype(float)
        obs = C.star_stat(y, Sx["L"], Sx["keep"])
        Dobs = C.observed_cell_density(Sx["L"], Sx["keep"], y)
        So = C.nc_setup(nuc_g, exp_g, Dobs, float(np.median(list(Dobs.values()))), float(y[Sx["L"]["exp"] == 0].mean()))
        U = np.random.default_rng(C.SEED).uniform(size=(5, len(y)))
        am = Sx["L"]["exp"] == 0
        IX = C.stat_index(Sx["L"], Sx["keep"])
        rhat, _ = C.invert_r_refined(obs, So, Sx["pos"], kern, d_img, c, U, IX, Sx["L"], Sx["keep"], Dobs, float(y[am].mean()), am)
        check(f"T7 nucleus/ring ratio inverted, r={r_true}", np.isfinite(rhat) and abs(rhat / r_true - 1) <= 0.03, f"r_hat={rhat:.4f} rel err={rhat / r_true - 1:+.4f}")
    # coerenza: il blur di riferimento e quello del modulo coincidono
    lam_ref = ref_blur_bins(C.nc_field(St, 2.0), Sx["pos"], 2.0, 0.05, 30.0, d_img)
    lam_mod = C.blur_at(C.nc_field(St, 2.0), Sx["pos"], 2.0, 2.0, 0.05, 30.0, d_img)
    inner = (Sx["pos"][:, 0] > 60) & (Sx["pos"][:, 0] < Sx["pos"][:, 0].max() - 60) & (Sx["pos"][:, 1] > 60) & (Sx["pos"][:, 1] < Sx["pos"][:, 1].max() - 60)
    rel = np.max(np.abs(lam_mod[inner] / lam_ref[inner] - 1))
    check("T7 module blur equals reference real-space blur", rel < 0.02, f"max rel diff={rel:.4f}")


# ---------------------------------------------------------------- T10, T11 (rinforzo dopo la revisione del codice, rilievi F03, F04, C10)
def t10():
    """expected_bins contro un'implementazione indipendente nello spazio reale (sottogriglia, convoluzione diretta, spostamento intero, media dei
    4 x 4 sottopixel di ogni bin), con sorgente asimmetrica e scarto diverso sulle due righe/colonne: verifica le convenzioni ASSOLUTE su entrambi gli assi."""
    rng = np.random.default_rng(10); nb = (50, 60); n = (nb[0] * C.SUB, nb[1] * C.SUB)
    T = np.zeros(n); T[60:120, 40:200] = 1.0; T[150:170, 100:110] = 3.0; T += 0.2 * ndimage.gaussian_filter(rng.random(n), 2)
    s, dr, dc = 1.5, 1.0, -1.5                       # multipli di 0.5 µm
    ps = C.pad_shape(n, 120); kr, kc = C.freq_grid(ps)
    g = C.expected_bins(C.source_ft(T, ps), kr, kc, ps, nb, s, s, 0.0, 30.0, dr, dc)
    R_ = int(np.ceil(5 * s / C.SUB_UM)); yy, xx = np.mgrid[-R_:R_ + 1, -R_:R_ + 1] * C.SUB_UM
    K = np.exp(-(yy ** 2 + xx ** 2) / (2 * s ** 2)); K /= K.sum()
    H = signal.fftconvolve(np.pad(T, R_), K, mode="same")[R_:-R_, R_:-R_]
    H = np.roll(H, (int(dr / C.SUB_UM), int(dc / C.SUB_UM)), axis=(0, 1))   # campo dei conteggi = campo H&E spostato di +delta
    Hb = H.reshape(nb[0], C.SUB, nb[1], C.SUB).mean(axis=(1, 3))             # media sul quadratino (integrazione esatta sui sottopixel)
    # il modulo integra il quadratino in modo continuo (sinc) su una sorgente a banda limitata: confronto lontano dai bordi, tolleranza relativa
    inner = (slice(8, -8), slice(8, -8))
    err = np.max(np.abs(g[inner] - Hb[inner])) / T.max()
    check("T10 expected_bins = independent real-space model (absolute, both axes, delta=(1,-1.5))", err < 0.01, f"max err/peak={err:.4f}")
    # mutante di convenzione su un solo asse: lo stesso confronto deve fallire
    g2 = C.expected_bins(C.source_ft(T, ps), kr, kc, ps, nb, s, s, 0.0, 30.0, dr + 0.5, dc)
    err2 = np.max(np.abs(g2[inner] - Hb[inner])) / T.max()
    check("T10b a 0.5 um row offset is detected by T10", err2 > 3 * err and err2 > 0.01, f"err with +0.5 um rows={err2:.4f}")

def t11():
    """Mappatura griglia -> immagine dei dati veri: la media delle coordinate in pixel dei 4 x 4 sottopixel di ogni bin (sub_points + affine locale,
    come in tools/R5_tiles.py) coincide con le coordinate del centro del bin di Space Ranger (tissue_positions), su una tessera per sezione."""
    import glob, R5_tiles as T5
    for key in ("A1", "A6"):
        f = sorted(glob.glob(f"/mnt/micron/geo_spatialtrans/R5/run1/tiles/{key}/t*.npz"))[0]
        z = np.load(f); Rr, Cr = int(z["R0"]), int(z["C0"]); M, t = z["M_local"], z["t_local"]
        rr, cc = T5.sub_points(Rr, Cr, 300); RR, CC = np.meshgrid(rr, cc, indexing="ij")
        PR = M[0, 0] * RR + M[0, 1] * CC + t[0]; PC = M[1, 0] * RR + M[1, 1] * CC + t[1]
        mr = PR.reshape(300, 4, 300, 4).mean(axis=(1, 3)); mc = PC.reshape(300, 4, 300, 4).mean(axis=(1, 3))
        PXR, PXC = T5.pixel_grids(T5.DS[key], 3350, 3350)
        a0, b0 = max(Rr, 0), max(Cr, 0); a1, b1 = min(Rr + 300, 3350), min(Cr + 300, 3350)
        er = np.nanmax(np.abs(mr[a0 - Rr:a1 - Rr, b0 - Cr:b1 - Cr] - PXR[a0:a1, b0:b1])); ec = np.nanmax(np.abs(mc[a0 - Rr:a1 - Rr, b0 - Cr:b1 - Cr] - PXC[a0:a1, b0:b1]))
        check(f"T11 sub-grid mapping centred on Space Ranger bin centres ({key})", max(er, ec) < 0.5, f"max |diff| px rows={er:.3f} cols={ec:.3f} (0.5 px = 0.14 um)")

# ---------------------------------------------------------------- T8 bootstrap riproducibile
def t8():
    v = np.random.default_rng(8).normal(size=40)
    a = C.boot_median(v, rng=np.random.default_rng(C.SEED)); b = C.boot_median(v, rng=np.random.default_rng(C.SEED))
    check("T8 bootstrap reproducible with seed", a[:3] == b[:3] and np.array_equal(a[3], b[3]), f"{a[:3]}")

if __name__ == "__main__":
    t0 = time.time()
    if os.environ.get("R5_ONLY_T7") == "1":
        t7(quick=False); print(f"T7 only: {sum(ok for _, ok, _ in RES)}/{len(RES)} PASS"); sys.exit(0)
    if os.environ.get("R5_CHEAP") == "1" or A.mutant:
        # test economici; un mutante e' rilevato se fallisce almeno un test che passa senza mutazioni (regola relativa, rilievo F02)
        t1(); t2(); t3(); t4(); t9(); t8(); t10(); t11(); t12_seed_reproducible(); t7(quick=True)
        base = {}
        if os.path.exists("results/R5/R5_test_cheap_baseline.json"):
            base = json.load(open("results/R5/R5_test_cheap_baseline.json"))
        if not A.mutant:
            json.dump({n: ok for n, ok, _ in RES}, open("results/R5/R5_test_cheap_baseline.json", "w"), indent=1)
            print(f"cheap: {sum(ok for _, ok, _ in RES)}/{len(RES)} PASS"); sys.exit(0)
        newly = [n for n, ok, _ in RES if not ok and base.get(n, False)]
        print(f"MUTANT {A.mutant} {'DETECTED' if newly else 'NOT DETECTED'} by: {newly[:4]}")
        sys.exit(0)
    t1(); t2(); t3(); t4(); t9(); t8(); t12_seed_reproducible()
    if A.mutant:
        detected = not all(ok for _, ok, _ in RES)
        if not detected:
            t7(quick=True)
            detected = not all(ok for _, ok, _ in RES)
        if not detected:
            D = t5(A.workers, 2); judge_recovery(D); detected = not all(ok for _, ok, _ in RES)
        print(f"MUTANT {A.mutant} {'DETECTED' if detected else 'NOT DETECTED'}")
        sys.exit(0)
    t7(quick=A.quick)
    D = t5(A.workers, 2 if A.quick else A.nrep)
    judge_recovery(D)
    os.makedirs("results/R5", exist_ok=True)
    if not A.quick and os.environ.get("R5_TAG") is None:
        D.to_csv("results/R5/R5_G1_recovery.csv", index=False)
    else:
        D.to_csv(f"/tmp/r5_quick_{os.environ.get('R5_TAG', 'x')}.csv", index=False)
    n_ok = sum(ok for _, ok, _ in RES)
    summary = f"R5_test: {n_ok}/{len(RES)} PASS ({time.time() - t0:.0f} s)"
    print(summary)
    if not A.quick:
        fo = "results/R5/R5_test_output.txt" if os.environ.get("R5_TAG") is None else f"/tmp/r5_test_output_{os.environ['R5_TAG']}.txt"
        with open(fo, "w") as fh:
            for n, ok, d in RES: fh.write(f"{'PASS' if ok else 'FAIL'}\t{n}\t{d}\n")
            fh.write(summary + "\n")
