"""tools/S2.0_test.py — S2.0, cancello G-1 (docs/steps/S2.0.md, cancelli afef06b, correzioni f2a23d7). Stampa PASS/FAIL per asserzione.
Uso: .venv/bin/python tools/S2.0_test.py [--workers N]          batteria completa (parti: ident64, ident32, recovery, stripes, proto, seeds, mutants)
     .venv/bin/python tools/S2.0_test.py --part <nome>          una parte sola (usata internamente nei sottoprocessi)"""
import os, sys, json, hashlib, subprocess, time
import numpy as np
from multiprocessing import Pool
from scipy import ndimage, fft as sfft
sys.path.insert(0, "tools")
import R5_common as C
import S2_0_common as S

RES = []
def check(name, ok, val=""):
    RES.append((name, bool(ok))); print(f"{'PASS' if ok else 'FAIL'} {name} {val}", flush=True)

NB = 300
TILE = (slice(75, 225), slice(75, 225))
PAR = np.array([0.6, 1.0, 0.9, 1.5])

def blob(nb=NB):
    n = nb * C.SUB
    yy, xx = (np.mgrid[0:n, 0:n] + 0.5) * C.SUB_UM; H = n * C.SUB_UM
    t = (yy - H / 2) ** 2 + (xx - H / 2) ** 2 < 130.0 ** 2
    t &= (yy - H / 2 - 40) ** 2 + (xx - H / 2 + 30) ** 2 > 25.0 ** 2
    t |= (yy - H / 2 + 60) ** 2 + (xx - H / 2 - 210) ** 2 < 15.0 ** 2
    return t.astype(float)

def grid_P(T, f, nb=NB, pad=900):
    ps = C.pad_shape(f.shape, pad); kr, kc = C.freq_grid(ps)
    Tp = np.zeros(ps, S.DT); Tp[:f.shape[0], :f.shape[1]] = T
    fp = np.zeros(ps, S.DT); fp[:f.shape[0], :f.shape[1]] = f
    dp = np.zeros(ps, S.DT); dp[:f.shape[0], :f.shape[1]] = ndimage.distance_transform_edt(f > 0.5) * C.SUB_UM
    P = {"Tft": sfft.rfft2(Tp, workers=1), "kr": kr, "kc": kc, "pshape": ps, "region_bins": (nb, nb)}
    return P, {"T": Tp, "f": fp, "din": dp}

def field_real(rule, prm, P, X):
    return sfft.irfft2(S.field_ft(rule, prm, P["Tft"], P["kr"], P["kc"], P["pshape"], X), s=P["pshape"], workers=1)

BASE = {"w": 0.08, "ell": 25.0, "dr": 0.5, "dc": -0.5}
def rel(a, b):
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))

# ------------------------------------------------------------------ parti
def part_ident64():
    """Doppia precisione: massa e annidamento (AM-1); C3 lascia invariato l'interno del tessuto."""
    assert S.DT is np.float64
    f = blob(); rng = np.random.default_rng([S.SEED, 1])
    T = f * (20.0 + 10.0 * ndimage.gaussian_filter(rng.normal(size=f.shape), 6))     # sorgente con struttura
    P, X = grid_P(T, f)
    tot = T.sum()
    prms = {"C0": dict(BASE, sigma=3.0), "C1": dict(BASE, q=0.4, sigma_m=3.5), "C2": dict(BASE, sigma_e=3.5, D=5.0), "C3": dict(BASE, sigma_g=4.0)}
    for r, prm in prms.items():
        m = field_real(r, prm, P, X).sum()
        check(f"a1_mass_{r}", abs(m - tot) / tot <= 1e-10, f"rel={abs(m - tot) / tot:.2e}")
    g0 = S.g_rule("C0", dict(BASE, sigma=3.5), P, X)
    check("a2_C1_q1_eq_C0", rel(S.g_rule("C1", dict(BASE, q=1.0, sigma_m=3.5), P, X), g0) <= 1e-10)
    d = rel(S.g_rule("C2", dict(BASE, sigma_e=3.5, D=1e15), P, X), g0); check("a2_C2_Dinf_eq_C0", d <= 1e-10, f"rel={d:.2e}")
    gn = S.g_rule("C0", dict(BASE, sigma=S.SIG_IN), P, X)
    d = rel(S.g_rule("C2", dict(BASE, sigma_e=3.5, D=1e-3), P, X), gn); check("a2_C2_D0_eq_narrow", d <= 1e-10, f"rel={d:.2e}")
    d = rel(S.g_rule("C3", dict(BASE, sigma_g=1e-7), P, X), gn); check("a2_C3_s0_eq_narrow", d <= 1e-10, f"rel={d:.2e}")
    # C3: il trasporto tocca solo cio' che atterra sul vetro: nell'interno del tessuto (>= 30 µm dal bordo) il campo resta quello stretto
    g3 = S.g_rule("C3", dict(BASE, sigma_g=8.0), P, X)
    db = C.edge_distance_bins(f); inner = db <= -30
    d = float(np.max(np.abs(g3[inner] - gn[inner]) / gn[inner])); check("a5_C3_interior_unchanged", d <= 1e-3, f"maxrel={d:.2e}")

def part_ident32():
    """Precisione singola (come il calcolo vero): C0 = modello di R5; C1(q=1) = C0; make_blur(C0) = R5 blur_at."""
    assert S.DT is np.float32
    f = blob(); T = f * 20.0
    P, X = grid_P(T, f)
    P["Tft"] = C.source_ft(T, P["pshape"])
    for sg, w in ((3.0, 0.08), (1.2, 0.0)):
        g_r5 = C.expected_bins(P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], sg, sg, w, 25.0, 0.5, -0.5)
        g_s = S.g_rule("C0", {"sigma": sg, "w": w, "ell": 25.0, "dr": 0.5, "dc": -0.5}, P, X)
        check(f"a3_C0_eq_R5_s{sg}", np.array_equal(g_r5, g_s), f"maxabs={np.max(np.abs(g_r5 - g_s)):.2e}")
    a = S.g_rule("C0", dict(BASE, sigma=3.5), P, X); b = S.g_rule("C1", dict(BASE, q=1.0, sigma_m=3.5), P, X)
    check("a3_C1_q1_eq_C0_f32", np.array_equal(a, b), f"maxabs={np.max(np.abs(a - b)):.2e}")
    rng = np.random.default_rng([S.SEED, 2]); F = rng.gamma(2.0, 2.0, size=(400, 400)); pos = rng.uniform(5, 195, size=(500, 2))
    kern = (2.0, 2.0, 0.1, 20.0); d_img = np.array([0.3, -0.7])
    v5 = C.blur_at(F, pos, *kern, d_img); vs = S.make_blur("C0", {"sigma": 2.0, "w": 0.1, "ell": 20.0}, pos, d_img)(F)
    check("a3_make_blur_C0_eq_R5", rel(vs, v5) <= 1e-10, f"rel={rel(vs, v5):.2e}")

def synth_edge(rule, prm, seed, mean_bin=60.0, b=0.3):
    f = blob(); T = f * mean_bin
    P0, X0 = grid_P(T, f)
    g = S.g_rule(rule, prm, P0, X0)
    rows, cols = np.mgrid[0:NB, 0:NB]
    c = (PAR / PAR.mean())[C.parity_class(rows, cols)]
    y = np.random.default_rng([S.SEED, 100 + seed]).poisson(c * (np.maximum(g, 0) + b)).astype(float)
    return y, f, rows, cols

TRUE = {"C1": dict(BASE, q=0.5, sigma_m=3.5), "C2": dict(BASE, sigma_e=3.5, D=5.0), "C3": dict(BASE, sigma_g=4.0)}

def recovery_job(args):
    rule, rep = args
    y, f, rows, cols = synth_edge(rule, TRUE[rule], rep)
    P = C.prepare_region(y, f, rows, cols); X = S.region_extra(P, f)
    sel = C.likelihood_mask(P["occ"], P["dist"], TILE)
    k2, _ = C.fit_k2_screened(P["Tft"], P["kr"], P["kc"], P["pshape"], P["region_bins"], y, P["c"], sel, None)
    fit = S.fit_rule(rule, P, X, y, sel, S.starts_for(rule, k2, None))
    return rule, rep, fit

def part_recovery(nw):
    jobs = [(r, k) for r in ("C1", "C2", "C3") for k in range(3)]
    with Pool(min(nw, len(jobs))) as p:
        out = p.map(recovery_job, jobs, chunksize=1)
    for rule in ("C1", "C2", "C3"):
        fits = [o[2] for o in out if o[0] == rule]
        tr = dict(TRUE[rule])
        keys = {"C1": ["p", "sigma_m"], "C2": ["sigma_e", "D"], "C3": ["sigma_g"]}[rule]
        if rule == "C1": tr["p"] = tr["q"] * (1 - tr["w"])
        for k in keys:
            est = float(np.median([fi[k] for fi in fits]))
            check(f"a4_recover_{rule}_{k}", abs(est - tr[k]) <= 0.10 * tr[k], f"est={est:.4g} true={tr[k]:.4g}")
        est = float(np.median([fi["w"] for fi in fits]))
        check(f"a4_recover_{rule}_w", abs(est - tr["w"]) <= 0.02, f"est={est:.4g} true={tr['w']:.4g}")
        json.dump(fits, open(f"/tmp/s20_recovery_{rule}.json", "w"))

def synth_section(sd, seed, n=800, mean=5.0, amp=0.1):
    rng = np.random.default_rng([S.SEED, 200 + seed])
    yy, xx = np.mgrid[0:n, 0:n]
    T = np.zeros((n, n), bool)
    for cy, cx, r in ((250, 250, 220), (560, 520, 260), (200, 620, 150), (640, 180, 140)):
        T |= (yy - cy) ** 2 + (xx - cx) ** 2 < r ** 2
    z = ndimage.gaussian_filter(rng.normal(size=(n, n)), 10); z /= z.std()
    m = mean * np.exp(amp * z)
    al = rng.normal(0, sd, n) if sd > 0 else np.zeros(n); be = rng.normal(0, sd, n) if sd > 0 else np.zeros(n)
    c = (PAR / PAR.mean())[C.parity_class(yy, xx)]
    lam = m * c * np.exp(al[:, None] + be[None, :]) * T
    return rng.poisson(lam).astype(float), T

def stripe_summary(U, T):
    al, be, mr, mc = S.stripe_factors(U, T)
    out = {}
    for ax, f, nm in ((0, al, "rows"), (1, be, "cols")):
        idx = np.arange(len(f))
        mu, _ = S.periodic_split(f, idx)
        h1, h2, n1, n2 = S.half_split_factors(U, T, ax)
        ok = (n1 >= 200) & (n2 >= 200)
        _, r1 = S.periodic_split(np.where(ok, h1, np.nan), idx); _, r2 = S.periodic_split(np.where(ok, h2, np.nan), idx)
        out[nm] = {"ratio": float(np.exp(mu[1] - mu[0])), "sd": S.repro_sd(r1, r2), "n_eligible": int(ok.sum())}
    tp = S.tile_parity(U, T)
    out["tile_sd_col"] = S.repro_sd(tp[:, 2], tp[:, 3]) if len(tp) else float("nan")
    out["tile_sd_row"] = S.repro_sd(tp[:, 4], tp[:, 5]) if len(tp) else float("nan")
    return out

def part_stripes():
    rr = (PAR[2] + PAR[3]) / (PAR[0] + PAR[1]); rc = (PAR[1] + PAR[3]) / (PAR[0] + PAR[2])
    for sd in (0.0, 0.1):
        U, T = synth_section(sd, int(sd * 10))
        o = stripe_summary(U, T)
        check(f"b1_parity_rows_sd{sd}", abs(o["rows"]["ratio"] / rr - 1) <= 0.01, f"est={o['rows']['ratio']:.4f} true={rr:.4f}")
        check(f"b1_parity_cols_sd{sd}", abs(o["cols"]["ratio"] / rc - 1) <= 0.01, f"est={o['cols']['ratio']:.4f} true={rc:.4f}")
        for nm in ("rows", "cols"):
            if sd == 0:
                check(f"b2_sd0_{nm}", o[nm]["sd"] < 0.01, f"est={o[nm]['sd']:.4f} n={o[nm]['n_eligible']}")
            else:
                check(f"b2_sd{sd}_{nm}", abs(o[nm]["sd"] - sd) <= 0.01, f"est={o[nm]['sd']:.4f} true={sd} n={o[nm]['n_eligible']}")
        if sd == 0:
            check("b3_tile_sd0_col", o["tile_sd_col"] < 0.01, f"est={o['tile_sd_col']:.4f}")
            check("b3_tile_sd0_row", o["tile_sd_row"] < 0.01, f"est={o['tile_sd_row']:.4f}")

def synth_roi_proto(seed):
    rng = np.random.default_rng([S.SEED, 300 + seed])
    R = C.synth_roi(rng, side_um=400.0, nuc_r_um=4.5, density=1 / 220.0)
    nuc_g = C.labels_to_grid(R["mask"], R["upp"]); exp_g = C.expand_labels_grid(nuc_g, 8.0)
    Dobs = {int(k): 5.0 for k in R["keep"]}
    St = C.nc_setup(nuc_g, exp_g, Dobs, 5.0, 0.5)
    c = (PAR / PAR.mean())[C.parity_class(R["rows"], R["cols"])]
    return R, St, c, rng

def part_proto():
    U2 = np.random.default_rng([S.SEED, 400]).poisson(3.0, size=(64, 48))
    manual = np.array([[U2[4*i:4*i+4, 4*j:4*j+4].sum() for j in range(12)] for i in range(16)])
    check("c1_bin8_sum", np.array_equal(S.bin8(U2), manual))
    R, St, c, rng = synth_roi_proto(0)
    F = C.nc_field(St, 1.0)
    inner = np.all((R["pos"] > 40) & (R["pos"] < 360), axis=1)
    for rule, prm in (("C0", {"sigma": 1.0, "w": 0.0, "ell": 30.0}), ("C1", {"q": 0.5, "sigma_m": 3.0, "w": 0.0, "ell": 30.0})):
        lam = S.make_blur(rule, prm, R["pos"], np.zeros(2))(F)
        # massa: somma dei bin del riquadro interno contro la media del campo vero sugli stessi bin (bin di 2 µm = 4 x 4 sottopixel)
        Fb = F[: (F.shape[0] // 4) * 4, : (F.shape[1] // 4) * 4].reshape(F.shape[0] // 4, 4, F.shape[1] // 4, 4).mean(axis=(1, 3))
        ib = (np.floor(R["pos"][:, 0] / 2).astype(int), np.floor(R["pos"][:, 1] / 2).astype(int))
        ok = inner & (ib[0] < Fb.shape[0]) & (ib[1] < Fb.shape[1])
        d = abs(lam[ok].sum() / Fb[ib[0][ok], ib[1][ok]].sum() - 1)
        check(f"c2_mass_{rule}", d <= 1e-3, f"rel={d:.2e}")
    lab = S.deep_labels(R["mask"], R["rr"], R["cc"], R["upp"], 1.0)
    lam = c * np.maximum(S.make_blur("C0", {"sigma": 0.05, "w": 0.0, "ell": 30.0}, R["pos"], np.zeros(2))(F), 1e-12)
    y = rng.poisson(lam).astype(float)
    psi, mb, nn = S.psi_index(y, lab, c, R["keep"])
    check("c3_psi_homogeneous", abs(psi - 1) <= 0.03, f"psi={psi:.4f} n_nuclei={nn} mbar={mb:.2f}")
    ps0 = np.mean([S.psi_index(np.random.default_rng([S.SEED, 500 + k]).poisson(lam).astype(float), lab, c, R["keep"])[0] for k in range(5)])
    for v in (0.1, 0.3):
        g = np.random.default_rng([S.SEED, 600 + int(v * 10)]).gamma(1 / v, v, size=len(lam))
        yv = np.random.default_rng([S.SEED, 700 + int(v * 10)]).poisson(lam * g).astype(float)
        psv, mbv, _ = S.psi_index(yv, lab, c, R["keep"])
        vh = (psv - ps0) / mbv
        check(f"c4_recover_v{v}", abs(vh - v) <= 0.2 * v, f"vhat={vh:.4f}")

def part_seedrun():
    """Pipeline piccola deterministica; stampa un'impronta (eseguita due volte in processi diversi)."""
    h = hashlib.md5()
    U, T = synth_section(0.05, 9, n=300)
    o = stripe_summary(U, T); h.update(json.dumps(o, sort_keys=True, default=str).encode())
    lo, hi, _ = S.block_boot_sd(np.random.default_rng(S.SEED).normal(size=300), np.random.default_rng(S.SEED + 1).normal(size=300), B=200)
    h.update(f"{lo:.12g},{hi:.12g}".encode())
    R, St, c, rng = synth_roi_proto(1)
    lam = c * np.maximum(S.make_blur("C1", {"q": 0.3, "sigma_m": 3.0, "w": 0.05, "ell": 20.0}, R["pos"], np.zeros(2))(C.nc_field(St, 1.3)), 1e-12)
    h.update(np.random.default_rng([S.SEED, 9]).poisson(lam).tobytes())
    print("SEEDHASH", h.hexdigest(), flush=True)

def part_seeds():
    hs = []
    for k in range(2):
        out = subprocess.run([sys.executable, "tools/S2.0_test.py", "--part", "seedrun"], capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED=str(k + 1)))
        hs.append([l.split()[1] for l in out.stdout.splitlines() if l.startswith("SEEDHASH")])
    check("d1_seeds_reproducible_across_processes", len(hs[0]) == 1 and hs[0] == hs[1], str(hs))

# ------------------------------------------------------------------ mutanti (applicati in un sottoprocesso con S20_MUTANT)
def apply_mutant(m):
    if m == "M1":      # nucleo non normalizzato: manca (1 - w) in C1
        orig = S.field_ft
        def ff(rule, prm, Tft, kr, kc, ps, X=None):
            if rule == "C1":
                q = prm["q"]; E = C.exp_ft(kr, kc, prm["ell"])
                return Tft * (((1 - q) * S._G(kr, kc, S.SIG_IN) + q * S._G(kr, kc, prm["sigma_m"])) + prm["w"] * E)
            return orig(rule, prm, Tft, kr, kc, ps, X)
        S.field_ft = ff
    elif m == "M2":    # pi applicato al vetro invece che alle sorgenti
        orig = S.field_ft
        def ff(rule, prm, Tft, kr, kc, ps, X=None):
            if rule == "C2":
                w = prm["w"]; E = C.exp_ft(kr, kc, prm["ell"])
                dout = ndimage.distance_transform_edt(X["f"] <= 0.5) * S.SUB_UM
                pi = np.where(X["f"] <= 0.5, np.exp(-dout / prm["D"]), 0.0).astype(S.DT)
                Te = S._rfft(X["T"] * pi, ps); Ti = Tft - Te
                return (1 - w) * (Te * S._G(kr, kc, prm["sigma_e"]) + Ti * S._G(kr, kc, S.SIG_IN)) + w * E * Tft
            return orig(rule, prm, Tft, kr, kc, ps, X)
        S.field_ft = ff
    elif m == "M3":    # trasporto applicato al tessuto invece che al vetro
        orig = S.field_ft
        def ff(rule, prm, Tft, kr, kc, ps, X=None):
            if rule == "C3":
                w = prm["w"]; E = C.exp_ft(kr, kc, prm["ell"])
                Fin = S._irfft(Tft * ((1 - w) * S._G(kr, kc, S.SIG_IN) + w * E), ps)
                Lf = Fin * (1 - X["f"]) + S._irfft(S._rfft(Fin * X["f"], ps) * S._G(kr, kc, prm["sigma_g"]), ps)
                return S._rfft(Lf, ps)
            return orig(rule, prm, Tft, kr, kc, ps, X)
        S.field_ft = ff
    elif m == "M4":    # parita' ignorata nelle strisce
        S.periodic_split = lambda f, idx: (np.zeros(2), np.asarray(f, float))
    elif m == "M5":    # meta' sovrapposte nella stima riproducibile
        orig = S.half_masks
        S.half_masks = lambda T, axis: (orig(T, axis)[0], orig(T, axis)[0] | orig(T, axis)[1])
    elif m == "M6":    # psi senza correzione di parita'
        orig = S.psi_index
        S.psi_index = lambda u, lab, c, keep, min_bins=4: orig(u, lab, np.ones(len(c)), keep, min_bins)

def part_cheap():
    part_ident64(); part_stripes(); part_proto()

def part_mutants():
    for m in ("M1", "M2", "M3", "M4", "M5", "M6"):
        out = subprocess.run([sys.executable, "tools/S2.0_test.py", "--part", "cheap"], capture_output=True, text=True,
                             env=dict(os.environ, S20_MUTANT=m, S20_F64="1"))
        nf = sum(1 for l in out.stdout.splitlines() if l.startswith("FAIL"))
        check(f"e_mutant_{m}_detected", nf >= 1, f"fails={nf} rc={out.returncode}")

if __name__ == "__main__":
    args = sys.argv[1:]
    if os.environ.get("S20_MUTANT"):
        apply_mutant(os.environ["S20_MUTANT"])
    if "--part" in args:
        part = args[args.index("--part") + 1]
        nw = int(args[args.index("--workers") + 1]) if "--workers" in args else 9
        {"ident64": part_ident64, "ident32": part_ident32, "recovery": lambda: part_recovery(nw), "stripes": part_stripes,
         "proto": part_proto, "seedrun": part_seedrun, "seeds": part_seeds, "mutants": part_mutants, "cheap": part_cheap}[part]()
        if part != "seedrun":
            print(f"PART {part}: {sum(o for _, o in RES)}/{len(RES)} PASS", flush=True)
        sys.exit(0)
    nw = int(args[args.index("--workers") + 1]) if "--workers" in args else 9
    t0 = time.time(); allres = []
    for part, env in (("ident64", {"S20_F64": "1"}), ("ident32", {}), ("stripes", {}), ("proto", {}), ("seeds", {}), ("mutants", {}), ("recovery", {})):
        out = subprocess.run([sys.executable, "tools/S2.0_test.py", "--part", part, "--workers", str(nw)], capture_output=True, text=True,
                             env=dict(os.environ, **env))
        lines = [l for l in out.stdout.splitlines() if l.startswith(("PASS", "FAIL"))]
        for l in lines: print(l, flush=True)
        if out.returncode != 0:
            print(f"FAIL part_{part}_crashed rc={out.returncode} {out.stderr[-800:]}", flush=True); lines.append("FAIL crash")
        allres += lines
    npass = sum(l.startswith("PASS") for l in allres)
    print(f"S2.0_test: {npass}/{len(allres)} PASS ({time.time() - t0:.0f} s)")
