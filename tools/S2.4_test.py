"""tools/S2.4_test.py — S2.4, G-1 (parte Python): strumenti di analisi dello step su dati sintetici con valori noti.
(1) parita' delle aree con l'indice assoluto della griglia (regressione dell'errore di E-1); (2) adattamento (b, w2) di R-3 a ell2 fisso;
(3) parte lenta di R-2 su fattori noti; (4) fasce dei lumi su un lumo sintetico. Uso: .venv/bin/python tools/S2.4_test.py"""
import sys, importlib.util
import numpy as np
sys.path.insert(0, "tools")
def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
EX = _load("ex", "tools/S2.4_explore.py"); PH = _load("ph5", "tools/S2.0_ph5_stripes_hp.py"); R3 = _load("r3", "tools/S2.4_R3_glass.py"); R2 = _load("r2", "tools/S2.4_R2_slow.py")
res = []
def ck(name, ok, val=""): res.append(bool(ok)); print(("PASS" if ok else "FAIL"), name, val)
rng = np.random.default_rng(1)
# (1) area che parte da una riga e da una colonna dispari della griglia: la parita' stimata deve essere quella della griglia
n = 520; i0, j0 = 1001, 333
I, J = np.meshgrid(np.arange(n) + i0, np.arange(n) + j0, indexing="ij")
par = np.where(I % 2 == 1, 1.5, 1.0) * np.where(J % 2 == 1, 1.7, 1.0)
al = rng.normal(0, 0.2, n); be = rng.normal(0, 0.2, n)
U = rng.poisson(20 * par * np.exp(al[:, None] + be[None, :])).astype(float)
Um, Tm, a0, b0 = EX.area_matrix(I.ravel(), J.ravel(), U.ravel(), np.ones(n * n, bool))
o = EX.stripe_summary(PH, Um, Tm)
ck("1 parita' con l'indice assoluto (righe 1.5, colonne 1.7)", abs(o["stripes_rows_odd_over_even"] - 1.5) < 0.03 and abs(o["stripes_cols_odd_over_even"] - 1.7) < 0.03,
   f"{o['stripes_rows_odd_over_even']:.3f} {o['stripes_cols_odd_over_even']:.3f}")
ck("1b deviazione standard della parte fine circa 0.2", abs(o["stripes_rows_sd_repro_hp"] - 0.2) < 0.02 and abs(o["stripes_cols_sd_repro_hp"] - 0.2) < 0.02,
   f"{o['stripes_rows_sd_repro_hp']:.3f} {o['stripes_cols_sd_repro_hp']:.3f}")
# (2) adattamento convesso di R-3
m = 400000; core = rng.gamma(2, 0.3, m); e2 = rng.gamma(2, 0.2, m); y = rng.poisson(0.15 + core + 0.04 * (e2 - core))
nll, x = R3.fit_bw(y.astype(float), core, e2, np.array([0.3, 0.01]))
ck("2 R-3: (b, w2) ritrovati", abs(x[0] - 0.15) < 0.01 and abs(x[1] - 0.04) < 0.01, f"b {x[0]:.4f} w2 {x[1]:.4f}")
nll0, x0 = R3.fit_bw(rng.poisson(core).astype(float), core, e2, np.array([0.3, 0.01]))
ck("2b R-3: b = 0 e w2 = 0 al bordo quando non ci sono", x0[0] < 0.01 and x0[1] < 0.01, f"b {x0[0]:.4f} w2 {x0[1]:.4f}")
# (3) parte lenta: media mobile su 21 di un'onda lenta piu' rumore indipendente
t = np.arange(3000); slow_true = 0.1 * np.sin(2 * np.pi * t / 400)
s1, _, _ = R2.slow(slow_true + rng.normal(0, 0.2, 3000)); s2, _, _ = R2.slow(slow_true + rng.normal(0, 0.2, 3000))
r, cv, nn = R2.cc(s1, s2)
# correlazione attesa fra le meta': 0.005 / (0.005 + 0.2^2 / 21) = 0.72 (il rumore mediato su 21 righe ha SD 0.044)
ck("3 R-2: parte lenta riproducibile ritrovata (SD circa 0.0707, correlazione circa 0.72)", abs(np.sqrt(cv) - 0.1 / np.sqrt(2)) < 0.01 and abs(r - 0.724) < 0.1, f"sd {np.sqrt(cv):.4f} r {r:.3f}")
# (4) fasce dei lumi: lumo quadrato di 100 um in un tessuto di 300 um, UMI 10 nel tessuto e 1 nel lumo
import shapely
RL = EX._load("s23real", "tools/S2.3_real.py"); RL.BANDS = EX.DEEP_BANDS
R = shapely.Polygon([(0, 0), (300, 0), (300, 300), (0, 300)], [[(100, 100), (200, 100), (200, 200), (100, 200)]])
yy, xx = np.mgrid[0:150, 0:150] * 2.0 + 1.0; xx = xx.ravel(); yy = yy.ravel()
inreg = shapely.contains_xy(R, xx, yy); u = np.where(inreg, 10.0, 1.0)
L = RL.lumen_stats(R, np.c_[yy, xx], u, inreg, 300.0)
ok = L["n_lumens"] == 1 and all(abs(L[f"band_{a:g}_{b:g}_sum"] / L[f"band_{a:g}_{b:g}_n"] / (L["ref_sum"] / L["ref_n"]) - 0.1) < 1e-12 for a, b in EX.DEEP_BANDS if L[f"band_{a:g}_{b:g}_n"] > 0)
ck("4 fasce dei lumi: rapporto 0.1 in ogni fascia", ok, f"lumi {L['n_lumens']}")
nf = sum(not x for x in res); print(f"\nTOTALE {len(res) - nf} PASS, {nf} FAIL"); sys.exit(1 if nf else 0)
