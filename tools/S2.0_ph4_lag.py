"""tools/S2.0_ph4_lag.py — S2.0, post hoc PH-4 (dal dubbio D-09 della revisione G-5 e): il passo 1 negativo dei fattori non periodici
delle strisce viene dallo stimatore? Sezioni sintetiche (generatore di G-1, tools/S2.0_test.py synth_section, struttura del tessuto amp)
con fattori di riga e colonna noti: indipendenti, oppure media mobile di ordine 1 con passo 1 vero rho. Stimatore e passo 1 riproducibile
come in M1 (tools/S2.0_stripes.py). Uscita: results/S2.0/S2.0_PH4_lag.csv. Uso: .venv/bin/python tools/S2.0_ph4_lag.py"""
import sys, importlib.util
import numpy as np, pandas as pd
sys.path.insert(0, "tools")
import S2_0_common as S, R5_common as C
sp = importlib.util.spec_from_file_location("t", "tools/S2.0_test.py"); T = importlib.util.module_from_spec(sp); sp.loader.exec_module(T)
sp2 = importlib.util.spec_from_file_location("st", "tools/S2.0_stripes.py"); ST = importlib.util.module_from_spec(sp2); sys.argv = ["x", "/tmp"]; sp2.loader.exec_module(ST)
from scipy import ndimage
def section(sd, rho, seed, amp, n=800, mean=5.0):
    rng = np.random.default_rng([S.SEED, 900 + seed])
    yy, xx = np.mgrid[0:n, 0:n]; Tm = np.zeros((n, n), bool)
    for cy, cx, r in ((250, 250, 220), (560, 520, 260), (200, 620, 150), (640, 180, 140)): Tm |= (yy - cy) ** 2 + (xx - cx) ** 2 < r ** 2
    z = ndimage.gaussian_filter(rng.normal(size=(n, n)), 10); z /= z.std(); m = mean * np.exp(amp * z)
    def fac():
        if rho == 0: return rng.normal(0, sd, n)
        th = (1 - np.sqrt(1 - 4 * rho ** 2)) / (2 * rho)            # rho = th / (1 + th^2), |th| < 1
        e = rng.normal(0, 1, n + 1); x = e[1:] + th * e[:-1]; return sd * x / np.sqrt(1 + th ** 2)
    al, be = fac(), fac()
    c = (T.PAR / T.PAR.mean())[C.parity_class(yy, xx)]
    U = rng.poisson(m * c * np.exp(al[:, None] + be[None, :]) * Tm).astype(float)
    return U, Tm, al, be
rows = []
for amp in (0.2, 0.4):
    for rho in (0.0, -0.2):
        for k in range(4):
            U, Tm, al, be = section(0.2, rho, k, amp)
            for ax, nm, f in ((0, "rows", al), (1, "cols", be)):
                idx = np.arange(U.shape[ax])
                h1, h2, n1, n2 = S.half_split_factors(U, Tm, ax); ok = (n1 >= 200) & (n2 >= 200)
                _, r1 = S.periodic_split(np.where(ok, h1, np.nan), idx); _, r2 = S.periodic_split(np.where(ok, h2, np.nan), idx)
                _, ft = S.periodic_split(np.where(ok, f, np.nan), idx)
                oo = np.isfinite(ft[:-1]) & np.isfinite(ft[1:])
                rows.append({"amp": amp, "rho_true": rho, "seed": k, "axis": nm, "n_eligible": int(ok.sum()), "sd_true_np": float(np.nanstd(ft)),
                             "sd_repro": S.repro_sd(r1, r2), "lag1_true_sample": float(np.corrcoef(ft[:-1][oo], ft[1:][oo])[0, 1]), "lag1_repro": ST.lag1_repro(r1, r2)})
                print(rows[-1], flush=True)
D = pd.DataFrame(rows); D.to_csv("results/S2.0/S2.0_PH4_lag.csv", index=False)
print(D.groupby(["amp", "rho_true"])[["lag1_true_sample", "lag1_repro", "sd_true_np", "sd_repro"]].mean().round(3).to_string())
