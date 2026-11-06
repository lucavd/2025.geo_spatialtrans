"""tools/S2.3_test.py — S2.3, cancello G-1 (d, i): campo e statistiche del banco (tools/S2_3_common.py). PASS/FAIL per asserzione.
Uso: .venv/bin/python tools/S2.3_test.py"""
import sys
import numpy as np
sys.path.insert(0, "tools")
import S2_3_common as S
res = []
def ck(name, ok, val=""):
    ok = bool(ok); res.append(ok); print(f"{'PASS' if ok else 'FAIL'} {name} {val}")
# (d) campo in Python
rng = np.random.default_rng([S.SEED, 1])
Z = S.gaussian_z((1200, 1200), 1.0, rng)
def cz(k): return np.mean([np.corrcoef(Z[:, :-k].ravel(), Z[:, k:].ravel())[0, 1], np.corrcoef(Z[:-k].ravel(), Z[k:].ravel())[0, 1]])
dz = max(abs(Z.var() - 1), *(abs(cz(k) - S.z_corr_theory(1.0, k)) for k in (1, 2, 4)))
ck("pd1 Z: varianza 1 e correlazione a 0.5, 1, 2 um come la teorica (entro 0.02)", dz < 0.02, f"{dz:.4f}")
F = S.gamma_field((1200, 1200), 0.5, 1.0, np.random.default_rng([S.SEED, 2]))
ck("pd2 F: media 1 (entro 0.02), varianza s2 (entro il 5 %)", abs(F.mean() - 1) < 0.02 and abs(F.var() / 0.5 - 1) < 0.05, f"{F.mean():.4f} {F.var():.4f}")
ck("pd3 campo spento = 1", np.all(S.gamma_field((5, 5), 0.0, 1.0, None) == 1))
ck("pd4 tabella finita", np.all(np.isfinite(S.gamma_table(0.5))) and np.all(np.isfinite(S.gamma_table(8.0))))
# (i) statistiche: griglia 400 x 400 quadratini, gruppi di 16 x 16
nb = 400; rows, cols = np.divmod(np.arange(nb * nb), nb); grp = (rows // 16) * 1000 + cols // 16; c = np.ones(nb * nb)
r = np.random.default_rng([S.SEED, 3])
y = r.poisson(10.0, nb * nb).astype(float)
st = S.group_stats(y, grp, c, rows, cols)
ck("pi1 psi = 1 su Poisson omogeneo (|psi - 1| <= 0.03)", abs(st["psi"] - 1) <= 0.03, f"{st['psi']:.4f}")
ck("pi2 rho1 = 0 su Poisson indipendente (|rho1| <= 0.02)", abs(st["rho1"]) <= 0.02, f"{st['rho1']:.4f}")
for v in (0.1, 0.3):
    g = r.gamma(1 / v, v, nb * nb); y = r.poisson(10.0 * g).astype(float)
    s = S.group_stats(y, grp, c, rows, cols); vh = (s["psi"] - 1) / s["mbar"]
    ck(f"pi3 recupero di v = {v} (entro il 20 %)", abs(vh / v - 1) <= 0.2, f"{vh:.4f}")
# fattori di cattura: psi corregge c
cc = np.exp(r.normal(0, 0.2, nb * nb)); y = r.poisson(10.0 * cc).astype(float)
s = S.group_stats(y, grp, cc, rows, cols)
ck("pi4 psi = 1 con fattori di cattura noti", abs(s["psi"] - 1) <= 0.03, f"{s['psi']:.4f}")
print(f"\nTOTALE: {sum(res)} PASS, {len(res) - sum(res)} FAIL su {len(res)}")
