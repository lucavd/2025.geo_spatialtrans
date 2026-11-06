"""tools/S2.4_R1_calib.py — S2.4, regola R-1 (docs/steps/S2.4.md): calibrazione estesa del passo 1 della parte fine delle strisce.
Sezioni sintetiche del generatore di S2.0 PH-5 (job_syn: 800 x 800, struttura del tessuto 0.2 e 0.4, fattori di SD 0.2, 4 semi),
passo 1 vero rho in {0, -0.1, -0.2, -0.3, -0.4, -0.45, -0.5}. Curva = media di lag1_repro_hp per rho (righe e colonne, ampiezze e
semi insieme), interpolata linearmente; passo 1 per il generatore = inverso della curva nel valore della sezione (S2.0 PH-5, 10
iterazioni), limitato a [-0.5, 0]. Controllo di continuita': per rho in {0, -0.2, -0.3} le righe coincidono con S2.0_PH5_calib.csv.
Uscite: results/S2.4/S2.4_R1_calib.csv, S2.4_R1_curve.csv, S2.4_R1_generator.csv. Uso: .venv/bin/python tools/S2.4_R1_calib.py <nproc>"""
import sys, importlib.util
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "tools")
_sp = importlib.util.spec_from_file_location("ph5", "tools/S2.0_ph5_stripes_hp.py"); PH = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(PH)
RHOS = [0.0, -0.1, -0.2, -0.3, -0.4, -0.45, -0.5]
def job(a): return PH.job_syn(a)      # involucro: il modulo caricato con importlib non si serializza
if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    jobs = [(amp, rho, s) for amp in (0.2, 0.4) for rho in RHOS for s in range(4)]
    with Pool(nproc) as p: rows = sum(p.map(job, jobs, chunksize=1), [])
    Cb = pd.DataFrame(rows); Cb.to_csv("results/S2.4/S2.4_R1_calib.csv", index=False)
    old = pd.read_csv("results/S2.0/S2.0_PH5_calib.csv"); key = ["amp", "rho_true", "seed", "axis"]
    m = old.merge(Cb, on=key, suffixes=("_s20", "_s24"))
    cont = float(np.max(np.abs(m.lag1_repro_hp_s20 - m.lag1_repro_hp_s24))) if len(m) else np.nan
    cur = Cb.groupby("rho_true").agg(lag1_meas=("lag1_repro_hp", "mean"), lag1_meas_sd=("lag1_repro_hp", "std"), n=("lag1_repro_hp", "size"),
                                     sd_meas=("sd_repro_hp", "mean")).reset_index().sort_values("rho_true")
    cur.to_csv("results/S2.4/S2.4_R1_curve.csv", index=False)
    mono = bool(np.all(np.diff(cur.lag1_meas.to_numpy()) > 0))
    sec = pd.read_csv("results/S2.0/S2.0_PH5_sections.csv"); sec = sec[sec.n_iter == 10]
    gen = []
    for _, r in sec.iterrows():
        x = float(r.lag1_repro_hp); rho = float(np.interp(x, cur.lag1_meas.to_numpy(), cur.rho_true.to_numpy()))
        gen.append({"key": r.key, "axis": r.axis, "lag1_measured": x, "rho_generator": float(np.clip(rho, -0.5, 0.0)),
                    "below_curve": bool(x < cur.lag1_meas.min()), "above_curve": bool(x > cur.lag1_meas.max())})
    G = pd.DataFrame(gen); G.to_csv("results/S2.4/S2.4_R1_generator.csv", index=False)
    print(f"continuita' con S2.0 (max |diff| lag1 su {len(m)} righe): {cont:.3g}; curva monotona: {mono}")
    print(cur.round(4).to_string(index=False)); print(G.round(4).to_string(index=False))
