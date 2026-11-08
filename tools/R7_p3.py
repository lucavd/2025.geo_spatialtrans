"""tools/R7_p3.py — R7: previsione P-3 (controprova CP-2) e misura M-e.
Valore del tessuto = media delle 5 aree (come R4) oppure dei 5 pezzi disgiunti. Intervallo della media sulle aree:
R4_refs.csv per le quantità di R4; bootstrap a uno stadio (seme R7) per le due quantità dei nuclei. Correzione: kappa_bin / kappa_nuc
(regola R-2). Esito: |scarto| <= 0.5 x semi-ampiezza corretta dal lato dello scarto -> confermata, altrimenti smentita.
-> results/R7/pieces/R7_P3.csv"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "tools")
import R7_common as C, R4_common as R4
V = pd.read_csv("results/R7/pieces/R7_values_area_piece.csv"); refs = pd.read_csv("results/R4/R4_refs.csv")
K = pd.read_csv("results/R7/coverage/R7_kappa_final_v2.csv").set_index("archetype")
AIDX = {a: i for i, a in enumerate(["A1", "A2", "A3", "A4", "A5", "A6"])}
QT = {"mean_umi2": ("bin", "mean_umi2"), "mean_umi8": ("bin", "mean_umi8"), "frac_empty2": ("bin", "frac_empty2"), "mean_genes8": ("bin", "mean_genes8"),
      "rho_obs": ("bin", "rho_obs nucleus/ring UMI per bin (uncorrected)"), "nuc_area_median": ("nuc", None), "nuc_density": ("nuc", None)}
rows = []
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    for qi, (q, (t, rname)) in enumerate(QT.items()):
        d = V[(V.archetype == A) & (V.quantity == q)]
        va, vp = d.area_value.to_numpy(), d.piece_value.to_numpy()
        ea, ep = va.mean(), vp.mean()
        if rname is not None:
            rr = refs[(refs.archetype == A) & (refs.quantity == rname)].iloc[0]; est, lo, hi = rr.estimate, rr.lo, rr.hi
        else:
            est, lo, hi, _ = R4.boot_one_stage(va, rng=np.random.default_rng([C.SEED, AIDX[A], 900 + qi]))
        k = K.loc[A, f"kappa_{t}"]; lo2, hi2 = C.widen(est, lo, hi, k)
        shift = ep - ea; half = (hi2 - est) if shift > 0 else (est - lo2)
        ok = abs(shift) <= 0.5 * half
        rows.append(dict(archetype=A, quantity=q, qtype=t, value_areas=ea, value_pieces=ep, shift=shift, rel_shift=shift / ea, est=est, lo=lo, hi=hi,
                         kappa=k, lo_corr=lo2, hi_corr=hi2, bound=0.5 * half, shift_over_bound=abs(shift) / (0.5 * half),
                         verdict=("confermata" if ok else "smentita") if A != "A6" else ("G-2: identico" if shift == 0 else "G-2 FAIL"),
                         sd_areas=va.std(ddof=1), sd_pieces=vp.std(ddof=1), sd_ratio=vp.std(ddof=1) / va.std(ddof=1)))
P = pd.DataFrame(rows); P.to_csv("results/R7/pieces/R7_P3.csv", index=False)
print(P[["archetype", "quantity", "rel_shift", "shift_over_bound", "verdict", "sd_ratio"]].round(4).to_string(index=False))
print(P[P.archetype != "A6"].groupby("archetype").verdict.apply(lambda s: (s == "confermata").sum()))
