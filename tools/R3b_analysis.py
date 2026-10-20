#!/usr/bin/env python
"""tools/R3b_analysis.py — R3b: previsioni P-1...P-4, controprova CP-1, regola di decisione e cancelli, dai CSV di results/R3b.
Regole di giudizio (docs/steps/R3b.md, protocollo v2): stima +/- 2 errori standard (P-1, P-3, P-4, CP-1) o intervallo bootstrap
a due stadi (P-2), estremi delle frazioni tagliati a [0, 1]; tutta dentro l'intervallo di rilevanza = confermata, tutta fuori =
smentita, a cavallo = indeterminata. Intervalli letti dai file: P-1 da results/R2b/R2b_archetype_summary.csv (IC Poisson del
totale manuale / totale, righe di B-046, B-048); P-2, P-3, CP-1 da results/R2c/R2c_posthoc_ci.csv (righe «due stadi», B-130, B-132).
Uscite: results/R3b/R3b_predictions.csv, R3b_decision.csv, R3b_gates.csv.  Uso: .venv/bin/python tools/R3b_analysis.py
"""
import glob
import numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; RES = ROOT / "results/R3b"; M = Path("/mnt/micron/geo_spatialtrans/R3b")
TIS = ["A4", "A6"]; NAME = {"consensus": "CO", "stardist_he": "SD", "spaceranger": "SR"}


def verdict(lo, hi, a, b):
    if lo >= a and hi <= b: return "confermata"
    if hi < a or lo > b: return "smentita"
    return "indeterminata"


def main():
    sm = pd.read_csv(ROOT / "results/R2b/R2b_archetype_summary.csv").set_index("archetype")
    ci = pd.read_csv(ROOT / "results/R2c/R2c_posthoc_ci.csv"); ci = ci[ci.ci == "due stadi"].set_index("archetype")
    I1 = {A: (sm.ci_lo[A] / sm.density_manual[A], sm.ci_hi[A] / sm.density_manual[A]) for A in TIS}
    I2 = {A: (ci.lo[A], ci.hi[A]) for A in TIS}
    rows = []
    p12 = pd.read_csv(RES / "R3b_p1p2_summary.csv")
    for r in p12.itertuples():
        a, b = I1[r.archetype]; lo, hi = r.p1_ratio_mean - 2 * r.p1_se, r.p1_ratio_mean + 2 * r.p1_se
        rows.append(dict(id="P-1", archetype=r.archetype, mask=NAME[r.mask], est=r.p1_ratio_mean, lo=lo, hi=hi, two_se=2 * r.p1_se, int_lo=a, int_hi=b,
                         unit="finestra (4)", verdict=verdict(lo, hi, a, b)))
        a, b = I2[r.archetype]
        rows.append(dict(id="P-2", archetype=r.archetype, mask=NAME[r.mask], est=r.p2_est, lo=r.p2_lo, hi=r.p2_hi, two_se=np.nan, int_lo=a, int_hi=b,
                         unit="bootstrap due stadi", verdict=verdict(r.p2_lo, r.p2_hi, a, b)))
    w = pd.read_csv(RES / "R3b_windows.csv")
    per_seed = w.groupby(["source", "archetype", "seed"]).s1.mean().reset_index()      # peso uguale per finestra
    for (src, A), g in per_seed.groupby(["source", "archetype"]):
        if src == "SD_sr": continue
        m, se = g.s1.mean(), g.s1.std(ddof=1) / np.sqrt(len(g)); lo, hi = max(0, m - 2 * se), min(1, m + 2 * se); a, b = I2[A]
        rows.append(dict(id="CP-1" if src == "CP1" else "P-3", archetype=A, mask={"CO_own": "CO", "SD_own": "SD", "CP1": "manuale+CO"}[src],
                         est=m, lo=lo, hi=hi, two_se=2 * se, int_lo=a, int_hi=b, unit=f"seme ({len(g)})", verdict=verdict(lo, hi, a, b)))
    gr = pd.read_csv(RES / "R3b_gen_rows.csv"); gr = gr[gr.phimode == "own"]
    per_roi = gr.groupby(["archetype", "method", "roi"]).trunc_all.mean().reset_index()
    for (A, meth), g in per_roi.groupby(["archetype", "method"]):
        m, se = g.trunc_all.mean(), g.trunc_all.std(ddof=1) / np.sqrt(len(g)); lo, hi = max(0, m - 2 * se), min(1, m + 2 * se)
        rows.append(dict(id="P-4", archetype=A, mask=NAME[meth], est=m, lo=lo, hi=hi, two_se=2 * se, int_lo=0.0, int_hi=0.05,
                         unit=f"ROI ({len(g)})", verdict=verdict(lo, hi, 0.0, 0.05)))
    P = pd.DataFrame(rows); P.to_csv(RES / "R3b_predictions.csv", index=False)
    # regola di decisione
    v = {(r.id, r.archetype, r.mask): r for r in P.itertuples()}
    dec = []
    for A in TIS:
        cp_sm = v[("CP-1", A, "manuale+CO")].verdict == "smentita"
        elig = [k for k in ("CO", "SD") if v[("P-1", A, k)].verdict != "smentita" and v[("P-2", A, k)].verdict != "smentita"]
        if len(elig) == 1: choice, why = elig[0], "una sola maschera eleggibile"
        elif not elig: choice, why = "nessuna", "nessuna eleggibile: decide Luca"
        else:
            p3 = [k for k in elig if v[("P-3", A, k)].verdict == "confermata"]
            if not cp_sm and len(p3) == 1: choice, why = p3[0], "P-3 confermata per una sola"
            else:
                d = {k: abs(np.log(v[("P-1", A, k)].est)) for k in elig}; choice = min(d, key=d.get)
                why = f"|log rapporto P-1| minimo ({', '.join(f'{k} {d[k]:.3f}' for k in elig)}); CP-1 smentita: {cp_sm}"
        dec.append(dict(archetype=A, eligible=";".join(elig), cp1_refuted=cp_sm, choice=choice, reason=why))
    D = pd.DataFrame(dec); D.to_csv(RES / "R3b_decision.csv", index=False)
    # cancelli
    g4 = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(M / "check/G4_*.csv")))])
    g4ok = bool((g4.labels_consecutive & g4.count_identity & g4.rerun_identical).all() and (g4[["n_empty", "sd_pixels_bad", "sd_pixels_bad_indep", "sr_pixels_bad", "sr_npx_bad"]] == 0).all().all())
    inv = pd.concat([pd.read_csv(RES / "R3b_gen_rows.csv"), pd.read_csv(RES / "R3b_cp1_rows.csv")])
    invok = bool((inv[["c3_geos", "c3_sf", "no_pairs", "n_invalid", "n_empty"]] == 0).all().all())
    g3 = [pd.read_csv(RES / f) for f in ("R3b_G3_p1p2.csv", "R3b_G3_windows.csv")]
    g3r = pd.read_csv(RES / "R3b_G3_R.csv"); nucr = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(M / "regress/nuclei_stardist_*.csv")))])
    phs = pd.read_csv(RES / "R3b_G3_phi_sr.csv")
    g3ok = bool(all(x.equal.all() for x in g3) and g3r[["gen", "cells", "nuc", "s14_cells"]].astype(bool).all().all()
                and g3r.s14_terr.dropna().astype(bool).all() and (nucr.n_cols_equal == nucr.n_cols).all() and (g4.cp2_sr_equal & g4.cp2_sd_equal).all()
                and np.allclose(np.round(phs.phi.to_numpy(), 3), [0.816, 0.180]))
    G = pd.DataFrame([dict(gate="G-3 regressione", passed=g3ok), dict(gate="G-4 invarianti CO", passed=g4ok),
                      dict(gate="G-4 invarianti nucleo simulato (contenimento, sovrapposizioni, validita')", passed=invok)])
    G.to_csv(RES / "R3b_gates.csv", index=False)
    pd.set_option("display.width", 200)
    print(P.round(4).to_string(index=False)); print(D.to_string(index=False)); print(G.to_string(index=False))


if __name__ == "__main__":
    main()
