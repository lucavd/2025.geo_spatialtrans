"""tools/F3.1_intervals.py - F3.1 fase A: misure vere per area, intervalli di rilevanza per tessuto (costruzione di S2.5: tools/S2_5_common.py,
DEC-060, DEC-061; quote limitate a [0, 1]), controprova CP-1 sul vero, ed esiti attesi dall'esplorazione (solo per scrivere le previsioni).
Uso: .venv/bin/python tools/F3.1_intervals.py <main_dir> <explore_csv> <out_dir>"""
import sys, os, json, glob, numpy as np, pandas as pd
sys.path.insert(0, "tools"); import S2_5_common as CM
main, explore_csv, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
MEAS = [("M-1", "share_assigned_exp", "bin", True), ("M-2", "umi_per_cell_exp", "nuc", False), ("M-3", "foreign_exp", "nuc", True),
        ("M-4", "rctd_frac_singlet", "bin", True), ("M-5", "rctd_agree_nucleus", "bin", True), ("M-6", "rctd_neighbour_same", "bin", True)]
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; SEED = 20261115
def load(d):
    r = {"tag": os.path.basename(d).replace("_k0", "").replace("_k1", "")}
    for lc, suf in [("label", ""), ("label_dimg", "_dimg"), ("label_shift", "_shift")]:
        f = f"{d}/b2c_{lc}.json"
        if os.path.exists(f): r.update({k + suf: v for k, v in json.load(open(f)).items() if isinstance(v, (int, float))})
    f = f"{d}/rctd.json"
    if os.path.exists(f): r.update({"rctd_" + k: v for k, v in json.load(open(f)).items() if isinstance(v, (int, float))})
    return r
real = pd.DataFrame([load(d) for d in sorted(glob.glob(f"{main}/real/*"))]); real["A"] = real.tag.str[:2]
real.to_csv(f"{out}/F3.1_real_areas.csv", index=False)
K = pd.read_csv(CM.KAP).set_index(["archetype", "unit"])
ex = pd.read_csv(explore_csv); ex["A"] = ex.tag.str[:2]
rows = []
for i, A in enumerate(ARCH):
    R = real[real.A == A]
    for j, (mid, col, q, frac) in enumerate(MEAS):
        v = R[col].to_numpy(float); rng = np.random.default_rng([SEED, i, j])
        ci = CM.corrected_interval(v, A, rng, side="real", q=q)
        if frac: ci["hi"] = min(1.0, ci["hi"])
        hw_boot = (ci["hi_boot"] - ci["lo_boot"]) / 2
        k_sim = float(K.loc[("all", "sim_5"), f"kappa_boot_{q}"])
        e = ex[ex.A == A].set_index("side")[col]
        shift = float(e["sim"] - e["real"]); exp_sim = ci["est"] + shift; hw_sim = k_sim * hw_boot
        v_exp = CM.verdict(exp_sim - hw_sim, exp_sim + hw_sim, ci["lo"], ci["hi"])
        row = {"id": mid, "measure": col, "archetype": A, "real_est": ci["est"], "lo": ci["lo"], "hi": ci["hi"], "half_width": (ci["hi"] - ci["lo"]) / 2,
               "cv": ci["cv"], "kappa": ci["kappa"], "cv_factor": ci["cv_factor"], "real_min": float(v.min()), "real_max": float(v.max()),
               "explore_sim_minus_real": shift, "expected_sim": exp_sim, "expected_sim_2se": hw_sim, "expected_verdict": v_exp}
        rows.append(row)
T = pd.DataFrame(rows); T.to_csv(f"{out}/F3.1_intervals.csv", index=False)
# CP-1 sul vero: nuclei spostati (M-3 sale), tipo del nucleo spostato (M-5 scende); differenze appaiate per area, 2 errori standard
cp = []
for A in ARCH:
    R = real[real.A == A]
    for mid, a, b, sgn in [("M-3", "foreign_exp_shift", "foreign_exp", 1), ("M-5", "rctd_agree_nucleus_shift", "rctd_agree_nucleus", -1)]:
        dlt = sgn * (R[a].to_numpy(float) - R[b].to_numpy(float)); m = dlt.mean(); se = dlt.std(ddof=1) / np.sqrt(len(dlt))
        cp.append({"archetype": A, "id": mid, "delta_mean": float(m), "se": float(se), "passes": bool(m - 2 * se > 0)})
pd.DataFrame(cp).to_csv(f"{out}/F3.1_CP1_real.csv", index=False)
pd.set_option("display.width", 250)
print(T[["id", "archetype", "real_est", "lo", "hi", "cv", "real_min", "real_max", "expected_sim", "expected_sim_2se", "expected_verdict"]].round(3).to_string(index=False))
print(pd.DataFrame(cp).round(4).to_string(index=False))
