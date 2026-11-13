"""tools/F3.1_PH3_analysis.py - F3.1, analisi post hoc PH-3: semina a compartimenti (Q6) contro semina uniforme (U6, k = 0, stessi semi) e vero.
Per tessuto e misura: valore vero, U6, Q6, differenza appaiata per area, frazione dello scarto chiusa (U6 - Q6) / (U6 - vero), e posizione dell'intervallo
di Q6 (costruzione del lato simulato) rispetto all'intervallo di rilevanza del vero (esito descrittivo, non pre-registrato).
Uso: .venv/bin/python tools/F3.1_PH3_analysis.py <main_dir> <results_dir>"""
import sys, os, json, glob, numpy as np, pandas as pd
sys.path.insert(0, "tools"); import S2_5_common as CM
main, res = sys.argv[1:3]
MEAS = [("M-1", "share_assigned_exp", "bin", True), ("M-2", "umi_per_cell_exp", "nuc", False), ("M-3", "foreign_exp", "nuc", True),
        ("M-4", "rctd_frac_singlet", "bin", True), ("M-5", "rctd_agree_nucleus", "bin", True), ("M-6", "rctd_neighbour_same", "bin", True)]
def load(d):
    b = os.path.basename(d); r = {"tag": b.split("_k")[0]}
    f = f"{d}/b2c_label.json"
    if os.path.exists(f): r.update({k: v for k, v in json.load(open(f)).items() if isinstance(v, (int, float))})
    f = f"{d}/rctd.json"
    if os.path.exists(f): r.update({"rctd_" + k: v for k, v in json.load(open(f)).items() if isinstance(v, (int, float))})
    r["umi_per_bin"] = r.get("umi_reg", np.nan) / r.get("n_reg", np.nan); return r
U = pd.DataFrame([load(d) for d in sorted(glob.glob(f"{main}/sim/*_k0"))]).set_index("tag")
Q = pd.DataFrame([load(d) for d in sorted(glob.glob(f"{main}/simq6/*_k0"))]).set_index("tag")
Q.reset_index().to_csv(f"{res}/F3.1_simq6_runs.csv", index=False)
IV = pd.read_csv(f"{res}/F3.1_intervals.csv").set_index(["id", "archetype"])
rows = []
for i, A in enumerate(["A1", "A2", "A3", "A4", "A5", "A6"]):
    tg = [t for t in Q.index if t.startswith(A) and t in U.index]
    for j, (mid, col, q, frac) in enumerate(MEAS):
        u, w = U.loc[tg, col].astype(float), Q.loc[tg, col].astype(float); real = IV.loc[(mid, A), "real_est"]
        lo, hi = IV.loc[(mid, A), "lo"], IV.loc[(mid, A), "hi"]
        if mid == "M-4": lo, hi = max(0.0, min(lo, real - 0.01)), min(1.0, max(hi, real + 0.01))
        ci = CM.corrected_interval(w.to_numpy(), A, np.random.default_rng([20261119, i, j]), side="simulated", q=q)
        dl = (u - w).to_numpy()
        rows.append({"id": mid, "archetype": A, "real": real, "rel_lo": lo, "rel_hi": hi, "U6": u.mean(), "Q6": w.mean(),
                     "fraction_closed": (u.mean() - w.mean()) / (u.mean() - real) if u.mean() != real else np.nan,
                     "paired_diff": dl.mean(), "paired_se": dl.std(ddof=1) / np.sqrt(len(dl)), "Q6_lo": ci["lo"], "Q6_hi": min(1.0, ci["hi"]) if frac else ci["hi"],
                     "Q6_vs_real": CM.verdict(ci["lo"], min(1.0, ci["hi"]) if frac else ci["hi"], lo, hi), "n": len(tg)})
    rows.append({"id": "UMI per quadratino", "archetype": A, "U6": U.loc[tg, "umi_per_bin"].mean(), "Q6": Q.loc[tg, "umi_per_bin"].mean(), "n": len(tg)})
    rows.append({"id": "cellule per area", "archetype": A, "U6": U.loc[tg, "n_cells"].mean(), "Q6": Q.loc[tg, "n_cells"].mean(), "n": len(tg)})
T = pd.DataFrame(rows); T.to_csv(f"{res}/F3.1_PH3.csv", index=False)
V = pd.read_csv(f"{res}/F3.1_verdicts.csv").set_index(["id", "archetype"])
T["U6_verdict"] = [V.verdict.get((r.id, r.archetype), "") for r in T.itertuples()]
pd.set_option("display.width", 250)
print(T[["id", "archetype", "real", "U6", "Q6", "fraction_closed", "paired_diff", "paired_se", "U6_verdict", "Q6_vs_real"]].round(4).to_string(index=False))
