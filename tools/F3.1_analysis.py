"""tools/F3.1_analysis.py - F3.1 fase B: esiti delle previsioni P-1...P-6, controprove CP-1...CP-4, regola R-1, riassunto della verita'.
Intervalli del vero: quelli fissati nella fase A (results/F3.1/F3.1_intervals.csv), con il margine minimo di DEC-086 per M-4.
Lato simulato: valore per area = media dei 2 semi; intervallo con la costruzione di S2.5 e il fattore sim_5 (B = 2000, seme 20261117).
Uso: .venv/bin/python tools/F3.1_analysis.py <main_dir> <out_dir>"""
import sys, os, json, glob, numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, "tools"); import S2_5_common as CM
main, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
MEAS = [("M-1", "share_assigned_exp", "bin", True), ("M-2", "umi_per_cell_exp", "nuc", False), ("M-3", "foreign_exp", "nuc", True),
        ("M-4", "rctd_frac_singlet", "bin", True), ("M-5", "rctd_agree_nucleus", "bin", True), ("M-6", "rctd_neighbour_same", "bin", True)]
ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; SEED = 20261117; FLOOR_M4 = 0.01
NAMES = {"A1": "intestino", "A2": "tumore del colon", "A3": "stroma", "A4": "linfonodo", "A5": "cervello", "A6": "cuore"}
def load(d):
    b = os.path.basename(d); r = {"run": b, "tag": b.split("_k")[0], "k": int(b.split("_k")[1]) if "_k" in b else -1}
    for lc, suf in [("label", ""), ("label_dimg", "_dimg"), ("label_shift", "_shift")]:
        f = f"{d}/b2c_{lc}.json"
        if os.path.exists(f): r.update({k + suf: v for k, v in json.load(open(f)).items() if isinstance(v, (int, float))})
    f = f"{d}/rctd.json"
    if os.path.exists(f): r.update({"rctd_" + k: v for k, v in json.load(open(f)).items() if isinstance(v, (int, float))})
    r["umi_per_bin"] = r.get("umi_reg", np.nan) / r.get("n_reg", np.nan)
    return r
def table(pattern):
    t = pd.DataFrame([load(d) for d in sorted(glob.glob(pattern)) if os.path.isdir(d)]); t["A"] = t.tag.str[:2]; return t
real = table(f"{main}/real/*"); sim = table(f"{main}/sim/*_k*"); nuc = table(f"{main}/simnuc/*_k0")
for nm, t in [("real", real), ("sim", sim), ("simnuc", nuc)]: t.to_csv(f"{out}/F3.1_{nm}_runs.csv", index=False)
num = sim.select_dtypes("number").columns.drop("k")
sim_a = sim.groupby(["A", "tag"])[list(num)].mean().reset_index()            # media dei 2 semi per area
sim_a.to_csv(f"{out}/F3.1_sim_areas.csv", index=False)
print("runs: real", len(real), "sim", len(sim), "simnuc", len(nuc), "| areas sim", len(sim_a), flush=True)
IV = pd.read_csv(f"{out}/F3.1_intervals.csv").set_index(["id", "archetype"])
rows = []
for i, A in enumerate(ARCH):
    for j, (mid, col, q, frac) in enumerate(MEAS):
        lo, hi, est = IV.loc[(mid, A), "lo"], IV.loc[(mid, A), "hi"], IV.loc[(mid, A), "real_est"]
        if mid == "M-4": lo, hi = max(0.0, min(lo, est - FLOOR_M4)), min(1.0, max(hi, est + FLOOR_M4))
        v = sim_a.loc[sim_a.A == A, col].to_numpy(float)
        if len(v) < 2 or not np.all(np.isfinite(v)):
            rows.append({"id": mid, "archetype": A, "verdict": "non valutabile", "n_sim": int(np.isfinite(v).sum())}); continue
        ci = CM.corrected_interval(v, A, np.random.default_rng([SEED, i, j]), side="simulated", q=q)
        s_lo, s_hi = ci["lo"], (min(1.0, ci["hi"]) if frac else ci["hi"])
        vd = CM.verdict(s_lo, s_hi, lo, hi)
        rows.append({"id": mid, "measure": col, "archetype": A, "tissue": NAMES[A], "real_est": est, "rel_lo": lo, "rel_hi": hi, "sim_est": ci["est"],
                     "sim_lo": s_lo, "sim_hi": s_hi, "verdict": vd, "direction": ("eccesso" if ci["est"] > est else "difetto"),
                     "gap": ci["est"] - est, "gap_in_halfwidths": (ci["est"] - est) / ((hi - lo) / 2), "n_sim": len(v)})
V = pd.DataFrame(rows); V.to_csv(f"{out}/F3.1_verdicts.csv", index=False)
# CP-1 sul simulato (differenze appaiate per area; media dei semi)
cp1 = []
for A in ARCH:
    S = sim_a[sim_a.A == A]
    for mid, a, b, sgn in [("M-3", "foreign_exp_shift", "foreign_exp", 1)]:
        dl = sgn * (S[a] - S[b]).to_numpy(float); m, se = dl.mean(), dl.std(ddof=1) / np.sqrt(len(dl))
        cp1.append({"archetype": A, "id": mid, "side": "sim", "delta_mean": m, "se": se, "passes": bool(m - 2 * se > 0)})
cpr = pd.read_csv(f"{out}/F3.1_CP1_real.csv"); cpr["side"] = "real"
CP1 = pd.concat([cpr, pd.DataFrame(cp1)], ignore_index=True); CP1.to_csv(f"{out}/F3.1_CP1.csv", index=False)
# nota: M-5 con il tipo spostato si calcola nello script di RCTD (agree_nucleus_shift) anche sul simulato
cp1b = []
for A in ARCH:
    S = sim_a[sim_a.A == A]; dl = (S.rctd_agree_nucleus - S.rctd_agree_nucleus_shift).to_numpy(float)
    m, se = dl.mean(), dl.std(ddof=1) / np.sqrt(len(dl)); cp1b.append({"archetype": A, "id": "M-5", "side": "sim", "delta_mean": m, "se": se, "passes": bool(m - 2 * se > 0)})
CP1 = pd.concat([CP1, pd.DataFrame(cp1b)], ignore_index=True); CP1.to_csv(f"{out}/F3.1_CP1.csv", index=False)
# CP-2: le misure senza verita' seguono la verita'? (10 esecuzioni per tessuto, Spearman)
cp2 = []
for A in ARCH:
    S = sim[sim.A == A]
    r3 = spearmanr(S.foreign_exp, S.truth_foreign_type_exp).statistic; r5 = spearmanr(S.rctd_agree_nucleus, S.rctd_truth_acc_first).statistic
    cp2.append({"archetype": A, "rho_M3_truth_foreign": r3, "rho_M5_truth_acc": r5, "n": len(S)})
CP2 = pd.DataFrame(cp2); CP2.to_csv(f"{out}/F3.1_CP2.csv", index=False)
# CP-3: profili dal solo nucleo (stessi semi di k = 0; confronto appaiato per area con il principale k = 0)
cp3 = []
s0 = sim[sim.k == 0].set_index("tag")
for A in ["A1", "A2", "A4", "A5", "A6"]:
    N = nuc[nuc.A == A].set_index("tag"); tags = N.index.intersection(s0.index)
    for mid, col, q, frac in MEAS:
        m0, mn = s0.loc[tags, col].astype(float), N.loc[tags, col].astype(float)
        real_est = IV.loc[(mid, A), "real_est"]; gap = m0.mean() - real_est
        closed = (m0.mean() - mn.mean()) / gap if gap != 0 else np.nan
        dl = (m0 - mn).to_numpy(); se = dl.std(ddof=1) / np.sqrt(len(dl))
        cp3.append({"archetype": A, "id": mid, "real": real_est, "sim_k0": m0.mean(), "sim_nuc": mn.mean(), "gap_k0": gap, "fraction_closed": closed,
                    "paired_diff": dl.mean(), "paired_se": se, "n": len(tags)})
    cp3.append({"archetype": A, "id": "UMI per quadratino", "real": np.nan, "sim_k0": s0.loc[tags, "umi_per_bin"].mean(), "sim_nuc": N.loc[tags, "umi_per_bin"].mean()})
CP3 = pd.DataFrame(cp3); CP3.to_csv(f"{out}/F3.1_CP3.csv", index=False)
cp3_pass = CP3[(CP3.id == "M-3") & CP3.archetype.isin(["A1", "A2", "A4", "A5"])]
cp3_ok = int((cp3_pass.fraction_closed >= 0.5).sum())
# R-1
cp1_ok = {(r.archetype, r.id): True for _, r in CP1.iterrows()}
for _, r in CP1.iterrows(): cp1_ok[(r.archetype, r.id)] = cp1_ok.get((r.archetype, r.id), True) and bool(r.passes)
ncell = (sim_a.groupby("A").n_cells.mean() / real.groupby("A").n_cells.mean())
att = []
for _, r in V[(V.verdict == "smentita") & (V.archetype != "A6")].iterrows():
    A, mid = r.archetype, r.id
    if mid in ("M-3", "M-5") and not cp1_ok.get((A, mid), True):
        att.append({**r.to_dict(), "limit": "escluso (CP-1)"}); continue
    lim = "non attribuito"
    c3 = CP3[(CP3.archetype == A) & (CP3.id == mid)]
    if A in ("A1", "A2", "A4", "A5") and len(c3) and np.isfinite(c3.fraction_closed.iloc[0]) and c3.fraction_closed.iloc[0] >= 0.5:
        lim = "profili dai nuclei con la corona"
    elif mid in ("M-1", "M-2") and np.sign(r.gap) == np.sign(ncell[A] - 1) and not (0.9 <= ncell[A] <= 1.1):
        lim = "densita' dei preset"
    elif mid in ("M-4", "M-6") and r.gap < 0:
        lim = "disposizione dei tipi"
    att.append({**r.to_dict(), "limit": lim, "cells_ratio": ncell[A]})
R1 = pd.DataFrame(att); R1.to_csv(f"{out}/F3.1_R1.csv", index=False)
summ = R1[~R1.limit.isin(["escluso (CP-1)"])].groupby("limit").agg(n=("id", "size"), mean_gap_hw=("gap_in_halfwidths", lambda x: float(np.mean(np.abs(x))))).sort_values(["n", "mean_gap_hw"], ascending=False)
summ.to_csv(f"{out}/F3.1_R1_summary.csv")
# verita' (descrittiva)
tcols = ["truth_purity_exp", "truth_recall_exp", "truth_foreign_type_exp", "truth_purity_nuc", "truth_recall_nuc", "rctd_truth_acc_first", "rctd_truth_acc_nucleus", "foreign_profile"]
TR = pd.concat([sim.groupby("A")[tcols].mean().assign(arm="main"), nuc.groupby("A")[tcols].mean().assign(arm="nuc")]).reset_index()
TR.to_csv(f"{out}/F3.1_truth.csv", index=False)
pd.set_option("display.width", 250)
print(V[["id", "archetype", "real_est", "rel_lo", "rel_hi", "sim_est", "sim_lo", "sim_hi", "verdict", "direction"]].round(4).to_string(index=False))
print(V.groupby(["id", "verdict"]).size().unstack(fill_value=0).to_string())
print(CP1.round(4).to_string(index=False)); print(CP2.round(3).to_string(index=False)); print(CP3.round(4).to_string(index=False))
print("CP-3 M-3 closed >= 1/2 in", cp3_ok, "of 4 (A1, A2, A4, A5)")
print(R1[["id", "archetype", "gap", "limit"]].round(4).to_string(index=False)); print(summ.to_string()); print(ncell.round(3).to_string())
print(TR.round(3).to_string(index=False))
