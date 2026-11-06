"""S2.2: punteggi della gara (regola R-1), vincitore per tessuto, previsioni P-1...P-5, controprove CP-1 e CP-2.
Uso: .venv/bin/python tools/S2.2_score.py   (legge /mnt/micron/geo_spatialtrans/S2.2/race/*; scrive results/S2.2/S2.2_*.csv)
Pre-registrazione: docs/steps/S2.2.md (regole R-1, R-2; previsioni; controprova)."""
import glob, json, os, numpy as np, pandas as pd
M = os.environ.get("S22_RACE_OUT", "/mnt/micron/geo_spatialtrans/S2.2/race"); RES = os.environ.get("S22_RES", "results/S2.2"); os.makedirs(RES, exist_ok=True)
METRICS = {"m1": "m1", "m2": "m2", "m3": "m3", "m4": "m4_ratio"}
MAIN = ["atlas_M", "atlas_G", "nuc_M", "nuc_G", "exp_M", "exp_G"]          # ordine di preferenza nei pareggi (R-1)
SEED_B = 20261107; B = 2000
rel = pd.read_csv(f"{RES}/S2.2_E1_relevance.csv").set_index("archetype")
race = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{M}/*_*.csv"))], ignore_index=True)
real = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{M}/*_real.json"))]).set_index(["archetype", "roi_id"])
err_rows = race[race.get("error").notna()] if "error" in race else race.iloc[0:0]
race = race[race.get("error").isna()] if "error" in race else race
# media delle repliche e errore Monte Carlo
agg = race.groupby(["archetype", "roi_id", "arm"])[list(METRICS.values())].agg(["mean", "std"])
rows = []
for (a, roi, arm), r in agg.iterrows():
    out = dict(archetype=a, roi_id=roi, arm=arm)
    for k, col in METRICS.items():
        rv = real.loc[(a, roi), col]
        out[f"sim_{k}"] = r[(col, "mean")]; out[f"real_{k}"] = rv; out[f"mc_sd_{k}"] = r[(col, "std")]
        out[f"e_{k}"] = (r[(col, "mean")] - rv) / rv
    out["S"] = np.mean([abs(out[f"e_{k}"]) for k in METRICS])
    rows.append(out)
per = pd.DataFrame(rows)
per.to_csv(f"{RES}/S2.2_race_per_roi.csv", index=False)


def verdict(est, se, half):
    lo, hi = est - 2 * se, est + 2 * se
    if lo >= -half and hi <= half: return "confermata"
    if hi < -half or lo > half: return "smentita"
    return "indeterminata"


summ, preds, cps, winners = [], [], [], []
for a, g in per.groupby("archetype"):
    arms = [x for x in MAIN if x in set(g.arm)]
    W = g.pivot(index="roi_id", columns="arm", values="S")
    Sa = W[arms].mean()
    best = Sa.idxmin()
    second = Sa.drop(best).idxmin()
    for x in arms:
        d = W[x] - W[best]; n = d.notna().sum()
        se = d.std(ddof=1) / np.sqrt(n) if x != best else 0.0
        summ.append(dict(archetype=a, arm=x, S=Sa[x], diff_vs_best=d.mean(), se_diff=se, n_roi=int(n), is_second=bool(x == second),
                         tie_with_best=bool(x == best or (x == second and d.mean() < 2 * se))))
    # R-1 come scritta: pareggio solo fra il secondo e il primo; allora si sceglie per ordine fra i due (correzione dopo la revisione G-5)
    tie = [x for x in (best, second) if summ[-len(arms) + arms.index(x)]["tie_with_best"]]
    win = [x for x in MAIN if x in tie][0]
    winners.append(dict(archetype=a, best_score=best, winner=win, tie_set=";".join(tie), rule="score" if win == best else "tie-break order"))
    gw = g[g.arm == win].set_index("roi_id")
    for i, k in enumerate(METRICS):
        e = gw[f"e_{k}"]; est = e.mean(); se = e.std(ddof=1) / np.sqrt(e.notna().sum())
        half = rel.loc[a, f"rel_{k}"]
        preds.append(dict(id=f"P-{i+1}", archetype=a, arm=win, quantity=k, estimate=est, se=se, lo=est - 2 * se, hi=est + 2 * se,
                          relevance=half, verdict=verdict(est, se, half), real_mean=gw[f"real_{k}"].mean(), sim_mean=gw[f"sim_{k}"].mean(),
                          mc_sd_max=gw[f"mc_sd_{k}"].max()))
    for nul in ("null_pooled", "null_other"):
        if nul not in W:
            continue
        d = W[nul] - W[win]; se = d.std(ddof=1) / np.sqrt(d.notna().sum())
        cps.append(dict(id="CP-1", archetype=a, null=nul, winner=win, S_null=W[nul].mean(), S_winner=W[win].mean(), diff=d.mean(), se=se,
                        null_loses=bool(d.mean() - 2 * se > 0)))
pd.DataFrame(summ).to_csv(f"{RES}/S2.2_race_scores.csv", index=False)
pd.DataFrame(winners).to_csv(f"{RES}/S2.2_race_winners.csv", index=False)
pd.DataFrame(cps).to_csv(f"{RES}/S2.2_CP1.csv", index=False)

# P-5 e CP-2 (intestino): mediana sulle coppie (area, tipo, gene), come R4 PH-2; bootstrap sulle aree
win_a1 = {w["archetype"]: w["winner"] for w in winners}.get("A1")
p5rows = []
if win_a1:
    files = sorted(glob.glob(f"{M}/A1_*_p5.npz"))
    def pairs(key_arm):
        per_roi = []
        for f in files:
            z = np.load(f); reps = {}
            for kk in z.files:
                parts = kk.split("|")
                if key_arm == "real" and parts[0] == "real":
                    reps.setdefault(0, []).append(z[kk])
                elif parts[0] == key_arm:
                    reps.setdefault(int(parts[1]), []).append(z[kk])
            per_roi.append({r: np.concatenate(v) for r, v in reps.items()})
        return per_roi
    real_p = pairs("real")
    rng = np.random.default_rng(SEED_B)
    for arm in (win_a1, win_a1 + "_nosub"):
        sim_p = pairs(arm)
        def stat(ix):
            rv = np.median(np.concatenate([real_p[i][0] for i in ix]))
            sv = np.mean([np.median(np.concatenate([sim_p[i][r] for i in ix])) for r in sim_p[0]])
            return sv, rv
        sv, rv = stat(range(len(files)))
        bs = np.array([stat(rng.integers(0, len(files), len(files))) for _ in range(B)])
        dd = bs[:, 0] - bs[:, 1]
        p5rows.append(dict(arm=arm, sim=sv, real=rv, diff=sv - rv, diff_lo=np.percentile(dd, 2.5), diff_hi=np.percentile(dd, 97.5),
                           sim_lo=np.percentile(bs[:, 0], 2.5), sim_hi=np.percentile(bs[:, 0], 97.5), n_pairs_real=int(sum(len(x[0]) for x in real_p))))
    pd.DataFrame(p5rows).to_csv(f"{RES}/S2.2_P5_CP2.csv", index=False)
pd.DataFrame(preds).to_csv(f"{RES}/S2.2_predictions.csv", index=False)
print(pd.DataFrame(winners).to_string(index=False))
print(pd.DataFrame(preds)[["id", "archetype", "arm", "estimate", "lo", "hi", "relevance", "verdict"]].round(4).to_string(index=False))
print(pd.DataFrame(cps).round(4).to_string(index=False))
if p5rows: print(pd.DataFrame(p5rows).round(4).to_string(index=False))
if len(err_rows): print("ERRORI:", err_rows[["archetype", "roi_id", "arm", "error"]].drop_duplicates().to_string(index=False))
