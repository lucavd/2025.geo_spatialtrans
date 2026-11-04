"""tools/S2.0_ph23.py — S2.0, analisi post hoc PH-2 (prototipo con i fattori di riga e colonna misurati) e PH-3 (correlazione dei
residui fra bin vicini dello stesso nucleo). Tabella results/S2.0/S2.0_PH23_disp.csv e figura S2.0_F8_disp_corr.png.
Uso: .venv/bin/python tools/S2.0_ph23.py <outdir>"""
import sys, os, json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = sys.argv[1]; RES = os.environ.get("S20_RES", "results/S2.0"); SEED, B = 20261103, 2000
def load(d): return pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/{d}/*.json"))])
a, n = load("disp"), load("disp_nbr")
m = a[["tag", "archetype", "ratio", "v"]].merge(n.drop(columns=["archetype"]), on="tag", suffixes=("_pre", ""))
for k in ("h1", "h2", "v1", "v2"): m[f"exc_{k}"] = m[f"rho_{k}_obs"] - m[f"rho_{k}_sim"]
m["exc_lag1"] = 0.5 * (m.exc_h1 + m.exc_v1); m["exc_lag2"] = 0.5 * (m.exc_h2 + m.exc_v2)
rng = np.random.default_rng(SEED); rows = []
for A_, t in m.groupby("archetype"):
    r = {"archetype": A_, "n_roi": len(t)}
    for q in ("ratio_pre", "v_pre", "ratio", "v", "exc_lag1", "exc_lag2"):
        x = t[q].to_numpy(float); bs = np.array([x[rng.integers(0, len(x), len(x))].mean() for _ in range(B)])
        r[q] = x.mean(); r[q + "_lo"], r[q + "_hi"] = np.percentile(bs, [2.5, 97.5])
    rows.append(r)
T = pd.DataFrame(rows); T.to_csv(f"{RES}/S2.0_PH23_disp.csv", index=False)
print(T[["archetype", "ratio_pre", "ratio", "ratio_lo", "ratio_hi", "v_pre", "v", "v_lo", "v_hi", "exc_lag1", "exc_lag1_lo", "exc_lag1_hi", "exc_lag2", "exc_lag2_lo", "exc_lag2_hi"]].round(3).to_string(index=False))
AN = {"A1": "intestine", "A2": "tumour", "A3": "stroma", "A4": "lymph node", "A5": "brain", "A6": "heart"}
plt.rcParams.update({"figure.facecolor": "white", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
for j, (A_, t) in enumerate(m.groupby("archetype")):
    col = plt.cm.tab10(j)
    for i, d in enumerate((1, 2)):
        ax[0].scatter(np.full(len(t), d) + (j - 2.5) * 0.06, 0.5 * (t[f"rho_h{d}_obs"] + t[f"rho_v{d}_obs"]), color=col, s=12, label=AN[A_] if d == 1 else None)
        ax[0].scatter(np.full(len(t), d) + (j - 2.5) * 0.06, 0.5 * (t[f"rho_h{d}_sim"] + t[f"rho_v{d}_sim"]), facecolors="none", edgecolors=col, s=12)
    r = T[T.archetype == A_].iloc[0]
    ax[1].errorbar([1 + (j - 2.5) * 0.06, 2 + (j - 2.5) * 0.06], [r.exc_lag1, r.exc_lag2], yerr=[[r.exc_lag1 - r.exc_lag1_lo, r.exc_lag2 - r.exc_lag2_lo], [r.exc_lag1_hi - r.exc_lag1, r.exc_lag2_hi - r.exc_lag2]], fmt="o-", color=col, capsize=2, lw=0.8)
ax[0].set_xticks([1, 2]); ax[0].set_xticklabels(["2 µm (adjacent)", "4 µm"]); ax[0].set_ylabel("residual correlation, same nucleus"); ax[0].axhline(0, color="k", lw=0.5)
ax[0].set_title("filled = observed, open = simulated (with row/column factors)"); ax[0].legend(frameon=False, fontsize=7, ncol=2)
ax[1].set_xticks([1, 2]); ax[1].set_xticklabels(["2 µm (adjacent)", "4 µm"]); ax[1].axhline(0, color="k", lw=0.5); ax[1].set_ylabel("observed - simulated")
ax[1].set_title("excess correlation: short range (mean over ROIs, 95% CI)")
fig.tight_layout(); fig.savefig("reports/figures/S2.0/S2.0_F8_disp_corr.png", dpi=130)
