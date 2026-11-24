"""tools/sonda/S5_diag_summary.py - S5 (docs/sonda/REGOLA_S5.md), diagnosi dell'arresto dopo il pilota, nessuna scelta di analisi: per ogni
esecuzione, blocchi di bin2cell con >= 20 UMI (UMI_min di tools/F3.1_rctd.R) e geni con >= 3 conteggi su quei blocchi sul pannello G_p
(results/sonda/S5/diag_rctd_panel.csv, tools/sonda/S5_diag_rctd.R), accanto alla quota di blocchi passati a RCTD sui geni completi nelle esecuzioni
originali (n_run / n_blocks di rctd.json di S3/S3b, sola lettura). Uscite: results/sonda/S5/diag_summary.csv, results/sonda/S5/fig_diag_rctd.png.
Uso: .venv/bin/python -B tools/sonda/S5_diag_summary.py"""
import os, sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.expanduser("~/2025.geo_spatialtrans/tools/sonda")); import S5_common as S
os.chdir(S.REPO); D = pd.read_csv("results/sonda/S5/diag_rctd_panel.csv"); rows = []
for r in D.itertuples():
    full = None
    if r.arm in S.SRC:
        p = f"{S.SRC[r.arm]}/{r.unit}/rctd.json"
        if os.path.exists(p): j = json.load(open(p)); full = (j["n_run"], j["n_blocks"])
    rows.append(dict(arm=r.arm, unit=r.unit, tissue=r.unit.split("_")[1], n_blocks=r.n_blocks, panel_blocks_umi20=r.n_blocks_umi20, panel_genes_ge3=r.n_genes_ge3,
                     threshold_genes=r.threshold_genes, panel_create_RCTD_ok=bool(r.create_ok), panel_message=r.message,
                     full_genes_rctd_n_run=None if full is None else full[0], full_genes_frac_blocks_run=None if full is None else full[0] / full[1]))
Q = pd.DataFrame(rows); Q.to_csv("results/sonda/S5/diag_summary.csv", index=False)
print(Q.groupby("arm").agg(n=("unit", "size"), ok=("panel_create_RCTD_ok", "sum"), med_blocks20=("panel_blocks_umi20", "median"),
                           med_full_frac=("full_genes_frac_blocks_run", "median")).to_string())
order = ["V", "T", "B", "P", "XB", "XV"]; fig, axs = plt.subplots(1, 2, figsize=(11, 4), facecolor="white")
for i, a in enumerate(order):
    q = Q[Q.arm == a]; x = i + np.linspace(-0.25, 0.25, len(q))
    axs[0].scatter(x, q.panel_blocks_umi20.fillna(-5), c=np.where(q.panel_create_RCTD_ok, "#2166ac", "#b2182b"), s=18)
    axs[1].scatter(x, q.panel_genes_ge3.fillna(-5), c=np.where(q.panel_create_RCTD_ok, "#2166ac", "#b2182b"), s=18)
axs[1].axhline(float(Q.threshold_genes.dropna().iloc[0]), color="k", ls="--", lw=1)
lab = {"V": "Vp", "T": "Tp", "B": "Bp", "P": "Pp", "XB": "XB", "XV": "XV"}
for ax, t in zip(axs, ("blocchi con >= 20 UMI sul pannello (UMI_min di RCTD)", "geni con >= 3 conteggi su quei blocchi (soglia di get_de_genes, tratteggio)")):
    ax.set_xticks(range(len(order))); ax.set_xticklabels([lab[a] for a in order]); ax.set_title(t, fontsize=9); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("S5, pilota: create.RCTD sul pannello G_p (blu = passa, rosso = errore; -5 = bin2cell fallito)", fontsize=10); fig.tight_layout()
fig.savefig("results/sonda/S5/fig_diag_rctd.png", dpi=150, facecolor="white"); print("fig ok")
