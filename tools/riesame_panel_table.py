"""tools/riesame_panel_table.py - riesame, Parte B: tabella per area di tutti i bracci (V, T, S, B) e insiemi di geni (panel, common), cancello
G-4 della regola (docs/riesame/REGOLA_B.md) e tempi della catena.
Ingressi: /mnt/micron/geo_spatialtrans/riesame/B/<arm>/measures/*.json (tools/riesame_panel_measures.py, etichette trasferite = labels tl);
results/F3.3/F3.3_area_values.csv (V e T sui geni comuni: righe IN e OUT/L2 del tumore e dello stroma, riusate come sono).
Uscite: results/riesame/B_area_values.csv (colonne di F3.3_area_values.csv piu' genes, arm, origin, n_genes), results/riesame/B_G4.csv,
results/riesame/B_timings.csv. Uso: .venv/bin/python -B tools/riesame_panel_table.py"""
import os, glob, json, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans"))
BR = "/mnt/micron/geo_spatialtrans/riesame/B"; OUT = "results/riesame"; os.makedirs(OUT, exist_ok=True)
F33 = ["unit", "run", "A", "src", "set", "level", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M2r", "M3r", "M3u", "u2"]
MEAS = ["M1", "M2", "M3", "M4", "M5", "M6", "M7"]
def load(arm):
    rows = [r for f in sorted(glob.glob(f"{BR}/{arm}/measures/*.json")) for r in json.load(open(f))]
    return pd.DataFrame(rows)
av = pd.read_csv("results/F3.3/F3.3_area_values.csv")
reuse = av[av.A.isin(["A2", "A3"]) & ((av.src == "IN") | ((av.src == "OUT") & (av.level == "L2")))].copy()
reuse["genes"] = "common"; reuse["arm"] = np.where(reuse.src == "IN", "V", "T"); reuse["origin"] = "F3.3"; reuse["n_genes"] = 18024
parts = [reuse]
for arm in ("V", "T", "S", "B"):
    d = load(arm)
    if not len(d): continue
    d = d[d.labels == "tl"].copy(); d["origin"] = "chain"; parts.append(d)
U = pd.concat(parts, ignore_index=True)[F33 + ["genes", "arm", "origin", "n_genes"]]
U = U.sort_values(["arm", "genes", "A", "unit", "run"]).reset_index(drop=True)
U.to_csv(f"{OUT}/B_area_values.csv", index=False)
exp = {("V", "panel"): 10, ("V", "common"): 10, ("T", "panel"): 19, ("T", "common"): 19, ("S", "panel"): 20, ("S", "common"): 20, ("B", "panel"): 20, ("B", "common"): 20}
cnt = U.groupby(["arm", "genes"]).size().to_dict()
print("rows per arm x genes (observed / expected):", {f"{k[0]}-{k[1]}": f"{cnt.get(k, 0)}/{v}" for k, v in exp.items()})
print("missing values in M1..M7:", int(U[MEAS].isna().sum().sum()))
# ---------- G-4: un'area per sezione rifatta con la catena, contro F3.3 ----------
g = load("G4"); g4 = []
if len(g):
    g = g[g.labels == "tl"]
    for _, r in g.iterrows():
        f = {"IN": av[(av.run == f"{r.unit}_cg") & (av.src == "IN")], "OUT": av[(av.run == r.run) & (av.src == "OUT")],
             "SIM": av[(av.run == r.run) & (av.src == "SIM")]}[r.src]
        assert len(f) == 1, r.run
        f = f.iloc[0]
        dev = {m: abs(r[m] - f[m]) / max(abs(f[m]), 1e-12) for m in MEAS}
        a = pd.read_csv(f"{BR}/G4/runs/{r.run}_{r.genes}/rctd_blocks.csv")
        src_dir = {"IN": f"/mnt/micron/geo_spatialtrans/F3.3/runs/{r.unit}_cg", "OUT": f"/mnt/micron/geo_spatialtrans/F3.3/runs/{r.run}",
                   "SIM": f"/mnt/micron/geo_spatialtrans/F3.1/main/sim/{r.run}"}[r.src]
        b = pd.read_csv(f"{src_dir}/rctd_blocks.csv")
        same_blocks = a.shape == b.shape and all(a[c].equals(b[c]) for c in ["block", "spot_class", "first", "second", "nuc_label", "nuc_bins"])
        g4.append({"run": r.run, "unit": r.unit, "A": r.A, "src": r.src, "genes": r.genes, "n_genes": r.n_genes, "f33_run": f.run,
                   **{f"rel_dev_{m}": v for m, v in dev.items()}, "max_rel_dev": max(dev.values()), "rctd_blocks_identical": same_blocks,
                   "pass_1e-9": max(dev.values()) <= 1e-9})
    G4 = pd.DataFrame(g4); G4.to_csv(f"{OUT}/B_G4.csv", index=False)
    print("G-4:", int(G4["pass_1e-9"].sum()), "of", len(G4), "runs within 1e-9; max rel dev", G4.max_rel_dev.max())
# ---------- tempi ----------
tm = []
for arm in ("G4", "V", "T", "S", "B"):
    for f in sorted(glob.glob(f"{BR}/{arm}/runs/*/chain_times.tsv")):
        t = pd.read_csv(f, sep="\t").iloc[0].to_dict(); tm.append({"arm": arm, "run_dir": os.path.basename(os.path.dirname(f)), **t})
for f in sorted(glob.glob(f"{BR}/sim_rds/*_check.json")):
    j = json.load(open(f)); tm.append({"arm": "B_sim", "run_dir": f"{j['tag']}_k{j['k']}", "simulate": j["elapsed"]})
T = pd.DataFrame(tm); T.to_csv(f"{OUT}/B_timings.csv", index=False)
if len(T): print(T.groupby("arm").sum(numeric_only=True).round(0).to_string())
