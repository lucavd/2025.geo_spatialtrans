"""tools/riesame_B_h.py - riesame, Parte B: margine h ricalcolato per insieme di geni (docs/riesame/REGOLA_B.md, «Tessuti, unita', geni»),
tumore (A2) e stroma (A3), M-1...M-6, con la costruzione di F3.1 sulle 5 aree di V: tools/S2_5_common.py corrected_interval() come in
tools/F3.1_intervals.py (M-1, M-2, M-4, M-6) e tools/F3.3_prereg_numbers.py (M-3, M-5 con le etichette trasferite): bootstrap della media
delle 5 aree (B = 2 000, generatore numpy default_rng([20261115, i, j]), i = posizione del tessuto in A1...A6, j = posizione della misura
in M-1...M-6), semi-ampiezze per kappa_boot del tessuto (results/R7/coverage/R7_kappa_final_v4.csv, unita' area; q = bin per M-1, M-4,
M-5, M-6, nuc per M-2, M-3) e per il fattore dell'asimmetria se CV > 0,45; limite inferiore >= 0; h = (min(1, hi) - lo) / 2 per le
quote, (hi - lo) / 2 per M-2; M-4 almeno 0,01 (DEC-086). Controllo: con gli stessi ingressi di F3.3 (righe IN_panel, geni completi)
il codice ritrova gli h di results/F3.3/F3.3_prereg_numbers.csv.
Uscita: results/riesame/B_h_panel.csv (genes = panel, common e la riga di controllo full_F33). Uso: .venv/bin/python -B tools/riesame_B_h.py"""
import os, sys, numpy as np, pandas as pd
os.chdir(os.path.expanduser("~/2025.geo_spatialtrans")); sys.path.insert(0, "tools"); import S2_5_common as CM
SEED = 20261115; ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]
MEAS = [("M-1", "M1", "bin", True), ("M-2", "M2", "nuc", False), ("M-3", "M3", "nuc", True), ("M-4", "M4", "bin", True), ("M-5", "M5", "bin", True), ("M-6", "M6", "bin", True)]
U = pd.read_csv("results/riesame/B_area_values.csv"); av = pd.read_csv("results/F3.3/F3.3_area_values.csv")
P = pd.read_csv("results/F3.3/F3.3_prereg_numbers.csv").set_index(["archetype", "id"])
sets = {"panel": U[(U.arm == "V") & (U.genes == "panel")], "common": U[(U.arm == "V") & (U.genes == "common")], "full_F33": av[av.src == "IN_panel"]}
rows = []
for gname, D in sets.items():
    for A in ("A2", "A3"):
        i = ARCH.index(A); R = D[D.A == A].sort_values("unit")
        assert list(R.unit) == [f"{A}_r{k}" for k in range(1, 6)], (gname, A, list(R.unit))
        for j, (mid, col, q, frac) in enumerate(MEAS):
            v = R[col].to_numpy(float); ci = CM.corrected_interval(v, A, np.random.default_rng([SEED, i, j]), side="real", q=q)
            hi = min(1.0, ci["hi"]) if frac else ci["hi"]; h_raw = (hi - ci["lo"]) / 2; h = max(h_raw, 0.01) if mid == "M-4" else h_raw
            h33 = float(P.loc[(A, mid), "h"])
            rows.append({"A": A, "measure": mid, "genes": gname, "n_areas": len(v), "V_mean": ci["est"], "lo": ci["lo"], "hi": hi, "h_raw": h_raw, "h": h,
                         "kappa": ci["kappa"], "cv": ci["cv"], "cv_factor": ci["cv_factor"], "h_F33_full_genes": h33, "h_over_h_F33": h / h33,
                         "check_rel_dev_vs_F33": abs(h - h33) / h33 if gname == "full_F33" else np.nan})
H = pd.DataFrame(rows); H.to_csv("results/riesame/B_h_panel.csv", index=False)
c = H[H.genes == "full_F33"]; print("control vs F3.3 prereg h: max rel dev", c.check_rel_dev_vs_F33.max(), "over", len(c), "values")
pd.set_option("display.width", 200); print(H[H.genes != "full_F33"][["A", "measure", "genes", "V_mean", "h", "h_F33_full_genes", "h_over_h_F33"]].round(5).to_string(index=False))
