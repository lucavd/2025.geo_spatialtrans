"""tools/S2.3_ph2_depth.py — S2.3, post hoc PH-2 (profondita dei nuclei simulati contro veri; v depurato dalla profondita) e tempi della funzione sul pannello."""

import pandas as pd, numpy as np
T = pd.read_csv("results/S2.3/S2.3_panel_roi.csv"); M = T[T.variant == "main"].copy()
M["v_sim"] = (M.panel_nucleus_psi_sim - 1) / M.nucleus_mbar_sim; M["v_real"] = (M.panel_nucleus_psi_real - 1) / M.nucleus_mbar_real; M["v_ratio"] = M.v_sim / M.v_real
M["mbar_ratio"] = M.nucleus_mbar_sim / M.nucleus_mbar_real
g = M.groupby("archetype").agg(nucleus_mbar_sim=("nucleus_mbar_sim","mean"), nucleus_mbar_real=("nucleus_mbar_real","mean"), mbar_ratio=("mbar_ratio","mean"),
    v_sim=("v_sim","mean"), v_real=("v_real","mean"), v_ratio=("v_ratio","mean"), v_ratio_sd=("v_ratio","std"), mean_umi_sim=("mean_umi_sim","mean"), mean_umi_real=("mean_umi_real","mean"))
g["v_ratio_lo2se"] = g.v_ratio - 2*g.v_ratio_sd/np.sqrt(5); g["v_ratio_hi2se"] = g.v_ratio + 2*g.v_ratio_sd/np.sqrt(5)
g.reset_index().to_csv("results/S2.3/S2.3_PH2_depth.csv", index=False); print(g.round(3).to_string())
G = pd.read_csv("/mnt/micron/geo_spatialtrans/S2.3/panel/main/G2.csv")
t = G.groupby("archetype").agg(n_tx_mean=("n_transcripts","mean"), sec_min=("elapsed_function_s","min"), sec_max=("elapsed_function_s","max"), lost=("n_lost","sum"), tx=("n_transcripts","sum"))
t["lost_frac"] = t.lost / t.tx; t.reset_index().to_csv("results/S2.3/S2.3_timing.csv", index=False); print(t.round(5).to_string())
