## Da dove vengono i numeri del riassunto

| frase del riassunto | fonte (file, riga) | valore letto |
|---|---|---|
| nuclei della maschera scelta nelle finestre di linfonodo / punti di Luca | R3b_p1p2_summary.csv A4 consensus n_obj, n_manual | 416.0 / 404.0 |
| StarDist da solo nel linfonodo | R3b_p1p2_summary.csv A4 stardist_he n_obj; p1_ratio_mean | 314.0; 0.768 |
| cuore: maschera scelta, StarDist, Luca | R3b_p1p2_summary.csv A6 n_obj, n_manual | 262.0, 202.0, 242.0 |
| poligoni Space Ranger rimasti nelle finestre di linfonodo, di cui con un punto di Luca | R3b_posthoc_srprec.csv A4 sr | 102, 65 |
| idem nel cuore | R3b_posthoc_srprec.csv A6 sr | 67, 47 |
| contatti simulati nel linfonodo: punti di Luca, maschera scelta, StarDist | R3b_predictions.csv CP-1, P-3 A4 | 0.715, 0.690, 0.408 |
| contatti simulati nel cuore (maschera scelta) | R3b_predictions.csv P-3 A6 CO | 0.0466 |
| valore manuale di Luca e margine (linfonodo, cuore) | results/R2c/R2c_posthoc_ci.csv, due stadi | 0.550 [0.35, 0.750]; 0.017 [0.00, 0.070] |
| frazione di nuclei con forma Space Ranger (celle interne) | R3b_posthoc_area.csv frac_sr_ev | 0.268, 0.267 |
| poligono Space Ranger / StarDist sullo stesso nucleo, linfonodo (mediana) | masks/A4_*_consensus_map.parquet, src = sd | 1.344 |
| variabilita' della taglia nel linfonodo: tutti, solo forma StarDist, oggi | R3b_posthoc_area.csv sigma_all, sigma_sd; R3b_refs.csv sigma old_recomputed | 0.288, 0.210, 0.212 |
| contatti simulati nel linfonodo con la taglia dai soli nuclei di forma StarDist | R3b_posthoc_sdarea_windows.csv A4 | 0.736 |
| nuclei focali del linfonodo: contatti nella maschera, di Luca, separati per Luca | R3b_posthoc_agreement.csv A4 consensus | 44.0, 32.0, 15.0 su 59.0 |
| contatti che passano da un poligono Space Ranger | R3b_posthoc_partners.csv A4, contact = True, any_sr_partner | 0.636 |
| cambiamento massimo delle misure dei territori | R3b_refs.csv rel_change, metriche dei territori | 0.0157 |
| rapporto nucleo/cellula mediano nel linfonodo: oggi, nuovo | R3b_refs.csv A4 nc | 0.606, 0.480 |
| tetto del nucleo nel cuore: oggi, nuovo (tutti / solo forma StarDist) | R3b_refs.csv A6 phi; R3b_posthoc_area.csv phi_sd | 0.180, 0.147 / 0.149 |
| nuclei con centro fuori dal tessuto valido: maschera scelta, Space Ranger, StarDist; dopo il filtro; generatori | R3b_posthoc_valid.csv | 59.0, 60.0, 107.0; 0; 0 |
| nuclei nella maschera scelta (10 ROI) | R3b_posthoc_valid.csv consensus n_labels | 145098.0 |
