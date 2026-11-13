| frase del riassunto | fonte (tabella o file, riga) | valore letto |
|---|---|---|
| 18 confronti su 36 | results/F3.1/F3.1_verdicts.csv, colonna verdict = smentita | 18 |
| tumore: 15 % contro meno del 5 % (marcatori estranei) | F3.1_verdicts.csv, M-3, A2: sim_est, real_est | 0.1515, 0.0452 |
| tumore: quasi un quadrato su cinque contro uno su venti (miscugli) | F3.1_verdicts.csv, M-4, A2 (1 − singoletti) | 0.1854, 0.0461 |
| stroma: 65 % contro 89 % (stesso tipo fra vicini) | F3.1_verdicts.csv, M-6, A3 | 0.6489, 0.8861 |
| trascritti per cellula tornano in tutti e sei i tessuti | F3.1_verdicts.csv, M-2 | confermata, confermata, confermata, confermata, confermata, confermata |
| intestino: quasi la metà contro un quinto dei vicini di altro tipo | F3.1_PH1_mixing.csv, mix_mean (media per area, poi per tessuto) | 0.472, 0.217 |
| tumore: un terzo contro uno su sedici | F3.1_PH1_mixing.csv, mix_mean | 0.345, 0.061 |
| quasi tutti (intestino) e due terzi (tumore) dei trascritti estranei in più | F3.1_PH1_mixing_reweight.csv, side = sim, fraction_of_gap_explained | 0.945, 0.651 |
| intestino: stroma in un quadrato su dodici contro uno su quattro | F3.1_PH1_mixing_rctd_composition.csv, A1, stromale e muscolo liscio | 0.078, 0.253 |
| stessa quota di cellule stromali, una su cinque | cells.parquet delle 5 aree simulate k = 0; ref_A1.json types_pool | 0.214; 8275/38681 = 0.214 |
| da un terzo a quasi tutta la differenza (miscugli, ordine, accordo) | F3.1_PH3.csv, fraction_closed, M-4...M-6 smentite fuori dal cervello | 0.37-0.95 |
| circa metà sulle cellule sporche (intestino, tumore) | F3.1_PH3.csv, M-3, A1 e A2 | 0.523, 0.595 |
| 5 delle 17 differenze spariscono | F3.1_PH3.csv con F3.1_verdicts.csv: smentite fuori dal cuore e Q6_vs_real diverso da smentita | 5 su 17 |
| con profili dal solo nucleo restano sporche, tranne il cervello | F3.1_CP3.csv, M-3, fraction_closed | A1 -0.06, A2 0.11, A4 0.08, A5 0.65 |
| linfonodo: plasmacellule 30 % contro 13 % | F3.1_PH1_mixing_rctd_composition.csv, A4, plasmacellule | 0.302, 0.130 |
| cuore: nessuna misura vede le cellule allungate | F3.1_verdicts.csv, M-6, A6 | confermata |
| una delle cinque aree del cervello si ferma | /mnt/micron/geo_spatialtrans/F3.1/main/logs_sim_q6.log | ERROR A5_r3 .s27_cdt() |
| simulazione di un'area 2-7 minuti; RCTD 4-20 minuti (mediana 11) | esplorazione E-7 (130-418 s); rctd_sec di F3.1_real_runs.csv, F3.1_sim_runs.csv, F3.1_simnuc_runs.csv | 232-1209 s, mediana 662 s |
