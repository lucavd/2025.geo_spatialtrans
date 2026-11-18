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
| linfonodo: geni degli anticorpi nei nuclei T quasi cinque volte | results/F3.1/F3.1_PH4c_ig_share.csv, linfociti T: sim_nucleus_bins, real_nucleus_bins | 0.0238, 0.0051 |
| un quinto di quello che si raccoglie sopra un nucleo T viene da B o plasmacellule | F3.1_PH4c_origin.csv, sim, linfociti T, B_pl_share | 0.196 |
| nel vero al massimo un ventesimo | F3.1_PH4c_ig_share.csv: 0.0051 / 0.0937 (immunoglobuline nei nuclei T veri / nei nuclei B veri) | 0.054 |
| profili dei linfociti B giusti | F3.1_PH4_lymph_calls.csv: B chiamati plasmacellule, sim_profile_thinned contro real_exp | 0.099, 0.117 |
| plasmacellule 30 % contro 13 %; tre quarti dell'eccesso da quadrati di altre cellule | results/F3.1/F3.1_PH4b_decomposition.csv: excess, rate_nonB | 0.302, 0.130; 0.127/0.172 = 0.74 |
| nuclei veri puliti quanto i profili, simulati no (intestino, tumore, stroma) | F3.1_PH4d_spill.csv, foreign_nuc (real, sim), foreign_profile | A1 0.210/0.324/0.224; A2 0.038/0.115/0.044; A3 0.181/0.291/0.222 |
| anche con i profili dal solo nucleo l'eccesso resta | F3.1_PH4d_spill_simnuc_areas.csv, foreign_nuc | A1 0.329, A2 0.103, A4 0.754 |
| i compartimenti ne tolgono circa metà nell'intestino e nel tumore | F3.1_PH4d_spill.csv, simq6 | A1 0.269 (eccesso da 0.114 a 0.059), A2 0.077 (da 0.077 a 0.039) |
| nell'intestino e nel linfonodo i vicini già mescolati come nel vero | F3.1_PH3.csv / PH-3, mescolanza Q6 contro vero | A1 0.23/0.22, A4 0.37/0.36 |
| la misura dell'orientamento vede il cuore e non gli altri tessuti | F3.1_PH5_aniso_diff.csv, aniso_umi | A6 0.0967 ± 0.0151; altri entro 2 ES |
| composizione da 30 a migliaia di volte, densità da 5 a 150 volte meno | F3.1_PH6_summary.csv, sim: median_cv_ratio_types, cv_ratio_density | 0.00024-0.033; 0.0065-0.21 |
| follicoli e paracorteccia simulati con il 49 % di B, veri dal 33 % al 64 % | F3.1_PH4_composition.csv, linfociti B | 0.495; 0.331-0.639 |
| errore di geometria capito e corretto; area di controllo identica | results/F3.1/F3.1_PH7.log; `tools/F3.1_PH7b_diag.R` | buco 1 tocca il bordo; A5_r3 ok (9 311 944 UMI); A5_r1 identica TRUE |
| nel cervello l'eccesso sta nei profili | F3.1_PH4d_spill.csv e simnuc, A5 | vero 0.203, simulato 0.261, profili dal nucleo 0.204 |
