| frase del riassunto | fonte (tabella o file, riga) | valore letto |
|---|---|---|
| 187 affermazioni, 124 coincidevano, 37 superate, 19 senza test, 5 scostamenti | `results/S2.8/S2.8_design_map.csv`, colonna `status` | 187; consistent 124, design_behind 37, untested 19, code_deviates 5, unclear 2 |
| pezzo solo da un 3 % in meno a un 23 % in più del tempo della funzione vecchia | `results/S2.8/S2.8_rule_R1.csv`, colonna `ratio` | 0.9709 (cuore) ... 1.2334 (cervello) |
| al più 6,2 GB in un pezzo solo | `results/S2.8/S2.8_predictions.csv`, P-3, colonna `est` | 6.2014 (tumore del colon) |
| pezzi da mezzo millimetro: da 85 secondi a 4 minuti e mezzo | `results/S2.8/S2.8_predictions.csv`, P-1 | 85.31 (cuore), 269.67 (tumore del colon) |
| meno di 3,5 GB con i pezzi da mezzo millimetro | `results/S2.8/S2.8_predictions.csv`, P-2 | massimo 3.2499 (tumore del colon) |
| la volta scorsa sembrava fuori bilancio | `docs/steps/S2.8.md`, «Stime di preparazione», E-0 | 385 s con 16 processi, carico 15-26; verità geometrica già calcolata per il motore di S2.5 |
| file identici bit per bit | `results/S2.8/S2.8_G2_identity.csv`; `results/S2.8/S2.8_G2bis_final_code.csv` | 1 080/1 080; G-2 bis 84/84 |
| i file cambiano in tutti i sei tessuti | `results/S2.8/S2.8_G2_identity.csv`, gruppo `d_cp2_vs_S2.7` | 6/6 diversi |
| i tipi di cellule cambiati, la ripresa restituiva i conteggi vecchi | `docs/steps/S2.8.md`, «Revisione avversariale», A15 | 1 474 cellule |
| meno del 2 % nel tessuto e nel vetro vicino; striscia da 0,4 a 1 UMI | `results/S2.8/S2.8_BL237.csv` | tessuto 0.9998-1.0000, fasce 0.997-1.002, vetro fra finestra e tessuto 1.018-1.020, striscia 0.3785 → 0.9819 |
| pezzo solo più veloce dei pezzi da mezzo millimetro da 1,1 a 1,2 volte | `results/S2.8/S2.8_summary.csv`, rapporto `time_mean` e504 / einf | 1.08 (linfonodo) ... 1.24 (tumore del colon) |
| 81 simulazioni, 3,9 ore-processo, 44 minuti, 86 °C | json delle esecuzioni; `/mnt/micron/geo_spatialtrans/S2.8/sched.log` | 81; 3.94 h; 19:59:23-20:43:25; Tctl massimo 85.6 |
