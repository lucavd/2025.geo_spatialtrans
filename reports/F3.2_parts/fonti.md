| frase del riassunto | fonte (tabella o file, riga) | valore letto |
|---|---|---|
| 380 dei 478 valori dalle stesse 30 aree | `results/F3.2/F3.2_census.csv`, area_set in R2_ROI o R2b_WINDOW | 368 + 12 = 380 su 478 |
| nessuno script legge aree non censite | `results/F3.2/F3.2_readers.csv` | 146 script, 0 non classificati |
| 48 aree in otto tessuti | `results/F3.2/F3.2_area_agreement.csv`, kept | 48 (L1 15, L2 23, L3 10) |
| pazienti 1 e 5: 7 di tumore e 4 di stroma | stesso file, sheet CRC_P1_SR300 e CRC_P5_SR300 | tumore 4 + 3, stroma 3 + 1 |
| secondo tumore: 5 di tumore e 3 di stroma; secondo cervello 4 | stesso file | 5, 3, 4 |
| rene e fegato 10 aree | stesso file | 5 + 5 |
| intestino 4, cuore 5, linfonodo 2 follicoli di 500 µm, paracorticale 0, stroma originale 0 | stesso file, livello L1 | 4, 5, 2, 0, 0 |
| file integri | `data/holdout/F3.2_download_status.csv` | 46 OK su 46 |
| circa l'1 % fra versioni | `results/F3.2/F3.2_bridge_predictions.csv` | da +0,0087 a +0,0094 (UMI e geni), −0,012 e −0,015 (vuoti) |
| due regole fallite | `results/F3.2/F3.2_rule_validation.csv`, `F3.2_rule2_validation.csv` | 12/30, 18/30 |
| la terza funziona solo sulle sezioni che conosce | `results/F3.2/F3.2_area_agreement.csv` | kappa con il lettore 0,78 (taratura), 0,26 (altre sezioni) |
| i due lettori concordano su 9 aree su 10 | stesso file | 0,907 (97/107), kappa 0,882 |
| 13 aree vere su 240 | `results/F3.2/F3.2_R4_dryrun.csv`, arm sim | 13 su 240 |
| compartimenti: più che raddoppia, circa un terzo | stesso file, rapporto mediano delle DS | 0,124 (uniforme), 0,299 (compartimenti) |
| un'ora e un quarto per scaricare | job di scaricamento, primo file 13:31, fine con estrazione 14:40 (2026-11-18) | circa 1 h 10 min |
| 3,5 simulazioni per processore all'ora | F3.1, riassunto: 85 simulazioni con i due metodi in circa 4 ore con 6 processi | 85 / (4 × 6) ≈ 3,5 |
