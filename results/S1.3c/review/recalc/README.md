# S1.3c — revisione avversariale «Ricalcolo» (recalc)

Revisore: sotto-agente «Ricalcolo», 2026-10-12. Codice indipendente: nessun import da `tools/S1.3c_ecm.py`, nessuna funzione copiata; il codice della pipeline è stato **letto solo dopo** il primo confronto, per trovare la causa delle discrepanze. Commit di riferimento `5fa0294`.

## File
- `recalc_umi_bg.py` — UMI per bin (geni della filtrata letti dalla grezza; tutti i geni), C-3, soglie di sfondo q95 (k = 3, 5) per i 5 dataset → `recalc_bg_<dataset>.json`.
- `recalc_roi.py` — per ROI: n bin, bin in tessuto valido, territori (generatore R3 più vicino, cKDTree), d_nuc (EDT nativa sui nuclei `keep`, scala 1), finestre, φ_obs, nullo (permutazione **mia**, seed 9001–9030, 30 permutazioni; 200 per A2 e A6), E, r*, BA, R, CP-2; `ctrl` = controllo positivo tipo A su A3 r1. Variabile `RECALC_MAP`: `rint_noclip` (prima versione del revisore), `floor_clip` (mappatura esatta della maschera R2), `px4` (mappatura della pipeline, per isolare la causa).
- `compare_bg.py`, `compare_roi.py`, `summarise.py` — confronti; `recalc_vs_reported.csv` = tabella completa (198 righe).
- Intermedi grandi (bin per dataset) in `/tmp/recalc_s13c/` sul server (fuori repo, non versionati).

## 1. UMI e sfondo (tutti e 5 i dataset; richiesto: Lymph Node)

| dataset | C-3 bin in tessuto | n filtrata | bin con UMI(geni comuni, grezza) ≠ filtrata | idem tutti i geni | soglia k=3 (rip./ric.) | n finestre k=3 (rip./ric.) | soglia k=5 (rip./ric.) | n finestre k=5 (rip./ric.) |
|---|---|---|---|---|---|---|---|---|
| Visium_HD_6p5mm_Human_Heart | 10615406 | 10615406 | 0 (rip. 0) | 2043438 (rip. 2043438) | 5 / 5 | 33390 / 33390 | 11 / 11 | 11780 / 11780 |
| Visium_HD_6p5mm_Mouse_Brain | 6343031 | 6343031 | 0 (rip. 0) | 196172 (rip. 196172) | 10 / 10 | 466521 / 466521 | 26 / 26 | 167325 / 167325 |
| Visium_HD_Human_Colon_Cancer | 8176117 | 8176117 | 0 (rip. 0) | 2707020 (rip. 2707020) | 15 / 15 | 228961 / 228961 | 42 / 42 | 81330 / 81330 |
| Visium_HD_Human_Lymph_Node_FFPE | 5482441 | 5482441 | 0 (rip. 0) | 408027 (rip. 408027) | 2 / 2 | 471420 / 471420 | 4 / 4 | 168102 / 168102 |
| Visium_HD_Mouse_Small_Intestine | 5479660 | 5479660 | 0 (rip. 0) | 1592446 (rip. 1592446) | 5 / 5 | 483562 / 483562 | 12 / 12 | 173182 / 173182 |

Anche media e frazione di zeri delle finestre di sfondo (com/all, k = 3/5) coincidono entro 1e-12 in tutti i dataset. Small Intestine (Space Ranger 3.0.0): 4111129 barcode di `tissue_positions` assenti dalla matrice grezza, di cui 0 con `in_tissue = 1`; vanno contati come 0 UMI (come fa la pipeline). La mia prima versione li lasciava NaN → soglia NaN: errore **mio**, corretto.

## 2. ROI (variante primaria, geni comuni)

Colonne: **rip.** = riportato; **ric. R2** = mio ricalcolo con la mappatura bin→maschera ds4 esatta di R2 (`floor(px/4)`, clip all’ultimo pixel); **ric. pipeline** = mio ricalcolo con la mappatura della pipeline (`px4 = side_um/913`); per φ_null ed E il ricalcolo usa le mie 30 permutazioni, quindi si riporta anche z = (ric. − rip.)/SE, SE = sd·√(1/20+1/30).


### k = 3 (s = 6 µm)

| ROI | grandezza | rip. | ric. R2 | ric. pipeline | ric.R2 − rip. | z (nullo) |
|---|---|---|---|---|---|---|
| A1_r1 | n_bins_roi (C-1) | 250152 | 250152 | 250152 | 0 |  |
| A1_r1 | n_tissue_bins (PH1) | 199707 | 199656 | 199707 | -51 |  |
| A1_r1 | n_gen | 8739 | 8739 | 8739 | 0 |  |
| A1_r1 | R [4,8)/[0,2) | 0.723031 | 0.723538 | 0.723031 | 0.000506319 |  |
| A1_r1 | n_lumen_bins CP-2 (PH1) | 20252 | 20328 | 20252 | 76 |  |
| A1_r1 | umi_lumen_mean (PH1) | 6.02024 | 5.99946 | 6.02024 | -0.020786 |  |
| A1_r1 | thr | 5 | 5 | 5 | 0 |  |
| A1_r1 | n_tissue_windows | 21066 | 21046 | 21066 | -20 |  |
| A1_r1 | n_empty | 2 | 2 | 2 | 0 |  |
| A1_r1 | phi_obs | 9.49397e-05 | 9.50299e-05 | 9.49397e-05 | 9.02211e-08 |  |
| A1_r1 | phi_null (media) | 0 | 0 | 3.16466e-06 | 0 |  |
| A1_r1 | phi_null_min | 0 | 0 | 0 | 0 |  |
| A1_r1 | phi_null_max | 0 | 0 | 4.74699e-05 | 0 |  |
| A1_r1 | E | 9.49397e-05 | 9.50299e-05 | 9.17751e-05 | 9.02211e-08 |  |
| A1_r1 | r_star | 39 | 39 | 39 | 0 |  |
| A1_r1 | BA | 0.499953 | 0.499952 | 0.499953 | -4.51191e-08 |  |
| A1_r1 | CP-2 n_windows | 1395 | 1388 | 1395 | -7 |  |
| A1_r1 | CP-2 n_empty | 30 | 30 | 30 | 0 |  |
| A2_r1 | n_bins_roi (C-1) | 250243 | 250243 | 250243 | 0 |  |
| A2_r1 | n_tissue_bins (PH1) | 186178 | 186252 | 186178 | 74 |  |
| A2_r1 | n_gen | 6546 | 6546 | 6546 | 0 |  |
| A2_r1 | R [4,8)/[0,2) | 0.923671 | 0.923554 | 0.923671 | -0.000116793 |  |
| A2_r1 | n_lumen_bins CP-2 (PH1) | 30443 | 30402 | 30443 | -41 |  |
| A2_r1 | umi_lumen_mean (PH1) | 5.32099 | 5.21926 | 5.32099 | -0.101731 |  |
| A2_r1 | thr | 15 | 15 | 15 | 0 |  |
| A2_r1 | n_tissue_windows | 19836 | 19839 | 19836 | 3 |  |
| A2_r1 | n_empty | 83 | 81 | 83 | -2 |  |
| A2_r1 | phi_obs | 0.00418431 | 0.00408287 | 0.00418431 | -0.000101444 |  |
| A2_r1 | phi_null (media) | 0.000259629 | 0.000257069 | 0.000238623 | -2.55955e-06 | -0.63 |
| A2_r1 | phi_null_min | 0 | 0.000100812 | 0.000100827 | 0.000100812 |  |
| A2_r1 | phi_null_max | 0.000403307 | 0.000604869 | 0.000403307 | 0.000201562 |  |
| A2_r1 | E | 0.00392468 | 0.0038258 | 0.00394569 | -9.88847e-05 | 0.63 |
| A2_r1 | r_star | 29.25 | 29.5 | 29.25 | 0.25 |  |
| A2_r1 | BA | 0.757998 | 0.752074 | 0.757998 | -0.00592416 |  |
| A2_r1 | CP-2 n_windows | 2845 | 2845 | 2845 | 0 |  |
| A2_r1 | CP-2 n_empty | 774 | 774 | 774 | 0 |  |
| A3_r1 | n_bins_roi (C-1) | 250229 | 250229 | 250229 | 0 |  |
| A3_r1 | n_tissue_bins (PH1) | 249353 | 249344 | 249353 | -9 |  |
| A3_r1 | n_gen | 4165 | 4165 | 4165 | 0 |  |
| A3_r1 | R [4,8)/[0,2) | 0.688517 | 0.68858 | 0.688517 | 6.27413e-05 |  |
| A3_r1 | n_lumen_bins CP-2 (PH1) | 876 | 885 | 876 | 9 |  |
| A3_r1 | umi_lumen_mean (PH1) | 7.77169 | 7.63503 | 7.77169 | -0.136661 |  |
| A3_r1 | thr | 15 | 15 | 15 | 0 |  |
| A3_r1 | n_tissue_windows | 27449 | 27447 | 27449 | -2 |  |
| A3_r1 | n_empty | 808 | 809 | 808 | 1 |  |
| A3_r1 | phi_obs | 0.0294364 | 0.029475 | 0.0294364 | 3.85788e-05 |  |
| A3_r1 | phi_null (media) | 0.0205017 | 0.0205523 | 0.0205011 | 5.06796e-05 | -0.01 |
| A3_r1 | phi_null_min | 0.0200372 | 0.0201115 | 0.0201829 | 7.43278e-05 |  |
| A3_r1 | phi_null_max | 0.0208751 | 0.0208766 | 0.0209115 | 1.52112e-06 |  |
| A3_r1 | E | 0.00893475 | 0.00892265 | 0.00893536 | -1.21008e-05 | 0.01 |
| A3_r1 | r_star | 13.5 | 13.5 | 13.5 | 0 |  |
| A3_r1 | BA | 0.544805 | 0.544751 | 0.544805 | -5.46363e-05 |  |
| A3_r1 | CP-2 n_windows | 62 | 62 | 62 | 0 |  |
| A3_r1 | CP-2 n_empty | 2 | 2 | 2 | 0 |  |
| A4_f1 | n_bins_roi (C-1) | 250082 | 250082 | 250082 | 0 |  |
| A4_f1 | n_tissue_bins (PH1) | 247689 | 247706 | 247689 | 17 |  |
| A4_f1 | n_gen | 23707 | 23707 | 23707 | 0 |  |
| A4_f1 | R [4,8)/[0,2) | 1.0169 | 1.01688 | 1.0169 | -1.5486e-05 |  |
| A4_f1 | n_lumen_bins CP-2 (PH1) | 2393 | 2376 | 2393 | -17 |  |
| A4_f1 | umi_lumen_mean (PH1) | 0.36356 | 0.366582 | 0.36356 | 0.00302211 |  |
| A4_f1 | thr | 2 | 2 | 2 | 0 |  |
| A4_f1 | n_tissue_windows | 27216 | 27217 | 27216 | 1 |  |
| A4_f1 | n_empty | 219 | 219 | 219 | 0 |  |
| A4_f1 | phi_obs | 0.00804674 | 0.00804644 | 0.00804674 | -2.95651e-07 |  |
| A4_f1 | phi_null (media) | 0.00520282 | 0.00524917 | 0.00525304 | 4.63484e-05 | 0.67 |
| A4_f1 | phi_null_min | 0.00492357 | 0.00473968 | 0.00459289 | -0.00018389 |  |
| A4_f1 | phi_null_max | 0.00554821 | 0.00595216 | 0.00562169 | 0.000403955 |  |
| A4_f1 | E | 0.00284392 | 0.00279727 | 0.0027937 | -4.6644e-05 | -0.67 |
| A4_f1 | r_star | 30.5 | 30 | 30.5 | -0.5 |  |
| A4_f1 | BA | 0.880353 | 0.8849 | 0.880353 | 0.00454772 |  |
| A4_f1 | CP-2 n_windows | 181 | 175 | 181 | -6 |  |
| A4_f1 | CP-2 n_empty | 117 | 109 | 117 | -8 |  |
| A5_r1 | n_bins_roi (C-1) | 250018 | 250018 | 250018 | 0 |  |
| A5_r1 | n_tissue_bins (PH1) | 186572 | 186618 | 186572 | 46 |  |
| A5_r1 | n_gen | 1267 | 1267 | 1267 | 0 |  |
| A5_r1 | R [4,8)/[0,2) | 0.56119 | 0.560893 | 0.56119 | -0.000296587 |  |
| A5_r1 | n_lumen_bins CP-2 (PH1) | 63446 | 63400 | 63446 | -46 |  |
| A5_r1 | umi_lumen_mean (PH1) | 19.338 | 19.338 | 19.338 | 4.01693e-05 |  |
| A5_r1 | thr | 10 | 10 | 10 | 0 |  |
| A5_r1 | n_tissue_windows | 18466 | 18470 | 18466 | 4 |  |
| A5_r1 | n_empty | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_obs | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_null (media) | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_null_min | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_null_max | 0 | 0 | 0 | 0 |  |
| A5_r1 | E | 0 | 0 | 0 | 0 |  |
| A5_r1 | r_star | 46.5 | 47.75 | 46.5 | 1.25 |  |
| A5_r1 | BA | — | — | — | — |  |
| A5_r1 | CP-2 n_windows | 5164 | 5140 | 5164 | -24 |  |
| A5_r1 | CP-2 n_empty | 0 | 0 | 0 | 0 |  |
| A6_r1 | n_bins_roi (C-1) | 249873 | 249873 | 249873 | 0 |  |
| A6_r1 | n_tissue_bins (PH1) | 230827 | 230827 | 230827 | 0 |  |
| A6_r1 | n_gen | 674 | 674 | 674 | 0 |  |
| A6_r1 | R [4,8)/[0,2) | 0.675913 | 0.675913 | 0.675913 | 0 |  |
| A6_r1 | n_lumen_bins CP-2 (PH1) | 15598 | 15598 | 15598 | 0 |  |
| A6_r1 | umi_lumen_mean (PH1) | 3.94326 | 3.94326 | 3.94326 | 0 |  |
| A6_r1 | thr | 5 | 5 | 5 | 0 |  |
| A6_r1 | n_tissue_windows | 24403 | 24403 | 24403 | 0 |  |
| A6_r1 | n_empty | 227 | 227 | 227 | 0 |  |
| A6_r1 | phi_obs | 0.00930213 | 0.00930213 | 0.00930213 | 0 |  |
| A6_r1 | phi_null (media) | 0.00189936 | 0.00177301 | 0.00177301 | -0.000126351 | -2.20 |
| A6_r1 | phi_null_min | 0.00147523 | 0.00139327 | 0.00139327 | -8.19571e-05 |  |
| A6_r1 | phi_null_max | 0.00241774 | 0.00217186 | 0.00217186 | -0.000245871 |  |
| A6_r1 | E | 0.00740278 | 0.00752913 | 0.00752913 | 0.000126351 | 2.20 |
| A6_r1 | r_star | 53.25 | 53.25 | 53.25 | 0 |  |
| A6_r1 | BA | 0.495388 | 0.495388 | 0.495388 | 0 |  |
| A6_r1 | CP-2 n_windows | 958 | 958 | 958 | 0 |  |
| A6_r1 | CP-2 n_empty | 49 | 49 | 49 | 0 |  |

### k = 5 (s = 10 µm)

| ROI | grandezza | rip. | ric. R2 | ric. pipeline | ric.R2 − rip. | z (nullo) |
|---|---|---|---|---|---|---|
| A1_r1 | thr | 12 | 12 | 12 | 0 |  |
| A1_r1 | n_tissue_windows | 7186 | 7183 | 7186 | -3 |  |
| A1_r1 | n_empty | 0 | 0 | 0 | 0 |  |
| A1_r1 | phi_obs | 0 | 0 | 0 | 0 |  |
| A1_r1 | phi_null (media) | 0 | 0 | 0 | 0 |  |
| A1_r1 | phi_null_min | 0 | 0 | 0 | 0 |  |
| A1_r1 | phi_null_max | 0 | 0 | 0 | 0 |  |
| A1_r1 | E | 0 | 0 | 0 | 0 |  |
| A1_r1 | r_star | 31 | 31 | 31 | 0 |  |
| A1_r1 | BA | — | — | — | — |  |
| A1_r1 | CP-2 n_windows | 317 | 318 | 317 | 1 |  |
| A1_r1 | CP-2 n_empty | 2 | 2 | 2 | 0 |  |
| A2_r1 | thr | 42 | 42 | 42 | 0 |  |
| A2_r1 | n_tissue_windows | 6853 | 6859 | 6853 | 6 |  |
| A2_r1 | n_empty | 9 | 8 | 9 | -1 |  |
| A2_r1 | phi_obs | 0.00131329 | 0.00116635 | 0.00131329 | -0.000146943 |  |
| A2_r1 | phi_null (media) | 7.29607e-06 | 1.45794e-05 | 4.86405e-06 | 7.28331e-06 | -0.19 |
| A2_r1 | phi_null_min | 0 | 0 | 0 | 0 |  |
| A2_r1 | phi_null_max | 0.000145921 | 0.000145794 | 0.000145921 | -1.27647e-07 |  |
| A2_r1 | E | 0.001306 | 0.00115177 | 0.00130843 | -0.000154226 | 0.19 |
| A2_r1 | r_star | 42.75 | 40 | 42.75 | -2.75 |  |
| A2_r1 | BA | 0.721857 | 0.687135 | 0.721857 | -0.0347218 |  |
| A2_r1 | CP-2 n_windows | 881 | 884 | 881 | 3 |  |
| A2_r1 | CP-2 n_empty | 232 | 232 | 232 | 0 |  |
| A3_r1 | thr | 42 | 42 | 42 | 0 |  |
| A3_r1 | n_tissue_windows | 9799 | 9798 | 9799 | -1 |  |
| A3_r1 | n_empty | 248 | 248 | 248 | 0 |  |
| A3_r1 | phi_obs | 0.0253087 | 0.0253113 | 0.0253087 | 2.58305e-06 |  |
| A3_r1 | phi_null (media) | 0.0183692 | 0.0184459 | 0.0184475 | 7.672e-05 | 0.80 |
| A3_r1 | phi_null_min | 0.017859 | 0.0178608 | 0.0177569 | 1.82272e-06 |  |
| A3_r1 | phi_null_max | 0.0190836 | 0.0191876 | 0.0189815 | 0.000104009 |  |
| A3_r1 | E | 0.00693948 | 0.00686535 | 0.00686124 | -7.4137e-05 | -0.80 |
| A3_r1 | r_star | 11.75 | 11.75 | 11.75 | 0 |  |
| A3_r1 | BA | 0.534854 | 0.534853 | 0.534854 | -1.20598e-06 |  |
| A3_r1 | CP-2 n_windows | 10 | 11 | 10 | 1 |  |
| A3_r1 | CP-2 n_empty | 0 | 0 | 0 | 0 |  |
| A4_f1 | thr | 4 | 4 | 4 | 0 |  |
| A4_f1 | n_tissue_windows | 9653 | 9654 | 9653 | 1 |  |
| A4_f1 | n_empty | 31 | 31 | 31 | 0 |  |
| A4_f1 | phi_obs | 0.00321144 | 0.0032111 | 0.00321144 | -3.32653e-07 |  |
| A4_f1 | phi_null (media) | 0.00129493 | 0.00138112 | 0.00136745 | 8.61859e-05 | 0.82 |
| A4_f1 | phi_null_min | 0.000725163 | 0.000621504 | 0.000828758 | -0.000103659 |  |
| A4_f1 | phi_null_max | 0.0019683 | 0.00207168 | 0.00238268 | 0.00010338 |  |
| A4_f1 | E | 0.0019165 | 0.00182998 | 0.00184399 | -8.65185e-05 | -0.82 |
| A4_f1 | r_star | 59 | 59.5 | 59 | 0.5 |  |
| A4_f1 | BA | 0.886733 | 0.870552 | 0.886733 | -0.016181 |  |
| A4_f1 | CP-2 n_windows | 39 | 38 | 39 | -1 |  |
| A4_f1 | CP-2 n_empty | 21 | 20 | 21 | -1 |  |
| A5_r1 | thr | 26 | 26 | 26 | 0 |  |
| A5_r1 | n_tissue_windows | 5873 | 5868 | 5873 | -5 |  |
| A5_r1 | n_empty | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_obs | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_null (media) | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_null_min | 0 | 0 | 0 | 0 |  |
| A5_r1 | phi_null_max | 0 | 0 | 0 | 0 |  |
| A5_r1 | E | 0 | 0 | 0 | 0 |  |
| A5_r1 | r_star | 34.5 | 34.5 | 34.5 | 0 |  |
| A5_r1 | BA | — | — | — | — |  |
| A5_r1 | CP-2 n_windows | 1286 | 1284 | 1286 | -2 |  |
| A5_r1 | CP-2 n_empty | 0 | 0 | 0 | 0 |  |
| A6_r1 | thr | 11 | 11 | 11 | 0 |  |
| A6_r1 | n_tissue_windows | 8367 | 8367 | 8367 | 0 |  |
| A6_r1 | n_empty | 13 | 13 | 13 | 0 |  |
| A6_r1 | phi_obs | 0.00155372 | 0.00155372 | 0.00155372 | 0 |  |
| A6_r1 | phi_null (media) | 4.1831e-05 | 1.59356e-05 | 1.59356e-05 | -2.58954e-05 | -2.17 |
| A6_r1 | phi_null_min | 0 | 0 | 0 | 0 |  |
| A6_r1 | phi_null_max | 0.000239034 | 0.000119517 | 0.000119517 | -0.000119517 |  |
| A6_r1 | E | 0.00151189 | 0.00153779 | 0.00153779 | 2.58954e-05 | 2.17 |
| A6_r1 | r_star | 60.75 | 60.75 | 60.75 | 0 |  |
| A6_r1 | BA | 0.499222 | 0.499222 | 0.499222 | 0 |  |
| A6_r1 | CP-2 n_windows | 150 | 150 | 150 | 0 |  |
| A6_r1 | CP-2 n_empty | 6 | 6 | 6 | 0 |  |

Intervallo del nullo: il mio min–max su 30 permutazioni si sovrappone a quello riportato (20 permutazioni) in tutti i ROI; con 200 permutazioni (A2, A6) il riportato cade sempre dentro il mio min–max (vedi rilievo RC-5).

## 3. Controllo positivo tipo A su A3 r1 (k = 3)

| modo | grandezza | rip. (C5 / PH-3 f=0.10) | ric. R2 | ric. pipeline | z |
|---|---|---|---|---|---|
| thr | frac_zeroed (C5 A_frac_zeroed) | 0.0999146 | 0.0999182 | 0.0999146 |  |
| thr | phi_obs (PH3 A f=0.1) | 0.0957776 | 0.0957846 | 0.0957776 |  |
| thr | phi_null (PH3 A f=0.1) | 0.0269281 | 0.0269671 | 0.0270975 | 1.51 |
| thr | E (C5 A_E) | 0.0688495 | 0.0688175 | 0.0686801 | -1.51 |
| thr | E (PH3 A f=0.1) | 0.0688495 | 0.0688175 | 0.0686801 | -1.51 |
| thr | BA (C5 A_BA) | 0.850557 | 0.850556 | 0.850557 |  |
| rank | frac_zeroed (C5 A_frac_zeroed) | 0.0999146 | 0.0999984 | 0.0999988 |  |
| rank | phi_obs (PH3 A f=0.1) | 0.0957776 | 0.0958575 | 0.0958141 |  |
| rank | phi_null (PH3 A f=0.1) | 0.0269281 | 0.0269708 | 0.0271097 | 1.63 |
| rank | E (C5 A_E) | 0.0688495 | 0.0688867 | 0.0687044 | -1.30 |
| rank | E (PH3 A f=0.1) | 0.0688495 | 0.0688867 | 0.0687044 | -1.30 |
| rank | BA (C5 A_BA) | 0.850557 | 0.850699 | 0.850629 |  |

`thr` = azzerati i bin di tessuto con d_nuc > r₉₀ (definizione della pipeline, nota 12); `rank` = esattamente il 10 % dei bin con d_nuc più alto (pareggi rotti a caso, seed 777), lettura letterale del compito. C5 `A_E` e PH-3 (A, f = 0.10) sono lo stesso numero (0.0688495) perché escono dallo stesso codice: non sono due conferme indipendenti.


## 4. Rilievi

### RC-1 — confermato
**Affermazione.** UMI per bin (geni della filtrata letti dalla grezza), C-3 e soglie di sfondo q95 coincidono esattamente con S1.3c_bg_*_bg.json in 5/5 dataset (richiesto solo Lymph Node).

**Evidenza.** recalc_bg_<dataset>.json vs S1.3c_bg_<dataset>_bg.json (compare_bg.py): LN soglia k=3 2/2, k=5 4/4, finestre di sfondo 471420/168102 identiche; bin in tessuto con UMI(geni comuni) ≠ filtrata = 0 in 5/5; tutti i geni: 408027 (LN) identico; media e frazione di zeri delle finestre entro 1e-12. Small Intestine: 4111129 barcode assenti dalla grezza (SR 3.0.0), tutti con in_tissue = 0, contati come 0 come nella pipeline.

**Azione suggerita.** Nessuna.

### RC-2 — bug
**Affermazione.** La pipeline mappa i bin sulle maschere ds4 di R2 con px4 = side_um / shape[1] invece che con 4 px nativi: la griglia ds4 viene allungata di un fattore 913/913.5 (side 3654), 913/913.25 (side 3653) o 912/912.25 (side 3649). Disallineamento da 0 all'origine fino a circa 2 px nativi (≈0.55 µm) sul bordo opposto. Effetto trascurabile su tutti i numeri chiave; nessun verdetto cambia.

**Evidenza.** tools/S1.3c_ecm.py r. 225–226 (load_roi_context) e r. 302–303 (CP-2); generators() r. 169 (solo il fallback sens_seg). R2 costruisce la maschera nativa come np.repeat(t4, 4) con estensione del bordo (tools/R2_common.py, masks_for/up), quindi la mappatura corretta è floor(px/4) con clip. Con la mia mappatura px4 tutti i 141 valori deterministici tornano identici; con floor_clip gli scarti sono: bin di tessuto da −51 a +74 su ~2·10^5 (A6, side 3652 = 4·913: 0), φ_obs |Δ| ≤ 1.5e-4, E |Δ| ≤ 1.6e-4 (max 1.54e-4, A2 r1 k=5, include seed), R |Δ| ≤ 5.1e-4, CP-2 finestre |Δ| ≤ 24 (A5 r1 k3: 5140 contro 5164), n vuote |Δ| ≤ 8 (A4 f1 k3: 109 contro 117).

**Azione suggerita.** Sostituire px4 con 4·um_per_px (cioè floor(px_nativo/4), clip all'ultimo pixel) in load_roi_context, CP-2 e generators; rieseguire per completezza e dichiararlo nel report come correzione senza effetto sui verdetti (Δ attesi della tabella del README).

### RC-3 — confermato
**Affermazione.** Con la stessa mappatura px4 il mio codice indipendente riproduce esattamente tutte le grandezze deterministiche dei 6 ROI (k = 3 e 5): n bin nel ROI (C-1), bin in tessuto valido, n generatori, R, n e UMI medio dei bin 'lume' (PH-1/CP-2), n finestre di tessuto, n vuote, φ_obs, r*, BA, finestre CP-2 e vuote.

**Evidenza.** recalc_vs_reported.csv, colonna recalc_px4: 141/141 righe deterministiche con |Δ| < 1e-9 (3 BA non definite in entrambi). Esempi: A3 r1 n_tissue_bins 249353, φ_obs(k=3) 0.0294364, R 0.688517; A4 f1 φ_obs 0.00804674, BA 0.880353; CP-2 A2 r1 k3 2845/774.

**Azione suggerita.** Nessuna (oltre a RC-2).

### RC-4 — confermato
**Affermazione.** φ_null ed E coincidono con il riportato entro la variabilità di seed; la mia permutazione dentro il territorio conserva i totali per territorio (assert) ed è equidistribuita come quella della pipeline (che è l'inversa della mia).

**Evidenza.** 30 permutazioni (seed 9001–9030) a parità di bin (px4): z = (ric − rip)/(sd·√(1/20+1/30)) tra −0.82 e 0.82 per A2, A3, A4; A6 r1 −2.20 (k=3) e −2.17 (k=5). Con 200 permutazioni A6 r1: z della media riportata = 1.93 (k=3), 0.85 (k=5); A2 r1: 0.67, 0.23. Su circa 14 confronti, |z| ≈ 2 non è significativo. Gli intervalli min–max riportati cadono dentro il mio min–max a 200 permutazioni (A6 k3: rip. 0.001475–0.002418, mio 0.001188–0.002541). A1 r1 e A5 r1: φ_null = 0 in entrambi (A1 k3 px4: 2 eventi su 30×21066 finestre, media 3.2e-6).

**Azione suggerita.** Nessuna. Opzionale: 100+ permutazioni nei ROI con φ_null ~1e-3 (A6), dove 20 permutazioni danno SE ≈ 5e-5.

### RC-5 — confermato
**Affermazione.** Controllo positivo tipo A su A3 r1 (10 % dei bin di tessuto con d_nuc più alto azzerati): E ricalcolato 0.0688 contro 0.0688495 riportato (C5 A_E = PH-3 A f=0.10); frazione azzerata, φ_obs e BA identici con la mappatura della pipeline.

**Evidenza.** recalc_ctrlA_A3r1_<map>.csv. Definizione d_nuc > r90: frac_zeroed 0.0999146 (identica), φ_obs 0.0957776 (identica), BA 0.850557 (identica), E 0.0686801 (z = −1.51, solo seed); mappatura R2 E 0.0688175. Lettura letterale «esattamente il 10 %» (pareggi casuali): frazione 0.0999988, E 0.0687044–0.0688867. Il verdetto A (E ≥ 0.05, BA ≥ 0.75) è confermato.

**Azione suggerita.** Nessuna. Nota: C5 A_E e PH-3 (A, f=0.10) sono lo stesso numero dallo stesso codice, non due repliche indipendenti; nel report non vanno presentati come conferme distinte.

### RC-6 — metodo
**Affermazione.** d_nuc e maschere vengono campionati nel pixel floor(pxl − c0); con la convenzione dei centroidi R2 (x_um = indice_colonna·µm/px, centri dei pixel sugli interi) il pixel che contiene il centro del bin è round(·). La scelta non è dichiarata e sposta il campionamento di ≤ 1 px (0.27 µm).

**Evidenza.** tools/S1.3c_ecm.py r. 177–178 (floor); R2_common.nuclei_table usa round per i centroidi. Il mio ricalcolo con rint (e bin oltre la copertura ds4 esclusi) cambia R di al massimo 0.0041 (A5 r1: 0.5653 contro 0.5612) e φ_obs di al massimo 1.5e-4. Effetto trascurabile.

**Azione suggerita.** Dichiarare la convenzione dei pixel (floor) nelle note di sviluppo; nessun rerun necessario.

### RC-7 — metodo
**Affermazione.** BA e r* sono instabili quando le finestre vuote sono poche: una sola finestra in più o in meno sposta BA di 0.035 e r* di 2.75 µm. La BA mediana per archetipo dà a questi ROI lo stesso peso di quelli con centinaia di finestre vuote.

**Evidenza.** A2 r1 k=5: n_empty 9 (riportato) contro 8 (mappatura R2): BA 0.7219 → 0.6871, r* 42.75 → 40.0. A1 r1 k=3: n_empty = 2, BA 0.49995. Nei dati osservati non ci sono effetti sui verdetti, perché E < 0.05 viene prima del criterio di forma; i controlli C-5 hanno centinaia di finestre vuote.

**Azione suggerita.** Nel report, riportare n_empty accanto a BA; valutare un minimo di finestre vuote (per esempio ≥ 30) sotto il quale la BA del ROI è 'non stimabile'. Sarebbe una regola post hoc: va dichiarata come tale.

## 5. Errori del revisore (trasparenza)
1. Prima versione di `recalc_umi_bg.py`: i barcode assenti dalla matrice grezza (Small Intestine, SR 3.0.0) restavano NaN, quindi soglia NaN e ROI A1 r1 fallito (l'assert di conservazione dei totali ha intercettato l'errore). Corretto: mancanti = 0.
2. Prima versione di `recalc_roi.py` (`rint_noclip`): i bin oltre la copertura della ds4 (pixel ≥ 4·913) erano esclusi dal tessuto e quindi contati come «lume» in CP-2 (A3 r1: 1209 contro 876 riportati, +333). Era un errore mio: R2 estende l'ultimo pixel (`up()` in `R2_common.py`). Corretto nella variante `floor_clip`.
3. Nella colonna z le differenze sono calcolate con la variante `px4` (stessi bin della pipeline), così riflettono solo la variabilità di seed.

## Sintesi
Ricalcolo indipendente di sfondo/UMI (5 dataset), dei 6 ROI r1 (k = 3, 5) e del controllo positivo A3 r1. Sfondo, C-3 e UMI: identici. Con la stessa convenzione bin→maschera della pipeline, 141/141 grandezze deterministiche sono identiche; nullo ed E coincidono entro la variabilità di seed (|z| ≤ 2.2 con 30 permutazioni, ≤ 1.93 con 200). Unico difetto: la mappatura sulle maschere ds4 usa px4 = side_um/913 invece di 4 px nativi (disallineamento ≤ 0.55 µm, |ΔE| ≤ 1.6e-4, nessun verdetto cambia). Controllo positivo A3 r1: E 0.0688 confermato. Rilievi di metodo: convenzione floor dei pixel non dichiarata; BA instabile con poche finestre vuote.
