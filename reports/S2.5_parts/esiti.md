Tabelle: `results/S2.5/S2.5_verdicts_<variante>.csv` (una riga per giudizio; varianti main, main_r2, depth_off), `S2.5_verdict_counts_<variante>.csv`, `S2.5_sensitivity_<variante>.csv`, `S2.5_measures_<variante>.csv`, `S2.5_R2_depth.csv`, `S2.5_real_refs.csv`, `S2.5_real_areas.csv`, `S2.5_sim_<variante>.csv`, `S2.5_G2.csv`, `S2.5_G2_truth8.csv`.

### Regole

- **R-1, blocchi**: 4 × 10⁶ trascritti per blocco; cellule per blocco da 93 (cervello) a 2 000 (linfonodo, tetto) nel pannello main; P-8 sotto. Dopo PH-1 il valore predefinito della funzione diventa 2 × 10⁶ (correzione C-6; il pannello resta quello con 4 × 10⁶).
- **R-2, fattore di profondità sulla catena completa** (`S2.5_R2_depth.csv`): f di S2.2 → f': intestino 0.867 → 1.013, tumore 0.893 → 0.954, stroma 1.584 → 1.561, linfonodo 0.787 → 0.819, cervello 1.252 → 1.428, cuore 2.625 → 2.939 (rapporto (U_vero − b)/(U_sim − b) da 0.986 a 1.168). I fattori ritarati diventano quelli del preset (`visium_hd_depth_factor()`, decisione tecnica [DEC-067](../docs/decisions/DEC-067.md)).

### Cancelli

| id | esito | dettaglio |
|---|---|---|
| G-1 | superato | batteria 40/40 (`S2.5_G1_tests.txt`; 35 asserzioni pre-registrate più 5 aggiunte dopo la revisione: e6, f5b, f6b, d5, d6); mutanti 16/16 rilevati (`S2.5_mutants.txt`; prima esecuzione 11/12, un frammento ambiguo nello strumento dei mutanti, correzione C-1; un mutante diventato equivalente dopo la correzione C-2, sostituito, C-5) |
| G-2 | superato | 150 esecuzioni (main 60, depth_off 30, main_r2 60): conservazione 150/150, nessun NA, catturati |z| massimo 3.54 e fondo 2.39 contro la soglia 4.59 (`S2.5_G2.txt`); verità a 8 µm per coppia (blocco, cellula) esatta in 150/150 per origine, sovrapposizione, frazioni e fondo (`S2.5_G2_truth8.txt`; controllo rafforzato dopo la revisione, C-4) |
| G-3 | superato | lato vero, pannello main k = 0 e misure rifatti in un'altra cartella con 2 e 3 processi invece di 3 e 6, con il codice dopo le correzioni: 241 file identici (`S2.5_G3.txt`) |
| G-4 | superato | (a) 230 file letti, md5 identici prima e dopo (`S2.5_md5_diff.txt` vuoto); (b) batterie S2.1 21/21, S2.2 36/36, S2.3 39/39, S2.4 37/37; (c) A2_r1 e A6_r1 con i semi e i blocchi di S2.4: totali per quadratino identici al pannello di S2.4 (`S2.5_G4_identity.txt`); (d) definizione sr = `R4_roi.csv`, 180 valori, scarto massimo 1.1e-13 (`S2.5_real_R4_reproduction.csv`) |
| G-5 | superato, con rilievi | tre revisori («Revisione avversariale»): ricalcoli entro la tolleranza; due difetti del codice (generatore dell'utente senza seme, strati senza nome) corretti senza effetto sui risultati; tre controlli deboli rafforzati; due esiti del cervello nel pannello main dipendono dal seme del bootstrap |

### Previsioni

Conteggi per esito (confermata / indeterminata / smentita):

| previsione | pannello main | dopo R-2 (main_r2) |
|---|---|---|
| P-1 UMI per quadratino di 2 µm | 4 / 1 / 1 | 6 / 0 / 0 |
| P-2 quadratini di 2 µm vuoti | 4 / 0 / 2 | 4 / 0 / 2 |
| P-3 geni per quadratino di 2 µm | 4 / 1 / 1 | 6 / 0 / 0 |
| P-4 UMI per blocco di 8 µm | 4 / 1 / 1 | 6 / 0 / 0 |
| P-5 blocchi di 8 µm vuoti | 6 / 0 / 0 | 6 / 0 / 0 |
| P-6 geni per blocco di 8 µm | 5 / 1 / 0 | 6 / 0 / 0 |
| P-7 rapporto nucleo/corona | 3 / 2 / 1 | 3 / 2 / 1 |
| P-8 tempo e memoria (12 giudizi) | 12 / 0 / 0 | 11 / 1 / 0 |
| P-9 orientamento | 2 / 3 / 1 | 2 / 3 / 1 |

- **P-1, P-3, P-4, P-6 (pannello main)**: intestino smentita come atteso (UMI a 2 µm 25.06 [24.72, 25.49] contro 26.21-32.34; geni a 2 µm 24.04 contro 24.46-29.83; UMI a 8 µm 398.9 contro 419.0-514.4; geni a 8 µm confermata, 328.3 contro 323.8-385.3); cervello indeterminata (attesa smentita): UMI a 2 µm 39.07 [37.77, 40.56] contro 40.43-47.77. Dopo R-2 tutte confermate in 6 su 6 (UMI a 2 µm simulato / vero 0.9999-1.0001: per costruzione, il fattore si ritara sulle stesse aree e si verifica con gli stessi semi; non per costruzione i geni e i blocchi di 8 µm: geni a 8 µm simulato / vero 1.005-1.064).
- **P-2**: smentita per difetto nell'intestino (dopo R-2 0.00040 contro 0.0047-0.0093, un diciassettesimo) e nel tumore (0.00042 contro 0.0041-0.0071, un tredicesimo), come attesa; confermata negli altri quattro (stroma 0.42 del vero, intervallo largo).
- **P-5**: confermata 6/6 (simulato 0 blocchi vuoti in tutte le esecuzioni).
- **P-7**: smentita per eccesso nel tumore del colon (dopo R-2 0.901 [0.888, 0.917] contro 0.818-0.883), come attesa; indeterminata nell'intestino (0.955 contro 0.922-0.951) e nel cervello (1.363 contro 1.357-1.411); confermata in stroma, linfonodo, cuore. Con la banda di ±10 % di S2.3 (sensibilità) confermata in 6 su 6.
- **P-8** (giudicata come pre-registrata: media delle 5 aree ± 2 errori standard per tessuto contro 0-300 s e 0-8 GB): pannello main 12/12 confermate (singole esecuzioni 19-237 s; memoria media del tumore 7.09 GB, ± 2 SE fino a 7.71 GB); main_r2 tempo del tumore indeterminato (276.8 [238.9, 314.8] s, server carico per la revisione e i mutanti). La regola non giudica le singole esecuzioni: oltre 8 GB ce ne sono 3 su 150, tutte nell'area A2_r4 (8 318 MB main, 8 749 MB main_r2, 9 858 MB depth_off). Rimisurato a macchina libera con i predefiniti finali (PH-2): 149-232 s e 6.0-7.2 GB nelle 5 aree del tumore.
- **P-9**: cuore smentita per difetto (dopo R-2 0.048 [0.035, 0.059] contro 0.086-0.41), come attesa (CP-2); indeterminata in intestino (0.203 contro 0.105-0.241), tumore (0.159 contro 0-0.177) e linfonodo (0.216 contro 0-0.184), dove il simulato supera o tocca il vero (vedi «Revisione avversariale» e BL-236); confermata in stroma e cervello. Con `kappa_boot_box` (sensibilità, C-3) nessun esito cambia.
- **Sensibilità con la tabella adottata di R7** (definizione sr, `S2.5_sensitivity_main_r2.csv`): dopo R-2 l'intestino supera i valori di Space Ranger (UMI a 2 e 8 µm indeterminate, geni a 2 e 8 µm smentiti per eccesso; rapporto nucleo/corona smentito per difetto, 0.943-0.965 contro 0.965-1.006), perché R-2 ritara sui quadratini delle regioni, che hanno più trascritti di quelli di Space Ranger (29.27 contro 27.14). Cambiano anche: vuoti a 2 µm del tumore confermati (l'intervallo di Space Ranger arriva a 0), geni a 2 e 8 µm del cervello indeterminati. Gli altri esiti come con la definizione reg. Mostra di quanto la definizione del tessuto pesa nel confronto ([DEC-064](../docs/decisions/DEC-064.md)).

### Controprove

- **CP-1** (`S2.5_verdicts_depth_off.csv`): con il fattore di profondità 1 l'UMI per quadratino di 2 µm è smentito nel linfonodo (14.08 contro 8.76-13.56), nel cervello (31.23 contro 40.43-47.77) e nel cuore (1.21 contro 2.06-5.05): superata, il confronto ha potere in quei tre tessuti. Intestino, tumore e stroma confermati, come atteso (nell'intestino senza fattore 28.80, più vicino al vero che con il fattore di S2.2).
- **CP-2**: P-9 smentita per difetto nel cuore nel pannello main e in main_r2: superata, la misura vede il limite noto delle cellule isotrope.

### Misure

- Fattore di sovra-generazione M per esecuzione: 2.93 (cuore) ... 5.94 (tumore); media per tessuto 3.20 (cuore) ... 5.00 (stroma). Valore atteso dei catturati / bersaglio di R4 0.997-1.004 per esecuzione (pannello main); catturati entro |z| ≤ 2.19 dall'atteso.
- BL-212: quadratini in tessuto senza territorio 155 su 6 616 090 (2.3e-5): 141 nell'area A5_r1 del cervello, 14 in tre aree del cuore, 0 altrove.
- BL-219: quota degli UMI veri dello stroma sui 522 geni con profilo nullo 0.99 % [0.45 %, 1.67 %] (E-3; limite dichiarato, [DEC-065](../docs/decisions/DEC-065.md)).
