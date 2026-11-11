
Tabelle in `results/S2.8/`: `S2.8_G1_summary.txt` e `S2.8_G1_test_S2.*.txt`, `S2.8_mutants.txt`, `S2.8_G2_identity.csv`, `S2.8_design_map.csv`, `S2.8_design_map_actions.csv`, `S2.8_exported_args.csv`, `S2.8_md5_before.txt`, `S2.8_md5_after.txt`, `S2.8_R_changed.txt`, `S2.8_runs.csv`, `S2.8_summary.csv`, `S2.8_predictions.csv`, `S2.8_rule_R1.csv`, `S2.8_rule_R1.json`, `S2.8_rule_R1_sens_k12.csv`, `S2.8_times_by_seed.csv`, `S2.8_BL237.csv`. Uscite grezze in `/mnt/micron/geo_spatialtrans/S2.8/` (81 esecuzioni, `sched.log`). Codice delle esecuzioni: commit `1ecf81c` (le modifiche successive, `21fc08f`, sono solo commenti e test).

### Cancelli

- **G-1** superato: sul codice delle esecuzioni `test_S2.8.R` 43/43 e mutanti 3/3; sul codice finale, dopo le correzioni della revisione, `test_S2.8.R` 46/46 (gruppi a, b, c 21; gruppo d 25) e mutanti 4/4 (dopo la correzione di due modelli di sostituzione rimasti alla riga di prima, D-8); batterie di S2.1...S2.7 con le stesse asserzioni di S2.7 (21, 36, 39, 37, 40, 33, 78), nessun FAIL, nessuna asserzione cambiata.
- **G-2** superato: 1 080 confronti di impronte su 1 080 come attesi (E504 contro S2.7: 324; E1008 ed EInf contro E504: 342 + 342; motore di S2.5 contro le uscite salvate di S2.5 e S2.7: 72). **G-2 bis** con il codice finale (`64cf92f`): le 6 esecuzioni E504 con k = 0 rifatte, 84 impronte su 84 uguali a S2.7 (`results/S2.8/S2.8_G2bis_final_code.csv`).
- **G-3** superato: mappa di 187 righe (124 coerenti, 37 documento indietro, 19 senza test, 5 codice diverso, 2 non chiare); le 63 righe non coerenti hanno tutte un'azione (39 annotazioni nella v0.6, 18 test aggiunti, 2 difetti corretti con test, 1 documentazione corretta con test, 1 documentata e misurata, 2 voci di BACKLOG); documento di progetto v0.6; controllo OKF 0 documenti in errore.
- **G-4** superato: 375 file letti identici prima e dopo; `R/` cambiato solo in `05b5`, `05b6` (perimetro), `05b2`, `05b3` (C-3, solo commenti) e nei test; `R/04b*.R` identici.
- **G-5** superato dopo le correzioni: due revisori, nessuna discrepanza numerica nei cancelli, nelle previsioni e nella regola; un difetto del codice (impronta senza il livello cellulare) e testi corretti (sotto «Revisione avversariale»).

### Previsioni

| id | giudizi | conf. | indet. | smentite | in breve |
|---|---|---|---|---|---|
| P-1 tempo della funzione a pezzi con 504 µm su 1 mm² | 6 | 6 | 0 | 0 | da 85 s (cuore) a 270 s (tumore del colon), 2 errori standard al più 5.4 s |
| P-2 picco di memoria con 504 µm | 6 | 6 | 0 | 0 | da 1.3 GB (cuore) a 3.2 GB (tumore del colon) |
| P-3 picco di memoria in un pezzo solo | 6 | 6 | 0 | 0 | da 1.7 GB (cuore) a 6.2 GB (tumore del colon); cervello 4.8 GB |

Attese sbagliate (le previsioni sono giudicate sugli intervalli, che restano quelli scritti): il tumore in P-1 e tumore e cervello in P-3 erano attesi fuori dal bilancio. P-1: i 385 s di S2.7 erano misurati con 16 processi sulla macchina; con 6 il tumore sta a 270 s. P-3: la stima di 9-10 GB scalava il picco dei quadrati di 4 mm² con il numero di quadratini della griglia (0.375), ma la memoria segue i trascritti, cioè il tessuto: 1 mm² contro 4 mm² (0.25), cioè circa 6.8 GB nel tumore, vicino ai 6.2 misurati.

### Regola dell'ingresso (R-1)

| tessuto | motore di S2.5 (s) | migliore a pezzi | tempo (s) | rapporto | con 504 µm |
|---|---|---|---|---|---|
| intestino | 160.9 | pezzo solo | 180.0 | 1.12 | 1.35 |
| tumore del colon | 212.1 | pezzo solo | 217.7 | 1.03 | 1.27 |
| stroma | 123.6 | pezzo solo | 132.3 | 1.07 | 1.28 |
| linfonodo | 173.2 | pezzo solo | 189.8 | 1.10 | 1.19 |
| cervello | 146.0 | pezzo solo | 180.1 | 1.23 | 1.51 |
| cuore | 78.8 | pezzo solo | 76.5 | 0.97 | 1.08 |

Rapporto ≤ 1.5 in 6 tessuti su 6: **un solo ingresso, `simulate_visium_hd_pieces()`**; il migliore non è 504 µm in nessun tessuto, quindi la scelta automatica del lato è una voce di BACKLOG (BL-252). Decisione tecnica: [S2.8: l'ingresso del simulatore dello Step 2 è la funzione a pezzi; la funzione di S2.5 resta come motore di riferimento (DEC-079)](../docs/decisions/DEC-079.md). Picco del motore di S2.5 con blocchi di 4 × 10⁶: da 3.4 a 7.4 GB, più del pezzo solo della funzione a pezzi.

### Controprova

- **CP-1** superata: la mappa, fatta senza sapere quali sezioni ci aspettavamo, segna come indietro le sezioni 1 (M-005, M-015, M-016), 12 (M-134...M-141) e 15 (M-163), più la 16 e, in parte, altre nove.
- **CP-2** superata: con `ceiling_q = 1e-4` l'impronta di `counts_2um` cambia in 6 tessuti su 6.

### Misure (senza intervallo)

- Tempi medi (s) per configurazione e tessuto, nell'ordine intestino, tumore, stroma, linfonodo, cervello, cuore: pezzi di 1 008 µm 211, 259, 154, 206, 221, 81; picchi (GB) motore di S2.5 5.8, 7.4, 5.3, 3.4, 5.9, 4.9; pezzi di 1 008 µm 3.7, 4.7, 2.8, 2.3, 4.8, 1.5 (`S2.8_summary.csv`).
- **Descrittiva di BL-237** (intestino A1_r1, 3 semi appaiati, finestra di R2 contro riquadro predefinito; distanze dal bordo sinistro del tessuto riga per riga, definizione corretta dopo la revisione dei numeri, vedi «Revisione avversariale»). Il bordo sinistro del tessuto è irregolare: il primo quadratino di tessuto di una riga sta da 6.5 a 184 µm dal lato della finestra (mediana 54 µm), e tocca il riquadro predefinito (x = 5.5 µm) solo nel 2.4 % delle righe. UMI medio per quadratino nel tessuto: 29.587 con il predefinito, 29.584 con la finestra (rapporto da 0.9998 a 1.0000 fra semi); nelle fasce di 0-8, 8-16, 16-32, 32-64 µm dal bordo (da 1 846 a 7 188 quadratini) rapporti da 0.997 a 1.002; vetro entro 8 µm dal bordo del tessuto 0.998; vetro fra il lato della finestra e il tessuto (17 231 quadratini) 1.02; striscia di vetro larga 5.5 µm lungo il lato della finestra, fuori dal riquadro predefinito (2 352 quadratini), da 0.38 a 0.98 UMI per quadratino (rapporto da 2.51 a 2.79): con il predefinito la seconda coda si ripiega a 5.5 µm e la striscia non la riceve. Vetro fuori dalla finestra 0.030 in entrambi. Conclusione: il riquadro conta solo su quella striscia; nel tessuto e nel vetro vicino al tessuto la differenza è sotto il 2 %. Il vero della striscia non è noto (D-3).

