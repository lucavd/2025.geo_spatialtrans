# S1.3b — revisione avversariale, revisore «code» (ricalcolo e codice)

Data: 2026-10-10. Repo `~/2025.geo_spatialtrans`, branch `step1-cell-layer`, commit `5ca9e36` (in sola lettura).
Uscite della revisione: questa cartella e `/mnt/micron/geo_spatialtrans/S1.3b/review_code/` (`regen_audit.rds`, `s13out/` vuota).
Gli stadi ufficiali non sono stati rilanciati.

## Sintesi
Ho ricalcolato in modo indipendente dai RDS grezzi tutti i verdetti PASS/FAIL/WARN di S1.3b, senza usare tools/S1.3b_analysis.R: C-0…C-6, le 10 coppie V0, mutanti, D-0, CP-1a/b/c, CP-2, B-1, B-2 e la decisione. Tutti i valori numerici del CSV dei verdetti coincidono e la decisione è V2, secondo la regola pre-registrata. Due rilievi cambiano un verdetto. (1) Il placebo CP-1b è nullo: aggiunge 8.9e-16 µm a coordinate vicine a 1000 µm, dove un ulp vale 1.1e-13, quindi nessuna coordinata cambia (0/20 seed) e il 20/20 è garantito. Con un placebo corretto, cioè lo stesso numero di ulp, il risultato regge: 20/20 su 766/1987 e su ciascuna delle 10 coppie. (2) Il PASS di CP-2 per V0 esiste solo grazie alla regola post hoc 'eccezione = FAIL'. Con la definizione pre-registrata V0 ha 0/180 FAIL, quindi la previsione non è confermata; la potenza era comunque 0.395. Controprove aggiunte dalla revisione rafforzano H1 nella forma 'vertice condiviso non identico sul bordo della regione'. Fondere i soli vertici non identici, senza giunzioni a T, corregge 15/15 unioni. Dallo stato annodato, 1 ulp sul vertice di bordo riproduce l'errore in 24/40 casi, 1 ulp su un vertice interno in 0/80. La condizione però non è sufficiente da sola: l'83 % dei nulli che passano ha vertici di questo tipo. 'La causa è isolata' va quindi riformulata come condizione necessaria e innesco riproducibile; il meccanismo interno di GEOS non è identificato. Su Step 2: C-6 comprendeva solo 1 dei 10 nulli con intersezione di coppia sbagliata. Esteso ai 10 dalla revisione, dà 0 violazioni sia in V0 sia in V2. 'Step 2 non esposto' è sostenibile solo come 'nessuna esposizione osservata', per la griglia da 2 µm e senza un controllo positivo di C-6. La scelta di V2 dipende soltanto da C-4: 2 celle invalide su 14 469 893. Senza C-4 sarebbe stata scelta V1.

## Script (tutti `R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla <script>` dalla root del repo)
| script | cosa fa | log |
|---|---|---|
| `recompute.R` | legge i 1 800 RDS di `null/` (9 000 righe variante×nullo → `null_long.csv`); C-0 contro `results/S1.3/S1.3_null_summary.csv`; C-1…C-5, B-1, B-2 (sommari), CP-3, tempi | `recompute.log` |
| `recompute2.R` | C-6 (`grid/`), mutanti, CP-1 (`cp1.rds`), CP-2 con e senza la regola post hoc, B-2 (verdetti, ingressi dei RDS modificati, margini), regola di decisione e sensibilità | `recompute2.log` |
| `regen_audit.R` | rigenera i 10 nulli con coppia V0 sbagliata: firma del vertice (coordinate, ulp, distanza dal bordo), CP-1a per coppia, placebo come quello ufficiale (verifica che le coordinate cambino), placebo corretto in ulp, C-6 indipendente V0/V2 | `regen_audit.log` |
| `union_audit.R` | 15 unioni sbagliate: sola fusione dei vertici non identici (senza giunzioni a T) e V3 | `union_audit.log`, `union_audit.csv` |
| `inverse_audit.R` | dallo stato annodato, ±1 e ±4 ulp sul vertice di bordo e su un vertice condiviso interno (80+80 prove); celle invalide P2 | `inverse_audit.log`, `inverse_audit.csv` |
| `invalid_cells.R` | celle invalide di V1/V3: vertici non adiacenti quasi coincidenti, con distinzione degli anelli | `invalid_cells.log` |

## Tabelle
- `findings.csv`: 14 rilievi (bug 1, discrepanza 1, metodo 7, confermato 5).
- `recomputed.csv`: 44 grandezze, valore del report accanto al valore ricalcolato. Non coincidono in 3: CP-1b, CP-2 V0 con la definizione pre-registrata, copertura di C-6.

## Rilievi principali
1. **RC-01 (bug) CP-1b, placebo nullo.** `x + 8.88e-16·cos θ` su coordinate di circa 1000 µm non cambia nulla: coordinate cambiate in 0/20 seed e aree identiche bit per bit 20/20. Placebo corretto (stesso numero di ulp): 20/20 errore intatto, con coordinate cambiate 20/20, su 766/1987 e su ciascuna delle 10 coppie.
2. **RC-02 (discrepanza) CP-2 V0.** Con la definizione pre-registrata la frequenza di FAIL è 0.000/0.000/0.000 → non confermata. Il PASS (0.050/0.083/0.100) esiste solo con «eccezione = FAIL», e le 14 eccezioni cadono tutte su input invalidi. La regola è stata committata alle 22:43:27, dopo che tutte le uscite null/grid/cp1 erano su disco.
3. **RC-03/RC-04 (causa).** H1, nella forma «vertice non identico sul bordo», è sostenuta meglio di quanto mostri il report: fusione senza T 15/15, inverso 24/40 sul bordo contro 0/80 all'interno, annodamento 10/10. Non è però sufficiente da sola (1 485/1 785 nulli che passano hanno vertici diversi) e l'esito è sensibile all'intera configurazione (placebo stile ufficiale: 15-17/20 in 3 coppie). Va riformulata «la causa è isolata».
4. **RC-05 (Step 2).** C-6 ufficiale: 1/10 nulli con coppia sbagliata; estensione: 0 violazioni in 10/10 (V0 e V2). C-6 non ha un mutante. Solo griglia da 2 µm.
5. **RC-06 (decisione).** V2 per la regola; l'unica discriminante è C-4 (2 celle su 14 469 893); senza C-4 → V1. Per V2 la validità è data per costruzione.

## Deviazioni / limiti di questa revisione
- C-unit non rilanciato: `R/testing/test_S1.3b.R` scrive in `results/S1.3b/S1.3b_test.csv`, fuori dal perimetro di scrittura. Ho letto il CSV (21/21, 21/21, 23/23).
- Il meccanismo dell'invalidità della cella 6 di A5 r1 CSR 12 non è identificato (RC-07, provvisorio).
- Il test inverso perturba un solo vertice interno per coppia (il primo vertice condiviso identico).
