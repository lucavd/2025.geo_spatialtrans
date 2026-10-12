# S1.3b — revisione avversariale, revisore 2 «exact» (aritmetica esatta e terzo motore)

Data: 2026-10-10 · repo `~/2025.geo_spatialtrans`, branch `step1-cell-layer`, commit `5ca9e36` (letto, non modificato).
Ambiente: R 4.6.1, sf 1.1.3, **GEOS 3.12.1** (motore del run); Python 3.12 + **shapely 2.2.0 / GEOS 3.14.1** (terzo motore, venv in
`/mnt/micron/geo_spatialtrans/S1.3b/review_exact/.venv`); aritmetica razionale `fractions.Fraction`.

## Esito in breve
Aritmetica razionale esatta (Fraction dai double esportati in %a, validata su 13 casi noti + 300 casuali) sulle 10 coppie P1: st_relate II=2 e' corretto in 10/10, con schegge reali di 2.1e-15–2.4e-14 µm², pari al triangolo (vertice condiviso, due piedi diversi di 1–4 ulp sul bordo della regione). st_intersection(i,j) restituisce invece l'intera cella i (errore 65.7–844.8 µm²). GEOS 3.14.1 ripete l'errore identico, mentre grid_size=1e-6 lo elimina. L'errore dipende dall'ordine degli argomenti (j,i corretto in 8/10). In A5 r5 RSA 15 la tassellazione e' una partizione entro 1.7e-14 µm². Il deficit di st_union (1358.5733 µm²) e' la cella 71 persa per intero, a causa di una coppia in configurazione a lacuna (intersezione esatta 0, relate 'touch') che C-2 non puo' vedere: le coppie sbagliate sono quindi almeno 11, non 10. Nessun verdetto cambia; servono correzioni di formulazione (X3, X4, X9).

## Metodo
1. `export_exact.R` rigenera gli 11 nulli (le 10 coppie P1 + A5 r5 RSA 15) con `regen()` di `tools/S1.3b_run.R`, nello stesso modo di
   `tools/S1.3b_posthoc.R`. Esporta in esadecimale (`sprintf('%a')`, bit-esatto) le celle i, j, la regione `w$poly[[1]]`, la geometria
   restituita da `st_intersection`, i metadati sf (relate, aree, WKB) e, per A5 r5 RSA 15, tutte le 929 celle.
   Uscite: `/mnt/micron/geo_spatialtrans/S1.3b/review_exact/{P1_pairs_coords_hex.csv, P1_pairs_sf_meta.csv, A5_r5_RSA_15_coords_hex.csv, A5_r5_RSA_15_sf_meta.csv}`.
   Controllo: hex == decimale `%.17g` per tutte le righe (0 discrepanze, in Python e in R).
2. `exactgeom.py` calcola l'area esatta di P∩Q con il teorema di Green sulla 1-catena del bordo: lati di ∂P spezzati in tutti i punti di
   contatto con ∂Q, tenuti se il punto medio è interno a Q; lo stesso per ∂Q; i lati collineari sono contati una volta, se le orientazioni
   concordano. Il punto-in-poligono è esatto, con lo stato «sul bordo». Controprova indipendente: Sutherland–Hodgman esatto, quando una
   delle due celle è convessa.
   `test_exactgeom.py` contiene 13 casi noti, tutti PASS: quadrati adiacenti, sovrapposti, identici e annidati; lati parziali comuni con
   interni opposti o concordi; scheggia di 1 ulp (area 2^-53) e lacuna di 1 ulp (area 0); L concava; buco; multipoligono. A questi si
   aggiungono 300 coppie casuali confrontate con shapely (differenza massima 1.7e-16).
3. `p1_exact.py`: per ogni coppia calcola aree esatte, intersezione esatta (in entrambi gli ordini e con S-H), vertici condivisi
   identici o diversi, ulp, distanza esatta dal bordo della regione, contenuto della geometria restituita da sf, e shapely in virgola
   mobile e con `grid_size=1e-6` / `set_precision`. `p1_triangle.py` confronta l'intersezione esatta con l'area del triangolo
   (vertice identico, piede i, piede j).
4. `sf_crosscheck.R` ricostruisce le geometrie dall'esadecimale in sf (senza rigenerare) e calcola `st_intersection` in entrambi gli
   ordini, `st_difference`, `st_union`, `st_set_precision(1e6)`. Per A5 r5 RSA 15 cerca le celle che mancano dall'unione.
5. `a5_exact.py` (A5 r5 RSA 15) calcola in modo esatto: somma delle aree contro area della regione (1 parte, 119 buchi, 9 418 vertici);
   intersezione esatta di tutte le 2 591 coppie che si toccano; contenimento esatto nella regione delle 325 celle di bordo.
   `a5_pair71.py` analizza la coppia 71/175; `a5_sweep.py` passa in rassegna con shapely tutte le coppie.

## Risultati per coppia
| coppia | i | j | area_i (µm²) | ∩ esatta (µm²) | ulp (x,y) | dist (µm) | su bordo regione (esatto) | sf ∩(i,j) | sf ∩(j,i) | shapely 3.14.1 ∩(i,j) | grid 1e-6 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 r3 RSArule 8 | 3780 | 5886 | 99.3963 | 3.190e-15 | (2,0) | 4.44e-16 | sì | 99.3963 | 3.55e-15 | 99.3963 | 0 |
| A2 r4 RSA 12 | 1884 | 5121 | 65.7719 | 2.192e-15 | (0,1) | 4.44e-16 | sì | 65.7719 | 2.19e-15 | 65.7719 | 0 |
| A2 r5 RSArule 19 | 1103 | 1714 | 65.7041 | 2.060e-15 | (1,0) | 4.44e-16 | sì | 65.7041 | **65.7041** | 65.7041 | 0 |
| A3 r3 CSR 3 | 766 | 1987 | 340.427 | 1.189e-14 | (0,2) | 8.88e-16 | sì | 340.427 | 4.44e-13 | 340.427 | 0 |
| A3 r3 CSR 5 | 3479 | 4200 | 140.716 | 3.957e-15 | (1,0) | 4.44e-16 | sì | 140.716 | 7.11e-15 | 140.716 | 0 |
| A3 r3 RSA 8 | 929 | 2686 | 223.841 | 5.061e-15 | (0,1) | 4.44e-16 | sì | 223.841 | 2.22e-13 | 223.841 | 0 |
| A3 r5 RSArule 2 | 705 | 2172 | 351.513 | 8.302e-15 | (0,1) | 8.88e-16 | sì | 351.513 | 8.30e-15 | 351.513 | 0 |
| A5 r2 CSR 16 | 197 | 335 | 635.480 | 2.149e-14 | (0,1) | 1.78e-15 | sì | 635.480 | 2.15e-14 | 635.480 | 0 |
| A6 r2 CSR 2 | 365 | 408 | 695.615 | 2.224e-14 | (4,0) | 1.78e-15 | sì | 695.615 | **695.615** | 695.615 | 0 |
| A6 r5 RSArule 19 | 640 | 917 | 844.819 | 2.396e-14 | (0,1) | 1.78e-15 | sì | 844.819 | 2.40e-14 | 844.819 | 0 |
| A5 r5 RSA 15 (non in P1) | 71 | 175 | 1358.57 | **0** (lacuna) | (1,0) | 8.88e-16 | sì | 1442.84 | 0 (shapely) | 1442.84 | 0 |

Il grassetto in «sf ∩(j,i)» indica gli ordini in cui anche l'intersezione inversa restituisce la cella intera. Lo stesso schema 8/2 si
ritrova in GEOS 3.14.1. Valori razionali esatti e tutte le colonne in [P1_exact_results.csv](P1_exact_results.csv).

## Rilievi
Vedi [findings.csv](findings.csv); le quantità ricalcolate sono in [recomputed.csv](recomputed.csv).
- **X1, confermato.** II=2 è corretto: la scheggia è reale ma vale 1e-15–1e-14 µm².
- **X2, confermato.** `st_intersection` sbaglia di un'intera cella.
- **X3, metodo.** C-2 non vede la configurazione a lacuna dello stesso difetto (71/175), quindi le coppie sbagliate sono ≥ 11, non 10.
- **X4, discrepanza.** L'errore dipende dall'ordine degli argomenti.
- **X5, confermato.** GEOS 3.14.1 ripete l'errore; con precisione 1e-6 l'errore sparisce.
- **X6, confermato.** area_v0 = area_i.
- **X7, confermato.** 1–4 ulp; vertice esattamente sul bordo della regione.
- **X8, confermato.** A5 r5 RSA 15 è una partizione entro 1.7e-14 µm²; il deficit è la cella 71 intera.
- **X9, discrepanza.** Il report scrive «−1 358.6», mentre il dato è un deficit di +1358.57.
- **X10, discrepanza.** La voce B/B di relate in GEOS 3.12.1 è sbagliata su 71/175.

## Limiti
- Il contenimento nella regione di A5 r5 RSA 15 è verificato in modo esatto solo sulle 325 celle a meno di 1e-3 µm dal bordo
  (distanza calcolata con shapely). Le altre sono lontane dal bordo almeno di 1e-3 µm.
- Le coppie candidate in A5 r5 RSA 15 sono scelte con il predicato `intersects` di shapely in virgola mobile (2 591 coppie).
- I nulli diversi dagli 11 elencati non sono stati esaminati, quindi «≥ 11» è un limite inferiore.
