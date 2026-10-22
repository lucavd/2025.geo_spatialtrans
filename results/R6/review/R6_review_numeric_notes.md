# R6 revisione numerica indipendente (sotto-agente revisore, 2026-10-22)

Script: /tmp/r6_review.py (pandas puro, nessun import da tools/), eseguito con ~/2025.geo_spatialtrans/.venv/bin/python.
Input: /mnt/micron/geo_spatialtrans/R6/labels/<A>_exp_umi20_seed20261022_nuclei{,_ph1}.parquet (A1-A6), cell_types.csv di Space Ranger (10x Cloud) per A5 e A6.
Confronto con results/R6/R6_refs_type{,_ph1}.csv (n, fraction, median_area_nuc, median_ecc_nuc per archetipo x tipo) e R6_predictions{,_ph1}.csv (P-1, P-2, P-3, P-6 informativo). Tolleranza 1e-6.
Tabella: R6_review_numeric.csv (316 confronti, 3 discrepanze).

## Discrepanza 1 e 2: P-6 in A6 (prereg e ph1), accordo con Space Ranger "matched, informative"
- mio: 0.7738028169014084 (prereg), 0.7943374858437147 (ph1); loro: 0.7502673796791444, 0.7765685019206146.
- Causa (verificata, differenza 0.0 riproducendo il loro numero): in tools/R6_analysis.py la riga
  `~d.sr_lowumi.fillna(True)` (funzioni agreement_all, agreement_informative, tabelle di confusione) non esclude nulla.
  Dopo il left merge la colonna sr_lowumi ha dtype object (NaN per i 33 nuclei di A6 senza sr_cell_id);
  `~` sui bool di Python restituisce -2/-1 e pandas 3.0.6, nell'operazione `bool_series & object_series`,
  converte entrambi a True. Risultato: i 421 nuclei di A6 appaiati a barcode con coarse 'Low UMI Barcode' (295) o 'cell' (126)
  restano nel calcolo (selezione = 4892 nuclei invece di 4471), e 213 di essi ricevono sr_type 'cardiomiociti' dal fine type.
  L'accordo fra questi 421 è basso (0.31 prereg, 0.49 ph1), da cui la sottostima.
- A5 non e' toccato: nelle ROI di A5 nessun barcode appaiato e' 'Low UMI Barcode' o 'cell' (2415 = 2415), e il mio numero coincide.
- Toccati anche: 'P-6 (all matched)' A6 (loro 0.5735895339329518 / 0.7438675388389206; con esclusione corretta 0.6144039364795347 / 0.7843882800268396) e R6_sr_confusion_A6{,_ph1}.csv (sommano a 4892).
- Fix suggerito: `d.sr_lowumi.fillna(True).astype(bool)` prima di `~`, oppure escludere i barcode in add_sr prima del merge. Gli intervalli bootstrap e i verdetti di P-6 A6 vanno ricalcolati.

## Discrepanza 3: P-2 in A2 (solo ph1), rapporto max/min delle aree mediane fra i tipi >= 5 %
- mio: 1.3876917906811153 (epiteliale normale 33.7855 / fibroblasti 24.3465); loro: 1.6223021582733814 (epiteliale normale / epiteliale tumorale 20.8256).
- Causa (verificata, differenza 0.0): ratio_max_min calcola la soglia del 5 % sulle frazioni rinormalizzate dopo aver tolto
  'non risolto' e 'non classificabile'. In A2 ph1 'epiteliale tumorale' e' il 4.954 % di tutti i nuclei (1621/32721) ma il 5.423 % dei nuclei risolti.
  La pre-registrazione (docs/steps/R6.md riga 81) dice "fra i tipi con almeno il 5 % dei nuclei", che ho letto come frazione su tutti i nuclei
  (la stessa riportata in R6_refs_type*.csv). E' l'unico caso, fra 12 archetipo x variante, in cui le due definizioni selezionano tipi diversi.
- Non e' un errore di codice ma una scelta di definizione da esplicitare nel report; con la definizione pre-registrata il valore e' 1.388, con quella del codice 1.622 (entrambi sopra la soglia 1.161).

## Tutto il resto
Tutti i 313 confronti rimanenti (n, frazioni, aree mediane, eccentricita' mediane per tipo in 12 archetipo x variante; P-1 A5; P-2 negli altri 11 casi; P-3 A6; P-6 A5) coincidono entro 1e-6.
