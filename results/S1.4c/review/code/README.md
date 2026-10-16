# S1.4c — revisione avversariale, ruolo CODICE (revisore: `code`)

Data: 2026-10-15. Commit rivisto: `d710cfb` (branch `step1-cell-layer`). Nessun file del repo modificato fuori da questa
cartella; nessuna scrittura nei grezzi ufficiali; scratch in `/mnt/micron/geo_spatialtrans/S1.4c_review/code/`. Nessuna
funzione di `tools/S1.4c_analysis.R` o `tools/S1.4c_decision.R` usata per ricalcolare (decision.R e' solo copiato in scratch
come bersaglio dei mutanti DM1-DM10 di 06). `tools/S1.4c_raster.py` e' importato solo come *oggetto* del confronto (04) o
per applicare lo stesso operatore a insiemi alternativi (05); le misure di confronto sono codice proprio.

## Script (rilanciabili dalla root del repo; Python = `.venv/bin/python`, R = `Rscript --vanilla` senza librerie di sito/utente)
| Script | Cosa fa | Uscite |
|---|---|---|
| `00_inspect.R` | struttura degli oggetti RDS di input | stdout |
| `01_export.R` | esporta in scratch celle reali S1.4, prep, loo, gen, regioni ROI (WKT), mutanti | `scratch/exp/` |
| `02_obs_recompute.py` | ricalcolo indipendente dello stadio obs su 900 file (Voronoi shapely, cerchi cKDTree, KS propria) | `scratch/02_obs_recomputed.csv`, `02_obs_recompute_summary.txt` |
| `02b_voronoi_diag.R` | diagnosi della differenza B-1 (frammenti dei territori multiparte) | scratch |
| `03_decision_recompute.py` | metri, verdetti, punteggi, CP-0/1/2, bootstrap (RNG proprio), phi, a_ref, B-4, robustezza | `03_*.csv`, `03_extra.json` (suffisso `_obs_rv` = stessa analisi sulle distanze ricalcolate da 02) |
| `04_raster_checks.py` | raster indipendente, contatto indipendente, fessura 2e-5 um, multiparte, solidita' su dischi | `04_*.csv` |
| `05_asym_A4A6.py` | effetto delle asimmetrie A4/A6 (SR contro StarDist) su B-7, B-8, S-1..S-3 e sul punteggio | `05_*.csv` |
| `06_mutant_power.py` | motivo di rilevamento dei mutanti, MU5 nel banco, 10 mutanti dello strato di decisione | `06_mutant_power.json` |
| `07_kclip_scan.py` | proprieta' geometriche di K-1/K-2 sui 300 insiemi (multiparte, buchi, K-1 in K-2, tetto) | `07_*.csv` |
| `08_at_vs_geometry.R` | a_t di `territory_df` contro area della geometria del territorio (30 ROI) | `08_at_vs_geometry.csv` |

## Esito in breve
Nessun bug che cambi verdetti o decisione. Tutti i numeri chiave ricalcolati con codice proprio coincidono (tabella
`03_recomputed.csv`). La scelta K-2 regge a margine 1-5, a p_min 0.8-0.99, a 10/10 leave-one-seed-out, al pannello completo a
7 metriche e a tutte le varianti A4/A6 provate; cade su K-c con margine 6 e su K-1 (parita' 7-7, regola della piu' semplice)
togliendo A2. K-2 batte K-1 di 1 solo punto (A2 x B-2, 4/5 contro 3/5 fold PASS). Con una lettura di «C-perf solo WARN» che
blocchi l'eleggibilita' la scelta sarebbe K-1 (D-6): la correzione e' coerente con il testo ma e' successiva ai risultati.

## Rilievi (id, gravita', cambia verdetto/decisione?)
- **CODE-1 conferma (no).** Stadio obs ricalcolato in modo indipendente su 900/900 file: n_eval, n_excluded, b7, b8, B-2, B-7,
  B-8, S-1, S-2, S-3 identici (max |diff| <= 1e-15). B-1 differisce in 57/900 file di al piu' una cella (1/n_ev, max 4.5e-3):
  `tessellate_voronoi()` assegna i frammenti staccati dei territori multiparte (regioni con buchi, A5 r4 ha 133 buchi) in modo
  diverso da un Voronoi ∩ regione semplice; i territori del pacchetto sono coerenti con se stessi (`08`: a_t = area della
  geometria, max rel 5e-14) e reali/simulati usano lo stesso operatore (C-OP). Con le distanze ricalcolate da me: verdetti
  126/126 identici, stessi punteggi, stessa scelta (`03_*_obs_rv`).
- **CODE-2 conferma (no).** Metri e fold: ROI di test = k-esimo ROI per i 30 fold, `yard$j` = altri ROI nello stesso ordine
  (60/60); `ks_nc`, `ks_an` di S1.4 ricalcolati dai reali (max 5e-16); phi = min(1, q99 tipo 7) solo sui 4 ROI di
  addestramento (max 4.9e-15); a_ref e n del fit M2 solo dai ROI di addestramento (30/30). Controprova: phi con il ROI di test
  incluso differirebbe (es. A6 fold 4 0.138 -> 0.180), quindi il valore usato non contiene il test.
- **CODE-3 conferma (no).** K-1/K-2 su 300 insiemi (2.4 M nuclei per variante): 0 vuoti, 0 buchi, 0 violazioni del tetto
  (area(K-2) <= phi a_t), 15 nuclei multiparte in totale (A5 10, A6 5 su K-1+K-2); dove il tetto non e' attivo K-1 ⊆ K-2 a
  meno del rumore dello snap (0 violazioni con tolleranza 3e-6 um; area di K-1 \ K-2 <= 3.5e-6 um²). Dove il tetto e' attivo
  K-2 e' piu' piccolo di K-1 per costruzione (K-1 non ha tetto).
- **CODE-4 conferma (no).** Fessura di 2e-5 um (rientro 1e-5): chiudendola (pixel non etichettati entro 5e-5 um da un nucleo
  assegnati al nucleo) su 24 insiemi K-1/K-2: |ΔS-1| <= 2.3e-4, Δmediana S-2 = 0, |Δmediana S-3| <= 7e-5, KS(S-2) <= 9e-4,
  KS(S-3) <= 7e-3 (A5; A4 <= 1.4e-3), Δarea media <= 3.5e-3 um². Trascurabile rispetto ai metri. Nota: in A4 ~490-660 centri di
  pixel per ROI cadono esattamente fra due territori SR (lati di Voronoi su multipli del passo), senza rientro sarebbero
  conflitti; con il rientro restano vuoti.
- **CODE-5 conferma (no).** Raster: rasterizzazione indipendente (shapely.contains_xy sui centri) = skimage.draw.polygon,
  0 pixel diversi su 30 insiemi; bordo e contatto δ=1,2 con codice proprio = sim_out (30 insiemi) e = real_shape (30 ROI),
  0 discrepanze; coppie di contatto di CP-0 coerenti con n_c2 (30/30); conflitti 0.
- **CODE-6 rilievo (no, S-3 gia' declassata).** Solidita' di skimage su oggetti piccoli: un disco perfetto ha solidita'
  0.958 ± 0.028 (60 px), 0.961 ± 0.011 (190 px), 0.975 ± 0.002 (800 px), 17-58 valori distinti. La mediana reale A4 (0.958, ~213
  px) coincide con quella di un disco: a questa taglia S-3 non distingue forme e dipende dalla taglia (confusa con B-2).
  Etichette reali sconnesse: cellpose 1-58 per ROI, StarDist 0-5; simulati 0-3 (coda bassa di S-3 asimmetrica).
- **CODE-7 rilievo (no).** Asimmetria A4/A6 quantificata (`05`). n_SR / n_StarDist = 1.30-1.40 (A4), 1.15-1.35 (A6);
  24-29 % dei territori SR in A4 (19-40 % in A6) senza nucleo StarDist. I riferimenti reali dipendono molto dalla definizione: A4 f1 B-8 reale
  0.335 (codice: a_n o a_own) -> 0.273 (a_own ricalibrata) -> 0.172 (sole celle appaiate) -> 0.148 (sola maschera StarDist);
  B-7 0.397 -> 0.210 (appaiate) -> 0.188 (StarDist). Punteggi (3, 7, 8) invariati in tutte le 8 combinazioni B-7 x B-8; il WARN
  di B-7 per K-1/K-2 in A4 diventa FAIL con entrambe le alternative a densita' coerente (non cambia il punteggio). S-1 in A4: la
  distanza simulato-reale di K-2 scende da 0.32 a 0.14 restringendo entrambi alle celle appaiate (meta' del divario e' l'asimmetria
  di densita'); in A6 S-1/S-2 di K-2 passano da FAIL a WARN. La pre-registrazione prescrive i territori R3 primari (SR) per i
  reali: il codice la rispetta; e' un limite dichiarato (dev notes), ora quantificato.
- **CODE-8 discrepanza (no).** CP-0: la pre-registrazione indica lo score r_i + r_j - d_ij, il codice usa (r_i + r_j)/d_ij (dev
  notes). Non sono trasformazioni monotone l'una dell'altra; con lo score pre-registrato AUC A1 0.959, A2 0.942, A3 0.965,
  A4 0.975 (codice 0.960/0.947/0.967/0.974): premessa non smentita in entrambi i casi.
- **CODE-9 rilievo (decisione: si', sotto lettura alternativa).** D-6: il testo «tutti i check C PASS (C-perf solo WARN)» e la
  colonna soglia di C-perf («solo WARN») si leggono naturalmente come «C-perf non blocca»: la correzione e' coerente. Ma e' stata
  fatta dopo aver visto i risultati e capovolge la scelta: con C-perf bloccante K-2 non e' eleggibile e la scelta e' K-1 (7 >= 3+3,
  P = 1.0). C-perf misurato 622.7 s (> 300 s) e' un WARN reale per S1.4d.
- **CODE-10 rilievo (no).** K-2 contro K-1: 8 contro 7, differenza dovuta a una sola cella (A2 x B-2, 4/5 contro 3/5 fold PASS).
  Bootstrap appaiato: P(K-2 > K-1) = 0.928 (RNG proprio, B = 2000) / 0.935 (file riportato, B = 200), parita' 7 %, mai K-1 > K-2.
  La regola 4 non chiede un criterio bootstrap fra K-1 e K-2. Togliendo A2 la regola sceglierebbe K-1 (parita' 7-7).
- **CODE-11 rilievo (no).** CP-1 come pre-registrata (7 metriche): guadagno A1-A4 6, A5+A6 3 -> smentita (3 >= 6/2), ricalcolato.
  Sul pannello dopo CP-2 (4 metriche) 4 contro 1: non smentita. Il guadagno in A5/A6 viene dalle metriche S declassate. La
  pre-registrazione non lega una conseguenza a CP-1: va riportata come riserva, non cambia la regola.
- **CODE-12 conferma (no).** CP-2: 18 celle metrica x archetipo ricalcolate uguali; S-1, S-2, S-3 declassate anche con due
  regole alternative (distanza fra segmentatori contro il massimo del metro del banco; mediana dei ROI contro mediana dei massimi).
- **CODE-13 rilievo (no).** Potere dei test di decisione: 7/10 mutanti costruiti rilevati da `test_S1.4c_decision.R`; non
  rilevati: DM1 bootstrap senza reinserimento, DM2 bootstrap non appaiato, DM3 verdetto sulla media del metro (D1 usa un metro
  con media = mediana; D6 usa distanze costanti fra semi). Sui dati: DM3 da' 2/6/8 (scelta K-2); bootstrap non appaiato
  P(K > K-c) = 1.0, P(K-2 > K-1) = 0.864. Nessun effetto sulla decisione, ma bootstrap e metro non sono coperti dai test.
- **CODE-14 conferma (no).** Mutanti C-mut rilevati per il motivo giusto: MU1 da U3/U4c (tetto), MU2 da U4/U4b (C3), MU5 da U5
  (Voronoi dei centroidi), MU3 da T1 (centroide 0.56 px), MU4 da T2/T3/T5 (contatto). Nel banco: MU1 C-area 0.06-0.85, MU2 C3
  su tutti i nuclei K-2 (in K-1 per lo piu' nuclei vuoti, C-empty), MU5 sposta B-7 simulato di K-c da 0.004 a 0.31 (A1).
  Nota: MU3 dipende da un'unica asserzione (centroide <= 0.25 px), MU5 da una sola (U5).
- **CODE-15 conferma (no).** Verdetti (126/126), punteggi 3/7/8, P(bootstrap) 1.0/1.0, scelta K-2, CP-0 AUC (max diff 4e-16, stesse
  coppie), B-4 K-2 6/6, C-0 330/330, C3/C-NO/vuoti 0, C-area 5.84e-7: ricalcolati e coincidenti.

## Robustezza (descrittiva, `03_robustness.csv`)
margine 1-5: K-2; margine 6: K-c. p_min 0.8/0.95/0.99: K-2. Senza un seme (10 casi): K-2 sempre (K-c 2-3, K-1 7, K-2 8-9).
Senza un archetipo: A1 K-2, A2 **K-1** (7-7), A3 K-2 (2/4/5), A4 K-2, A5 K-2, A6 K-2. Pannello a 7 metriche: 3/10/12 K-2.
PASS stretto (d < mediana): invariato. Bootstrap: quantili 2.5-97.5 % K-c [2, 3], K-1 [6, 8], K-2 [7, 10].

## Non fatto / ridotto
- `frac_out` (D-1, descrittiva) non reimplementata in modo indipendente (resta il test analitico T6 del banco).
- Bootstrap: non replicata la sequenza RNG di R; usato un bootstrap appaiato proprio (B = 2000) e letto il file riportato (B = 200).
- Prova della fessura e del raster indipendente su 30 insiemi (fold 1 s = 1 per archetipo x 3 varianti, fold 3 s = 7 x K-1/K-2),
  non sui 900.
- C-area non ricalcolata dalla geometria (area target di K-c non esportata): verificato il tetto e K-1 ⊆ K-2.
