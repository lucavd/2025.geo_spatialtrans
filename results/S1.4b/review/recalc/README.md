# S1.4b — revisione avversariale: ricalcolo indipendente (revisore `recalc`)

Data 2026-10-15; repo commit `1c13c4c` (branch step1-cell-layer), in sola lettura. Codice scritto da zero in questa cartella
(non e' stato fatto il source di `tools/S1.4b_analysis.R`). Dipendenze in sola lettura: `R/04b4_derive_cell_geometry.R` tramite `dcg_env()`
(`tools/S1.4_mutants.R`), `R/04b1_extract_regions.R` e `tools/R3_voronoi_metrics.R` (`r3_roi_window`), solo per il punto 4.

## Script (dalla root del repo, `R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla <script>`)
| Script | Cosa fa | Output |
|---|---|---|
| `01_recalc_verdicts.R` | punti 1-3: distanze di fold (mediana semi 1..10), PASS/WARN/FAIL contro mediana/max del metro, regola d'archetipo, B-4/B-8/B-9, CP-1, P1-P6, regola di decisione | `recalc_folds.csv`, `recalc_archetype_B.csv`, `recalc_B4_B8_B9.csv`, `recalc_cp1_roi.csv`, `recalc_cp1_archetype.csv`, `recalc_ab_identity.csv` |
| `02_compare_csv.R` | confronto con `S1.4b_verdicts.csv` | stdout |
| `03_compare_folds.R` | confronto a piena precisione con `S1.4b_folds.csv`; fold piu' vicini alle soglie | stdout |
| `04_rebuild_vb.R` | punto 4: V-b ricostruita da S1.4 (A1, A4, A6 fold 1, s = 1) | `recalc_vb_rebuild_A1f1_A4f1_A6f1.csv` |
| `05_b7real_b8.R` | punto 5: b7_real contro PH-1 `circ_gen_1` e B-8 (tutti i 30 ROI) | `recalc_b7real_b8.csv` |
| `06_seed_sensitivity.R` | sensibilita' della regola ai semi (leave-one-seed-out, bootstrap B = 2000) | `recalc_seed_bootstrap_*.csv` |
| `07_cl2_recount.py` (`.venv/bin/python`) | riconteggio di C-L2 e confronto con PH-B | stdout |

## Sintesi
Tutti i numeri chiave sono IDENTICI a S1.4b_verdicts.csv / S1.4b_folds.csv: 60/60 righe di fold (diff 0 a piena precisione), 48/48 righe B archetipo x metrica x variante (stringhe di fold ed esiti). Esiti V-a/b | V-c: B-1 PASS solo A3 | solo A3; B-2 PASS A5, A6 | A5, A6; B-3 nessun PASS | A6; B-7 PASS A5 | nessuno (V-c 0 per costruzione). B-4 V-a/b 0.380 0.557 0.357 0.525 0.185 0.087 (tutti FAIL); V-c 0.164 0.327 0.179 0.478 0.068 0.0030 (PASS solo A6). CP-1: 30/30 ROI, 6/6 archetipi -> (a) non smentita. P1 FAIL, P2 FAIL, P3 FAIL (0/6; B-1 1, B-2 2), P4 FAIL (3/6; V-c FAIL 5/5), P5 PASS, P6 PASS. Conteggio PASS {B-1,B-2,B-3,B-7}: {a,b} 4 vs c 4 -> (c). C-L2 22/30 (68 vs 15.28) ricontato identico. Verifiche da S1.4: V-b ricostruita per A1/A4/A6 fold 1 s=1 identica (ks_an, trunc_all diff 0); b7_real = PH-1 circ_gen_1 e B-8 identici in 30/30 ROI. Nessuna discrepanza numerica; l'unica differenza di definizione e' P1 (rapporto delle medie vs media dei rapporti: A5 0.66 vs 0.75), che non cambia l'esito. Rilievi di metodo: il pareggio 4-4 che decide (c) e' fragile rispetto ai semi (leave-one-seed-out: 5/10 -> {a,b} 4 vs 3; bootstrap P(parita')=0.33); C-L2 usa una soglia gaussiana mal calibrata (PH-B); il placebo di CP-1 e' debole.

## Rilievi

### RC-1 — osservazione — cambia il verdetto: no
Ricalcolo indipendente (codice da zero, dai grezzi loo/*.rds e cp1/*.csv) dei verdetti B-1, B-2, B-3, B-7, B-4, B-8, B-9, CP-1, P1-P6 e regola di decisione: nessuna discrepanza con S1.4b_verdicts.csv e S1.4b_folds.csv.

*Evidenza.* 03_compare_folds.R: 60/60 righe di fold, max |diff| = 0 su ks_nc, ks_an, d_sp, b7, trunc_all, loss_all, b7_sim, b7_real, y_b7_med, y_b7_max; verdetti di fold 60/60. 02_compare_csv.R: 48/48 righe archetipo x metrica x variante con stringa di fold ed esito identici, distanze a 2 cifre 36/36, B-7 sim|reale diff 0, B-4 e perdita diff 0. CP-1: 30/30 ROI battono, 6/6 archetipi, (a) non smentita; P1 F, P2 F, P3 F (0; 1|2), P4 F (3/6; 5/5), P5 P, P6 P; PASS {B-1,B-2,B-3,B-7}: V-a 4, V-b 4, V-c 4 -> (c). File: recalc_folds.csv, recalc_archetype_B.csv, recalc_B4_B8_B9.csv, recalc_cp1_roi.csv, recalc_cp1_archetype.csv.

### RC-2 — difetto di metodo — cambia il verdetto: no
La decisione finale (c) poggia su un pareggio 4 vs 4 risolto dalla regola di parita' aggiunta nelle dev notes (non nella pre-registrazione, 'da confermare con Luca'), e il pareggio non e' robusto alla scelta dei 10 semi Monte Carlo: le celle decisive sono A6 V-c B-2 e B-3 (ROI di test A6 con n_ev 170-688).

*Evidenza.* 06_seed_sensitivity.R: leave-one-seed-out (stesso seme tolto in tutti i fold) -> V-c scende a 3 PASS in 5/10 casi (semi 2,3,5,8,10), V-a/V-b restano 4 => la regola sceglierebbe {a,b}. Bootstrap dei semi per fold (B=2000, seme 20261015): P({a,b}>c)=0.609, P(parita')=0.334, P(c>{a,b})=0.058. Stabilita' dell'esito base: A6 V-c B-3 0.554, A6 V-c B-2 0.572, A6 V-a B-7 0.826. A6 fold 3 V-c ks_an: mediana 0.24539 vs mediana del metro 0.24586 (margine 0.19%); valori per seme 0.2462 0.2077 0.2044 0.2513 0.2446 0.2529 0.2831 0.2245 0.2546 0.2245. Il verdetto con i semi pre-registrati e' applicato correttamente; e' la sua stabilita' a mancare.

### RC-3 — difetto di metodo — cambia il verdetto: no
L'ineleggibilita' di (a) discende da C-L2 FAIL (22/30), ricontato identico, ma la soglia 'celle con |z|>4 <= E+3sqrt(E)' usa l'attesa gaussiana, mal calibrata con pochi punti per cella; con la nulla esatta (PH-B) 7 degli 8 insiemi falliti sono compatibili con il caso. Con i semi ufficiali non cambia la decisione (parita' -> c), ma negli scenari di RC-2 in cui {a,b} vince, la scelta (a) vs (b) dipende proprio da questo.

*Evidenza.* 07_cl2_recount.py: 22/30 PASS, 68 osservati vs E=15.28 (identico al CSV); falliti A1f3, A3f2, A3f4, A4f1, A4f3, A4f4, A4f5, A5f3. S1.4b_posthoc_arbiter_PHB.csv: in A4 (6.4-8.4 punti/cella) media nulla esatta 8.2-12.5 celle vs E gaussiana 1.50-1.97; osservati 8-11 (p 0.59-0.93). Solo A1 fold3 supera q99 della nulla (8 vs 6; p<1/500) ed e' l'unico sotto Bonferroni 0.05/30; A5f3 p=0.044. PH-A (4e6 punti): 8/8 PASS.

### RC-4 — osservazione — cambia il verdetto: no
V-a e V-b sono identiche su tutte le metriche B e su B-4 per costruzione (stessi nuclei; operatore di osservazione = Voronoi dei generatori): il pannello B e il conteggio {a,b} vs (c) non possono discriminare (a) da (b); lo fanno solo CP-1 e i check C.

*Evidenza.* recalc_ab_identity.csv: 30/30 fold, semi 1..10, max |V-a - V-b| = 0 su ks_nc, ks_an, d_sp, trunc_all, trunc_ev, loss_all, b7.

### RC-5 — osservazione — cambia il verdetto: no
B-7 di V-c vale 0 per costruzione (clamp 0.85 dal POI): V-c non puo' mai fare PASS su B-7, quindi nel conteggio a 4 metriche V-c concorre di fatto su 3. Asimmetria pre-registrata (P4), che favorisce {a,b}: nonostante cio' si ha parita'.

*Evidenza.* recalc_folds.csv: b7_sim V-c = 0 in 30/30 fold x 10 semi; verdetti di fold B-7 V-c 29 FAIL + 1 WARN (A3 fold 2).

### RC-6 — difetto di metodo — cambia il verdetto: no
CP-1: il placebo a pesi permutati e' peggiore del Voronoi stesso, quindi 'battere 20/20 permutazioni' e' un criterio debole: e' soddisfatto anche dove Laguerre con pesi reali peggiora rispetto al Voronoi. La riduzione rispetto al Voronoi e' piccola in A1-A4. '(a) non smentita' e' percio' evidenza debole a favore di (a).

*Evidenza.* recalc_cp1_roi.csv: min_perm > Voronoi in 30/30 ROI su m1 e 29/30 su m2; Laguerre < Voronoi in 30/30 (m1) e 29/30 (m2); un ROI di A6 ha riduzione di m2 = -0.025 ma batte comunque 20/20 permutazioni. Riduzione di m2 per ROI (min-max): A1 0.150-0.219, A2 0.126-0.238, A3 0.128-0.362, A4 0.016-0.082, A5 0.440-1.000, A6 -0.025-0.310.

### RC-7 — osservazione — cambia il verdetto: no
P1: la riduzione e' calcolata come rapporto delle medie sui ROI (1 - media m2_lag / media m2_vor), mentre le dev notes dicono '(media sui ROI)', leggibile come media dei rapporti. Con la media dei rapporti cambiano A5 e A6, ma l'esito di P1 resta FAIL.

*Evidenza.* recalc_cp1_archetype.csv: rapporto delle medie A1 0.184, A2 0.172, A3 0.226, A4 0.043, A5 0.656, A6 0.145 (= CSV: 0.18 0.17 0.23 0.04 0.66 0.14); media dei rapporti A1 0.186, A2 0.172, A3 0.226, A4 0.045, A5 0.755, A6 0.167. In A1-A4 entrambe < 0.5.

### RC-8 — osservazione — cambia il verdetto: no
Verifica indipendente di V-b dai dati S1.4 (punto 4): partendo da derive_cell_geometry() (dcg_env(), parametri M2 del fold, seme seed+100000) e con la mia regola di non sovrapposizione (coppie da spatstat closepairs, s_i = min(1, min_j 0.99 d_ij/(R_i+R_j), 0.99 d_reg/R_i), k32 = sqrt(2pi/(32 sin(2pi/32)))), i grezzi di A1, A4 e A6 (fold 1, s = 1) si riproducono esattamente.

*Evidenza.* recalc_vb_rebuild_A1f1_A4f1_A6f1.csv: diff = 0 su ks_an (0.1084004 | 0.1208197 | 0.3000000), trunc_all (0.41251859 | 0.51022905 | 0.08160237), ks_nc, d_sp, trunc_ev, b7, trunc_pair, trunc_reg e C3 V-b (2095 | 6815 | 34); loss_all diff 2.4e-10 (area analitica s^2 A contro area del poligono sulla griglia). V-c ricostruita = S1.4 seeds[1, ks_nc] (diff 0). Coppie in sovrapposizione prima: 2362 | 7845 | 27 (= n_pairs dei grezzi); dopo lo scaling 0 violazioni.

### RC-9 — osservazione — cambia il verdetto: no
Verifica di b7_real e di B-8 (punto 5) col mio codice sui 30 ROI (ne erano richiesti 2): identici.

*Evidenza.* recalc_b7real_b8.csv: b7_real = PH-1 circ_gen_1 in 30/30 ROI (max |diff| 5.0e-16) e = loo$b7_real (diff 0); n_ev = n di PH-1 in 30/30; B-8 'circle' e 'rule' = loo$b8 (diff 0) in 30/30.

### RC-10 — osservazione — cambia il verdetto: no
B-8 in A4/A6 e' dominato da vicini non valutabili (celle Space Ranger non appaiate, raggio da a_own): limitando i vicini all'insieme ev, B-8 circle si dimezza in A4 e quasi si annulla in A6. Metrica descrittiva; a_own (maschera Space Ranger) potrebbe non essere lo stesso oggetto del nucleo StarDist.

*Evidenza.* recalc_b7real_b8.csv: A4 0.335-0.383 (tutti i vicini) vs 0.171-0.205 (vicini ev); A6 0.037-0.059 vs 0.000-0.032; A1 0.325-0.347 vs 0.303-0.330.

### RC-11 — osservazione — cambia il verdetto: no
Fold vicini alla soglia: A4 fold1 V-a/b B-2 e' FAIL per lo 0.1% (non cambia A4, che resta FAIL); A6 fold3 V-c B-2 e' PASS per lo 0.19% (decisivo, vedi RC-2). Il KS su A_N e' identico tra V-c e V-a/b in A5 fold3 (0.19366): il sup del KS cade nella parte non troncata della distribuzione, quindi B-2 e' poco sensibile al troncamento dei nuclei grandi.

*Evidenza.* 03_compare_folds.R: A4 f1 V-a ks_an 0.12106 vs metro max 0.12092; A6 f3 V-c ks_an 0.24539 vs metro mediano 0.24586; A5 f3 ks_an V-c = V-a = 0.19366.

## Note su definizioni scelte dal revisore
- Distanza B-7 di fold = mediana sui semi di |b7_s - b7_reale| (come da compito). Controllo: |mediana(b7_s) - b7_reale| produce verdetti di fold identici (90/90).
- B-9: mediana sui fold delle mediane sui semi di `dat_med` (riproduce il CSV); la media sui fold darebbe A1 0.0513, A2 0.0770, A3 0.0509, A4 0.0491 (P5 resta PASS).
- B-8 (punto 5): vicini cercati su tutte le celle del ROI, raggi da a_n se appaiato, altrimenti a_own; frazione su ev. Riproduce i grezzi in 30/30 ROI.
- Nessun file del repo o dei grezzi e' stato modificato; nessun commit.
