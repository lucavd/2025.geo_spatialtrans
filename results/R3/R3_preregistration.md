# R3 — Pre-registrazione: Voronoi sui nuclei reali

**Data**: 2026-10-08 · **Branch**: `step1-cell-layer` · **Base**: `acc301f`
**Approvazione**: Luca, in chat, 2026-10-08 («approvo»), prima di qualsiasi codice e di qualsiasi calcolo sui territori reali. L'approvazione copre il testo qui sotto e le tre decisioni proposte in apertura (D-R3.1…D-R3.3) così come formulate.
**Obiettivo unico**: tassellazione di Voronoi dei centroidi nucleari reali dei 30 ROI R2 (5 per archetipo A1–A6), ritagliata sulla maschera valida (tessuto − bolle) → distribuzioni empiriche per archetipo di area / eq_radius del territorio, rapporto N/C, eccentricità e orientazione del territorio, contenimento dei nuclei; quanto il Voronoi «puro» spiega la geometria reale.
**Fuori scope**: codice del pacchetto (`tessellate_voronoi()` è S1.3), preset, diagramma di potenza/Laguerre (eventualmente BACKLOG), regioni graphclust (BL-049 non toccato: il ritaglio è la maschera H&E), area nucleare manuale (BL-039).

## Decisioni prese in apertura (2026-10-08, Luca «approvo»)
- **D-R3.1** Segmentatore: posizioni come decisione R2b — Cellpose-SAM RGB in A1, A2, A3, A5; Space Ranger 4 in A4, A6. N/C con lo stesso segmentatore delle posizioni; **sensibilità** in A4 e A6 con StarDist `2D_versatile_he` nativo (posizioni e aree StarDist), perché i poligoni SR4 sono quantizzati a 4 µm² (BL-027).
- **D-R3.2** BL-035: riferimento primario = segmentatore sui ROI 1 mm² (n grande, distorto dai nuclei persi); calibrazione nelle 24 finestre R2b (Voronoi sui punti di Luca vs Voronoi sul segmentatore primario nella stessa finestra, stesso ritaglio) → riferimento come intervallo [corretto, segmentatore].
- **D-R3.3** Roadmap §4: **R3 prima di S1.3** (la tabella §4 lo richiede per i check B di S1.3; il testo «R3 prima di S1.4» è superato).

## Definizioni operative (fissate qui; gli script che le implementano sono committati prima del run sui dati reali, BL-059)
- **Ingressi**: `results/R2/R2_nuclei_all.parquet`, righe `scale == 1`, `keep == TRUE`, `method` = segmentatore di D-R3.1; maschera `results/R2/tissue_masks/<A>_<roi>_valid_ds4.png` (tessuto ∧ ¬bolla, 913×913, pixel = side_um/913). Coordinate in µm del ROI, y verso il basso; centroide = centroide regionprops × µm/px (centro del pixel 0 = 0).
- **Finestra / ritaglio**: poligono della maschera da `extract_regions(maschera, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1)` (contorno sui lati dei pixel, area esatta), come in S1.2. Un nucleo entra se il suo centroide cade in un pixel valido.
- **Territorio**: cella `deldir` (riquadro = ROI allargato di 50 µm) del centroide, intersecata con il poligono della maschera (`sf::st_intersection`). Ogni territorio è associato al proprio generatore per indice (deldir scarta i duplicati: controllato in C-R3.2).
- **Cellula interna**: territorio non toccato dal ritaglio, cioè area ritagliata / area della cella non ritagliata > 0.999, e cella non ritagliata interamente dentro il ROI. Tutte le metriche di forma e di area sono calcolate solo sulle interne; la frazione di interne è riportata.
- **Metriche per cellula (interne)**: area A (µm²); eq_r = √(A/π); area normalizzata a = A / media(A) del ROI; eccentricità e_T e orientazione θ_T dai momenti del secondo ordine del poligono (ecc = √(1 − λ₂/λ₁), stessa definizione di regionprops usata in R2); numero di lati (vicini di Voronoi); N/C = area nucleo / A; rapporto dei raggi = √(N/C) (= `nucleus_to_eq_ratio` del design).
- **CV locale**: area normalizzata per l'intensità locale, a_loc = A · λ̂(x), λ̂ = intensità di kernel gaussiano dei centroidi con banda σ = 5/√λ_ROI (≈ 5 spaziature medie) e correzione di bordo (spatstat `density.ppp`, `edge = TRUE`), valutata al generatore; CV_loc = sd(a_loc)/media(a_loc).
- **Δθ**: angolo fra asse maggiore del territorio e asse maggiore del proprio nucleo, in [0°, 90°], per cellule interne con ecc nucleo ≥ 0.8 ed e_T ≥ 0.5. Nullo uniforme: mediana 45°. Permutazione: orientazioni nucleari permutate fra le cellule eleggibili del ROI (2 000 permutazioni), p = frazione di mediane permutate ≤ osservata.
- **Contenimento** (CP-3): per ogni nucleo, frazione dei suoi pixel nativi (maschera di label R2) il cui generatore più vicino (KD-tree su tutti i generatori del ROI) non è il proprio; «tagliato» se > 5 %. Solo nuclei con territorio interno.
- **Nulli** (CP-1, CP-2): CSR (punti uniformi nei pixel validi della maschera, n = n osservato) e RSA con d* (`seed_centroids()` di S1.2, un tipo, `min_dist_um` = d* B-049…B-054: 2.5, 3.5, 3.5, 4, 7, 2.5 µm per A1…A6; densità scelta perché n_target = n osservato, a meno dell'arrotondamento stocastico); 20 repliche per (ROI, modello), stesso ritaglio e stesse metriche. Seed = 20261008 + 1000·indice_ROI (1…30, ordine di `R2_rois_checked.csv`) + 100·modello (1 CSR, 2 RSA) + replica.
- **Unità**: il ROI. Valori di archetipo = distribuzione in pool delle cellule interne dei 5 ROI, con min–max delle mediane per ROI. Decisioni per ROI («≥ 4/5 ROI»), esito «fuori dal nullo» = valore reale fuori dall'intervallo [min, max] delle 20 repliche.
- **Calibrazione (CP-4)**: finestre R2b (`results/R2b/R2b_windows.csv`), punti di Luca (`R2b_points.csv`, rater Luca), stessi centroidi del segmentatore primario dentro la finestra; finestra = quadrato della finestra ∩ maschera valida del ROI; stesse definizioni di cellula interna e metriche.

## Check C (correttezza software) — per ogni check l'errore plausibile che lo farebbe fallire (BL-048)
| id | Asserzione | Soglia | Errore che lo farebbe fallire |
|---|---|---|---|
| C-R3.1 | Σ aree dei territori = area del poligono maschera = n pixel validi × px², per ROI (30/30) | rel. 1e-9 (poligono), 1e-6 (pixel) | ritaglio errato, celle perse o sovrapposte |
| C-R3.2 | n territori = n nuclei; 100 % dei generatori dentro il proprio territorio | esatto | disallineamento tile ↔ punti (duplicati scartati da deldir) |
| C-R3.3a | reticolo quadrato passo a: interne con area a², e_T = 0, 4 lati; esagonale: √3/2·a², 6 lati | 1e-9 rel.; e_T < 1e-6 | formula dell'area / dei lati |
| C-R3.3b | ellisse poligonale (1 024 vertici) con semiassi e angolo noti: e_T e θ_T | \|Δe\| < 1e-3, \|Δθ\| < 0.5° | formula dei momenti, convenzione dell'angolo |
| C-R3.3c | ellisse rasterizzata: e_T e θ_T dai momenti del poligono vs regionprops (`eccentricity`, `orientation` convertita nelle coordinate del ROI) | \|Δe\| < 0.02, \|Δθ\| < 2° | convenzione diversa da R2 (segno o assi dell'orientazione) |
| C-R3.4 | CSR, 20 000 punti in un quadrato, interne: media area normalizzata 1, varianza ≈ 2/7 (gamma a = b = 7/2, Ferenc & Néda 2007, doi:10.1016/j.physa.2007.07.063), n lati medio 6 (Euler) | 1 ± 0.01; 0.286 ± 0.015; 6.00 ± 0.02 | bias del meno-campionamento, aree errate |
| C-R3.5 | contenimento raster: due dischi tangenti di raggio diverso con generatori ai centroidi → frazione del disco grande oltre il bisettore analitica; dischi separati → 0 | ±0.01 assoluto; 0 esatto | assi x/y invertiti, etichette sfasate, convenzione del centro del pixel |
| C-R3.6 | riproducibilità: secondo run → md5 identico delle tabelle per cellula (reale e nulli) | identico | seed non per unità, ordine non deterministico |
| C-R3.7 | mutanti M1 tile permutati, M2 niente ritaglio, M3 x↔y nel contenimento, M4 flag interne ignorato, M5 varianza al posto della deviazione standard in e_T | 5/5 rilevati da C-R3.1–C-R3.5 | — |
| C-R3.8 | ingressi: 30/30 ROI; n nuclei usati = conteggio `keep` con centroide in pixel valido; md5 degli ingressi registrati | esatto | filtro di ingresso errato |

## Check B (plausibilità, per archetipo; fonti: righe BIO_REFERENCES esistenti)
Valori derivati prima del run (eq_r della media = √(10⁶/(πρ)), densità manuale B-043…B-048 [totale, evidenti]):
| | A1 | A2 | A3 | A4 | A5 | A6 |
|---|---|---|---|---|---|---|
| area media (µm²) [tot, evid] | 80.5–97.2 | 112.1–135.6 | 322.6–432.7 | 35.6–37.3 | 843.9–1049.3 | 1012.1–1103.8 |
| eq_r della media (µm) | 5.06–5.56 | 5.97–6.57 | 10.13–11.74 | 3.37–3.44 | 16.39–18.28 | 17.95–18.74 |
| √(A_nuc·ρ_seg) da R2 (B-002…B-037; A4, A6 StarDist) | 0.37 | 0.50 | 0.37 | 0.59 | 0.31 | 0.15 |
| catalogo design §4.4 (`nucleus_to_eq_ratio`) | 0.45 epithelial_dense | 0.45 (epiteliale) | 0.45 stromal_loose | 0.55 immune_T | — | — |

| id | Asserzione | Previsione |
|---|---|---|
| B-R3.1 | eq_radius mediano (interne, pool) ∈ intervallo C8 provvisorio [5, 25] µm | PASS A1, A2, A3, A5, A6; **FAIL A4** (≈ 3.3–3.6 µm) → intervalli C8 per archetipo (BL-006) |
| B-R3.2 | mediana del rapporto dei raggi entro ±20 % del catalogo del design | A1 FAIL (≈ 0.37), A2 PASS, A3 FAIL (≈ 0.37), A4 PASS; A5, A6 registrati senza verdetto |
| B-R3.3a | roadmap §3, A4 «N/C ≈ 0.8–0.9» (N/C mediano ∈ [0.8, 0.9]) | **FAIL** (≈ 0.35–0.55) |
| B-R3.3b | roadmap §3, A3 «N/C basso»: N/C mediano A3 < A1 e < A2 | PASS |
| B-R3.3c | roadmap §3, A5 «aree bimodali»: test dip sulle aree dei territori interni, p < 0.05 in ≥ 4/5 ROI | **FAIL** (p > 0.05 in ≥ 4/5 ROI: il Voronoi non conosce la taglia della cellula) |
| B-R3.4 | design §4.4 «nucleo ∝ territorio»: Spearman ρ(area nucleo, area territorio) ≥ 0.3 (pool, interne) | sostenuta in ≤ 2/6 archetipi (se confermato: S1.4 con raggio nucleare assoluto per tipo) |

## Controprove
| id | Test | Previsione |
|---|---|---|
| CP-1a | CV globale e CV locale delle aree interne: reale vs CSR (20 repliche per ROI) | CV_loc reale < CSR in A4, A5 (≥ 4/5 ROI fuori dal nullo, sotto); CV_loc reale > CSR in A6 (≥ 4/5 ROI, sopra); A1–A3 nessuna previsione direzionale |
| CP-1b | RSA(d*) riproduce il CV_loc reale (\|CV_loc,RSA mediano − CV_loc,reale\| / CV_loc,reale ≤ 0.10) | sì in A4, A5, A6 (≥ 4/5 ROI); no in A1, A2, A3 (≥ 3/5 ROI), coerente con D(d*)/D_CSR di S1.2 |
| CP-2a | A6: mediana e_T reale − mediana e_T RSA(d*) ≥ 0.05 | sì in ≥ 4/5 ROI |
| CP-2b | A6: mediana Δθ ≤ 35° e permutazione p < 0.01 | sì in ≥ 4/5 ROI; A1 (palizzata epiteliale) idem, **bassa confidenza** |
| CP-2c | A4: \|mediana e_T reale − mediana e_T RSA(d*)\| < 0.05 | sì in ≥ 4/5 ROI |
| CP-3 | frazione di nuclei tagliati (> 5 % dei pixel fuori dal proprio territorio), interne | A1 ≥ 10 %; A5 ≤ 5 %; A2, A3, A4, A6 registrati. Limite inferiore: i nuclei persi allargano i territori |
| CP-4 | finestre R2b: CV globale delle aree interne, segmentatore vs punti di Luca (pool delle 4 finestre per archetipo) | CV_seg > CV_manuale in ≥ 5/6 archetipi. Il rapporto delle mediane di eq_r (seg/manuale) è stimato (fattore di correzione), senza verdetto |

## Descrittivi (senza verdetto)
- D-1 oggetti `cellpose_hed` (nucleo + alone) vs territorio di Voronoi, appaiati per nucleo (il generatore cade nell'oggetto hed): rapporto delle aree (BL-028).
- D-2 distanza al primo vicino (con correzione di bordo per meno-campionamento) e Clark–Evans nella forma standard (spatstat `clarkevans`, correzione `cdf`) (BL-053/054).
- D-3 distribuzione del numero di lati, reale vs nulli.
- D-4 sensibilità: A3 con e senza r3 (BL-036); A4 follicolo (f1, f2) vs paracorticale (p1–p3); A4 e A6 con StarDist (D-R3.1).
- Limite dichiarato: bolle sfocate in A6 non mascherate (BL-040); le aree dei territori adiacenti sono gonfiate.

## Regole
Le attese e le soglie qui sopra non si modificano dopo il risultato. Ogni divergenza fra definizioni operative implementate e questo testo è una deviazione e va dichiarata come tale nel report (BL-059).
