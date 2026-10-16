# R2c — Pre-registrazione: verità manuale sul contatto fra nuclei (BL-120)

**Approvata da Luca il 2026-10-16** («si»), prima di qualsiasi annotazione, calcolo di candidati sui segmentatori o risultato.
Scelte di disegno di Luca (2026-10-16): compito «sì/no su coppie proposte»; copertura «nuclei focali + Federica a metà».

## Obiettivo unico
Verità di riferimento manuale per il contatto fra nuclei sulle 24 finestre di R2b, per decidere per archetipo quale segmentatore
(cellpose_rgb, stardist_he) dà una frazione di nuclei a contatto (S-1) credibile, e dove cade K-2 (S1.4c). Dipendenze: S1.4c chiusa
con riserve, finestre e punti di R2b presenti — soddisfatte.

## Disegno
- **Punti**: punti manuali di Luca R2b (`results/R2b/R2b_points.csv`, rater = Luca), in µm del ROI.
- **Candidati**: lati della triangolazione di Delaunay dei punti di Luca per finestra con distanza d ≤ d_max(archetipo).
  Regola per d_max (fissata ora): q99.5 della distanza fra i due punti manuali appaiati alle due etichette di ogni coppia a contatto
  (Chebyshev ≤ δ = 2 px, operatore di S1.4c) in cellpose_rgb ∪ stardist_he, sulle coppie con entrambi i punti appaiati, + 1 µm.
- **Nuclei focali**: 15 per finestra (60 per archetipo), campionati con seme fisso 20261016 fra i punti di Luca a distanza ≥ d_max
  dal bordo della finestra e dalle zone escluse. Si giudicano tutte le coppie candidate che contengono almeno un nucleo focale.
- **Compito**: ritaglio ingrandito centrato sulla coppia, due marcatori piccoli, nessun contorno; risposte «si toccano» /
  «separati» / «non giudicabile». Ordine casuale (seme), alla cieca rispetto ai segmentatori. ~10 % delle coppie ripresentate
  (intra-rater); coppie-trappola con d > 20 µm (fra punti della stessa finestra), ~5 % del totale.
- **Annotatori**: Luca tutte le coppie; Federica le coppie delle sue 2 finestre per archetipo di R2b (flag `rater2`), HTML autonomo.
- **S-1 manuale** per nucleo focale = almeno una coppia «si toccano». Primaria: «non giudicabile» escluse (nucleo escluso se
  tutte le sue coppie sono non giudicabili e nessuna è «si toccano»); intervallo con le varianti «non giudicabile = contatto» /
  «= separati». IC 95 %: bootstrap a grappolo per finestra (B = 2000, finestre ricampionate con reinserimento entro archetipo).
- **S-1 dei segmentatori** sugli stessi nuclei focali appaiati in quel segmentatore, sulle maschere native del ROI con il codice
  raster di S1.4c (contatto con qualunque etichetta, δ = 2 px). Space Ranger 4 come terzo confronto descrittivo in A4/A6 (BL-121).

## Check C — correttezza software
| id | attesa | soglia |
|---|---|---|
| C-1 copertura | frazione delle coppie a contatto di ciascun segmentatore (entrambi i punti appaiati) che cade fra i candidati | PASS ≥ 0.99, WARN 0.97–0.99, FAIL < 0.97, per archetipo × segmentatore |
| C-2 annotatore | ritagli pixel-identici alla finestra R2b; ordine riproducibile dal seme; HTML senza informazioni dei segmentatori (grep); test del core JS | tutti PASS |
| C-3 operatore | S-1 dei segmentatori calcolato dal codice R2c = valori per nucleo di S1.4c (`cp2`) sugli stessi nuclei | identico |
| C-4 mutanti | bootstrap senza reinserimento; bootstrap non a grappolo; κ con marginali scambiate | 3/3 rilevati |
| C-5 eleggibilità (BL-125a) | regole di candidati/focali eseguite e ispezionate su una finestra per archetipo prima di generare l'annotatore | eseguito |

## Check B — previsioni (di Claude, non misure)
- **P-1**: S-1 manuale ≥ 0.30 in A1–A4 e ≤ 0.15 in A5, A6 (soglie di P1 di S1.4c).
- **P-2**: S-1 manuale strettamente fra StarDist e cellpose in ≥ 3 di A1–A4 (cellpose più largo, StarDist più convesso; S1.4c).
- **P-3**: sulle coppie con entrambi i punti appaiati in entrambi i segmentatori: sensibilità cellpose > StarDist, specificità StarDist > cellpose.
- **B-dec** (per archetipo): un segmentatore è «credibile» se il suo S-1 cade nell'IC 95 % manuale; se entrambi, quello con |Δ|
  minore; se nessuno, riferimento = intervallo manuale. Esiti: cellpose / StarDist / nessuno / indeciso (entrambi credibili e
  |Δ| che differiscono < 0.02).
- **B-K2**: S-1 di K-2 dai raster di S1.4c (10 semi, mediana), nuclei simulati con centroide nelle finestre (stessa regola di
  bordo dei focali). PASS dentro l'IC manuale; WARN dentro l'IC allargato di |S-1 Luca − S-1 Federica| dell'archetipo; FAIL
  altrimenti. Limite dichiarato: insieme di nuclei diverso (generatori dai segmentatori).

## Controprove
| id | test | soglia |
|---|---|---|
| CP-1 | κ di Cohen Luca–Federica sulle coppie comuni (non giudicabile esclusa) | ≥ 0.60 PASS; 0.40–0.60 WARN; < 0.40 FAIL → nessuna verità, nessuna decisione. Convenzione Landis & Koch 1977 (doi:10.2307/2529310) |
| CP-2 | κ intra-annotatore sulle ripetute | ≥ 0.70 |
| CP-3 | coppie-trappola giudicate «separati» | ≥ 95 % |
| CP-4 | AUC di −d per «si toccano» in A1–A4 | ≥ 0.80 per archetipo |
| CP-5 | B-dec con le sole risposte di Federica sul suo sottoinsieme | stesso esito in ≥ 3/4 di A1–A4 |

## Potenza (dichiarata prima)
60 focali, S-1 ≈ 0.3: semi-ampiezza IC ~0.12–0.14. Divario cellpose–StarDist ~0.18 (A1), ~0.28 (A2): risolvibile; A3 ~0.13:
rischio «indeciso». S-2/S-3 restano declassate (non misurate da un giudizio sì/no).

## Perimetro
Solo `tools/R2c_*`, `results/R2c/`, `R/testing/test_R2c.sh`; pacchetto invariato. Ogni ridefinizione rispetto a questo testo
va dichiarata come deviazione (BL-124).


## Addendum 2026-10-16 — deviazioni dichiarate prima di qualsiasi annotazione (approvate da Luca in chat: «ok»)

Primo run di `tools/R2c_candidates.py` con la regola pre-registrata (output conservati in `results/R2c/prereg_rule/`):
- C-3 PASS: 0 differenze fra le etichette a contatto di R2c e lo stadio cp2 di S1.4c su 40/40 ROI × segmentatore.
- C-1 (copertura dei candidati di Delaunay): A1 cellpose 0.970 **FAIL** (3/99 coppie non di Delaunay); A2 StarDist 0.967 **FAIL**
  (1/30); A2 cellpose 0.972 WARN; A3 cellpose 0.971 WARN; A1 StarDist, A3 StarDist, A4, A5, A6 PASS.
- C-5 (eleggibilità): in A6 d_max = 5.09 µm stimato da **una sola** coppia a contatto → 0 coppie giudicabili, S-1 manuale = 0 per
  costruzione; A5 da 14 coppie.

Deviazioni (le soglie di tutti i check restano quelle pre-registrate):
- **D-1** candidati = tutte le coppie di punti di Luca della finestra con d ≤ d_max (tolta la restrizione ai lati di Delaunay).
- **D-2** d_max(A) = max(regola pre-registrata, q99.5 delle distanze fra centroidi delle coppie a contatto (δ = 2 px) di
  cellpose_rgb ∪ stardist_he su tutti i ROI dell'archetipo che contengono finestre + 1 µm). La stima su più coppie può solo
  allargare d_max.
- **D-3** verifica della fascia (proposta di Luca: le coppie lontane si classificano «separate» per regola). d_geo(A) =
  max sui due segmentatori di [2 · q99.5(r_max) + 2 · upp], r_max = distanza massima dal centroide del nucleo ai suoi pixel
  (+ 0.5 px), sugli stessi ROI. Coppie con d_max < d ≤ d_geo e almeno un estremo focale: campione di 15 per archetipo (tutte se
  meno), seme fisso, mischiate alla cieca. Le coppie oltre d_max non campionate sono «separate» per regola.
  Regola (risposte di Luca, «non giudicabile» escluse, totale della fascia): limite superiore IC 95 % Clopper–Pearson della
  frazione «si toccano» ≤ 5 % → PASS; esattamente 1 contatto e limite > 5 % → WARN, S-1 manuale dell'archetipo riportata anche
  con le coppie di fascia giudicate; ≥ 2 contatti → FAIL, la classificazione per regola è smentita, limite dichiarato e BACKLOG.
- **D-4** trappole: max(20 µm, d_geo) < d ≤ lato del ritaglio − 6 µm (in A6 i candidati arrivano a 27 µm: una trappola a
  20–27 µm non sarebbe certamente separata).
- Dettaglio d'implementazione (non pre-registrato): lato del ritaglio = max(2.5 · d_max, 30 µm), fisso per archetipo, per tutti
  gli item (coppie, fascia, ripetute, trappole), così l'ingrandimento non rivela il tipo di item.

Carico atteso: ~1 309 coppie + ≤ 90 di fascia + 10 % ripetute + 5 % trappole.
