# R2b — Pre-registrazione (scritta il 2026-09-21, approvata da Luca il 2026-09-21: «procedi»)

**Obiettivo unico**: verità manuale (punti + contorni) su finestre campionate dai 30 ROI R2, per archetipo A1–A6, per (a) fissare valore e intervallo delle righe di densità e area nucleare di `docs/BIO_REFERENCES.md` (chiudere il "provvisorio", BL-024/BL-025/BL-030) e (b) misurare recall/precisione di Cellpose-SAM RGB, StarDist HE e Space Ranger 4 per archetipo e scegliere il segmentatore per R3 (BL-032, BL-033). Nessun codice del cell layer, nessun Voronoi, nessuna nuova segmentazione.

**Dipendenze**: R2 chiusa (2026-09-21, `d4c8a9e`) — soddisfatte.

## Disegno (deviazione dichiarata dalla ROADMAP §7 «finestre 100×100 µm»)
A densità R2 una finestra 100×100 µm contiene ~109 nuclei in A1 ma ~10 in A5/A6 (errore di Poisson su 3 finestre ±18 %, peggiore del disaccordo da arbitrare). Finestre quadrate a **lato adattato alla densità** (bersaglio 70–110 nuclei): A1 100 µm, A2 120, A3 150, A4 60, A5 250, A6 250. **4 finestre per archetipo, una per ROI su 4 ROI distinti** scelti a caso fra i 5 di R2 (seed 20260921), posizione uniforme dentro il ROI con margine 20 µm, accettata se tessuto valido ≥ 90 % e 0 % pixel nella maschera bolle R2. **40 contorni per archetipo**: 10 bersagli (croci) seedati per finestra; l'annotatore contorna il nucleo più vicino a ciascuna croce (selezione non soggettiva).

Annotazione **alla cieca** (nessuna segmentazione mostrata) con annotatore HTML autonomo (`tools/R2b_annotator.py`, nucleo geometrico `tools/R2b_annotator_core.js` testato con node). Annotatori: **Luca** su tutte le 24 finestre; **Federica Pezzuto** (anatomopatologa) su 2 finestre per archetipo (12, scelte con seed) per l'accordo inter-annotatore.

Aggiunta dichiarata il 2026-09-21 **prima di qualsiasi annotazione**, dopo aver visto la tavola delle finestre (bolle sfocate non mascherate in A6, margine di ghiandola tumorale in A3_w2): l'annotatore può tracciare **zone di esclusione** di due tipi — *illeggibile* (bolla/artefatto) e *fuori archetipo*. La densità primaria è quella su **area netta** (finestra − esclusioni); la densità lorda è riportata a confronto. Gli oggetti dei segmentatori dentro le esclusioni sono conteggiati a parte (misura per BL-026).

## Check C (correttezza software)
| id | Attesa | Soglia |
|---|---|---|
| C-R2b.1 | Finestre dentro il ROI e la maschera tessuto; coordinate finestra → ROI → immagine intera ricostruite | tessuto ≥ 0.90; bolle = 0; ritaglio dal ROI pixel-identico al ritaglio dall'immagine intera |
| C-R2b.2 | Round-trip dello strumento: schermo↔immagine, zoom ancorato, export→JSON→import | errore < 0.5 px su punti/poligoni sintetici; conversione px→µm con `um_per_px` di `results/R1/R1_scale.csv` |
| C-R2b.3 | Appaiamento punto manuale ↔ oggetto (punto dentro la maschera, altrimenti centroide entro 3 µm, assegnazione 1:1) su dati sintetici con recall e precisione noti | recall e precisione recuperati entro 1 % |
| C-R2b.4 | Riproducibilità del campionamento delle finestre con lo stesso seed | md5 degli array identici |

## Check B (plausibilità biologica, contro R2) — fonti: `results/R2/R2_consensus_density.csv`, `R2_archetype_summary.csv`
| id | Attesa | Soglia |
|---|---|---|
| B-R2b.1 | La densità manuale per archetipo cade fra la stima Cellpose e quella StarDist (allargate del 10 %) | 6/6 archetipi; se fuori da entrambi → entrambi di parte nello stesso verso, il consenso R2 non è un limite |
| B-R2b.2 | Il consenso R2 è vicino alla densità manuale | ±15 % in A1–A4; ±20 % in A5–A6 |
| B-R2b.3 | Direzione dei richiami (R2 §4, BL-032/033): recall Cellpose > StarDist in A1, A2, A3; StarDist > Cellpose in A5, A6; A4 entrambi ≥ 0.85; precisione Cellpose < 0.85 in A5 e < 0.70 in A6 | come scritto |
| B-R2b.4 | Ordinamento delle densità manuali | A4 > A1 > A2 > A3 > A5 > A6 |
| B-R2b.5 | Area nucleare da contorni manuali. Alternativa 1 (BL-025 vera): mediana manuale A4 ≥ 1.20 × mediana Cellpose (13.9 µm²). Alternativa 2 (contorni automatici corretti): rapporto manuale/Cellpose in [0.85, 1.15]. Attesa: alt. 1 in A4, alt. 2 in A2/A3/A5; nessuna attesa per A1/A6 | come scritto |

## Controprove
| id | Test | Attesa / soglia |
|---|---|---|
| CP-R2b.1 | Accordo inter-annotatore (Luca vs Federica) su 12 finestre | differenza relativa di conteggio < 10 %; F1 punto-punto (3 µm) > 0.90. Se l'accordo umano fosse peggiore del disaccordo fra segmentatori (11–55 %), la verità manuale non può fare da arbitro |
| CP-R2b.2 | Falsi negativi comuni: punti manuali non appaiati da nessun segmentatore | < 5 % per archetipo; se > 5 % il consenso R2 sottostima strutturalmente |
| CP-R2b.3 | Oggetti fantasma: oggetti Cellpose senza punto manuale entro 3 µm in A5 e A6 | ≥ 15 % in A5, ≥ 30 % in A6 (coerente con BL-032/033); se < 5 % BL-032/033 vanno chiusi come infondati |

Le attese non si modificano dopo il risultato: una discrepanza si scrive nel report e va in BACKLOG o in revisione del design.
