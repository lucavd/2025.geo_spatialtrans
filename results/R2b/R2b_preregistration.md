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

## Modifiche di metodo dichiarate il 2026-09-21 prima di qualsiasi annotazione (richiesta di Luca: «la complessità introduce errore»)
Annotatore **v1.1 — solo conteggio** (documento `R2b_main_v11` / `R2b_rater2_v11`): un solo compito per l'annotatore (click sul centro del nucleo);
esclusione «illeggibile» disponibile solo nelle finestre A6 (bolle); nessuna croce, nessun contorno, nessun giudizio istologico.
1. Il giudizio «fuori archetipo» è tolto all'annotatore: Federica Pezzuto risponde sì/no per finestra («tessuto omogeneo dell'archetipo?») sulla tavola delle 24 miniature (`R2b_windows_sheet.png`); le finestre non omogenee (candidata: A3_w2, margine di ghiandola tumorale) si riportano **con e senza** nel sommario.
2. Area nucleare (B-R2b.5, BL-025): non più poligoni. Eventuale **fase 2** con calibri (2 click asse maggiore + 2 click asse minore) su 5 nuclei per finestra scelti dai bersagli seedati (`R2b_targets.csv`, i primi 5); area manuale = π·a·b/4, confrontata con l'ellisse equivalente dei segmentatori (`major_um`, `minor_um` di regionprops, stessa formula). Attese e soglie di B-R2b.5 immutate. Se la fase 2 non si fa, B-R2b.5 = «non valutata in R2b» e BL-025 resta aperta.
3. Le attese B-R2b.1–B-R2b.4 e le controprove CP-R2b.1–CP-R2b.3 non cambiano.

## Aggiunta dichiarata il 2026-09-21 (prima di qualsiasi dato): annotatori automatici
Luca fa valutare i 24 ritagli (`R2b_windows_for_models.zip`, solo area di conteggio, nessun contesto, nessun conteggio atteso) a **uno o due modelli visione-linguaggio** con il prompt fissato in `results/R2b/R2b_model_prompt.md` (una immagine per richiesta, contesto pulito). I loro output (`n_nuclei`, `points` in frazioni dell'immagine, `unreadable_fraction`) sono convertiti nel formato dell'annotatore da `tools/R2b_import_model.py` ed entrano nella stessa catena come `rater = model:<nome>`.
Ruolo: **annotatori automatici aggiuntivi**, sullo stesso piano di Cellpose/StarDist/Space Ranger, valutati contro l'annotatore umano primario (Luca). Non sostituiscono la verità manuale; la strategia resta due umani (Luca 24, Federica 12). Attese per i modelli (nessuna fonte in letteratura per questi tessuti a questa scala; registrate come previsione da smentire): errore relativo di conteggio |n_model − n_Luca| / n_Luca mediano per archetipo **< 20 %**; F1 punto-punto (3 µm, se i punti sono forniti) atteso **basso (< 0.6)** perché i VLM localizzano male. Se i modelli battessero i segmentatori dedicati su recall e precisione, sarebbe un risultato da riportare, non da nascondere.

## Controprova aggiunta il 2026-09-21 dopo l'import di GPT-6 Astra, prima di qualsiasi dato umano
| id | Test | Attesa / soglia |
|---|---|---|
| CP-R2b.4 | Modello nullo per la localizzazione: per ogni finestra e annotatore (umano o modello), N punti uniformi nel tessuto valido (N = punti dell'annotatore, 20 repliche) appaiati agli stessi oggetti con le stesse regole → F1 nullo. | F1 dell'annotatore > F1 nullo + 0.20 in ogni archetipo; per gli umani atteso F1 ≫ nullo; per i VLM nessuna attesa (misura) |
Stato al momento della dichiarazione: import GPT-6 Astra (24/24 finestre) eseguito; il modello conta 1.07–2.83 volte i segmentatori; nessun conteggio umano ancora disponibile.

## Controprova aggiunta il 2026-09-21 dopo il conteggio di Luca, prima di calcolarla (obiezione di Luca: «i persi da CP e SD sono le cellule dubbie su cui anche due umani discuterebbero»)
| id | Test | Attesa / soglia |
|---|---|---|
| CP-R2b.5 | Stratificazione del disaccordo per evidenza di ematossilina: per ogni punto manuale, OD del canale ematossilina (Ruifrok) in un disco di raggio 1.5 µm attorno al punto meno OD dell'anello 3–5 µm (stesso arbitro di R2, `R2_exclusive_check.py`). Confronto fra punti appaiati da almeno un segmentatore e punti persi da entrambi (CP e SD). | Attesa (tesi di Luca): OD mediano dei «persi» < OD mediano degli «appaiati» in ≥ 5/6 archetipi, con differenza ≥ 25 % della mediana degli appaiati. Se confermata, il disaccordo sta sugli oggetti pallidi/dubbi e la densità di riferimento va data come intervallo [nuclei evidenti, nuclei totali]. Se smentita (OD simile), i segmentatori perdono nuclei con evidenza di ematossilina pari agli altri. |
Conseguenza dichiarata per BIO_REFERENCES: la densità per archetipo sarà riportata come **intervallo** — estremo inferiore = nuclei "evidenti" (punti manuali appaiati da almeno un segmentatore su Cellpose/StarDist/SR), estremo superiore = conteggio manuale totale — con il valore di Federica, quando disponibile, come terzo punto.
