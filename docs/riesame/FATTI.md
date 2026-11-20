---
type: Facts
title: "Riesame, Parte A: i fatti riletti dai file"
description: "Solo numeri riletti dai file del branch f3-methods (40d3b68) e dai calcoli del riesame; ogni riga ha la sua fonte."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T13:00:00+01:00 }
sources:
  - id: roadmap
    resource: "docs/ROADMAP.md"
    title: "ROADMAP, sezione 7 (stato delle sessioni)"
  - id: f31
    resource: "/steps/F3.1.md"
    title: "I primi metodi di analisi si comportano sui dati simulati come sui dati veri?"
  - id: f32
    resource: "/steps/F3.2.md"
    title: "Come separiamo i dati che tarano il simulatore da quelli che lo verificano?"
  - id: f33
    resource: "/steps/F3.3.md"
    title: "Il simulatore di oggi regge fuori dalle aree su cui è stato tarato?"
  - id: tetto
    resource: "tools/riesame_tetto.py"
    title: "Calcolo del tetto dalle tabelle di F3.3 (questo riesame)"
  - id: xen
    resource: "/mnt/micron/geo_spatialtrans/riesame/xenium/"
    title: "Verifica dei dati Xenium del tumore del colon (sotto-agente, 2026-11-20)"
---
# Fatti della consegna, riletti

Tutti confermati sul file; nessuno smentito.

| fatto | valore riletto | fonte |
|---|---|---|
| sessioni dal 2026-09-18 | 43 righe di sessione, da S0 (2026-09-18) a F3.3 (2026-11-19) | `docs/ROADMAP.md` sezione 7[^roadmap] |
| prima prova d'uso, 36 giudizi | 14 confermati, 4 indeterminati, 18 smentiti | `docs/steps/F3.1.md`, riga 207[^f31] |
| valori del simulatore dalle 30 aree di taratura | 380 su 478 | `docs/steps/F3.2.md`, righe 159 e 175[^f32] |
| per tessuto | da 40 a 80 valori scalari e 5 tabelle da 5 aree di una sola sezione | `results/F3.3/F3.3_D1.csv` |
| fuori campione, stesse sezioni | spostamento −0,20 h (−0,81; 1,12): non rilevato; uno scarto sotto circa 1 h non si vedrebbe | `results/F3.3/F3.3_summary.csv`, P-1[^f33] |
| fuori campione, altri pazienti e sezioni | vero: 2,38 h (2,07; 2,73); degrado del simulato 2,94 h (2,33; 3,43) | `F3.3_summary.csv`, P-2 e P-3 L2 |
| trascritti per cellula, tumore | P2 (taratura) 1 120, P1 640, P5 1 997, tumore da 6,5 mm 270; simulato 1 096 | `results/riesame/tetto_by_measure.csv` (da `F3.3_area_values.csv`) |
| valori veri dentro la variabilità simulata | 18 su 228 | `F3.3_summary.csv`, P-7 |
| vicini dello stesso tipo, stroma | P2 0,886; P1 0,752; P5 0,716; tumore da 6,5 mm 0,885; simulato 0,649 | `tetto_by_measure.csv` |
| backlog | 273 voci, 155 aperte | `docs/BACKLOG.md` (conteggio della colonna Stato) |
| confronti con simulatori esistenti | uno: `SIMPSPACE.md` sul branch `simspace`, datato 28 luglio 2025 | `git show origin/simspace:SIMPSPACE.md` |

Il margine h è la semi-ampiezza dell'intervallo della media vera delle 5 aree di taratura (F3.1). È la precisione di una media, non la variabilità biologica fra sezioni.

# Sessioni per componente (lettura del riesame)

| componente | sessioni | n |
|---|---|---|
| generatore dello strato cellulare (regioni, semina, tessellazione, nuclei) | S1.1, S1.2, S1.3, S1.3b, S1.3c, S1.3d, S1.4, S1.4b, S1.4c, S1.4d, S1.4e, S1.4f, S1.5, S1.5b | 14 |
| disposizione dei tipi e densità (raggruppamento, vuoti) | S1.5c, S1.5d, S1.5e | 3 |
| profili di espressione | S2.2 | 1 |
| spostamento dei trascritti | S2.3 | 1 |
| cattura e griglia | S2.0, S2.1, S2.4, S2.5 | 4 |
| riferimenti sui Visium HD veri | R1, R2, R2b, R2c, R3, R3b, R4, R5, R6, R7 | 10 |
| banco di prova e dati di verifica | F3.1, F3.2, F3.3 | 3 |
| infrastruttura, scala, consolidamento | S0, S1.6, S1.6b, S2.6, S2.7, S2.8, S2.8b | 7 |
| totale | | 43 |

La consegna conta 15 sessioni sulla geometria (S1.3 e S1.4 con le loro ripetizioni, più R2, R2b, R2c, R3, R3b): il conteggio torna con questa tabella.

# Il tetto, calcolato in questo riesame

Calcolo: `tools/riesame_tetto.py`[^tetto]. Legge solo `F3.3_cells.csv` e `F3.3_area_values.csv` (committati); non apre file sigillati, non giudica un cambiamento, quindi non consuma letture dei dati di verifica (DEC-093, R-6). Stimatore quadratico corretto di DEC-098, radice con segno, in unità di h; bootstrap parametrico B = 2 000, seme 20261121. Controllo: lo script ritrova P-1 L1 −0,1977, P-2 L2 2,3752, P-3 L2 2,9368 di F3.3. md5 degli ingressi uguali prima e dopo.

Distanze medie, misure M-1...M-6 (V = sezione di taratura, T = altra sezione vera, S = simulato):

| tessuto | vero T da V | vero T da vero T' | simulato da V | simulato da T |
|---|---|---|---|---|
| tumore del colon | 2,25 (1,63; 2,95) | 2,76 (2,26; 3,32) | 10,27 (9,70; 10,83) | 10,00 (9,61; 10,41) |
| stroma | 1,98 (1,67; 2,31) | 2,06 (1,69; 2,49) | 4,79 (4,62; 4,98) | 4,34 (4,18; 4,52) |
| cervello di topo | 3,53 (2,93; 4,20) | una sola sezione T | 1,32 (0,98; 1,72) | 4,25 (3,75; 4,74) |
| tumore e stroma | 2,12 (1,78; 2,55) | 2,44 (2,11; 2,80) | 8,01 (7,64; 8,38) | 7,71 (7,45; 7,96) |
| tutti (18 celle) | 2,38 (2,07; 2,73) | 2,44 (2,08; 2,78) | 6,59 (6,30; 6,88) | 7,32 (7,08; 7,56) |

Solo i pazienti P1 e P5 (senza il tumore da 6,5 mm, paziente incerto): tumore e stroma, vero T da V 2,24 (1,84; 2,70). Due pazienti veri fra loro: D(P1) 2,04 e D(P5) 2,43, differenza −0,39 ± 0,67.

Per misura, simulato contro tetto (`results/riesame/tetto_by_measure.csv`), 18 celle: 7 oltre il tetto (tumore: marcatori estranei, singoletti, concordanza con il nucleo, vicini dello stesso tipo; stroma: marcatori estranei, singoletti, vicini dello stesso tipo), 5 dentro (trascritti per cellula nei tre tessuti; cervello: concordanza con il nucleo, vicini dello stesso tipo), 6 indeterminate.

Tessuti con una sola sezione (`results/riesame/single_section_tissues.csv`), scarto dentro campione di F3.3 contro il tetto più grande visto altrove (3,53 h, cervello): intestino 3 misure oltre (UMI assegnati 8,0 h; marcatori estranei 5,3; vicini dello stesso tipo −9,0); linfonodo nessuna (marcatori estranei 3,4; vicini −2,4); cuore nessuna.

# Dati Xenium per la Parte B

Fonte: sotto-agente del 2026-11-20[^xen], GEO GSE280314, articolo PMC12165841. Non verificati: pagina 10x del dataset (non aperta), colonne di `transcripts.parquet`, qualità della registrazione Xenium sul Visium HD.

| fatto | valore |
|---|---|
| pazienti con Xenium | P1, P2, P5 (tumore); P2 è la sezione di taratura |
| stessi blocchi | sezioni seriali di 5 µm degli stessi blocchi FFPE, «entro 5-10 µm» (schede GEO) |
| pannello | Xenium Human Colon con 100 geni aggiunti: 422 geni |
| nei 18 024 geni comuni | 407 per simbolo, 408 per identificativo Ensembl |
| contorni | nuclei e cellule per P1, P2, P5 (P2: 326 946 nuclei su 340 837 cellule); trascritti con coordinate |
| marcatori sul pannello, su 50 per tipo | tumore: epiteliale 11, fibroblasti 13; stroma: endoteliale 7, epiteliale 12, fibroblasti 8, immunitario 18 |
| quota degli UMI sui geni del pannello | 3,5 % (area di tumore A2_r1), 4,1 % (stroma A3_r1) |
| RCTD con i soli geni del pannello, un'area per tessuto | tipo principale uguale nel 98,3 % dei blocchi (tumore) e 95,6 % (stroma); singoletti 99,9 % contro 94,7 % (tumore), 99,98 % contro 98,4 % (stroma); blocchi elaborati 8 970 su 10 674 e 5 947 su 12 391 |

[^roadmap]: `docs/ROADMAP.md`, sezione 7.
[^f31]: [I primi metodi di analisi si comportano sui dati simulati come sui dati veri? (F3.1)](/steps/F3.1.md).
[^f32]: [Come separiamo i dati che tarano il simulatore da quelli che lo verificano? (F3.2)](/steps/F3.2.md).
[^f33]: [Il simulatore di oggi regge fuori dalle aree su cui è stato tarato? (F3.3)](/steps/F3.3.md).
[^tetto]: `tools/riesame_tetto.py`, uscite in `results/riesame/`.
[^xen]: `/mnt/micron/geo_spatialtrans/riesame/xenium/xenium_check.md`.
