---
type: Decisions
title: "Riesame, Parte B: scelte di attuazione della catena del pannello e del braccio banale B"
description: "Una riga per scelta, con la data, scritta e committata prima del calcolo dei bracci V, T, S e B sul pannello; nessuna scelta presa guardando D."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T13:45:00+01:00 }
sources:
  - id: regola
    resource: "/riesame/REGOLA_B.md"
    title: "Riesame, Parte B: regola di esito scritta prima del calcolo"
  - id: f33
    resource: "/steps/F3.3.md"
    title: "Il simulatore di oggi regge fuori dalle aree su cui è stato tarato?"
---
# Scelte di attuazione della Parte B (braccio «Catena»)

Regola di esito: `docs/riesame/REGOLA_B.md` (commit `45fdbbe`)[^regola]. Misure e catena di partenza: F3.1 e F3.3[^f33]. Ogni riga è scritta prima del calcolo del braccio che tocca; il commit che la contiene è il primo commit di questo file, prima di ogni esecuzione dei bracci (una prova della catena su un'area di taratura, A2_r1, sui geni comuni, è stata fatta prima per controllare che il codice giri: è un'esecuzione del cancello G-4, non di un braccio). Nessuna scelta dipende da D, che questo passo non calcola.

## Catena condivisa (bracci V, T, S, B e, dopo, Xenium)

| id | data | scelta | motivo |
|---|---|---|---|
| C-1 | 2026-11-20 | I conteggi si restringono nell'esportazione dell'area (`X.parquet`, `genes.tsv`), quadratino per quadratino, tenendo l'ordine dei geni e tutti i quadratini; i geni fuori dall'insieme si tolgono, senza rinormalizzare. Insiemi: `panel` = 408 identificativi Ensembl di `riesame/xenium/panel_ids_in_common.txt`; `common` = `results/F3.3/F3.3_crc_common_genes.txt` (18 024); `all` = nessuna restrizione. Script `tools/riesame_panel_restrict.py`. | È quello che vede un metodo su un dato a pannello: esistono solo i conteggi dei geni del pannello. |
| C-2 | 2026-11-20 | RCTD con `tools/F3.1_rctd.R` invariato e con il riferimento di F3.1 (`F3.1/main/ref`): lo script usa le sole righe del riferimento dei geni presenti nell'area, cioè il riferimento ristretto agli stessi geni; nUMI delle cellule del riferimento = totale su tutti i geni del riferimento (come F3.1, F3.3 e la prova di `riesame/xenium/rctd_panel.R`); nUMI dei blocchi = somma dei conteggi ristretti. | Identità con la catena di F3.1 e F3.3 (cancello G-4); è l'uso standard di spacexr (nUMI del riferimento = totale della cellula); restringere anche il nUMI del riferimento toglierebbe con UMI_min = 20 molte cellule del riferimento sui 408 geni. Differenza dichiarata dalla prova di Parte A: là il nUMI dei blocchi era il totale su tutti i geni; qui solo i geni del pannello esistono, quindi i singoletti sul pannello non sono confrontabili con lo 0,9994 di Parte A. |
| C-3 | 2026-11-20 | bin2cell con `tools/F3.1_b2c.py` invariato sui conteggi ristretti (destripe e pareggi con la PCA sui soli geni dell'insieme; nuclei invariati). | Come F3.1. |
| C-4 | 2026-11-20 | Etichette trasferite di F3.3 (`tools/F3_3_common.py`, `transfer_labels`): profili dal riferimento di F3.1 sui soli geni dell'area (quindi ristretti), nuclei con almeno 20 UMI sui geni dell'insieme; per V profili fuori area (si esclude l'area stessa), per T, S e B tutte le cellule del riferimento. Misure M-1...M-7 con `measures()` e varianti M2r, M3r, M3u con `measures_variants()`. Script `tools/riesame_panel_measures.py`. | Come F3.3 (righe IN, OUT e SIM). |
| C-5 | 2026-11-20 | Simulato ristretto (S e B su `panel` e `common`): `origin`, `profiles` e `background` non si copiano e bin2cell e RCTD girano con il lato `simr`, che salta le sole misure di verità. | `origin` non ha i geni e non si può restringere; M-1...M-7 non usano la verità (stesso codice). |
| C-6 | 2026-11-20 | V e T sui geni comuni: si riusano le righe di `results/F3.3/F3.3_area_values.csv` (IN e OUT di livello L2, tumore e stroma), e la catena le rifà su un'area per sezione (V: A2_r1 e A3_r1 dall'esportazione di F3.1 ristretta; T: un'area di P1, P5 e tumore da 6,5 mm). S: F3.3 l'ha misurato sui 18 132 geni del simulato, non sui 18 024 comuni (il ponte di DEC-097 punto 4 era sotto 0,1 h), quindi G-4 per S si fa con `all` su A2_r1_k0 e A3_r1_k0; S sui geni comuni si calcola di nuovo sulle 20 esecuzioni, per confrontare B e S sugli stessi geni. | Cancello G-4 della regola; la lettera di G-4 per S non è applicabile così com'è scritta e lo si dichiara. |
| C-7 | 2026-11-20 | Tabella per area `results/riesame/B_area_values.csv`: colonne di `F3.3_area_values.csv` più `genes`, `arm`, `origin` (chain o F3.3), `n_genes`; una riga per braccio × geni × area × seme; unità come F3.3 (unit = area, u2 = area per la media dei 2 semi). | Formato richiesto dal passo successivo. |
| C-8 | 2026-11-20 | Margine h per insieme di geni (`tools/riesame_B_h.py`): costruzione di F3.1 come `tools/F3.1_intervals.py` e `tools/F3.3_prereg_numbers.py` (funzione `corrected_interval` di `tools/S2_5_common.py`): bootstrap della media delle 5 aree di V, B = 2 000, generatore `default_rng([20261115, i, j])` con i = posizione del tessuto in A1...A6 e j = posizione della misura in M-1...M-6, kappa_boot del tessuto (`R7_kappa_final_v4.csv`, q = bin per M-1, M-4, M-5, M-6 e nuc per M-2, M-3), fattore dell'asimmetria se CV > 0,45; h = (min(1, hi) − lo)/2 per le quote, (hi − lo)/2 per M-2; M-4 almeno 0,01 (DEC-086). Controllo: lo stesso codice sulle righe IN_panel di F3.3 ritrova gli h di `F3.3_prereg_numbers.csv`. | Costruzione già approvata da Luca per queste misure (F3.1, F3.3). |
| C-9 | 2026-11-20 | Al più 2 processi, un thread ciascuno (OMP, OPENBLAS, MKL, NUMBA = 1); esecuzioni in `/mnt/micron/geo_spatialtrans/riesame/B/<braccio>/`, un file per esecuzione, rilanciabili. md5 dei file letti prima e dopo con `tools/riesame_B_md5.sh`. | Server condiviso (regola del compito). |

## Braccio B (banale), `tools/riesame_B_sim.R` e `tools/riesame_B_run.sh`

| id | data | scelta | motivo |
|---|---|---|---|
| B-1 | 2026-11-20 | Geometria: le 5 aree per tessuto di S (A2_r1...A2_r5, A3_r1...A3_r5), con le regioni (`region_df`, `region_polygons`) dei livelli U6 di S1.6; lato dell'esportazione uguale a quello di S. | «Stesse 5 aree (geometria)» della regola. |
| B-2 | 2026-11-20 | Densità del preset: `.TP_TABLE$density` (tumore 8 924, stroma 3 100 nuclei/mm²), la stessa che il preset realizza in S (8 924,0 e 3 099,8 sui 5 livelli U6). In ogni regione il numero di centri è Poisson(densità × area) e i centri sono uniformi (rigetto nel riquadro della regione), senza distanza minima. | Poisson omogeneo della regola. |
| B-3 | 2026-11-20 | Territori: `tessellate_voronoi()` del pacchetto, per regione, senza smussatura (`corner_smoothing = 0`, come S). | Voronoi della regola, stessa funzione di S. |
| B-4 | 2026-11-20 | Area nucleare media del preset: media di `nucleus_area` dei 5 livelli U6 del tessuto (tumore 34,62 µm² su 37 808 nuclei; stroma 36,18 µm² su 15 429). Nucleo = poligono regolare di 64 lati con quell'area esatta (fattore `.dcg_poly_factor` del pacchetto), centrato nel centro, intersecato con il territorio e non riportato all'area. | La regola chiede il cerchio con l'area media del preset, ritagliato. Il preset non ha un'area media unica (modello per tipo e territorio): si usa la media che il preset realizza sulle stesse aree. |
| B-5 | 2026-11-20 | Tipi estratti iid dalla composizione del preset (`tissue_preset_spec()$composition`: tumore epiteliale 0,776428, fibroblasti 0,223572; stroma endoteliale 0,0814, epiteliale 0,4071, fibroblasti 0,3251, immunitario 0,1864). | Regola; S assegna invece le quote esatte. |
| B-6 | 2026-11-20 | Trascritti per cellula: `simulate_cell_profiles()` con un profilo per tipo (`use_subtypes = FALSE`, media dei sottotipi pesata con le quote), conteggi multinomiali (`noise = "M"`, senza fattore gamma per gene), totale della cellula dalla distribuzione dei totali del tipo per il fattore di profondità del preset `visium_hd_depth_factor()` (tumore 0,954, stroma 1,561) per la sovra-generazione M delle strisce. | «Profilo medio del tipo con la profondità del preset»; il pacchetto non ha un totale fisso, e la distribuzione dei totali è quella del preset. |
| B-7 | 2026-11-20 | Nascita dei trascritti: `place_and_move_transcripts()` con r = 1 (densità uguale in nucleo e citoplasma, cioè uniforme nel territorio) e campo di variabilità fine spento (`field = FALSE`). | «Messi a caso dentro la cellula». |
| B-8 | 2026-11-20 | Spostamento spento: w = 0 (niente coda esponenziale), w2 = 0 (niente seconda coda), sigma_in = sigma_e = 1e-9 µm (la funzione chiede sigma > 0; lo spostamento resta sotto 1e-8 µm, sotto la precisione degli overlay di 1e-6 µm); `frame_open = TRUE` come S. | Regola: «cattura dello Step 2 con lo spostamento spento». |
| B-9 | 2026-11-20 | Cattura: `simulate_stripes()` con `capture_params()` del tessuto e `capture_and_bin()` con il fondo b del tessuto (`background = TRUE`), nell'ordine dell'orchestratore `simulate_visium_hd()`. Differenza da S (che usa `simulate_visium_hd_pieces()`): sovra-generazione M = massimo dei fattori di striscia invece del tetto fisso `stripe_ceiling()`; l'assottigliamento c/M dà la stessa profondità attesa. | Usare le funzioni del pacchetto senza modificarlo: l'orchestratore non permette di spegnere lo spostamento, quindi i passi si chiamano uno per uno. |
| B-10 | 2026-11-20 | Semi: seme base 20261124; tabella Mersenne-Twister sulle 10 aree ordinate (A2_r1...A3_r5) × k = 0, 1 (`riesame/B/sim_rds/seed_table.csv`); per ogni esecuzione 6 semi dei passi (poisson, types, stripes, profiles, place, capture) con `sample.int(1e8, 6)`. 2 semi per area. | Semi dichiarati prima del calcolo. |
| B-11 | 2026-11-20 | Esportazione con `tools/F3.1_export_sim.R` e `tools/F3.1_reg_sim.py` invariati (etichetta del quadratino = nucleo che ne copre almeno metà; in_reg come S); misura sul pannello e sui geni comuni con la catena C-1...C-5. | Come S. |
| B-12 | 2026-11-20 | Controllo per esecuzione: quota dei trascritti catturati in un quadratino che tocca il territorio della propria cellula = 1; trascritti persi = 0; densità, composizione e area nucleare realizzate in `riesame/B/sim_rds/<area>_k<k>_check.json`. | Verifica che lo spostamento sia davvero spento. |

[^regola]: [Riesame, Parte B: regola di esito scritta prima del calcolo](/riesame/REGOLA_B.md).
[^f33]: [Il simulatore di oggi regge fuori dalle aree su cui è stato tarato? (F3.3)](/steps/F3.3.md).

## Scelte prese durante il calcolo, prima di vedere D di S e di B (2026-11-20)

| id | data | scelta | motivo |
|---|---|---|---|
| B-13 | 2026-11-20 | La simulazione banale A2_r2 con il secondo seme si ferma in `.s27_cdt()` (triangolazione fallita anche dopo la semplificazione, geometria 6635; probabile coppia di centri di Poisson quasi coincidenti, B-2 non ha distanza minima). A2_r2 entra in B con il solo primo seme; non si estrae un seme sostitutivo. 19 esecuzioni su 20. | Nessun seme scelto dopo il fallimento; la media per area resta definita e l'incertezza fra aree la assorbe. |
| C-10 | 2026-11-20 | Con lo Xenium non usabile (XD-23) i bracci V, T, S, B sul pannello non si calcolano: servivano solo alla domanda sulla cattura (G-1 sul pannello, G-6). Si calcolano S e B sui geni comuni. `tools/riesame_B_h.py` salta gli insiemi di geni senza righe. | REGOLA_B, ultimo capoverso di «Esiti e conseguenze»: il calcolo serve solo a cio' che decide. |
