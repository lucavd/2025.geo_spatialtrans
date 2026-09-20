# geo_spatialtrans — Roadmap v4.1 (2026-09-18, in repo da S0)

## 0. Principi (fissati da Luca, 2026-09-18)

1. **Obiettivo**: dati simulati controllabili in modo granulare **e** biologicamente solidi. La solidità biologica è un'affermazione che va **verificata**, non dichiarata.
2. **Non abbiamo fretta.** Siamo scienziati, non developer: nessun prodotto da consegnare, solo buona scienza.
3. **Nessuno dei due è infallibile.** Ogni affermazione richiede dati, grafici, conferme **e controprove**.
4. **Si sviluppa solo sul server** (lesexp-server, R 4.6.1, 128 core, RTX 4090).
5. **Generalità**: la prima obiezione attesa è "come si comporta con tessuti diversi?". Dobbiamo essere pronti: ogni test gira su un pannello fisso di archetipi tissutali, non su un solo tessuto.

## 1. Stato (sintesi; dettagli nella v2)

Dopo S0 (2026-09-18): server `~/2025.geo_spatialtrans` su branch `step1-cell-layer` da `dev@7836829`, libreria R di progetto costruita (PPM binari), baseline `full_test.R` riproducibile (seed 42) e verificato; design Step 1 in **v1.1** (decisioni 1–11, check C1–C8, 6 preset A1–A6); nessun `04b*` in `R/`. Clone 2024 archiviato in `~/Projects/_archive/`; clone del Mac rimosso. Dopo R1 (2026-09-19): 5 dataset Visium HD reali su `/mnt/micron/geo_spatialtrans/data_real` (A1, A2/A3, A4, A5, A6; 79.8 GB; scala 0.274 µm/px verificata), inventario in `data/real/`.

## 2. Impianto della verifica

Ogni step è chiuso da un **report HTML** (Quarto, generato sul server, consegnato come artifact in chat) con quattro sezioni fisse:

| Sezione | Contenuto | Chi la verifica |
|---|---|---|
| **Pre-registrazione** | Valori attesi e soglie scritti **prima** di eseguire, con fonte (letteratura o dato reale) | Luca approva le attese prima del run |
| **Check C (correttezza software)** | C1-C8 del design: geometria valida, conservazione aree, nuclei contenuti, riproducibilità con seed | Claude esegue; Luca rilancia il comando |
| **Check B (plausibilità biologica)** | Distribuzioni simulate vs distribuzioni empiriche di riferimento (QQ-plot, KS/energy distance, intervalli): densità cellulare, area/raggio cellulare, area nucleare, rapporto N/C, spaziatura (funzione g(r) / Ripley K), per ogni archetipo del pannello | Entrambi, sul report |
| **Controprova** | Il test che, se fallisse, smentirebbe il passo: modello nullo (es. CSR vs Poisson-disk vs nuclei reali), archetipo in cui prevediamo il fallimento (anisotropia), confronto con un'alternativa | Entrambi |

Regole:
- Un valore di riferimento entra in `docs/BIO_REFERENCES.md` solo con **fonte tracciabile** (DOI, oppure dataset reale + script che lo ha stimato). Nessun numero "a memoria" (né mio né tuo).
- Le attese pre-registrate non si aggiustano dopo il risultato: se il risultato le smentisce, la discrepanza si scrive nel report e diventa una voce di `docs/BACKLOG.md` o una revisione del design (sez. 9).
- Ogni report è sottoposto a una **revisione avversariale indipendente** (sotto-agenti che ricalcolano i numeri dal codice e dai dati) prima di essere dichiarato chiuso.
- Uno step per sessione; ogni sessione finisce con commit su branch dedicato dal server e push su origin.

## 3. Pannello di archetipi tissutali (6, fissato il 2026-09-18)

Le densità e i raggi qui sotto sono **ipotesi di partenza** da confermare in `BIO_REFERENCES.md` (binario R, sessione R2-R3), non valori accertati. Dataset reali assegnati in R1 (2026-09-18): A1 `Visium_HD_Mouse_Small_Intestine`; A2 e A3 `Visium_HD_Human_Colon_Cancer` (A3 = regioni stromali); A4 `Visium_HD_Human_Lymph_Node_FFPE`; A5 `Visium_HD_6p5mm_Mouse_Brain`; A6 `Visium_HD_6p5mm_Human_Heart`. Dettagli in `data/real/README.md`.

| # | Archetipo | Perché è nel pannello | Attesa qualitativa (da verificare) | Dataset reale candidato |
|---|---|---|---|---|
| A1 | Epitelio semplice denso (mucosa intestinale/colon) | Caso d'uso primario; cellule addossate | densità alta, poligoni compatti e regolari | Visium HD mouse small intestine (H&E microscopica `.btf`), human CRC FFPE (mucosa normale adiacente) |
| A2 | Tumore solido con microambiente (CRC / mammella) | Eterogeneità intra-regione (tumore + TIL + CAF) | composizione mista, densità intermedia-alta, alta variabilità di area | Visium HD human CRC FFPE |
| A3 | Stroma lasso / connettivo | Cellule rade, molta ECM: il territorio Voronoi ≫ cellula | densità bassa, poligoni grandi, N/C basso | stroma nei dataset A2 / polmone |
| A4 | Tessuto linfoide (follicoli, tonsilla/linfonodo) | Estremo di densità: cellule piccole, nucleo ≈ cellula | densità molto alta, N/C ≈ 0.8-0.9, poligoni minuscoli | da individuare (tonsilla / linfonodo Visium HD, oppure aggregati linfoidi in A2) |
| A5 | Cervello (corteccia: neuroni + glia) | Cellule rade di taglia molto diversa in una matrice non cellulare (neuropilo) | densità bassa-media, distribuzione bimodale delle aree | Visium HD mouse brain |
| A6 | Muscolo / epitelio colonnare (anisotropo) | **Controllo negativo previsto**: il Voronoi isotropo non può riprodurre cellule allungate | ci aspettiamo che i check B su forma (eccentricità) **falliscano**; documentare quanto e come → motiva Step 1.5 (Voronoi pesato/anisotropo) | muscolo liscio/mucosa in A1; cuore se disponibile |

Il pannello entra in ogni test dello Step 1 tramite 6 `tissue_preset` (i 3 del design + 3 nuovi) e 6 tessuti sintetici di riferimento generati con seed fisso. Un preset è "biologicamente solido" solo quando i check B passano contro il dataset reale della sua riga.

## 4. Due binari paralleli

**Binario S — Simulatore** (Step 1 → Step 2): implementa il cell layer e poi la griglia HD.
**Binario R — Riferimenti reali**: costruisce le distribuzioni empiriche contro cui il binario S viene verificato. Parte subito, perché senza riferimenti i check B dello Step 1 non esistono.

Il binario R deve stare **davanti** al binario S: R1-R2 prima di S1.2, R3 prima di S1.4, R4-R5 prima di Step 2.

### Binario R — sessioni

| Sessione | Obiettivo | Output verificabile |
|---|---|---|
| **R1** | Inventario dataset Visium HD pubblici 10x con H&E ad alta risoluzione; download sul server; verifica integrità; mappa dataset → archetipi A1-A6 (e archetipi scoperti: A4, A6 potrebbero non avere un dataset dedicato) | `data/real/README.md` con URL, checksum, risoluzione immagine (µm/px), tessuto; report R1 |
| **R2** | Segmentazione nucleare su H&E (Cellpose/StarDist su GPU; confronto con la segmentazione di Space Ranger v4 dove disponibile) → per archetipo: densità nucleare (cell/mm²), area nucleare, distanza al primo vicino, g(r) | `docs/BIO_REFERENCES.md` v1 (sezione "geometria"); figure per archetipo; **controprova**: due segmentatori concordano? (Dice, differenza di densità) |
| **R3** | Voronoi sui nuclei reali → distribuzione empirica di area/eq_radius del territorio, rapporto N/C, eccentricità dei territori (misura l'anisotropia per A6) | riferimento diretto per i check B di S1.3-S1.4; quantifica quanto il Voronoi "puro" spiega la geometria reale (controprova del modello geometrico stesso) |
| **R4** | Library size e sparsità a 2 µm e 8 µm per archetipo; frazione di bin vuoti; UMI nucleo vs citoplasma | riferimento per Step 2 (ricalibrazione 06e/06f) |
| **R5** | Spot di confine e diffusione: profilo degli spot a cavallo fra nuclei; stima del range di diffusione (C-B di Daniele) | `diffusion_sigma_um` empirico con intervallo; riferimento per Step 2.3 |

### Binario S — sessioni

| Sessione | Obiettivo | Check C | Check B (contro Binario R) | Controprova |
|---|---|---|---|---|
| **S0** | Server: `git pull` su `~/2025.geo_spatialtrans`, archivio del checkout 2024, libreria R (sf, deldir, polylabelr, imager, testthat, quarto) da PPM, `renv.lock` aggiornato; baseline `full_test.R`; design doc v1.1 (C-A, C-B, C8, pannello 6 preset); branch `step1-cell-layer`; `docs/ROADMAP.md`, `CHANGELOG.md`, `docs/BACKLOG.md`, `docs/BIO_REFERENCES.md` (vuoto, con schema); template report Quarto | pacchetti caricano; `full_test.R` riproduce l'RDS esistente | — | — |
| **S1.1** | `extract_regions()` | Σ aree = n_px·px² (2%); `st_is_valid`; n_regioni ≥ k | — (macro-regioni, non biologia) | tessuto con un solo cluster; regione sotto soglia esclusa |
| **S1.2** | `seed_centroids()` (catalogo, composizione, Poisson-disk) sui 6 archetipi | C1 densità ±20%, C4 frazioni ±10%, centroidi dentro la regione, d_min ≥ r_min | densità e g(r) simulati vs R2 per archetipo | Poisson-disk vs CSR vs nuclei reali: quale g(r) è più vicina al reale? Se CSR è altrettanto buona, Poisson-disk non è giustificata |
| **S1.3** | `tessellate_voronoi()` | C2 conservazione area (1%); nessun poligono degenere con smoothing | distribuzione di area/eq_radius vs R3 per archetipo; CV dell'area | A6: eccentricità simulata vs reale → fallimento atteso documentato |
| **S1.4** | `derive_cell_geometry()` | C3 nuclei contenuti, C5, C7 | area nucleare e N/C vs R2-R3 per archetipo | il clamp del nucleo si attiva in quale frazione di cellule? se > 5% il ratio del preset è sbagliato |
| **S1.5** | Orchestrator `place_cells()`, `04b5` (C1-C8), 6 preset, `estimate_pixel_size()`, `plot_cell_layer()`, edge cases, testthat | tutte le asserzioni del design sez. 6.2; 4000×4000 px < 5 min, < 8 GB; riproducibilità | pannello completo: tabella 6 archetipi × metriche B con PASS/WARN/FAIL | `estimate_pixel_size()` su immagini con scala nota ma diversa (1, 2, 5, 10 µm/px): recupera la scala? |
| **S1.6** | Consolidamento: `CLAUDE.md`, `README.md`, revisione avversariale del codice contro il design, merge in `dev`, dossier per Daniele (report S1.1-S1.5 + R1-R3) | — | — | revisori indipendenti |
| **S2.0** | Design doc Step 2: mapping spot↔cellula, profili a livello cellula, aggregazione **con diffusione (C-B, σ da R5)**, binning 2→8 µm, griglia proporzionale (C-A), ricalibrazione da R4. Inviato a Daniele prima del codice | — | — | — |
| **S2.1-2.5** | Implementazione Step 2, una funzione per sessione, stesso impianto (C / B contro R4-R5 / controprova) | | | |
| **Fase 3** | Validazione end-to-end: distribuzioni simulate vs reali per archetipo; test di deconvoluzione; benchmark metodi | | | |

## 5. Workflow su server

- Lavoro in `~/2025.geo_spatialtrans` su branch `step1-cell-layer` (poi `step2-*`); commit e push dal server (remote https con credenziale già funzionante).
- Libreria R del progetto: `~/2025.geo_spatialtrans/renv/library/.../R-4.6/...` da Posit PPM binari (procedura collaudata su MMM: `--vanilla` + `.libPaths()` esplicito).
- (S0) `.Rprofile` del repo imposta `RENV_PATHS_LIBRARY=renv/library` e attiva renv: una sessione R interattiva o `quarto render` vedono la stessa libreria degli script `--vanilla`. I test girano con `R_LIBS_SITE`/`R_LIBS_USER` inesistenti (vedi `R/testing/run_S0.sh`) così l'isolamento non dipende dal contenuto delle site-library di sistema. `tools/setup_r_library.R` è idempotente e ricostruisce la libreria (chiusura ricorsiva, 186 pacchetti).
- Dati reali in `~/2025.geo_spatialtrans/data/real/` (fuori da git, `.gitignore`), con `README.md` e checksum tracciati.
- Segmentazione: ambiente Python dedicato sul server (Cellpose/StarDist, CUDA) — creato in R1/R2.
- Report: `reports/<sessione>.qmd` → HTML sul server → artifact in chat. Luca legge il report e, se vuole, rilancia `Rscript --vanilla R/testing/test_<step>.R` via ssh.
- Il server ha già mostrato cadute di rete a metà job (note MMM): i job lunghi scrivono su percorsi stabili e sono rilanciabili.

## 6. Decisioni prese

| Data | Decisione |
|---|---|
| 2026-09-18 | Sviluppo **solo su server** |
| 2026-09-18 | C-A → `pixel_size_um` obbligatorio + `estimate_pixel_size()` + check C8 |
| 2026-09-18 | Roadmap, changelog, backlog e riferimenti biologici in repo (`docs/`) |
| 2026-09-18 | Riferimenti empirici **solo da dataset pubblici 10x** Visium HD con H&E |
| 2026-09-18 | Pannello di **6 archetipi** (A1-A6), incluso un controllo negativo previsto (A6) |
| 2026-09-18 | Verifica di Luca tramite **report HTML per step** (Quarto, artifact in chat) |
| 2026-09-18 | Ogni step: pre-registrazione → check C → check B → controprova → revisione avversariale |
| 2026-09-18 | (S0) Il check S0 "full_test.R riproduce l'RDS esistente" è riformulato: gli RDS legacy non sono riproducibili (prodotti da `full_test_visHD.R` e da un run a seed casuale) → conservati in `R/testing/legacy/`, nuovo baseline `full_test_result_S0.rds` |
| 2026-09-18 | (R1) Dati reali sul volume `/mnt/micron/geo_spatialtrans/data_real` (symlink `data/real/datasets`); budget R1 ridefinito a 79.8 GB prima del download |
| 2026-09-18 | (R1) Mappa dataset → archetipi: A1 Mouse_Small_Intestine (SR 3.0.0), A2/A3 Human_Colon_Cancer (4.0.1), A4 Human_Lymph_Node_FFPE (4.0.1), A5 6p5mm_Mouse_Brain (4.1.0), A6 6p5mm_Human_Heart (4.1.0); A3 = regione stromale interna al CRC |
| 2026-09-18 | (R1) Per ogni dataset si scaricano: immagine microscopica, immagine CytAssist, binned_outputs, spatial, segmented_outputs + barcode_mappings (SR ≥ 4), metrics, web_summary. Mai FASTQ/cloupe/molecule_info |
| 2026-09-19 | (R1) La scala µm/px di un dataset reale viene da `scalefactors_json.json` verificata con la griglia dei bin; l'OME-XML `PhysicalSizeX` è ammesso solo se verificato contro la griglia; il tag TIFF di risoluzione non è mai una fonte |
| 2026-09-19 | (R1) I riferimenti biologici R2–R5 si stimano su regioni annotate per archetipo, non sull'intera maschera `in_tissue` |
| 2026-09-18 | (S0) Design Step 1 v1.1 = contratto per S1.1–S1.6: decisioni 10 (C-A, `pixel_size_um` obbligatorio, C8) e 11 (C-B, diffusione in Step 2); pannello A1–A6 = 6 preset |
| 2026-09-18 | (S0) Clone del Mac rimosso; il solo repo di lavoro è `~/2025.geo_spatialtrans` sul server |
| 2026-09-18 | **Una sessione = una chat.** Ogni step del protocollo si svolge in una chat nuova; a chiusura ci si ferma. Nessuna eccezione |

## 7. Stato delle sessioni

| Sessione | Data | Stato | Report | Note |
|---|---|---|---|---|
| S0 | 2026-09-18 | **chiusa con riserve** | `reports/S0.html` | Check C 17/17 + libreria 35/35; controprove PASS; revisione avversariale 10/10 confermate. Riserve: BL-007 (preset A4–A6 senza fonte), BL-008 (restore da zero non testato), BL-010/011 (pipeline attuale: n_cells non rispettato, filtro >50k UMI) |
| R1 | 2026-09-18/19 | **chiusa con riserve** | `reports/R1.html` | 46 pagine inventariate; 5 dataset / 79.8 GB su micron, md5 38/38; check C 30/30; controprova scala 5/5; istologia confermata (F. Pezzuto). Riserve: BL-014 (in_tissue ≠ archetipo), BL-015 (soglia µm/px senza fonte), BL-018 (A3 senza dataset dedicato), BL-019 (A1 solo SR 3.0.0); revisione avversariale 30/30 confermate, BL-020…BL-023 |
| R2 | — | prossima | — | segmentazione nucleare su H&E (Cellpose/StarDist, GPU) per archetipo; dipende da R1 (chiusa) |

**Prossima sessione: R2** (segmentazione nucleare; ambiente Python GPU da creare in `.venv` con `tools/setup_python_env.sh`) oppure **S1.1** (`extract_regions()`), senza dipendenze reciproche. Da S1.2 in avanti nessuno step del binario S parte se il riferimento R corrispondente non è pronto; S1.1 (`extract_regions()`) non ha dipendenze dal binario R e può alternarsi con R1–R2.

## 8. Coordinamento fra sessioni

**Principio**: la memoria del progetto è il repo sul server. Ciò che non è scritto in `docs/` non esiste per la sessione successiva.

### File di coordinamento (tutti versionati)
- `docs/ROADMAP.md` — questo documento: stato di ogni sessione (aperta / chiusa / chiusa con riserve), tabella decisioni datata.
- `CHANGELOG.md` — una voce per sessione: cosa è cambiato nel codice e nei documenti.
- `docs/BACKLOG.md` — idee e problemi emersi ma non affrontati, con la sessione di origine.
- `docs/BIO_REFERENCES.md` — valori di riferimento biologici, ognuno con fonte (DOI o script + dataset).
- `docs/sessions/YYYY-MM-DD_<id>.md` — nota di passaggio: obiettivo, fatto, verificato (check C/B/controprova con esito), verifiche di Luca, aperto, prossimo step, comandi esatti per riprodurre ogni numero del report.
- `reports/<id>.html` — report Quarto dello step (anche come artifact in chat).

### Rituale di apertura (Claude)
1. `git pull`; lettura di ROADMAP, CHANGELOG e ultima nota di sessione.
2. Stato del server: branch, working tree, ultimo commit, libreria R carica, dati reali presenti (checksum).
3. Dichiarazione in chat: obiettivo unico, dipendenze soddisfatte (sì/no), attese pre-registrate da approvare.
4. Avvio solo dopo l'approvazione di Luca.

### Rituale di chiusura (Claude)
1. Report generato; revisione avversariale eseguita e riportata.
2. CHANGELOG, ROADMAP (stato sessione), BACKLOG aggiornati.
3. Nota di sessione scritta.
4. Commit (messaggio con id sessione) + push dal server; report salvato come artifact.
5. Consegna in chat: link al report, esito in tre righe, la questione aperta più importante, prossimo step proposto.

### Regole per Luca
- Una decisione vale solo quando è nella tabella "Decisioni" con data.
- Le verifiche indipendenti di Luca (report letto, comando rilanciato, obiezioni) vengono scritte in chat e trascritte nella nota di sessione.

### Identificativi di sessione
`S0`, `S1.1`…`S1.6`, `S2.0`…, `R1`…`R5`, `F3.x`. Sessioni ripetute per lo stesso step: suffisso lettera (`S1.2b`).
