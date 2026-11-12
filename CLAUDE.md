# Istruzioni per chi lavora sul codice di geo.spatialtrans

Simulatore di dati Visium HD con un livello cellulare (territori di Voronoi, griglia a 2 µm, binning a 8 µm), costruito per mettere alla prova i
metodi di analisi. PI: Luca Vedovelli; revisore: Daniele Sabbatini. Questo file vale per ogni agente o persona che modifica il repository; il
quadro del progetto per chi lo usa è in `README.md` (in inglese).

## Stato

- Step 1 (livello cellulare, `R/04b*.R`): consolidato in S1.6 e unito in `dev` in S1.6b.
- Step 2 (griglia di Visium HD, `R/05b*.R`): consolidato in S2.8 e unito in `dev` in S2.8b. Ingresso unico `simulate_visium_hd_pieces()`
  ([DEC-079](docs/decisions/DEC-079.md)); `simulate_visium_hd()` resta come motore di riferimento delle batterie di S2.5 e S2.6. Guida
  d'uso in inglese: `vignettes/visium_hd_simulation.qmd`; documento di progetto `docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md` (v0.6).
- Il prossimo step lo sceglie Luca: proposta e alternative in `docs/ROADMAP.md`, sezione 7.
- `main` è fermo al 2025 e non si tocca finché Luca non decide (S2.8b): si lavora su rami di step che avanzano `dev`.
- Il codice del 2025 (`R/01*.R`-`R/03*.R`, `R/04_clustering.R`, `R/05_grid_sampling.R`, `R/06*.R`, `R/test functions/`) è la pipeline di
  espressione precedente: non è validata con il livello cellulare e non va estesa. `generate_synthetic_tissue()` e `cluster_image()` restano in uso per gli input
  sintetici dei test (`tools/S1.1_inputs.R`).

## Principi (fissati da Luca il 2026-09-18)

1. Dati simulati controllabili in modo granulare e biologicamente solidi; la solidità biologica si verifica, non si dichiara.
2. Non abbiamo fretta: scienza, non prodotto.
3. Nessuno è infallibile: ogni affermazione richiede dati, grafici, conferme e controprove.
4. Sviluppo solo sul server del laboratorio (`~/2025.geo_spatialtrans`), mai su un portatile.
5. Ogni test gira sul pannello dei sei tessuti (intestino di topo, tumore del colon, stroma del tumore, linfonodo, cervello di topo, cuore),
   mai su uno solo.
6. Nessun valore biologico a memoria: entra in `docs/BIO_REFERENCES.md` solo con un DOI o con lo script e il dataset che lo hanno stimato.
7. Uno step per sessione; la sessione successiva riparte solo da `docs/`.

## Regole scientifiche

- **Numeri dai file**: ogni numero di un report o di un documento si rilegge dal file che l'ha prodotto.
- **Dati reali e baseline intatti**: `data/real/datasets` (link al volume dati del server) è in sola lettura; md5 dei file letti prima e dopo.
- **Semi dichiarati e invariati** fra rerun dello stesso confronto; ogni deroga è una decisione registrata.
- **Incertezza sull'unità indipendente** (ROI, finestra, seme), non sulla cellula: bootstrap per blocchi, a due stadi quando servono due livelli.
- **Pre-registrazione**: cancelli (controlli binari del software) e previsioni con intervallo di rilevanza scritti e committati prima del
  calcolo; l'ok di Luca sulle previsioni si registra nel file dello step. Dopo il risultato non si cambiano.
- **Controprova**: ogni step ne ha almeno una, il test che, se fallisse, smentirebbe il passo.
- Il cuore è un controllo negativo previsto per il Voronoi isotropo: il suo fallimento è una previsione, non un'eccezione.

## Ambiente

- R 4.6 sul server (Ubuntu 24.04; 4.6.2 dal 2026-11-05, aggiornata da apt, prima 4.6.1), libreria del progetto `renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu`, ricostruibile da zero
  con `Rscript --vanilla tools/restore_library.R` (versioni di `renv.lock`, binari Posit PPM). `tools/setup_r_library.R` aggiunge i pacchetti
  mancanti all'ultima versione disponibile, senza garantire quelle del lockfile.
- Ogni script si lancia con `Rscript --vanilla` e chiama `.libPaths()` sulla libreria del progetto: `.Rprofile` attiva renv e altrimenti
  sposterebbe la libreria. `quarto render reports/<id>.qmd` dalla radice del repository funziona così com'è.
- Test isolati: `R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent`, perché `.libPaths()` da solo non toglie le librerie dell'utente.
- La batteria gira solo sul server: legge i dati reali del volume dati, file derivati ignorati da git (maschere e nuclei di R2, finestre di R2b)
  e `.venv`. In un clone nuovo li collega `tools/link_server_inputs.sh` (S1.6b). Un test nuovo che legge un file ignorato da git lo dichiara
  nella sua intestazione.
- Lo Step 2 legge i profili di espressione dal volume dati (`/mnt/micron/geo_spatialtrans/S2.2/profiles`, opzione
  `geo.spatialtrans.profiles_dir`); gli script di taratura e di prova dello Step 2 (`tools/S2.*`) leggono anche i livelli cellulari delle
  aree vere da `/mnt/micron/geo_spatialtrans/S1.6/layers`. Fuori dal server non gira (BL-181).
- La guida d'uso si rende dalla radice con `R_LIBS=$PWD/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu quarto render
  vignettes/visium_hd_simulation.qmd` (circa 3 minuti): nella cartella `vignettes/` non c'è `.Rprofile`, quindi la libreria va data con `R_LIBS`.
- Python: `.venv` (Cellpose-SAM, torch con GPU); StarDist con `tools/py_stardist.sh`.
- Calcolo: al più 50 core in tutto fra i processi R concorrenti; per run lunghi il limite termico pratico è di 20-28 worker a un thread.
  Lanciare con `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` e `data.table::setDTthreads(1L)` nei worker; numero di worker letto da un file
  `WORKERS` del run; uscite su percorsi stabili, un file per unità, così che il run si possa riprendere.

## Comandi

```bash
Rscript --vanilla tools/restore_library.R        # libreria dal lockfile
bash tools/link_server_inputs.sh                 # in un clone nuovo sul server: input dei test fuori da git
bash R/testing/run_S1.6.sh test                  # batteria completa dello Step 1 (GATE G-1 PASS con 4 FAIL noti)
bash R/testing/run_S1.6.sh mutants               # copie alterate del codice che i test devono rilevare
bash R/testing/run_S1.6.sh env                   # DESCRIPTION contro il codice, restore in una libreria vuota
for s in 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8; do Rscript --vanilla R/testing/test_S$s.R; done   # batterie dello Step 2 (contare i FAIL)
python3 -B tools/okf.py check --bundle docs      # controllo dei documenti OKF, prima di ogni commit
quarto render reports/<id>.qmd                   # report dello step
```

## Stile del codice

- R in snake_case; un file numerato per stadio (`R/04b0_tissue_presets.R` ... `R/04b6_estimate_pixel_size.R`, orchestratore
  `R/04b_cell_placement.R`; `R/05b1_build_hd_grid.R` ... `R/05b6_simulate_visium_hd_pieces.R` per lo Step 2); funzioni piccole a
  responsabilità singola con intestazione roxygen; funzioni interne con il prefisso del file o dello step che le ha introdotte (`.sc_`
  nella semina, `.tv_` nella tassellazione, `.dcg_` nella geometria; `.hd_`, `.s22_`, `.s23_`, `.s24_`, `.s25_`, `.p27_` nello Step 2).
- Lo Step 2 ha una sola legge e due motori: una modifica del modello va fatta in entrambi, oppure si dichiara che il motore di riferimento
  resta indietro. La funzione a pezzi deve dare lo stesso risultato per ogni lato dei pezzi e numero di processi (test di S2.7, gruppo a).
- Ogni modifica del codice raggiunto da `simulate_visium_hd_pieces()` (le sue funzioni interne `.p27_*`, quelle dei passi che chiama e le
  costanti che leggono) cambia l'impronta della ripresa (`.geo_code_hash()`): le cartelle `out_dir` scritte prima non si riprendono più,
  ed è voluto ([DEC-080](docs/decisions/DEC-080.md)). Le modifiche al solo motore di riferimento o alle funzioni dei passi che la funzione a
  pezzi non chiama (`build_hd_grid()`, `place_and_move_transcripts()`, `capture_and_bin()` e altre) non la cambiano: se una correzione del
  modello tocca una di queste, va riportata nella funzione a pezzi. Una funzione chiamata con `do.call("nome")` o `match.fun()` sfugge
  all'impronta: non usarle. Le versioni dei pacchetti non entrano nell'impronta; la versione di R sì, indirettamente: le impronte dei profili e del
  livello cellulare serializzano gli oggetti, e l'intestazione della serializzazione contiene la versione di R, quindi dopo un aggiornamento di
  R una cartella `out_dir` vecchia non si riprende. Per lo stesso motivo un file salvato con `saveRDS()` porta la versione di R che l'ha
  scritto: confronti di md5 di file fra versioni di R diverse falliscono anche con oggetti identici (S2.8b: R 4.6.1 → 4.6.2 il 2026-11-05).
- Nessun interruttore di prova nel codice di produzione: i mutanti sono copie modificate del codice (`tools/mutate.py`,
  `R/testing/mutants_<step>.py`).
- Test in `R/testing/test_<step>.R`, eseguibili con un solo comando che stampa PASS o FAIL per asserzione e `RESULT n/m PASS`; runner dei
  cancelli in `R/testing/run_<step>.sh` (alcuni test vecchi usano altri formati, per esempio `RISULTATO` in S1.5b). I test non scrivono file
  versionati (log nella cartella dei log del run).
- Un'opzione nuova che cambia il risultato predefinito è una decisione; per riprodurre gli step passati si aggiunge un argomento o
  un'opzione di sessione (`options(geo.spatialtrans.*)`), non si tocca il predefinito.

## Documenti

- Report e documenti del progetto in italiano; codice, tabelle e figure in inglese.
- Coordinamento: `docs/ROADMAP.md` (sezione 7 = stato delle sessioni e prossima proposta; sezione 6 = decisioni fino al 2026-10-20),
  `CHANGELOG.md`, `docs/BACKLOG.md`, `docs/BIO_REFERENCES.md`; design dello Step 1 in `docs/superpowers/specs/`.
- Dal 2026-10-20 i documenti nuovi seguono l'Open Knowledge Format (OKF v0.2): `docs/steps/<id>.md` (uno per step), `docs/decisions/DEC-<nnn>.md`
  (una decisione per file, con chi l'ha presa e con quali parole), `docs/index.md`, `docs/log.md`. Frontmatter YAML con `type`, `title`,
  `description`, `status`, `generated`, `sources`; ogni codice di decisione è un link al suo file. Indici generati da `tools/okf.py`, mai a mano.
- Una decisione vale quando è in `docs/decisions/` o nella tabella della sezione 6 di ROADMAP. Le scelte che cambiano il modello generativo,
  i riferimenti biologici, il pannello, il perimetro degli step o gli intervalli di rilevanza sono di Luca.
- Le note di sessione fino a R2c restano in `docs/sessions/` come sono.

## Git

- Remote https con credenziale del server (`origin`); il push via ssh verso GitHub non funziona dal server.
- Un insieme di regole di GitHub vale per tutti i rami salvo il ramo della tesi «thesis/expression-sources»: creazione, aggiornamento, cancellazione e push non
  fast-forward sono riservati a chi può aggirare le regole. La credenziale del server può: il push passa e stampa «Bypassed rule violations»
  o «Cannot update this protected ref»; verificare con `git ls-remote origin <ramo>`. Mai push forzati.
- Messaggi di commit `<id dello step>: <titolo>`; lavoro non committato sul server vuol dire che un'altra sessione è attiva: fermarsi e chiedere.

## Trappole note

- `git log --oneline | head -1` in uno script con `pipefail` esce con 141: usare `git log --oneline -n 1`.
- `grep -c` esce con 1 quando conta zero righe: aggiungere `|| true`.
- `pkill -f <schema>` uccide anche la shell che lo lancia se lo schema compare nella riga di comando.
- `sf::st_as_sfc()` su geometrie di `geos::geos_set_precision()` eredita la precisione come fattore di scala: passare dal WKB (`geos::geos_write_wkb()`).
- `sf::st_within(<molti poligoni piccoli>, <un poligono grande>)` non finisce: mettere il poligono grande per primo (`st_contains`, `st_intersects`).
- `deldir` scala male oltre 1e5 punti; `sf::st_voronoi()` è decine di volte più veloce sugli stessi punti.
- Le liste R indicizzate per nome sono lineari: oltre 1e5 elementi indicizzare per posizione.
- `metadata$package_version` viene da `.geo_version()`, che legge la `DESCRIPTION` della cartella di lavoro (opzione `geo.spatialtrans.root`):
  lanciare dalla radice del repository, altrimenti la versione risulta sconosciuta.
- Le sezioni intere si simulano a pezzi con `assemble = "light"`, ma il processo principale che cuce i pezzi arriva a circa 40 GB nel tumore
  del colon (S2.7): contare la sua memoria, non solo quella di un pezzo.
- La seconda coda dello spostamento arriva oltre il margine di 280 µm: con `frame_open = FALSE` i trascritti che escono dalla griglia si
  perdono (`truth$cells$n_lost`), con `frame_open = TRUE` la seconda coda si ripiega sul riquadro `frame` (qualche trascritto della prima
  coda può ancora uscire: `n_lost` piccolo ma non nullo).
- `R/testing/test_S2.3.R` esce con stato 0 anche quando un'asserzione fallisce: contare le righe `FAIL` (BL-256).
- Una cartella `out_dir` con i file dei pezzi ma senza *run.rds* si riusa senza controlli: non cancellare quel file da solo (BL-254).
