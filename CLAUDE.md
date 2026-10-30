# Istruzioni per chi lavora sul codice di geo.spatialtrans

Simulatore di dati Visium HD con un livello cellulare (territori di Voronoi, griglia a 2 µm, binning a 8 µm), costruito per mettere alla prova i
metodi di analisi. PI: Luca Vedovelli; revisore: Daniele Sabbatini. Questo file vale per ogni agente o persona che modifica il repository; il
quadro del progetto per chi lo usa è in `README.md` (in inglese).

## Stato

- Step 1 (livello cellulare, `R/04b*.R`): chiuso e consolidato, unito in `dev` (S1.6, S1.6b).
- Step 2 (griglia di Visium HD): da progettare. Il prossimo step lo sceglie Luca: proposta e alternative in `docs/ROADMAP.md`, sezione 7.
- Il codice del 2025 (`R/01`-`R/03`, `R/04_clustering.R`, `R/05`-`R/06k`, `R/test functions/`) è la pipeline di espressione precedente:
  non è validata con il livello cellulare e non va estesa. `generate_synthetic_tissue()` e `cluster_image()` restano in uso per gli input
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

- R 4.6.1 sul server (Ubuntu 24.04), libreria del progetto `renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu`, ricostruibile da zero
  con `Rscript --vanilla tools/restore_library.R` (versioni di `renv.lock`, binari Posit PPM). `tools/setup_r_library.R` aggiunge i pacchetti
  mancanti all'ultima versione disponibile, senza garantire quelle del lockfile.
- Ogni script si lancia con `Rscript --vanilla` e chiama `.libPaths()` sulla libreria del progetto: `.Rprofile` attiva renv e altrimenti
  sposterebbe la libreria. `quarto render reports/<id>.qmd` dalla radice del repository funziona così com'è.
- Test isolati: `R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent`, perché `.libPaths()` da solo non toglie le librerie dell'utente.
- La batteria gira solo sul server: legge i dati reali del volume dati, file derivati ignorati da git (maschere e nuclei di R2, finestre di R2b)
  e `.venv`. In un clone nuovo li collega `tools/link_server_inputs.sh` (S1.6b). Un test nuovo che legge un file ignorato da git lo dichiara
  nella sua intestazione.
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
python3 -B tools/okf.py check --bundle docs      # controllo dei documenti OKF, prima di ogni commit
quarto render reports/<id>.qmd                   # report dello step
```

## Stile del codice

- R in snake_case; un file numerato per stadio (`R/04b0_tissue_presets.R` ... `R/04b6_estimate_pixel_size.R`, orchestratore
  `R/04b_cell_placement.R`); funzioni piccole a responsabilità singola con intestazione roxygen; funzioni interne con il prefisso del file
  (`.sc_` nella semina, `.tv_` nella tassellazione, `.dcg_` nella geometria).
- Nessun interruttore di prova nel codice di produzione: i mutanti sono copie modificate del codice (`tools/mutate.py`,
  `R/testing/mutants_<step>.py`).
- Test in `R/testing/test_<step>.R`, eseguibili con un solo comando che stampa PASS o FAIL per asserzione e `RESULT n/m PASS`; runner dei
  cancelli in `R/testing/run_<step>.sh`. I test non scrivono file versionati (log nella cartella dei log del run).
- Un'opzione nuova che cambia il risultato predefinito è una decisione; per riprodurre gli step passati si aggiunge un argomento o
  un'opzione di sessione (`options(geo.spatialtrans.*)`), non si tocca il predefinito.

## Documenti

- Report e documenti del progetto in italiano; codice, tabelle e figure in inglese.
- Coordinamento: `docs/ROADMAP.md` (sezione 7 = stato delle sessioni e prossima proposta; sezione 6 = decisioni fino al 2026-10-20),
  `CHANGELOG.md`, `docs/BACKLOG.md`, `docs/BIO_REFERENCES.md`; design dello Step 1 in `docs/superpowers/specs/`.
- Dal 2026-10-20 i documenti nuovi seguono l'Open Knowledge Format (OKF v0.2): `docs/steps/<id>.md` (uno per step), `docs/decisions/DEC-nnn.md`
  (una decisione per file, con chi l'ha presa e con quali parole), `docs/index.md`, `docs/log.md`. Frontmatter YAML con `type`, `title`,
  `description`, `status`, `generated`, `sources`; ogni codice di decisione è un link al suo file. Indici generati da `tools/okf.py`, mai a mano.
- Una decisione vale quando è in `docs/decisions/` o nella tabella della sezione 6 di ROADMAP. Le scelte che cambiano il modello generativo,
  i riferimenti biologici, il pannello, il perimetro degli step o gli intervalli di rilevanza sono di Luca.
- Le note di sessione fino a R2c restano in `docs/sessions/` come sono.

## Git

- Remote https con credenziale del server (`origin`); il push via ssh verso GitHub non funziona dal server.
- I rami `main`, `dev` e `step1-cell-layer` hanno regole di GitHub (niente cancellazioni né push forzati). Il push della credenziale del server
  passa e stampa «Bypassed rule violations» o «Cannot update this protected ref»: verificare con `git ls-remote origin <ramo>`.
- Messaggi di commit `<id dello step>: <titolo>`; lavoro non committato sul server vuol dire che un'altra sessione è attiva: fermarsi e chiedere.

## Trappole note

- `git log --oneline | head -1` in uno script con `pipefail` esce con 141: usare `git log --oneline -n 1`.
- `grep -c` esce con 1 quando conta zero righe: aggiungere `|| true`.
- `pkill -f <schema>` uccide anche la shell che lo lancia se lo schema compare nella riga di comando.
- `sf::st_as_sfc()` su geometrie di `geos_set_precision()` eredita la precisione come fattore di scala: passare dal WKB (`geos_write_wkb()`).
- `sf::st_within(<molti poligoni piccoli>, <un poligono grande>)` non finisce: mettere il poligono grande per primo (`st_contains`, `st_intersects`).
- `deldir` scala male oltre 1e5 punti; `sf::st_voronoi()` è decine di volte più veloce sugli stessi punti.
- Le liste R indicizzate per nome sono lineari: oltre 1e5 elementi indicizzare per posizione.
