# CHANGELOG — geo_spatialtrans

Una voce per sessione (id sessione nel titolo). Lo stato di ogni sessione è in `docs/ROADMAP.md`; le idee rinviate in `docs/BACKLOG.md`.

## S0 — 2026-09-18 — Infrastruttura server, baseline riproducibile, design v1.1 (branch `step1-cell-layer`)

### Server / ambiente
- Checkout `~/2025.geo_spatialtrans` portato da `dev@53103b0` (17 commit indietro) a `dev@7836829`; creato branch `step1-cell-layer`.
- Clone stantio `~/Projects/2025.geo_spatialtrans` (main, dic. 2024, remote git@ non funzionante) rinominato in `~/Projects/_archive/2025.geo_spatialtrans_2024-12_stale` (nessuna cancellazione).
- Clone locale sul Mac di Luca rimosso (sviluppo solo su server).
- Nuova libreria R di progetto `renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu` costruita da binari Posit PPM con `tools/setup_r_library.R` in tre passaggi (127 → 164 → 186 pacchetti: richiesti, namespace risolti dal sistema, chiusura ricorsiva Depends/Imports/LinkingTo); 0 compilazioni da sorgente. Include sf, deldir, polylabelr, imager, testthat, quarto e le dipendenze della pipeline attuale. Autosufficiente: nessun namespace dalla site-library di sistema.
- `.Rprofile`: imposta `RENV_PATHS_LIBRARY` su `renv/library` prima di `renv/activate.R` (il progetto ha un `DESCRIPTION`, quindi renv avrebbe usato `~/.cache/R/renv/library/...`); così `quarto render` e le sessioni R interattive vedono la stessa libreria degli script `--vanilla`.
- `renv.lock` rigenerato per R 4.6.1 (era R 4.5.0), 186 pacchetti.

### Codice
- `R/testing/full_test.R` e `full_test_visHD.R`: header con `.libPaths()` esplicito sulla libreria di progetto; nuove opzioni `--seed=<int>` e `--tag=<str>`; gli output non condividono più lo stesso path (`R/testing/full_test_result_<tag>.rds`, `results/full_simulation_data_<tag>.rds`). Nessun cambiamento alla configurazione o alla pipeline.
- Nuovi `R/testing/test_S0_library.R`, `R/testing/test_S0_baseline.R`, `R/testing/run_S0.sh` (un comando, PASS/FAIL per asserzione; `R_LIBS_SITE`/`R_LIBS_USER` inesistenti per isolamento esplicito; scrive `results/S0_checksums.md5`).
- `R/testing/full_tissue_complex.png` rigenerato dal run baseline seed 42 (effetto collaterale di `full_test.R`, BL-013).
- Tabelle di verifica tracciate in `results/S0_*.csv` (i `.rds` restano fuori git).
- RDS legacy spostati in `R/testing/legacy/` (vedi nota di provenienza `R/testing/legacy/README.md`): non erano riproducibili da `full_test.R` (prodotti da `full_test_visHD.R` e da un run a seed casuale 1200).

### Documenti
- `docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md` → **v1.1**: C-A (scala px→µm obbligatoria, `estimate_pixel_size()`, `pixel_size_source`, check C8), C-B (diffusione mRNA come requisito di Step 2), pannello A1–A6 → 6 preset, refuso C6 corretto, storico §9.
- Nuovi `docs/ROADMAP.md` (v4 + stato sessioni), `CHANGELOG.md`, `docs/BACKLOG.md`, `docs/BIO_REFERENCES.md` (solo schema, zero valori), `docs/sessions/2026-09-18_S0.md`, `reports/_template.qmd`, `reports/S0.qmd` + `S0.html`.
