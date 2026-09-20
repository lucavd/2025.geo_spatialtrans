# CHANGELOG — geo_spatialtrans

Una voce per sessione (id sessione nel titolo). Lo stato di ogni sessione è in `docs/ROADMAP.md`; le idee rinviate in `docs/BACKLOG.md`.

## R1 — 2026-09-18/19 — Inventario e download dei dataset Visium HD di riferimento (branch `step1-cell-layer`)

### Dati
- Inventario di 46 pagine dataset Visium HD 10x (582 file con md5 pubblicato), letto con Claude in Chrome il 2026-09-18: `data/real/inventory/datasets_10x_pages_2026-09-18.csv`, `files_10x_pages_2026-09-18.csv`.
- Scaricati 5 dataset (38 file, 79.84 GB) su `/mnt/micron/geo_spatialtrans/data_real` (symlink `data/real/datasets`), md5 = pubblicato 38/38 (download + ricalcolo indipendente): A1 `Visium_HD_Mouse_Small_Intestine` (SR 3.0.0), A2/A3 `Visium_HD_Human_Colon_Cancer` (4.0.1), A4 `Visium_HD_Human_Lymph_Node_FFPE` (4.0.1), A5 `Visium_HD_6p5mm_Mouse_Brain` (4.1.0), A6 `Visium_HD_6p5mm_Human_Heart` (4.1.0). Niente FASTQ/cloupe/molecule_info.
- `data/real/README.md` (mappa dataset → archetipi, scala, verifica istologica), `data/real/checksums.md5`, `data/real/inventory/R1_download_manifest.csv`, `R1_archetype_map.csv`, `R1_download_status.tsv`.

### Codice / ambiente
- `tools/R1_download.sh` (download riprendibile con verifica md5), `tools/R1_verify.py` (metadati BigTIFF, scala, controprova griglia, ritagli), `R/testing/test_R1.sh` (un comando, PASS/FAIL; `R1_FAST=1` salta il md5).
- Ambiente Python di progetto: `tools/setup_python_env.sh` (uv, Python 3.12, `.venv/` ignorato) + `tools/requirements-r1.txt` / `tools/requirements.lock` (tifffile, imagecodecs, zarr, pyarrow, h5py, scikit-image).
- Tabelle di verifica tracciate in `results/R1/*.csv`; ritagli in `results/R1/crops/`.

### Esiti
- Check C 30/30 PASS; controprova 1 (scala dalla griglia dei bin) 5/5 entro 0.05 %; controprova 2 (istologia) confermata da Federica Pezzuto (anatomopatologa) per A1, A2, A4, A5, A6 con due rilievi (A4 ritaglio non rappresentativo; A6 bolle di montaggio).
- Deviazioni registrate: budget disco ridefinito (79.8 GB su micron); previsione B-R1.2 smentita (esistono dataset A4 e A6); soglia B-R1.3 senza fonte (il paper di Cellpose non fissa un minimo di µm/px).
- Reperto: il tag TIFF di risoluzione è placeholder (96 dpi) in tutte le immagini, mentre l'OME-XML `PhysicalSizeX` è corretto (entro 0.08 % da `scalefactors_json`) → `estimate_pixel_size()` può usare l'OME-XML, mai il tag, e deve verificare contro la griglia dei bin.
- Revisione avversariale (2026-09-20, sotto-agente indipendente): 30/30 affermazioni confermate, 0 discrepanze, 8 rilievi non bloccanti (5 accolti, 4 in BACKLOG BL-020…BL-023).

### Documenti
- `reports/R1.qmd` + `R1.html`; `docs/sessions/2026-09-18_R1.md`; ROADMAP §1, §3, §6, §7; BACKLOG BL-014…BL-023; BIO_REFERENCES (stato); `results/R1/R1_cp2_histology.csv` (verifica istologica), `results/R1/R1_adversarial_review.csv`.

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
