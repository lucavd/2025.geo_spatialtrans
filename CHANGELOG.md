# CHANGELOG — geo_spatialtrans

Una voce per sessione (id sessione nel titolo). Lo stato di ogni sessione è in `docs/ROADMAP.md`; le idee rinviate in `docs/BACKLOG.md`.

## R2 — 2026-09-20/21 — Segmentazione nucleare su H&E e BIO_REFERENCES v1 (geometria) (branch `step1-cell-layer`)

### Ambiente
- Due venv Python (uv): `.venv` += torch 2.14 (cu130) + Cellpose 4.2.1.1 (`tools/requirements-r2-cellpose.txt` → `requirements-cellpose.lock`); `.venv-stardist` = TensorFlow 2.21 + StarDist 0.9.2 (`requirements-r2-stardist.txt` → `requirements-stardist.lock`). `tools/py_stardist.sh` mette le librerie CUDA del venv nel loader path (senza, TensorFlow non vede la GPU). `tools/setup_python_env.sh` esteso (idempotente).
- Libreria R: aggiunti `diptest`, `arrow`, `spatstat.random` (`tools/setup_r_library.R`).

### Codice (tutto in `tools/`, dati derivati fuori git su `/mnt/micron/geo_spatialtrans/R2/`)
- `R2_overview.py` (panoramiche H&E 1/16 + cluster graphclust 8 µm, griglia mm), `R2_rois.py` (30 ROI 1×1 mm da `results/R2/R2_rois.csv`, check C-R2.4, ritagli nativi, miniature), `R2_common.py` (maschera tessuto, maschera bolle BL-017, canale ematossilina, tabella nuclei, cache maschere), `R2_segment_cellpose.py` (Cellpose-SAM `cpsam_v2`, varianti rgb/hed, scale 1/2/4, riprendibile), `R2_segment_stardist.py` (2D_versatile_he, riprendibile), `R2_spaceranger_rois.py` (rasterizzazione poligoni `nucleus_segmentations.geojson` nei ROI + C-R2.5), `R2_metrics.py` (sommari, appaiamento F1@IoU0.5, effetto scala, riproducibilità), `R2_exclusive_check.py` (arbitro: OD ematossilina degli oggetti esclusivi), `R2_consensus.py` (densità di consenso), `R2_export_masks.py`, `R2_spatial.R` (g(r) con envelope CSR, dip test, Clark–Evans), `R2_figures.py` (9 figure), `R2_variant_eval.py` (rivalutazione varianti Cellpose). `R/testing/test_R2.sh` (un comando, PASS/FAIL).
- Risultati tracciati: `results/R2/*.csv`, `results/R2/spatial/`, `results/R2/figures/`, `results/R2/rois_thumbs/`, `results/R2/variant_eval/`, `results/R2/R2_preregistration.md`. Fuori git: ritagli, label image, `R2_nuclei_all.parquet` (in `results/R2/nuclei/` ignorato se > 50 MB — vedi `.gitignore`).

### Esiti
- Check C: 4 PASS, 3 WARN (C-R2.5 criterio pre-registrato mal posto, orientamento confermato dal criterio indipendente; C-R2.6 maschera bolle parziale; C-R2.7 14 % di frammenti in A5_r1). Riproducibilità bit-identica.
- Check B: B-R2.1, B-R2.2, B-R2.3, B-R2.4 **FAIL** (attese non modificate: A4 ≈ 24 500/mm² vs 1 900–6 300; area linfociti 14 µm² vs 20–45; A5 ≈ 1 300 vs 1 900–3 500; ordinamento con A3 > A5); B-R2.5 PASS; B-R2.6 WARN.
- Controprova: i due segmentatori disaccordano del 12–55 % sulla densità con direzione tessuto-dipendente (CP-R2.1 WARN); risoluzione nativa necessaria (CP-R2.3, BL-015 → BL-029); inibizione a corto raggio in A2–A5 ma non in A1/A6 (CP-R2.4); A3 distinguibile da A2 (CP-R2.5 PASS).
- Deviazione metodologica documentata: variante Cellpose-ematossilina prima liquidata come «tassellatura», poi rivalutata su richiesta di Luca e riclassificata come segmentazione nucleo+alone (BL-028).
- `docs/BIO_REFERENCES.md` v1: B-001…B-042 (6 archetipi × 7 metriche), tutte **provvisorie** (BL-024). BACKLOG BL-024…BL-031; chiusi BL-015, BL-020.

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
