# Strategia per simulazione full‑size Visium HD

Questo documento descrive l'analisi e il piano operativo per simulare un esperimento Visium HD reale
su immagine `granuloma.png` (300 px = 3 mm) con:
- dimensione griglia fissa 6.5 × 6.5 mm
- risoluzione bin 2 × 2 µm
- 10 tipi cellulari
- 1000 geni

## 1. Calcolo dei parametri di scala
- L'immagine `granuloma.png` ha dimensione ~1430×2400 px.
- Scala spaziale: 300 px ↔ 3 mm ⇒ pixel_size_um = 10 µm/px.

## 2. Centrare la griglia fissa
- Vogliamo una griglia di 6.5 mm × 6.5 mm (6500 µm × 6500 µm) centrata sull'area di tessuto.
- In `create_sampling_grid()`, modificare il calcolo di `offset_x` e `offset_y` in:
  ```r
  offset_x <- max((img_width_um  - grid_width_um )/2, 0)
  offset_y <- max((img_height_um - grid_height_um)/2, 0)
  ```
- In questo modo la griglia è sempre centrata, indipendentemente dalle dimensioni relative.

## 3. Disabilitare sampling casuale
- Assicurarsi di eseguire sempre la simulazione con:
  ```r
  grid_mode      = TRUE
  use_fixed_grid = TRUE
  ```
- In questo modo NON verrà mai eseguito il ramo di sampling casuale in `create_sampling_grid()`;
  i punti saranno generati da `expand.grid()` e filtrati solo per threshold.

## 4. Impatto computazionale e chunked processing
- Grid fissa: 6500 µm / 2 µm = 3250 bin per lato ⇒ ~10.6 milioni di spot.
- Simulare ~10 M spot in un'unica matrice densa è proibitivo.
- Soluzione:
  1. Estrarre **solo** `cell_df` (coordinate + cluster) con `run_simulation_pipeline()` fino al
     passaggio **prima** di `generate_expression_profiles()`.
  2. Suddividere `cell_df` in blocchi (chunk) di ~5k–10k spot.
  3. In ciascun chunk chiamare `generate_expression_profiles()`, convertire in matrice sparsa.
  4. Usare `do.call(rbind, sparse_blocks)` per assemblare la matrice finale.

## 5. Parametri di simulazione full‑size
- `image_path`          = "images/granuloma.png"
- `pixel_size_um`       = 10
- `grid_mode`           = TRUE
- `use_fixed_grid`      = TRUE
- `grid_resolution`     = 2
- `grid_spacing`        = 0
- `fixed_grid_width_mm` = 6.5
- `fixed_grid_height_mm`= 6.5
- `n_genes`             = 1000
- `k_cell_types`        = 10
- `difficulty_level`    = "medium"  # o "hard"
- `use_spatial_correlation` = TRUE
- `correlation_method`      = "grf"
- `random_seed`         = 42

## 6. Script R di esecuzione
1. **Caricamento funzioni**:
   ```r
   files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
   for (f in sort(files)) source(f)
   ```
2. **Calcolo pixel_size_um**:
   ```r
   pixel_size_um <- 3e3 / 300  # = 10
   ```
3. **Initialize config e difficulty**:
   ```r
   cfg <- initialize_simulation_config(
     image_path = "images/granuloma.png",
     output_path = "results/visiumHD_full.rds",
     output_plot = "results/visiumHD_full.png",
     n_cells = NULL, # ignorato in grid_mode
     n_genes = 1000,
     k_cell_types = 10,
     threshold_value = 0.7,
     random_seed = 42,
     pixel_size_um = pixel_size_um,
     grid_mode = TRUE,
     grid_resolution = 2,
     grid_spacing = 0,
     use_fixed_grid = TRUE,
     fixed_grid_width_mm = 6.5,
     fixed_grid_height_mm = 6.5
   )
   diff_cfg <- configure_difficulty_level(
     "medium"
   )
   ```
4. **Estrapolazione `cell_df`**:
   ```r
   img_dat <- prepare_image(cfg$image_path, cfg$threshold_value)
   clust  <- cluster_image(img_dat$img_df_thresh, cfg$k_cell_types, cfg$random_seed)
   cell_df <- create_sampling_grid(
     img_df_thresh = clust,
     img_array      = img_dat$img_array,
     img_width      = img_dat$width,
     img_height     = img_dat$height,
     grid_mode      = cfg$grid_mode,
     grid_resolution = cfg$grid_resolution,
     use_fixed_grid = cfg$use_fixed_grid,
     fixed_grid_width_mm = cfg$fixed_grid_width_mm,
     fixed_grid_height_mm= cfg$fixed_grid_height_mm,
     pixel_size_um  = cfg$pixel_size_um,
     threshold_value= cfg$threshold_value,
     random_seed    = cfg$random_seed
   )
   ```
5. **Chunked expression generation**:
   ```r
   chunk_size   <- 5000
   idx_chunks   <- split(seq_len(nrow(cell_df)), ceiling(seq_along(cell_df$x)/chunk_size))
   sparse_list  <- vector("list", length(idx_chunks))
   for (i in seq_along(idx_chunks)) {
     sub_df <- cell_df[idx_chunks[[i]], ]
     expr_r <- generate_expression_profiles(
       cell_df = sub_df,
       n_genes = cfg$n_genes,
       k_cell_types = cfg$k_cell_types,
       marker_params = diff_cfg$marker_params,
       spatial_params= diff_cfg$spatial_params,
       dropout_params= diff_cfg$dropout_params,
       library_size_params = diff_cfg$library_size_params,
       hybrid_params       = diff_cfg$hybrid_params,
       cell_specific_params= diff_cfg$cell_specific_params,
       use_spatial_correlation = TRUE,
       correlation_method = "grf",
       random_seed       = cfg$random_seed
     )
     sparse_list[[i]] <- Matrix::Matrix(expr_r$expression, sparse = TRUE)
   }
   full_expr <- do.call(rbind, sparse_list)
   ```
6. **Salvataggio e plotting**:
   ```r
   risultato <- list(
     expression = full_expr,
     coordinates= as.matrix(cell_df[, c("x","y")]),
     intensity_cluster = cell_df$intensity_cluster
   )
   saveRDS(risultato, cfg$output_path)
   # plot:
   generate_and_save_plots(
     cell_df, cfg, diff_cfg, cfg$output_plot, full_expr
   )
   ```

## 7. Validazione in scala ridotta
- Usare `grid_resolution = 10` µm ⇒ ~422k spot
- `n_genes = 200`, `k_cell_types = 5`
- Controllare end‑to‑end prima di lanciare i ~10M spot.

---
*Fine della strategia completa.*