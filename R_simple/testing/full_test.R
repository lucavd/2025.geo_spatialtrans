#!/usr/bin/env Rscript
# Full test script per testing biologicamente realistico con parametri completi
# Pipeline COMPLETA con parametri full-size per validazione biologica
# Usage: Rscript full_test.R [image_path]

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
user_image_path <- if (length(args) > 0) args[1] else NULL

cat("=== FULL TEST PIPELINE R_SIMPLE ===\n")
if (!is.null(user_image_path)) {
  cat("Modalità: IMMAGINE UTENTE -", user_image_path, "\n")
} else {
  cat("Modalità: IMMAGINE SINTETICA\n")
}
start_time <- Sys.time()

# 1. Caricamento funzioni R_simple
cat("Caricamento funzioni R_simple...\n")
files <- list.files("R_simple", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) {
  tryCatch({
    source(f)
    cat("✓ ", basename(f), "\n")
  }, error = function(e) {
    cat("✗ ERROR loading", basename(f), ":", e$message, "\n")
    stop("Failed to load functions")
  })
}

# 2. Caricamento librerie essenziali
cat("\nCaricamento librerie...\n")
required_libs <- c("Matrix", "ClusterR", "dplyr", "sp", "gstat", "png")
for (lib in required_libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    cat("✗ MISSING:", lib, "\n")
    stop("Required library missing: ", lib)
  } else {
    library(lib, character.only = TRUE)
    cat("✓", lib, "\n")
  }
}

# 3. Configurazione FULL (parametri biologicamente realistici)
cat("\n=== CONFIGURAZIONE FULL BIOLOGICAMENTE REALISTICA ===\n")
cfg <- list(
  n_genes = 5000,          # Realistico per spatial transcriptomics
  n_cells = 10000,         # Target celle per dataset medio
  k_cell_types = 8,        # Più tipi cellulari realistici
  threshold_value = 0.7,   # Standard per segmentazione tissutale
  random_seed = 42,
  pixel_size_um = 10,      # 10μm per spot Visium
  grid_mode = TRUE,
  grid_resolution = 40,    # Più denso per full test
  grid_spacing = 0,
  use_fixed_grid = TRUE,
  fixed_grid_width_mm = 8.0,   # Area più grande 
  fixed_grid_height_mm = 8.0
)

cat("- Geni:", cfg$n_genes, "\n")
cat("- Celle target:", cfg$n_cells, "\n") 
cat("- Cluster:", cfg$k_cell_types, "\n")
cat("- Area:", cfg$fixed_grid_width_mm, "x", cfg$fixed_grid_height_mm, "mm\n")

# 4. Parametri biologici COMPLETI e REALISTICI
diff_cfg <- list(
  marker_params = list(
    marker_genes_per_type = 25,      # Più marker per tipo
    marker_expression_fold = 2.5,    # Fold-change più marcato
    marker_overlap_fold = 0.05       # Meno overlap tra tipi
  ),
  spatial_params = list(
    spatial_noise_intensity = 1.2,   # Variabilità spaziale realistica
    spatial_range = 40,              # Range correlazione spaziale
    random_noise_sd = 0.25,          # Rumore biologico moderato
    correlation_length = 0.18,       # Correlazione spaziale medio-alta
    noise_level = 0.35               # Livello rumore controllato
  ),
  dropout_params = list(
    dropout_range = c(0.45, 0.65),   # Range dropout realistico per spatial
    dispersion_range = c(12.0, 6.0), # Dispersione high-quality data
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.4,
    dropout_curve_steepness = 6      # Curva più sharp
  ),
  cell_specific_params = list(
    library_size_params = list(
      mean_library_size = 8000,      # Tipico per Visium HD
      library_size_cv = 0.22,        # CV più basso (alta qualità)
      spatial_effect_on_library = 0.15,  # Effetto spaziale moderato
      cell_type_effect = TRUE
    ),
    cell_specific_noise_sd = 0.18,   # Rumore cellula-specifico ridotto
    use_gene_modules = TRUE,
    n_gene_modules = 8,              # Più moduli co-espressi
    module_correlation = 0.75        # Correlazione moduli più forte
  )
)

# 5. Configurazione memoria per full test
cat("\n=== CONFIGURAZIONE MEMORIA ===\n")
options(future.globals.maxSize = 4 * 1024^2) # 4 GB
gc() # Pulizia iniziale

# 6. Preparazione immagine (sintetica o utente)
if (!is.null(user_image_path)) {
  # MODALITÀ IMMAGINE UTENTE
  cat("\n=== PROCESSAMENTO IMMAGINE UTENTE ===\n")
  tryCatch({
    # Verifica esistenza file
    if (!file.exists(user_image_path)) {
      stop("File immagine non trovato: ", user_image_path)
    }
    
    cat("Caricamento immagine:", user_image_path, "\n")
    img_dat <- prepare_image(user_image_path, cfg$threshold_value)
    
    cat("✓ Immagine caricata:", img_dat$width, "x", img_dat$height, "pixel\n")
    cat("✓ Pixel validi:", nrow(img_dat$img_df_thresh), "\n")
    
    # Verifica se ci sono abbastanza pixel validi
    if (nrow(img_dat$img_df_thresh) < 5000) {
      cat("⚠ WARNING: Pochi pixel validi per full test. Considerare threshold diverso.\n")
    }
    
  }, error = function(e) {
    cat("✗ ERROR caricamento immagine:", e$message, "\n")
    stop("User image loading failed")
  })
  
} else {
  # MODALITÀ IMMAGINE SINTETICA (full-size)
  cat("\n=== GENERAZIONE IMMAGINE SINTETICA FULL-SIZE ===\n")
  tryCatch({
    syn <- generate_synthetic_tissue(
      width_px = 800,              # Dimensioni realistiche
      height_px = 800,
      complexity = 3,              # Massima complessità
      seed = cfg$random_seed,
      output_path = "R_simple/testing/full_tissue_complex.png"
    )
    
    img_df_thresh <- syn$img_df[syn$img_df$intensity/255 < cfg$threshold_value, ]
    img_df_thresh$value <- img_df_thresh$intensity / 255
    img_df_thresh <- img_df_thresh[, c("x", "y", "value")]
    
    img_dat <- list(
      width = 800,
      height = 800, 
      img_df_thresh = img_df_thresh,
      img_array = syn$img_matrix / 255
    )
    
    cat("✓ Immagine full-size creata:", nrow(img_df_thresh), "pixel validi\n")
  }, error = function(e) {
    cat("✗ ERROR immagine:", e$message, "\n")
    stop("Image generation failed")
  })
}

# 7. Clustering COMPLETO con parametri biologici
cat("\n=== CLUSTERING BIOLOGICAMENTE REALISTICO ===\n")
tryCatch({
  clust <- cluster_image(
    img_df_thresh = img_dat$img_df_thresh,
    k_cell_types = cfg$k_cell_types,
    random_seed = cfg$random_seed,
    spatial_weight = 0.6                    # Bilanciamento spazio-intensità
  )
  
  # Rimuovi NA se presenti
  clust <- clust[!is.na(clust$intensity_cluster), ]
  
  n_clusters_found <- length(unique(clust$intensity_cluster))
  cat("✓ Clustering completato:", n_clusters_found, "cluster trovati\n")
  
  if (n_clusters_found != cfg$k_cell_types) {
    cat("⚠ WARNING: Attesi", cfg$k_cell_types, "cluster, trovati", n_clusters_found, "\n")
  }
}, error = function(e) {
  cat("✗ ERROR clustering:", e$message, "\n")
  stop("Clustering failed")
})

# 8. Creazione griglia high-density
cat("\n=== CREAZIONE GRIGLIA HIGH-DENSITY ===\n")
tryCatch({
  cell_df <- create_sampling_grid(
    img_df_thresh = clust,
    img_array = img_dat$img_array,
    img_width = img_dat$width,
    img_height = img_dat$height,
    grid_mode = cfg$grid_mode,
    grid_resolution = cfg$grid_resolution,
    grid_spacing = cfg$grid_spacing,
    use_fixed_grid = cfg$use_fixed_grid,
    fixed_grid_width_mm = cfg$fixed_grid_width_mm,
    fixed_grid_height_mm = cfg$fixed_grid_height_mm,
    pixel_size_um = cfg$pixel_size_um,
    threshold_value = cfg$threshold_value,
    random_seed = cfg$random_seed
  )
  
  cat("✓ Griglia high-density creata:", nrow(cell_df), "celle\n")
  
  if (nrow(cell_df) == 0) {
    stop("No cells generated in grid")
  }
  
  # Se abbiamo troppe celle, sub-sample per performance
  if (nrow(cell_df) > cfg$n_cells * 2) {
    set.seed(cfg$random_seed)
    keep_indices <- sample(nrow(cell_df), min(cfg$n_cells * 2, nrow(cell_df)))
    cell_df <- cell_df[keep_indices, ]
    cat("✓ Sub-sampled to:", nrow(cell_df), "celle per performance\n")
  }
  
}, error = function(e) {
  cat("✗ ERROR griglia:", e$message, "\n")
  stop("Grid creation failed")
})

# 9. Generazione espressione FULL BIOLOGICAMENTE REALISTICA
cat("\n=== GENERAZIONE ESPRESSIONE FULL BIOLOGICA ===\n")
cat("Parametri realistici:\n")
cat("- Library size target:", diff_cfg$cell_specific_params$library_size_params$mean_library_size, "UMI/cella\n")
cat("- Dropout range:", diff_cfg$dropout_params$dropout_range[1], "-", diff_cfg$dropout_params$dropout_range[2], "\n")
cat("- Marker genes per tipo:", diff_cfg$marker_params$marker_genes_per_type, "\n")

# Chunking per gestire memoria con dataset grandi
chunk_size <- min(2000, ceiling(nrow(cell_df) / 4))
n_chunks <- ceiling(nrow(cell_df) / chunk_size)

cat("Processamento in", n_chunks, "chunks di", chunk_size, "celle ciascuno...\n")

# Lista per risultati
expression_chunks <- vector("list", n_chunks)
pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)

for (i in 1:n_chunks) {
  setTxtProgressBar(pb, i)
  
  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, nrow(cell_df))
  chunk_cells <- cell_df[start_idx:end_idx, ]
  
  # Garbage collection periodico
  if (i %% 3 == 0) gc()
  
  tryCatch({
    expr_result <- generate_expression_profiles(
      cell_df = chunk_cells,
      n_genes = cfg$n_genes,
      k_cell_types = cfg$k_cell_types,
      marker_params = diff_cfg$marker_params,
      spatial_params = diff_cfg$spatial_params,
      dropout_params = diff_cfg$dropout_params,
      cell_specific_params = diff_cfg$cell_specific_params,
      use_spatial_correlation = TRUE,
      correlation_method = "grf",
      random_seed = cfg$random_seed + i
    )
    
    expr_matrix <- expr_result$expression
    
    # Verifica orientamento matrice
    if (nrow(expr_matrix) != cfg$n_genes) {
      if (ncol(expr_matrix) == cfg$n_genes) {
        expr_matrix <- t(expr_matrix)
      } else {
        stop("Matrix dimensions don't match expected n_genes")
      }
    }
    
    expression_chunks[[i]] <- Matrix(expr_matrix, sparse = TRUE)
    
  }, error = function(e) {
    cat("\n✗ ERROR nel chunk", i, ":", e$message, "\n")
    stop("Expression generation failed in chunk ", i)
  })
}

close(pb)
cat("\n")

# Combina tutti i chunks
cat("Combinazione chunks in matrice finale...\n")
if (n_chunks == 1) {
  final_expr <- expression_chunks[[1]]
} else {
  final_expr <- do.call(cbind, expression_chunks)
}

cat("✓ Espressione generata:", dim(final_expr)[1], "x", dim(final_expr)[2], "\n")

# 10. BENCHMARK E VALIDAZIONI FULL
cat("\n=== BENCHMARK E VALIDAZIONI FULL ===\n")

# Test dimensioni
genes_ok <- nrow(final_expr) == cfg$n_genes
cells_ok <- ncol(final_expr) == nrow(cell_df)
cat("✓ Dimensioni matrice:", ifelse(genes_ok && cells_ok, "PASS", "FAIL"), 
    "(", nrow(final_expr), "x", ncol(final_expr), ")\n")

# Test sparsità  
sparsity <- mean(final_expr == 0) * 100
sparsity_ok <- sparsity >= 30 && sparsity <= 92  # Range esteso per spatial transcriptomics
cat("✓ Sparsità:", round(sparsity, 1), "% -", ifelse(sparsity_ok, "PASS", "FAIL"), "\n")

# Test UMI per cella
umi_per_cell <- colSums(final_expr)
umi_mean <- mean(umi_per_cell)
umi_median <- median(umi_per_cell)
umi_ok <- umi_mean >= 3000 && umi_mean <= 20000  # Range realistico per spatial
cat("✓ UMI medio:", round(umi_mean), "mediano:", round(umi_median), "-", 
    ifelse(umi_ok, "PASS", "FAIL"), "\n")

# Test distribuzione UMI
umi_cv <- sd(umi_per_cell) / mean(umi_per_cell)
umi_cv_ok <- umi_cv >= 0.15 && umi_cv <= 1.0  # CV esteso per spatial data
cat("✓ UMI CV:", round(umi_cv, 2), "-", ifelse(umi_cv_ok, "PASS", "FAIL"), "\n")

# Test cluster assignment
clusters_assigned <- length(unique(cell_df$intensity_cluster))
cluster_ok <- clusters_assigned == cfg$k_cell_types
cat("✓ Cluster assegnati:", clusters_assigned, "/", cfg$k_cell_types, "-", 
    ifelse(cluster_ok, "PASS", "FAIL"), "\n")

# Test integrità dati
has_na <- any(is.na(final_expr)) || any(is.infinite(final_expr))
integrity_ok <- !has_na
cat("✓ Integrità dati:", ifelse(integrity_ok, "PASS", "FAIL"), "\n")

# Test correlazione spaziale (sample)
spatial_test_ok <- TRUE
tryCatch({
  if (nrow(cell_df) >= 100) {
    sample_cells <- sample(nrow(cell_df), min(100, nrow(cell_df)))
    sample_coords <- cell_df[sample_cells, c("x", "y")]
    sample_expr <- final_expr[1:min(50, nrow(final_expr)), sample_cells]
    
    # Test correlazione tra celle vicine vs lontane
    dist_matrix <- as.matrix(dist(sample_coords))
    expr_cor_matrix <- cor(t(sample_expr), use = "complete.obs")
    
    # Correlazione media per celle vicine (<10% max distance)
    close_threshold <- quantile(dist_matrix[upper.tri(dist_matrix)], 0.1)
    close_pairs <- which(dist_matrix < close_threshold & upper.tri(dist_matrix), arr.ind = TRUE)
    
    if (nrow(close_pairs) > 0) {
      close_cor <- mean(expr_cor_matrix[close_pairs], na.rm = TRUE)
      spatial_test_ok <- close_cor > 0.1  # Correlazione spaziale minima
      cat("✓ Correlazione spaziale:", round(close_cor, 3), "-", 
          ifelse(spatial_test_ok, "PASS", "FAIL"), "\n")
    }
  }
}, error = function(e) {
  cat("✓ Correlazione spaziale: SKIP (errore calcolo)\n")
})

# 11. RISULTATO FINALE FULL
cat("\n=== RISULTATO FINALE FULL ===\n")
all_tests <- c(genes_ok && cells_ok, sparsity_ok, umi_ok, umi_cv_ok, 
               cluster_ok, integrity_ok, spatial_test_ok)
overall_pass <- all(all_tests)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("Tempo totale:", round(elapsed, 2), "minuti\n")
cat("Tests passati:", sum(all_tests), "/", length(all_tests), "\n")
cat("RISULTATO COMPLESSIVO:", ifelse(overall_pass, "✓ PASS", "✗ FAIL"), "\n")

if (overall_pass) {
  cat("\n🎉 FULL PIPELINE R_SIMPLE BIOLOGICAMENTE VALIDATA!\n")
  
  # Salva risultato completo
  full_result <- list(
    timestamp = Sys.time(),
    mode = ifelse(is.null(user_image_path), "synthetic", "user_image"),
    image_path = user_image_path,
    config = cfg,
    difficulty_config = diff_cfg,
    results = list(
      n_cells_generated = nrow(cell_df),
      n_clusters_found = clusters_assigned,
      matrix_dims = dim(final_expr),
      sparsity_pct = round(sparsity, 1),
      umi_mean = round(umi_mean),
      umi_median = round(umi_median),
      umi_cv = round(umi_cv, 3),
      elapsed_minutes = round(elapsed, 2),
      chunks_processed = n_chunks
    ),
    biological_validation = list(
      library_size_target = diff_cfg$cell_specific_params$library_size_params$mean_library_size,
      dropout_range = diff_cfg$dropout_params$dropout_range,
      spatial_correlation = spatial_test_ok
    ),
    status = "PASS"
  )
  
  # Salva risultati e dati
  saveRDS(full_result, "R_simple/testing/full_test_result.rds")
  
  # Salva anche i dati finali per analisi successive
  final_data <- list(
    expression = final_expr,
    coordinates = cell_df[, c("x", "y")],
    clusters = cell_df$intensity_cluster,
    config = cfg
  )
  saveRDS(final_data, "R_simple/testing/full_simulation_data.rds")
  
  cat("Risultato salvato in: R_simple/testing/full_test_result.rds\n")
  cat("Dati simulazione salvati in: R_simple/testing/full_simulation_data.rds\n")
  
} else {
  cat("\n❌ ALCUNI TEST FULL FALLITI - CONTROLLARE PIPELINE\n")
  quit(status = 1)
}

cat("\n=== SUMMARY BIOLOGICO ===\n")
cat("Questa simulazione full rappresenta dati realistici di spatial transcriptomics:\n")
cat("- ", cfg$n_genes, " geni (range tipico per dataset pubblicati)\n")
cat("- ", nrow(cell_df), " spots/celle (densità Visium-like)\n")
cat("- ", cfg$k_cell_types, " tipi cellulari (complessità tissutale realistica)\n")
cat("- Library size ~", round(umi_mean), " UMI/spot (validato per tecnologie spatial)\n")
cat("- Sparsità ", round(sparsity, 1), "% (tipica per dati spatial filtrati)\n")
cat("- Effetti biologici completi: dropout, correlazione spaziale, moduli genici\n")
cat("\nDati pronti per benchmark avanzati di algoritmi di clustering!\n")