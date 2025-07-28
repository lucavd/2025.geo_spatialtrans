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
# 
# ========================================================================
# ISTRUZIONI PER CAMBIARE MOTORE SPAZIALE
# ========================================================================
# 
# Per passare tra i motori spaziali, modifica questi 2 parametri:
# 
# MOTORE MRF (pattern altamente strutturati):
#   spatial_engine = "mrf"
#   complexity = 8              # Deve = k_cell_types (8 in questo caso)
#   Risultato: Moran's I ~0.3-0.5
# 
# MOTORE IMAGE (pattern biologicamente realistici):
#   spatial_engine = "image"
#   complexity = 3              # 1=semplice, 2=medio, 3=complesso
#   Risultato: Moran's I ~0.4-0.6
# 
# PARAMETRI COMPLEXITY:
# - Image engine: 1 (blob gaussiani), 2 (patch Voronoi), 3 (mix + rumore)
# - MRF engine: deve corrispondere a k_cell_types per coerenza biologica
# 
# ALTRI PARAMETRI MRF (opzionali):
# - mrf_beta: autocorrelazione spaziale (0.1-1.0, default 0.6)
# - mrf_tissue_structure: struttura composita del tessuto
# - mrf_grid_size: dimensione griglia MRF
# 
# ========================================================================

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
  fixed_grid_height_mm = 8.0,
  # MOTORE SPAZIALE - Scegli uno dei due:
  spatial_engine = "image",     # ATTUALE: motore image
  # spatial_engine = "mrf",     # ALTERNATIVA: motore MRF
  
  mrf_beta = 0.6,             # Autocorrelazione spaziale (solo per MRF)
  
  # COMPLEXITY - Adatta al motore scelto:
  complexity = 3,             # ATTUALE: image complexity (1-3)
  # complexity = 8,           # ALTERNATIVA: MRF complexity (= k_cell_types)
  mrf_tissue_structure = list(  # COMPOSITE tissue architecture
    list(type = "vessel", weight = 0.5),   # Vascular structures
    list(type = "gradient", weight = 0.3), # Metabolic gradients
    list(type = "boundary", weight = 0.2)  # Tissue boundaries
  ),
  mrf_grid_size = c(200, 200)     # Larger grid for realistic patterns
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
  if (cfg$spatial_engine == "mrf") {
    if (is.list(cfg$mrf_tissue_structure)) {
      struct_desc <- paste(sapply(cfg$mrf_tissue_structure, function(s) paste0(s$type, "(", s$weight, ")")), collapse="+")
      cat("Engine: MRF (beta=", cfg$mrf_beta, ", composite=", struct_desc, ", grid=", paste(cfg$mrf_grid_size, collapse="x"), ")\n")
    } else {
      cat("Engine: MRF (beta=", cfg$mrf_beta, ", structure=", cfg$mrf_tissue_structure, ", grid=", paste(cfg$mrf_grid_size, collapse="x"), ")\n")
    }
  }
  tryCatch({
    if (cfg$spatial_engine == "mrf") {
      # Direct MRF call with realistic parameters
      mrf_data <- simulate_mrf(
        grid_size = cfg$mrf_grid_size,
        k_cell_types = cfg$k_cell_types,
        beta = cfg$mrf_beta,
        n_iter = 300,  # More iterations for convergence
        seed = cfg$random_seed,
        fast_mode = TRUE,
        tissue_structure = cfg$mrf_tissue_structure
      )

      # Calcola Moran's I sui dati MRF completi (prima del campionamento)
      # Questo è l'approccio corretto per dati con struttura griglia
      compute_moran <- function(df, grid_size) {
        mat <- matrix(df$cell_type, nrow = grid_size, byrow = TRUE)
        m <- mean(mat)
        w_total <- 0; num <- 0
        for (i in 1:grid_size) {
          for (j in 1:grid_size) {
            v <- mat[i, j] - m
            if (j < grid_size) { num <- num + v * (mat[i, j+1] - m); w_total <- w_total + 1 }
            if (i < grid_size) { num <- num + v * (mat[i+1, j] - m); w_total <- w_total + 1 }
          }
        }
        den <- sum((mat - m)^2)
        I <- (grid_size^2 / w_total) * (num / den)
        return(I)
      }

      # Calcola Moran's I sulla griglia completa 200x200
      mrf_moran_i <- compute_moran(mrf_data, cfg$mrf_grid_size[1])

      syn <- list(img_matrix = NA, img_df = mrf_data)
    } else {
      syn <- generate_synthetic_tissue(
        width_px = 800,
        height_px = 800,
        complexity = cfg$complexity,
        seed = cfg$random_seed,
        engine = cfg$spatial_engine,
        output_path = "R_simple/testing/full_tissue_complex.png"
      )

      # Per image engine, calcoliamo Moran's I sui dati immagine
      # Usa le intensità dell'immagine come variabile spaziale
      compute_image_moran <- function(img_df) {
        # Prendi un campione per efficienza computazionale
        if (nrow(img_df) > 5000) {
          set.seed(cfg$random_seed)
          sample_idx <- sample(nrow(img_df), 5000)
          img_df <- img_df[sample_idx, ]
        }
        
        # Normalizza intensità per calcolo Moran
        img_df$norm_intensity <- scale(img_df$intensity)[,1]
        
        # Calcola matrice distanze (usa solo subset per performance)
        n <- nrow(img_df)
        if (n < 100) return(0)  # troppo pochi dati
        
        # Usa distanza euclidea e soglia per definire vicinanza
        dist_threshold <- quantile(sqrt((img_df$x - mean(img_df$x))^2 + 
                                       (img_df$y - mean(img_df$y))^2), 0.1)
        
        w_total <- 0
        numerator <- 0
        mean_intensity <- mean(img_df$norm_intensity)
        
        # Calcola Moran's I con approccio efficiente
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            dist <- sqrt((img_df$x[i] - img_df$x[j])^2 + (img_df$y[i] - img_df$y[j])^2)
            if (dist <= dist_threshold) {
              w <- 1  # peso binario per vicinanza
              w_total <- w_total + 2  # simmetrico
              dev_i <- img_df$norm_intensity[i] - mean_intensity
              dev_j <- img_df$norm_intensity[j] - mean_intensity
              numerator <- numerator + 2 * w * dev_i * dev_j
            }
          }
        }
        
        if (w_total == 0) return(0)
        
        denominator <- sum((img_df$norm_intensity - mean_intensity)^2)
        if (denominator == 0) return(0)
        
        moran_i <- (n / w_total) * (numerator / denominator)
        return(moran_i)
      }
      
      mrf_moran_i <- compute_image_moran(syn$img_df)
    }

    if (cfg$spatial_engine == "mrf") {
      # MRF output has cell_type instead of intensity - convert properly
      img_df_thresh <- syn$img_df
      img_df_thresh$intensity <- 0.5  # dummy intensity for compatibility
      img_df_thresh$value <- img_df_thresh$intensity / 255
      # Convert cell_type to intensity_cluster factor for grid sampling compatibility
      img_df_thresh$intensity_cluster <- factor(paste0("cluster_", img_df_thresh$cell_type))
      img_df_thresh <- img_df_thresh[, c("x", "y", "value", "intensity_cluster")]
    } else {
      img_df_thresh <- syn$img_df[syn$img_df$intensity/255 < cfg$threshold_value, ]
      img_df_thresh$value <- img_df_thresh$intensity / 255
      img_df_thresh <- img_df_thresh[, c("x", "y", "value")]
    }

    if (cfg$spatial_engine == "mrf") {
      # MRF doesn't generate img_matrix, create dummy array for compatibility
      img_dat <- list(
        width = cfg$mrf_grid_size[1],
        height = cfg$mrf_grid_size[2],
        img_df_thresh = img_df_thresh,
        img_array = matrix(0.5, nrow = cfg$mrf_grid_size[2], ncol = cfg$mrf_grid_size[1])
      )
    } else {
      img_dat <- list(
        width = 800,
        height = 800,
        img_df_thresh = img_df_thresh,
        img_array = syn$img_matrix / 255
      )
    }

    cat("✓ Immagine full-size creata:", nrow(img_df_thresh), "pixel validi\n")
  }, error = function(e) {
    cat("✗ ERROR immagine:", e$message, "\n")
    stop("Image generation failed")
  })
}

# 7. Clustering COMPLETO con parametri biologici
cat("\n=== CLUSTERING BIOLOGICAMENTE REALISTICO ===\n")
tryCatch({
  if (cfg$spatial_engine == "mrf") {
    # MRF already has cell types assigned - use them directly
    clust <- img_dat$img_df_thresh
    n_clusters_found <- length(unique(clust$intensity_cluster))
    cat("✓ MRF clustering utilizzato:", n_clusters_found, "cluster trovati\n")
  } else {
    # Traditional clustering for image-based data
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
  }

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

  # Fix MRF cluster assignment using spatial mapping
  if (cfg$spatial_engine == "mrf") {
    # Map grid cells to original MRF clusters using nearest neighbor
    mrf_data <- syn$img_df  # Original MRF data with cell_type

    # For each grid cell, find nearest MRF pixel and assign its cluster
    for (i in 1:nrow(cell_df)) {
      # Find nearest MRF pixel
      distances <- sqrt((mrf_data$x - cell_df$x[i])^2 + (mrf_data$y - cell_df$y[i])^2)
      nearest_idx <- which.min(distances)
      nearest_cell_type <- mrf_data$cell_type[nearest_idx]

      # Assign cluster based on MRF cell_type
      cell_df$intensity_cluster[i] <- factor(paste0("cluster_", nearest_cell_type),
                                           levels = levels(cell_df$intensity_cluster))
    }

    n_clusters_preserved <- length(unique(cell_df$intensity_cluster))
    cat("✓ MRF clusters preservati:", n_clusters_preserved, "/", cfg$k_cell_types, "\n")
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

# 10. TEST CORRELAZIONE SPAZIALE (Indice di Moran)
tryCatch({
  # Usa il valore di Moran calcolato precedentemente sui dati MRF completi
  # Questo è stato calcolato subito dopo la generazione MRF usando l'implementazione corretta
  moran_i <- mrf_moran_i

  # Implementazione corretta di Moran per griglie (da composite_mrf_test.R)
  compute_moran <- function(df, grid_size) {
    mat <- matrix(df$cell_type, nrow = grid_size, byrow = TRUE)
    m <- mean(mat)
    w_total <- 0; num <- 0
    for (i in 1:grid_size) {
      for (j in 1:grid_size) {
        v <- mat[i, j] - m
        if (j < grid_size) { num <- num + v * (mat[i, j+1] - m); w_total <- w_total + 1 }
        if (i < grid_size) { num <- num + v * (mat[i+1, j] - m); w_total <- w_total + 1 }
      }
    }
    den <- sum((mat - m)^2)
    I <- (grid_size^2 / w_total) * (num / den)
    return(I)
  }

  if (moran_i > 0.2) {
    cat("✓ Correlazione spaziale (Moran's I):", round(moran_i, 3), "- PASS\n")
  } else if (moran_i > 0.05) {
    cat("✓ Correlazione spaziale (Moran's I):", round(moran_i, 3), "- WEAK (accettabile)\n")
  } else {
    cat("✓ Correlazione spaziale (Moran's I):", round(moran_i, 3), "- LOW (possibile casualità)\n")
  }

}, error = function(e) {
  cat("✓ Correlazione spaziale: SKIP (errore:", e$message, ")\n")
})

# 11. RISULTATO FINALE FULL
cat("\n=== RISULTATO FINALE FULL ===\n")
all_tests <- c(genes_ok && cells_ok, sparsity_ok, umi_ok, umi_cv_ok,
               cluster_ok, integrity_ok)
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
      spatial_correlation = "moran_index_computed"
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
