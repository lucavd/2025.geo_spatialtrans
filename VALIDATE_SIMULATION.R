#!/usr/bin/env Rscript
# Script per generare una simulazione a dimensione ridotta e creare plot di validazione

# 1. Caricamento delle funzioni
cat("Caricamento delle funzioni...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# 2. Caricamento librerie necessarie
cat("Caricamento librerie...\n")
library(png)
library(ggplot2)
library(dplyr)
library(Matrix)
library(ClusterR)
library(sp)
library(gstat)
library(tictoc)

# Configurazione per la memoria
options(future.globals.maxSize = 4 * 1024^2) # 4 GB
options(future.rng.onMisuse = "ignore")

# 3. Parametri di simulazione (dimensione ridotta)
tic("Configurazione")
pixel_size_um <- 3000 / 300  # 300 px = 3 mm ⇒ 10 µm/px
cfg <- initialize_simulation_config(
  image_path          = "images/synthetic_tissue.png", # Immagine più semplice
  output_path         = "results/test_simulation.rds",
  output_plot         = "results/test_simulation_plot.png",
  n_genes             = 50,                           # Ridotto per test
  k_cell_types        = 4,                            # Ridotto per test
  threshold_value     = 0.7,
  random_seed         = 42,
  pixel_size_um       = pixel_size_um,
  grid_mode           = TRUE,
  grid_resolution     = 50,  # Risoluzione MOLTO più bassa per test veloce
  grid_spacing        = 0,
  use_fixed_grid      = FALSE                         # Griglia adattiva anziché fissa
)
diff_cfg <- configure_difficulty_level("easy")        # Difficoltà ridotta
toc()

# 4. Preparazione immagine e clustering
tic("Preparazione immagine")
img_dat <- prepare_image(cfg$image_path, cfg$threshold_value)
cat("Dimensioni immagine:", img_dat$width, "x", img_dat$height, "pixel\n")
clust <- cluster_image(
  img_df_thresh = img_dat$img_df_thresh,
  k_cell_types  = cfg$k_cell_types,
  random_seed   = cfg$random_seed
)
toc()

# 5. Creazione griglia di campionamento
tic("Creazione griglia")
cell_df <- create_sampling_grid(
  img_df_thresh        = clust,
  img_array            = img_dat$img_array,
  img_width            = img_dat$width,
  img_height           = img_dat$height,
  grid_mode            = cfg$grid_mode,
  n_cells              = 500,                         # Limitato per test MOLTO veloce
  grid_resolution      = cfg$grid_resolution,
  grid_spacing         = cfg$grid_spacing,
  use_fixed_grid       = cfg$use_fixed_grid,
  fixed_grid_width_mm  = cfg$fixed_grid_width_mm,
  fixed_grid_height_mm = cfg$fixed_grid_height_mm,
  pixel_size_um        = cfg$pixel_size_um,
  threshold_value      = cfg$threshold_value,
  random_seed          = cfg$random_seed
)
cat("Punti generati:", nrow(cell_df), "\n")
toc()

# 6. Generazione chunked dei profili di espressione
tic("Generazione espressione")
# Dimensione del chunk - adattala in base alla memoria disponibile
chunk_size <- 200  # Dimensione ridotta per evitare problemi di memoria
spots <- seq_len(nrow(cell_df))
chunks <- split(spots, ceiling(seq_along(spots) / chunk_size))
cat("Divisione in", length(chunks), "chunks\n")

# Prepara lista per risultati a matrice sparsa
sparse_list <- vector("list", length(chunks))

# Elabora ogni chunk
for (i in seq_along(chunks)) {
  cat(sprintf("Elaborazione chunk %d/%d...\n", i, length(chunks)))
  idx <- chunks[[i]]
  sub_df <- cell_df[idx, ]
  
  # Forza garbage collection prima di ogni chunk
  gc()
  
  # Genera espressione per questo subset di celle
  expr_r <- generate_expression_profiles(
    cell_df = sub_df,
    n_genes = cfg$n_genes,
    k_cell_types = cfg$k_cell_types,
    marker_params = diff_cfg$marker_params,
    spatial_params = diff_cfg$spatial_params,
    dropout_params = diff_cfg$dropout_params,
    cell_specific_params = diff_cfg$cell_specific_params,
    use_spatial_correlation = FALSE,  # Disattivato per velocità
    correlation_method = "grf",
    random_seed = cfg$random_seed + i  # Diverso per ogni chunk
  )
  
  # Converti in matrice sparsa e salva
  # Verifica che la matrice sia orientata come geni x celle
  # nel chunk i le celle sono nel range chunks[[i]]
  if (nrow(expr_r$expression) == cfg$n_genes) {
    # Matrice già nel formato corretto (geni x celle)
    sparse_list[[i]] <- Matrix(expr_r$expression, sparse = TRUE)
  } else if (ncol(expr_r$expression) == cfg$n_genes) {
    # Trasponi per avere geni x celle
    cat("Trasposizione della matrice chunk", i, "da celle x geni a geni x celle...\n")
    sparse_list[[i]] <- Matrix(t(expr_r$expression), sparse = TRUE)
  } else {
    stop(paste("Dimensioni incompatibili nel chunk", i, ":", 
               dim(expr_r$expression)[1], "x", dim(expr_r$expression)[2],
               ", atteso", cfg$n_genes, "geni"))
  }
  
  # Forza garbage collection dopo ogni chunk
  gc()
}

# Ora tutte le matrici hanno lo stesso formato: geni x celle
# Dobbiamo combinarle per colonna, non per riga
# Controlla dimensioni
cat("Dimensioni di alcune matrici sparse:\n")
for (i in 1:min(3, length(sparse_list))) {
  cat("Chunk", i, "dimensioni:", dim(sparse_list[[i]]), "\n")
}

# Crea una matrice per tutti i geni
n_genes <- nrow(sparse_list[[1]])
n_cells <- ncol(sparse_list[[1]]) * length(sparse_list)

# Verifica i rownames
if (is.null(rownames(sparse_list[[1]]))) {
  # Crea nomi generici
  rownames_genes <- paste0("gene_", 1:n_genes)
} else {
  rownames_genes <- rownames(sparse_list[[1]])
}

# Prepara una nuova matrice sparsa per tutti i geni su tutte le celle
# Nota: cbind per matrici sparse può richiedere molta memoria
cat("Combinazione di", length(sparse_list), "matrici sparse (geni x celle)...\n")

# Crea una matrice vuota con le dimensioni finali
full_expr <- Matrix(0, nrow = n_genes, ncol = nrow(cell_df), sparse = TRUE)
rownames(full_expr) <- rownames_genes

# Ora riempi la matrice chunk per chunk
col_offset <- 0
for (i in seq_along(sparse_list)) {
  chunk_size <- ncol(sparse_list[[i]])
  col_indices <- col_offset + (1:chunk_size)
  
  # Assegna il chunk alla matrice completa
  full_expr[, col_indices] <- sparse_list[[i]]
  
  # Aggiorna l'offset per il prossimo chunk
  col_offset <- col_offset + chunk_size
}

# Imposta i rownames
rownames(full_expr) <- rownames_genes

cat("Dimensione matrice espressione finale:", dim(full_expr), "\n")
toc()

# 7. Salvataggio dei risultati
tic("Salvataggio risultati")
# Conversione della matrice sparsa in matrice densa per i plot di validazione
cat("Conversione matrice da sparsa a densa...\n")
dense_expr <- as.matrix(full_expr)
cat("Dimensione matrice densa:", dim(dense_expr), "\n")

# Verifica dell'orientazione della matrice e correzione se necessario
if (ncol(dense_expr) == nrow(cell_df)) {
  # La matrice è già orientata correttamente (geni x celle)
  cat("Matrice già nel formato corretto (geni x celle)\n")
} else if (nrow(dense_expr) == nrow(cell_df)) {
  # La matrice è orientata come celle x geni, dobbiamo trasporla
  cat("Trasposizione della matrice da celle x geni a geni x celle...\n")
  dense_expr <- t(dense_expr)
  cat("Nuova dimensione matrice dopo trasposizione:", dim(dense_expr), "\n")
} else {
  cat("ATTENZIONE: Dimensioni incongruenti tra matrice espressione e celle\n")
  cat("Dimensione matrice:", dim(dense_expr), "\n")
  cat("Numero celle:", nrow(cell_df), "\n")
}

risultato <- list(
  expression = dense_expr,  # Usiamo la matrice densa per compatibilità
  coordinates = cell_df[, c("x", "y")],
  intensity_cluster = factor(cell_df$intensity_cluster),  # Conversione esplicita a factor
  parameters = list(
    marker_params = diff_cfg$marker_params,
    spatial_params = diff_cfg$spatial_params,
    dropout_params = diff_cfg$dropout_params,
    cell_specific_params = diff_cfg$cell_specific_params
  )
)
dir.create(dirname(cfg$output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(risultato, cfg$output_path)
toc()

# 8. Visualizzazione e salvataggio del plot
tic("Generazione plot")
generate_and_save_plots(
  cell_df = cell_df,
  config = cfg,
  difficulty_config = diff_cfg,
  output_plot = cfg$output_plot
)
toc()

# 9. Creazione dei plot di validazione
tic("Plot di validazione")
dir.create("plots/validazione", recursive = TRUE, showWarnings = FALSE)

# Utilizza il modulo di validazione creato
plots <- tryCatch({
  # Per evitare problemi con umap/Rtsne non installati
  cat("Generazione plot di validazione...\n")
  generate_validation_plots(
    sim_results = risultato,
    marker_genes = NULL,  # Rileva automaticamente
    n_markers = 3,
    output_dir = "plots/validazione",
    file_prefix = "test_sim",
    file_format = "png",
    width = 8,
    height = 6,
    skip_dim_reduction = TRUE  # Salta UMAP/tSNE se non disponibili
  )
  cat("Plot di validazione generati in plots/validazione/\n")
}, error = function(e) {
  cat("Errore nella generazione dei plot di validazione:", conditionMessage(e), "\n")
  NULL
})
toc()

cat("\nSimulazione completata con successo!\n")
cat("Risultati salvati in:", cfg$output_path, "\n")
cat("Plot salvati in:", cfg$output_plot, "e nella cartella plots/validazione/\n")