#!/usr/bin/env Rscript
# Script per simulazione full‑size Visium HD su granuloma.png

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
library(ClusterR)  # Fast kmeans++ implementation (provides KMeans_rcpp)
library(sp)
library(gstat)
library(tictoc)

# Configurazione per la memoria
options(future.globals.maxSize = 100 * 1024^2)  # 100 GB per simulazioni full-size
options(future.rng.onMisuse = "ignore")

## 3. Parametri di simulazione
cat("Configurazione simulazione...\n")
pixel_size_um <- 3000 / 300  # 300 px = 3 mm ⇒ 10 µm/px
cfg <- initialize_simulation_config(
  image_path          = "images/granuloma.png",
  output_path         = "results/visiumHD_full.rds",
  output_plot         = "results/visiumHD_full.png",
  n_genes             = 200,
  k_cell_types        = 10,
  threshold_value     = 0.7,
  random_seed         = 42,
  pixel_size_um       = pixel_size_um,
  grid_mode           = TRUE,
  grid_resolution     = 10,
  grid_spacing        = 0,
  use_fixed_grid      = TRUE,
  fixed_grid_width_mm = 6.5,
  fixed_grid_height_mm= 6.5
)
diff_cfg <- configure_difficulty_level("medium")
validate <- TRUE  # Imposta FALSE per simulazioni molto grandi per evitare problemi di memoria

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

# 5. Creazione griglia fissa centrata
tic("Creazione griglia")
cell_df <- create_sampling_grid(
  img_df_thresh         = clust,
  img_array             = img_dat$img_array,
  img_width             = img_dat$width,
  img_height            = img_dat$height,
  grid_mode             = cfg$grid_mode,
  grid_resolution       = cfg$grid_resolution,
  grid_spacing          = cfg$grid_spacing,
  use_fixed_grid        = cfg$use_fixed_grid,
  fixed_grid_width_mm   = cfg$fixed_grid_width_mm,
  fixed_grid_height_mm  = cfg$fixed_grid_height_mm,
  pixel_size_um         = cfg$pixel_size_um,
  threshold_value       = cfg$threshold_value,
  random_seed           = cfg$random_seed
)
cat("Punti generati:", nrow(cell_df), "\n")
toc()

# 6. Generazione chunked dei profili di espressione
tic("Generazione espressione")
chunk_size <- 5000
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
    cell_df               = sub_df,
    n_genes               = cfg$n_genes,
    k_cell_types          = cfg$k_cell_types,
    marker_params         = diff_cfg$marker_params,
    spatial_params        = diff_cfg$spatial_params,
    dropout_params        = diff_cfg$dropout_params,
    cell_specific_params  = diff_cfg$cell_specific_params,
    use_spatial_correlation = TRUE,
    correlation_method    = "grf",
    random_seed           = cfg$random_seed + i  # Diverso per ogni chunk
  )
  
  # Converti in matrice sparsa e salva
  # Verifica che la matrice sia orientata come geni x celle
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

# Mostra dimensioni di alcuni chunk per verifica
cat("Dimensioni di alcune matrici sparse:\n")
for (i in 1:min(3, length(sparse_list))) {
  cat("Chunk", i, "dimensioni:", dim(sparse_list[[i]]), "\n")
}

# Combina i chunk in una singola matrice
# La direzione di combinazione dipende dall'orientazione delle matrici
if (nrow(sparse_list[[1]]) == cfg$n_genes) {
  # Se la matrice è orientata come geni x celle, combina per colonna
  cat("Combinazione di", length(sparse_list), "matrici sparse (geni x celle)...\n")
  
  # Crea una matrice vuota con le dimensioni finali
  full_expr <- Matrix(0, nrow = cfg$n_genes, ncol = nrow(cell_df), sparse = TRUE)
  
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
} else {
  # Se la matrice è orientata diversamente, usa rbind standard
  cat("Usando rbind per combinare le matrici...\n")
  full_expr <- do.call(rbind, sparse_list)
}
cat("Dimensione matrice espressione finale:", dim(full_expr), "\n")
toc()

# 7. Salvataggio dei risultati
tic("Salvataggio risultati")
# Converti matrice a seconda della necessità di validazione
if (validate) {
  cat("Conversione matrice da sparsa a densa per validazione...\n")
  dense_expr <- as.matrix(full_expr)
  cat("Dimensione matrice densa:", dim(dense_expr), "\n")
  
  # Verifica dell'orientazione della matrice finale
  if (ncol(dense_expr) == nrow(cell_df)) {
    cat("Matrice finale nel formato corretto (geni x celle)\n")
  } else if (nrow(dense_expr) == nrow(cell_df)) {
    cat("Trasposizione della matrice finale...\n")
    dense_expr <- t(dense_expr)
    cat("Nuova dimensione matrice dopo trasposizione:", dim(dense_expr), "\n")
  } else {
    cat("ATTENZIONE: Dimensioni incongruenti nella matrice finale\n")
    cat("Dimensione matrice:", dim(dense_expr), "\n")
    cat("Numero celle:", nrow(cell_df), "\n")
  }
  
  # Assicura nomi di riga corretti per la matrice
  if (is.null(rownames(dense_expr))) {
    cat("Aggiunta nomi di riga mancanti...\n")
    rownames(dense_expr) <- paste0("gene_", seq_len(nrow(dense_expr)))
  }
  
  # Creazione del risultato per validazione
  risultato <- list(
    expression = dense_expr,
    coordinates = cell_df[, c("x", "y")],
    intensity_cluster = factor(cell_df$intensity_cluster, 
                             levels = sort(unique(as.numeric(cell_df$intensity_cluster)))),  # Conversione a factor con livelli corretti
    parameters = list(
      marker_params = diff_cfg$marker_params,
      spatial_params = diff_cfg$spatial_params,
      dropout_params = diff_cfg$dropout_params,
      cell_specific_params = diff_cfg$cell_specific_params
    )
  )
} else {
  # Per simulazioni molto grandi, mantieni la matrice sparsa
  risultato <- list(
    expression = full_expr,
    coordinates = cell_df[, c("x", "y")],
    intensity_cluster = factor(cell_df$intensity_cluster, 
                              levels = sort(unique(as.numeric(cell_df$intensity_cluster)))),
    parameters = list(
      marker_params = diff_cfg$marker_params,
      spatial_params = diff_cfg$spatial_params,
      dropout_params = diff_cfg$dropout_params,
      cell_specific_params = diff_cfg$cell_specific_params
    )
  )
}
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

# 9. Creazione dei plot di validazione (se richiesto)
if (validate) {
  tic("Plot di validazione")
  dir.create("plots/validazione", recursive = TRUE, showWarnings = FALSE)
  
  # Funzione di generazione plot robusta
  generate_fixed_validation_plots <- function(sim_results, output_dir) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1. Plot di base con coordinate spaziali colorate per cluster
    cat("Generazione plot spaziale base...\n")
    plot_df <- data.frame(
      x = sim_results$coordinates$x,
      y = sim_results$coordinates$y,
      cluster = as.factor(sim_results$intensity_cluster)
    )
    
    p_spatial <- ggplot(plot_df, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 0.5, alpha = 0.7) +
      theme_minimal() +
      labs(title = "Distribuzione spaziale dei cluster",
           x = "Coordinata X", 
           y = "Coordinata Y") +
      theme(aspect.ratio = 1)
    
    ggsave(file.path(output_dir, "spatial_clusters.png"), p_spatial, width = 10, height = 8)
    
    # 2. Plot mean-variance (versione semplificata)
    cat("Generazione plot mean-variance...\n")
    gene_means <- rowMeans(sim_results$expression)
    gene_vars <- apply(sim_results$expression, 1, var)
    
    plot_df <- data.frame(
      gene = seq_along(gene_means),
      mean = gene_means,
      variance = gene_vars
    )
    
    p_meanvar <- ggplot(plot_df, aes(x = mean, y = variance)) +
      geom_point(alpha = 0.7) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      scale_x_log10() + 
      scale_y_log10() +
      theme_minimal() +
      labs(title = "Relazione Media-Varianza dell'Espressione Genica", 
           x = "Media dell'Espressione", 
           y = "Varianza dell'Espressione")
    
    ggsave(file.path(output_dir, "mean_variance.png"), p_meanvar, width = 10, height = 8)
    
    # 3. Plot dropout rate
    cat("Generazione plot dropout rate...\n")
    zero_fraction <- rowMeans(sim_results$expression == 0)
    
    plot_df <- data.frame(
      gene = seq_along(gene_means),
      mean_expression = log1p(gene_means),
      zero_fraction = zero_fraction
    )
    
    p_dropout <- ggplot(plot_df, aes(x = mean_expression, y = zero_fraction)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.2) +
      theme_minimal() +
      labs(title = "Relazione Dropout-Espressione Media", 
           x = "Log(Media dell'Espressione + 1)", 
           y = "Frazione di Zeri (Dropout Rate)")
    
    ggsave(file.path(output_dir, "dropout.png"), p_dropout, width = 10, height = 8)
    
    # 4. Plot di alcuni marker genes (primi 4 geni più variabili)
    cat("Generazione plot marker genes...\n")
    top_gene_idx <- order(gene_vars, decreasing = TRUE)[1:4]
    
    for (i in seq_along(top_gene_idx)) {
      gene_idx <- top_gene_idx[i]
      gene_expr <- sim_results$expression[gene_idx, ]
      gene_name <- ifelse(!is.null(rownames(sim_results$expression)), 
                          rownames(sim_results$expression)[gene_idx],
                          paste0("gene_", gene_idx))
      
      plot_df <- data.frame(
        x = sim_results$coordinates$x,
        y = sim_results$coordinates$y,
        expression = log1p(gene_expr)
      )
      
      p_marker <- ggplot(plot_df, aes(x = x, y = y, color = expression)) +
        geom_point(size = 0.5) +
        scale_color_viridis_c() +
        theme_minimal() +
        labs(title = paste0("Espressione spaziale del gene ", gene_name),
             x = "Coordinata X", 
             y = "Coordinata Y", 
             color = "Log Expr") +
        theme(aspect.ratio = 1)
      
      ggsave(file.path(output_dir, paste0("marker_gene_", i, ".png")), p_marker, width = 10, height = 8)
    }
    
    cat("Plot di validazione generati con successo in:", output_dir, "\n")
  }
  
  # Prova prima con il modulo standard
  plots <- tryCatch({
    cat("Tentativo con modulo standard di validazione...\n")
    generate_validation_plots(
      sim_results = risultato,
      marker_genes = NULL,  # Rileva automaticamente
      n_markers = 3,
      output_dir = "plots/validazione",
      file_prefix = "visiumHD",
      file_format = "png",
      width = 10,
      height = 8,
      skip_dim_reduction = TRUE  # Salta UMAP/tSNE che richiedono molto calcolo
    )
    cat("Plot di validazione generati con metodo standard in plots/validazione/\n")
    TRUE
  }, error = function(e) {
    cat("Errore nel generatore standard:", conditionMessage(e), "\n")
    cat("Utilizzo generatore alternativo robusto...\n")
    # In caso di errore, usa la versione robusta
    generate_fixed_validation_plots(risultato, "plots/validazione")
    TRUE
  })
  toc()
}

cat("\nSimulazione full‑size completata con successo!\n")
cat("Risultati salvati in:", cfg$output_path, "\n")
cat("Plot salvati in:", cfg$output_plot, "\n")
if (validate) {
  cat("Plot di validazione salvati in: plots/validazione/\n")
}