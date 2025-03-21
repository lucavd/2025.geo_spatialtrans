#!/usr/bin/env Rscript

# Script per la creazione di visualizzazioni diagnostiche a partire da simulazioni preesistenti
# Prende un file RDS generato dalla funzione simulate_spatial_transcriptomics() in image_simulations.R
# e produce una serie di grafici diagnostici e di visualizzazione

library(tidyverse)
library(reshape2)
library(scales)
library(MASS)
library(sp)
library(here)

generate_diagnostic_plots <- function(
  input_rds_path,           # Path del file RDS da analizzare
  output_dir = here::here("images/simulation_plots"), # Directory di output per i grafici
  base_name = NULL,         # Prefisso per i nomi dei file (se NULL, viene estratto dal nome del file RDS)
  plot_width = 8,           # Larghezza dei grafici in pollici
  plot_height = 6,          # Altezza dei grafici in pollici
  dpi = 300                 # Risoluzione dei grafici
) {
  # Verifica che il file di input esista
  if (!file.exists(input_rds_path)) {
    stop(sprintf("Il file %s non esiste", input_rds_path))
  }
  
  # Crea la directory di output se non esiste
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Estrai il nome base del file se non specificato
  if (is.null(base_name)) {
    base_name <- tools::file_path_sans_ext(basename(input_rds_path))
  }
  
  # Carica i dati simulati
  cat(sprintf("Caricamento dei dati da %s\n", input_rds_path))
  result <- readRDS(input_rds_path)
  
  # Estrai i componenti principali dai dati
  coordinates <- result$coordinates
  intensity_cluster <- result$intensity_cluster
  expression_data <- result$expression
  
  # Estrai eventuali parametri salvati (se disponibili)
  parameters <- if ("parameters" %in% names(result)) result$parameters else NULL
  
  # Ricalcola le metriche principali necessarie per i grafici
  cat("Calcolo delle metriche di distanza e densità...\n")
  
  # Calcolo delle distanze tra cellule
  cell_df <- cbind(coordinates, cluster = intensity_cluster)
  N <- nrow(cell_df)
  cluster_labels <- cell_df$cluster
  dist_mat <- as.matrix(dist(cell_df[, c("x", "y")]))
  
  # Calcola la distanza media di ciascuna cellula rispetto alle altre del proprio cluster
  mean_dist <- numeric(N)
  for (i in 1:N) {
    mean_dist[i] <- mean(dist_mat[i, cluster_labels == cluster_labels[i]])
  }
  
  # Calcola la densità locale
  local_density <- apply(dist_mat, 1, function(row) mean(row < quantile(row, 0.1)))
  
  # ============================
  # 1) Distribuzione spaziale dei cluster
  # ============================
  cat("Generazione del grafico di distribuzione spaziale...\n")
  
  p1 <- ggplot(cell_df, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    labs(title = "Distribuzione spaziale delle cellule",
         subtitle = "Colorazione per tipo cellulare",
         color = "Cluster")
  
  ggsave(file.path(output_dir, paste0(base_name, "_spatial_distribution.png")), 
         plot = p1, width = plot_width, height = plot_height, dpi = dpi)
  
  # ============================
  # 2) Profilo medio di espressione per cluster
  # ============================
  cat("Generazione del profilo di espressione per cluster...\n")
  
  # Calcola il profilo medio di espressione per ogni cluster
  n_genes <- ncol(expression_data)
  k_cell_types <- length(levels(intensity_cluster))
  mean_expression_list <- list()
  
  for (k in 1:k_cell_types) {
    cluster_cells <- which(as.integer(intensity_cluster) == k)
    mean_expr <- colMeans(expression_data[cluster_cells, , drop = FALSE])
    mean_expression_list[[k]] <- mean_expr
  }
  
  # Combina i profili medi in un unico data frame
  mean_expr_df <- do.call(rbind, lapply(1:length(mean_expression_list), function(k) {
    data.frame(Cluster = factor(k),
               Gene = 1:n_genes,
               Expression = mean_expression_list[[k]])
  }))
  
  # Crea una heatmap dei profili medi di espressione per ogni cluster
  p2 <- ggplot(mean_expr_df, aes(x = Cluster, y = Gene, fill = log1p(Expression))) +
    geom_tile() +
    scale_fill_gradientn(colors = c("steelblue4", "lightskyblue", "yellow", "red"), 
                         name = "log(expr+1)") +
    labs(title = "Profilo medio di espressione per cluster",
         x = "Cluster",
         y = "Gene") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylim(c(0, min(n_genes, 100)))  # Limita la visualizzazione ai primi 100 geni
  
  ggsave(file.path(output_dir, paste0(base_name, "_expression_profile.png")), 
         plot = p2, width = plot_width, height = plot_height, dpi = dpi)
  
  # ============================
  # 3) Distribuzione delle distanze per cluster
  # ============================
  cat("Generazione del grafico di distanze per cluster...\n")
  
  # Aggiunge a cell_df le informazioni sulla distanza media e densità locale
  scatter_df <- cell_df %>%
    mutate(mean_distance = mean_dist,
           density = local_density)
  
  # Distribuzione spaziale delle cellule colorata in base alla distanza media
  p3 <- ggplot(scatter_df, aes(x = x, y = y, color = mean_distance)) +
    geom_point(alpha = 0.7, size = 0.5) +
    facet_wrap(~ cluster, ncol = 2) +
    scale_color_gradient(low = "blue", high = "red", name = "Distanza media") +
    labs(title = "Distanza media intra-cluster",
         subtitle = "Distribuzione spaziale per cluster",
         x = "x",
         y = "y") +
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "white", size = 0.3),
          panel.grid.minor = element_line(color = "white", size = 0.3))
  
  ggsave(file.path(output_dir, paste0(base_name, "_distances_by_cluster.png")), 
         plot = p3, width = plot_width + 2, height = plot_height + 1, dpi = dpi)
  
  # ============================
  # 4) Relazione tra distanza e dropout
  # ============================
  cat("Analisi della relazione tra distanza e dropout...\n")
  
  # Calcola la proporzione di zeri per cella
  zero_prop <- rowMeans(expression_data == 0)
  
  # Crea dataframe per il grafico
  dropout_df <- data.frame(
    mean_distance = mean_dist,
    zero_proportion = zero_prop,
    cluster = intensity_cluster
  )
  
  # Crea il grafico per la relazione tra distanza e dropout
  p4 <- ggplot(dropout_df, aes(x = mean_distance, y = zero_proportion, color = cluster)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "loess", se = TRUE, aes(group = 1), color = "black") +
    labs(title = "Relazione tra distanza e dropout",
         subtitle = "Proporzione di geni non espressi (zero counts)",
         x = "Distanza media intra-cluster",
         y = "Proporzione di zeri",
         color = "Cluster") +
    theme_minimal()
  
  ggsave(file.path(output_dir, paste0(base_name, "_dropout_distance.png")), 
         plot = p4, width = plot_width, height = plot_height, dpi = dpi)
  
  # ============================
  # 5) Distribuzione dell'espressione genica
  # ============================
  cat("Generazione del grafico di distribuzione dell'espressione genica...\n")
  
  # Identifica potenziali geni stabili (con bassa varianza)
  gene_var <- apply(expression_data, 2, var)
  gene_mean <- colMeans(expression_data)
  cv2 <- gene_var / (gene_mean^2)  # Coefficiente di variazione quadratico
  
  # Definisci geni a bassa e alta variabilità
  low_var_genes <- which(cv2 < quantile(cv2, 0.25))
  high_var_genes <- which(cv2 > quantile(cv2, 0.75))
  
  # Crea dataframe per i grafici
  df_expr_low <- data.frame(
    Expression = as.vector(expression_data[, low_var_genes]),
    GeneType = "Geni a bassa variabilità"
  )
  
  df_expr_high <- data.frame(
    Expression = as.vector(expression_data[, high_var_genes]),
    GeneType = "Geni ad alta variabilità"
  )
  
  df_expr_combined <- rbind(df_expr_low, df_expr_high)
  
  # Crea istogrammi separati per tipo di gene
  p5 <- ggplot(df_expr_combined, aes(x = Expression, fill = GeneType)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    labs(title = "Distribuzione dell'espressione genica",
         subtitle = "Confronto tra geni a bassa e alta variabilità",
         x = "Conta di espressione",
         y = "Frequenza",
         fill = "Tipo di gene") +
    facet_wrap(~GeneType, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Geni a bassa variabilità" = "palegreen3",
                                "Geni ad alta variabilità" = "steelblue")) +
    theme_minimal() +
    xlim(c(0, min(200, max(df_expr_combined$Expression)))) +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, paste0(base_name, "_expression_distribution.png")), 
         plot = p5, width = plot_width, height = plot_height, dpi = dpi)
  
  # ============================
  # 6) Heatmap di espressione per un sottoinsieme di geni marker
  # ============================
  cat("Generazione della heatmap di espressione per geni marker...\n")
  
  # Identifica i geni più discriminanti tra i cluster
  cluster_means <- matrix(0, nrow = k_cell_types, ncol = n_genes)
  for (k in 1:k_cell_types) {
    cluster_cells <- which(as.integer(intensity_cluster) == k)
    cluster_means[k, ] <- colMeans(expression_data[cluster_cells, , drop = FALSE])
  }
  
  # Calcola quanto ogni gene è specifico per un cluster (max - media degli altri)
  gene_specificity <- apply(cluster_means, 2, function(x) {
    max_cluster <- which.max(x)
    max_value <- x[max_cluster]
    other_mean <- mean(x[-max_cluster])
    return(c(max_cluster, max_value - other_mean))
  })
  
  # Seleziona i top N geni più specifici
  n_top_genes <- min(50, n_genes)
  top_genes <- order(gene_specificity[2, ], decreasing = TRUE)[1:n_top_genes]
  
  # Prepara i dati per la heatmap
  set.seed(42)  # per riproducibilità nel campionamento
  cells_per_cluster <- min(500, N / k_cell_types)  # al massimo 500 cellule per cluster
  sampled_cells <- c()
  
  for (k in 1:k_cell_types) {
    cluster_cells <- which(as.integer(intensity_cluster) == k)
    if (length(cluster_cells) > cells_per_cluster) {
      sampled_cells <- c(sampled_cells, 
                         sample(cluster_cells, cells_per_cluster))
    } else {
      sampled_cells <- c(sampled_cells, cluster_cells)
    }
  }
  
  # Crea matrice di espressione per le cellule campionate e i geni top
  expr_subset <- expression_data[sampled_cells, top_genes]
  
  # Normalizza per riga (z-score)
  expr_z <- t(scale(t(expr_subset)))
  
  # Converti in dataframe lungo per ggplot
  heatmap_data <- as.data.frame(expr_z) %>%
    mutate(Cell = 1:n(), Cluster = intensity_cluster[sampled_cells]) %>%
    pivot_longer(cols = -c(Cell, Cluster), 
                names_to = "Gene", values_to = "Expression_Z") %>%
    mutate(Gene = factor(Gene))
  
  # Ordina le cellule per cluster e poi per espressione media
  cell_order <- heatmap_data %>%
    group_by(Cell) %>%
    summarise(Cluster = first(Cluster),
              Mean_Expr = mean(Expression_Z)) %>%
    arrange(Cluster, Mean_Expr) %>%
    pull(Cell)
  
  heatmap_data$Cell <- factor(heatmap_data$Cell, levels = cell_order)
  
  # Crea la heatmap
  p6 <- ggplot(heatmap_data, aes(x = Gene, y = Cell, fill = Expression_Z)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("navy", "white", "firebrick"), 
                         limits = c(-2, 2), 
                         name = "Z-score") +
    labs(title = "Heatmap di espressione per i geni marker top",
         subtitle = paste("Mostra", n_top_genes, "geni più discriminanti"),
         x = "Gene",
         y = "Cellula") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(Cluster ~ ., scales = "free_y", space = "free_y")
  
  ggsave(file.path(output_dir, paste0(base_name, "_marker_heatmap.png")), 
         plot = p6, width = plot_width, height = plot_height + 2, dpi = dpi)
  
  # ============================
  # 7) Distribuzione spaziale dell'espressione per geni rappresentativi
  # ============================
  cat("Generazione di mappe di espressione spaziale per geni selezionati...\n")
  
  # Seleziona un gene rappresentativo per ogni cluster
  representative_genes <- sapply(1:k_cell_types, function(k) {
    # Trova il gene con la massima specificità per questo cluster
    which.max(gene_specificity[2, ] * (gene_specificity[1, ] == k))
  })
  
  # Combina tutti i dati spaziali con i valori di espressione
  spatial_expr_data <- cbind(
    cell_df[, c("x", "y", "cluster")],
    expression_data[, representative_genes]
  )
  
  # Rinomina le colonne con nomi significativi
  colnames(spatial_expr_data)[(ncol(spatial_expr_data) - k_cell_types + 1):ncol(spatial_expr_data)] <- 
    paste0("marker_cluster_", 1:k_cell_types)
  
  # Converti in formato lungo per ggplot
  spatial_expr_long <- spatial_expr_data %>%
    pivot_longer(cols = starts_with("marker_"), 
                names_to = "marker_gene", 
                values_to = "expression")
  
  # Crea una mappa spaziale per ogni gene rappresentativo
  p7 <- ggplot(spatial_expr_long, aes(x = x, y = y, color = expression)) +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_color_gradientn(colors = c("gray90", "blue", "red"), 
                         name = "Espressione") +
    facet_wrap(~ marker_gene, ncol = 2) +
    scale_y_reverse() +
    labs(title = "Distribuzione spaziale dell'espressione genica",
         subtitle = "Un gene marker rappresentativo per ogni cluster",
         x = "Coordinata X", 
         y = "Coordinata Y") +
    theme_minimal()
  
  ggsave(file.path(output_dir, paste0(base_name, "_spatial_expression.png")), 
         plot = p7, width = plot_width + 2, height = plot_height + 2, dpi = dpi)
  
  # Restituisci l'elenco dei file creati
  output_files <- list.files(output_dir, pattern = paste0(base_name, "_"), full.names = TRUE)
  cat("Generazione dei grafici completata.\n")
  cat("File generati:\n")
  cat(paste(" -", basename(output_files)), sep = "\n")
  
  return(invisible(output_files))
}

# ========================
# Esempi d'uso
# ========================

# Esempio: genera grafici diagnostici per un dataset simulato
# generate_diagnostic_plots(
#   input_rds_path = "data/simulated_image_correlation.rds",
#   output_dir = "images/simulation_plots",
#   base_name = "sim_analysis"  # opzionale, altrimenti usa il nome del file
# )