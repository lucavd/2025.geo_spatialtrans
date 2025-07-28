#!/usr/bin/env Rscript
# Script per visualizzare i cluster MRF colorati sulla griglia originale
# Simile a visualize_clusters.R ma specifico per dati MRF

cat("=== VISUALIZZAZIONE CLUSTER MRF ===\n")

# Carica librerie necessarie
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)

# Carica librerie opzionali
grid_extra_available <- requireNamespace("gridExtra", quietly = TRUE)
if (grid_extra_available) {
  library(gridExtra)
  library(grid)
}

# Carica i dati della simulazione
if (!file.exists("R_simple/testing/full_simulation_data.rds")) {
  stop("File dati simulazione non trovato. Eseguire prima full_test.R")
}

cat("Caricamento dati simulazione...\n")
sim_data <- readRDS("R_simple/testing/full_simulation_data.rds")

# Verifica che sia modalità MRF
if (sim_data$config$spatial_engine != "mrf") {
  stop("Questo script è specifico per la modalità MRF. Usare visualize_clusters.R per modalità image.")
}

# Estrai componenti
coordinates <- sim_data$coordinates
clusters <- sim_data$clusters
config <- sim_data$config

cat("✓ Dati caricati:\n")
cat("  - Coordinate:", nrow(coordinates), "punti\n")
cat("  - Cluster unici:", length(unique(clusters)), "\n")
cat("  - Griglia MRF:", config$mrf_grid_size[1], "x", config$mrf_grid_size[2], "\n")

# Carica i dati MRF originali se disponibili
mrf_result_path <- "R_simple/testing/full_test_result.rds"
if (file.exists(mrf_result_path)) {
  cat("Caricamento dati MRF originali...\n")
  full_result <- readRDS(mrf_result_path)

  # Estrai la griglia MRF originale
  if (!is.null(full_result$syn_tissue) && !is.null(full_result$syn_tissue$img_df)) {
    mrf_grid <- full_result$syn_tissue$img_df
    cat("✓ Griglia MRF caricata:", nrow(mrf_grid), "pixel\n")
  } else {
    cat("⚠ Dati MRF originali non trovati nel risultato\n")
    mrf_grid <- NULL
  }
} else {
  cat("⚠ File risultato MRF non trovato, creo solo plot cluster campionati\n")
  mrf_grid <- NULL
}

# Prepara dati per visualizzazione dei cluster campionati
cluster_data <- data.frame(
  x = coordinates$x,
  y = coordinates$y,
  cluster = as.factor(clusters)
)

# Colori per i cluster (simili a quelli dell'immagine di riferimento)
n_clusters <- length(unique(clusters))
if (n_clusters <= 8) {
  # Usa palette specifica per replicare i colori dell'immagine
  colors <- c("#FF1493", "#00CED1", "#32CD32", "#FFD700",
              "#FF4500", "#9370DB", "#40E0D0", "#FF69B4")[1:n_clusters]
} else {
  # Per più cluster, usa palette automatica
  colors <- rainbow(n_clusters)
}

cat("Creazione visualizzazioni...\n")

# ==================== PLOT 1: GRIGLIA MRF ORIGINALE ====================
if (!is.null(mrf_grid)) {
  # Prepara dati griglia MRF per visualizzazione
  mrf_plot_data <- data.frame(
    x = mrf_grid$x,
    y = mrf_grid$y,
    cell_type = as.factor(mrf_grid$cell_type)
  )

  # Colori per i tipi cellulari MRF originali
  n_cell_types <- length(unique(mrf_plot_data$cell_type))
  mrf_colors <- rainbow(n_cell_types)

  p1 <- ggplot(mrf_plot_data, aes(x = x, y = y, color = cell_type)) +
    geom_point(size = 0.1, alpha = 0.7) +
    scale_color_manual(values = mrf_colors, name = "MRF Cell Type") +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(
      title = paste("Original MRF Grid (", config$mrf_grid_size[1], "x", config$mrf_grid_size[2], ")"),
      x = "X",
      y = "Y"
    )

  # Salva plot griglia MRF originale
  ggsave("R_simple/testing/mrf_original_grid.png", p1,
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("✓ Plot griglia MRF originale salvato\n")
}

# ==================== PLOT 2: CLUSTER CAMPIONATI ====================
p2 <- ggplot(cluster_data, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 1.0, alpha = 0.8) +
  scale_color_manual(values = colors, name = "Cluster") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = paste("Sampled Spatial Transcriptomics Clusters (", nrow(cluster_data), " cells)"),
    x = "X",
    y = "Y"
  )

# Salva plot cluster campionati
ggsave("R_simple/testing/mrf_cluster_visualization.png", p2,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("✓ Plot cluster MRF campionati salvato\n")

# ==================== PLOT 3: PANNELLO COMBINATO ====================
if (!is.null(mrf_grid) && grid_extra_available) {
  # Crea plot combinato che mostra sia la griglia originale che i cluster campionati
  p1_panel <- p1 + theme(plot.margin = margin(10, 10, 5, 10))
  p2_panel <- p2 + theme(plot.margin = margin(5, 10, 10, 10))

  # Combina i plot
  combined_plot <- grid.arrange(
    p1_panel, p2_panel,
    nrow = 2,
    heights = c(1, 1.2)  # Più spazio per il plot con legenda
  )

  # Salva plot combinato
  ggsave("R_simple/testing/mrf_combined_visualization.png", combined_plot,
         width = 10, height = 12, dpi = 300, bg = "white")
  cat("✓ Plot MRF combinato salvato\n")
} else if (!is.null(mrf_grid) && !grid_extra_available) {
  cat("⚠ gridExtra non disponibile, salto plot combinato\n")
}

# ==================== PLOT 4: OVERLAY COMPARISON ====================
if (!is.null(mrf_grid)) {
  # Crea un plot che sovrappone i cluster campionati sulla griglia MRF
  p4 <- ggplot() +
    # Griglia MRF di sfondo (più trasparente)
    geom_point(data = mrf_plot_data, aes(x = x, y = y, color = cell_type),
               size = 0.05, alpha = 0.3) +
    scale_color_manual(values = mrf_colors, name = "MRF Background") +
    # Cluster campionati in primo piano
    new_scale_color() +
    geom_point(data = cluster_data, aes(x = x, y = y, color = cluster),
               size = 1.5, alpha = 0.9) +
    scale_color_manual(values = colors, name = "Sampled Clusters") +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(
      title = "MRF Grid with Sampled Clusters Overlay",
      x = "X",
      y = "Y"
    )

  # Nota: il pacchetto ggnewscale è necessario per multiple scale_color
  # Se non disponibile, creiamo una versione semplificata
  tryCatch({
    library(ggnewscale)
    ggsave("R_simple/testing/mrf_overlay_visualization.png", p4,
           width = 12, height = 8, dpi = 300, bg = "white")
    cat("✓ Plot overlay MRF salvato\n")
  }, error = function(e) {
    cat("⚠ ggnewscale non disponibile, salto plot overlay\n")
  })
}

# ==================== STATISTICHE CLUSTER ====================
cluster_stats <- cluster_data %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_x = round(mean(x), 1),
    mean_y = round(mean(y), 1),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_cells))

cat("\n=== STATISTICHE CLUSTER MRF ===\n")
print(cluster_stats)

# Statistiche spaziali
cat("\n=== STATISTICHE SPAZIALI ===\n")
cat("Range X:", min(cluster_data$x), "-", max(cluster_data$x), "\n")
cat("Range Y:", min(cluster_data$y), "-", max(cluster_data$y), "\n")
cat("Densità media:", round(nrow(cluster_data) / ((max(cluster_data$x) - min(cluster_data$x)) *
                                                   (max(cluster_data$y) - min(cluster_data$y))), 4), "celle/unità²\n")

if (!is.null(mrf_grid)) {
  cat("Rapporto campionamento:", round(nrow(cluster_data) / nrow(mrf_grid) * 100, 2), "%\n")
  cat("Griglia originale:", nrow(mrf_grid), "pixel\n")
  cat("Celle campionate:", nrow(cluster_data), "celle\n")
}

cat("\n✓ Visualizzazione MRF completata!\n")
cat("File generati:\n")
if (!is.null(mrf_grid)) {
  cat("  - mrf_original_grid.png: Griglia MRF originale\n")
  if (grid_extra_available && file.exists("R_simple/testing/mrf_combined_visualization.png")) {
    cat("  - mrf_combined_visualization.png: Pannello combinato\n")
  }
  if (file.exists("R_simple/testing/mrf_overlay_visualization.png")) {
    cat("  - mrf_overlay_visualization.png: Overlay comparison\n")
  }
}
cat("  - mrf_cluster_visualization.png: Cluster campionati\n")
