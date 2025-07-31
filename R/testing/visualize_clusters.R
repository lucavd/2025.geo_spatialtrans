#!/usr/bin/env Rscript
# Script per visualizzare i cluster colorati sull'immagine originale
# Simile alla visualizzazione mostrata nello screenshot

cat("=== VISUALIZZAZIONE CLUSTER COLORATI ===\n")

# Carica librerie necessarie
library(ggplot2)
library(png)
library(dplyr)
library(RColorBrewer)

# Carica i dati della simulazione
if (!file.exists("R/testing/full_simulation_data.rds")) {
  stop("File dati simulazione non trovato. Eseguire prima full_test.R")
}

cat("Caricamento dati simulazione...\n")
sim_data <- readRDS("R/testing/full_simulation_data.rds")

# Estrai componenti
coordinates <- sim_data$coordinates
clusters <- sim_data$clusters
config <- sim_data$config

cat("✓ Dati caricati:\n")
cat("  - Coordinate:", nrow(coordinates), "punti\n")
cat("  - Cluster unici:", length(unique(clusters)), "\n")

# Carica immagine originale se disponibile
img_path <- "R/testing/full_tissue_complex.png"
if (file.exists(img_path)) {
  cat("Caricamento immagine originale...\n")
  img_array <- readPNG(img_path)

  # Converti in dataframe per ggplot
  img_df <- expand.grid(x = 1:dim(img_array)[2], y = 1:dim(img_array)[1])
  img_df$intensity <- as.vector(img_array[nrow(img_array):1, ])

  cat("✓ Immagine caricata:", dim(img_array)[2], "x", dim(img_array)[1], "pixel\n")
} else {
  cat("⚠ Immagine originale non trovata, creo solo plot cluster\n")
  img_df <- NULL
}

# Prepara dati per visualizzazione
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

# ==================== PLOT 1: IMMAGINE ORIGINALE ====================
if (!is.null(img_df)) {
  p1 <- ggplot(img_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = intensity), alpha = 0.8) +
    scale_fill_gradient(low = "white", high = "darkblue", guide = "none") +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(title = "Original Synthetic Tissue")

  # Salva plot immagine originale
  ggsave("R/testing/original_tissue_plot.png", p1,
         width = 8, height = 6, dpi = 300, bg = "white")
  cat("✓ Plot immagine originale salvato\n")
}

# ==================== PLOT 2: CLUSTER COLORATI ====================
# Flippa le coordinate X per matchare l'orientamento dell'immagine originale
cluster_data$x_flipped <- max(cluster_data$x) - cluster_data$x

p2 <- ggplot(cluster_data, aes(x = x_flipped, y = y, color = cluster)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_manual(values = colors, name = "Cell Type") +
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
    title = "Spatial Transcriptomics Clusters",
    x = "X",
    y = "Y"
  )

# Salva plot cluster
ggsave("R/testing/cluster_visualization.png", p2,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("✓ Plot cluster colorati salvato\n")

# ==================== PLOT 3: PANNELLO COMBINATO ====================
if (!is.null(img_df)) {
  # Crea plot combinato simile all'immagine di riferimento
  library(gridExtra)
  library(grid)

  # Aggiusta dimensioni per il pannello
  p1_panel <- p1 + theme(plot.margin = margin(10, 10, 5, 10))
  p2_panel <- p2 + theme(plot.margin = margin(5, 10, 10, 10))

  # Combina i plot
  combined_plot <- grid.arrange(
    p1_panel, p2_panel,
    nrow = 2,
    heights = c(1, 1.2)  # Più spazio per il plot con legenda
  )

  # Salva plot combinato
  ggsave("R/testing/combined_visualization.png", combined_plot,
         width = 10, height = 12, dpi = 300, bg = "white")
  cat("✓ Plot combinato salvato\n")
}

# ==================== STATISTICHE CLUSTER ====================
cat("\n=== STATISTICHE CLUSTER ===\n")
cluster_stats <- cluster_data %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_x = round(mean(x), 1),
    mean_y = round(mean(y), 1),
    .groups = 'drop'
  ) %>%
  arrange(cluster)

print(cluster_stats)

cat("\n✓ Visualizzazioni create con successo!\n")
cat("Files generati:\n")
if (!is.null(img_df)) {
  cat("  - original_tissue_plot.png (immagine originale)\n")
  cat("  - combined_visualization.png (pannello combinato)\n")
}
cat("  - cluster_visualization.png (cluster colorati)\n")

cat("\nLa visualizzazione replica lo stile dell'immagine di riferimento\n")
cat("con", n_clusters, "cluster distinti colorati spazialmente.\n")
