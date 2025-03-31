#!/usr/bin/env Rscript

# Carico la funzione analyze_and_compare_clusters
source("R/analyze_and_compare_clusters.R")

# Impostazioni per la parallelizzazione
library(future)
future::plan(multisession, workers=4)

# Eseguo l'analisi sui dati simulati
result <- analyze_and_compare_clusters(
  rds_path = "data/simulated_visiumhd.rds",
  output_path = "results/visiumhd_metrics.csv",
  k_cell_types = 7,               # Numero di tipi cellulari nella simulazione
  seurat_resolution = 0.4,        # Risoluzione per il clustering Seurat (aumentata per ottenere più cluster)
  hdbscan_minPts = 10             # Parametro minPts per HDBSCAN
)

print("Analisi completata. Risultati salvati in results/visiumhd_metrics.csv")
print("Un oggetto Seurat con i cluster è stato salvato in results/visiumhd_seurat_hdbscan_combined.rds")