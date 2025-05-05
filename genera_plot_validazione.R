#!/usr/bin/env Rscript
# Script per generare solo i plot di validazione dai risultati salvati

# 1. Caricamento delle funzioni
cat("Caricamento delle funzioni...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) {
  cat("Caricamento:", f, "\n")
  source(f)
}

# 2. Caricamento librerie necessarie
cat("Caricamento librerie...\n")
library(png)
library(ggplot2)
library(dplyr)
library(Matrix)
library(sp)
library(tictoc)

# 3. Caricamento dei risultati salvati
cat("Caricamento risultati salvati...\n")
risultato <- readRDS("results/visiumHD_full.rds")

# 4. Visualizza informazioni sui dati caricati
cat("Dimensioni matrice espressione:", dim(risultato$expression), "\n")
cat("Numero di punti spaziali:", nrow(risultato$coordinates), "\n")
cat("Livelli di cluster:", nlevels(risultato$intensity_cluster), "\n")

# 5. Creazione dei plot di validazione
cat("Generazione plot di validazione...\n")
dir.create("plots/validazione", recursive = TRUE, showWarnings = FALSE)

tic("Plot di validazione")
plots <- tryCatch({
  cat("Chiamata alla funzione generate_validation_plots...\n")
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
  cat("Plot di validazione generati in plots/validazione/\n")
}, error = function(e) {
  cat("ERRORE nella generazione dei plot di validazione:", conditionMessage(e), "\n")
  print(e)
  NULL
})
toc()

# Verifica dei file generati
cat("\nFile generati nella directory plots/validazione:\n")
system("ls -la plots/validazione/")