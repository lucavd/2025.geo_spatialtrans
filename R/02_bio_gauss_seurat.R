#!/usr/bin/env Rscript

# Script per l'analisi con Seurat (metodo bioinformatico con errori gaussiani)

library(tidyverse)
library(Seurat)
library(telegram.bot)
library(future)
library(future.apply)

# Configurazione del bot Telegram
bot <- Bot(token = Sys.getenv("R_TELEGRAM_BOT_lucavdRbot"))
chat_id <- Sys.getenv("R_TELEGRAM_BOT_lucavdRbot_chatID")

# Funzione per inviare messaggi su Telegram
send_telegram_message <- function(message) {
  tryCatch({
    bot$send_message(chat_id = chat_id, text = message)
  }, error = function(e) {
    warning("Errore nell'invio del messaggio Telegram: ", e$message)
  })
}

# Funzione per analizzare un dataset con Seurat
analyze_seurat <- function(data_path) {
  tryCatch({
    # Notifica inizio analisi
    send_telegram_message(paste("Iniziata analisi Seurat per:", data_path))

    # Carica il dataset
    data <- readRDS(data_path)

    # Crea oggetto Seurat
    seu <- CreateSeuratObject(counts = t(data$expression))

    # Aggiungi coordinate spaziali
    seu@meta.data$x <- data$coordinates$x
    seu@meta.data$y <- data$coordinates$y

    # Normalizzazione e selezione delle feature variabili
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

    # Scaling
    seu <- ScaleData(seu, verbose = FALSE)

    # PCA (aggiunta di un controllo dinamico sul numero di dimensioni computate)
    seu <- RunPCA(seu, verbose = FALSE)

    # Calcolo dinamico delle dimensioni disponibili
    dims_available <- ncol(seu@reductions$pca@cell.embeddings)
    dims_to_use <- min(10, dims_available)

    # Clustering
    seu <- FindNeighbors(seu, dims = 1:dims_to_use, verbose = FALSE)
    seu <- FindClusters(seu, resolution = 0.9, verbose = FALSE)

    # Genera il vettore degli hotspots
    hotspots <- as.integer(Idents(seu))

    # Salva risultati
    saveRDS(list(
      hotspots = hotspots,
      seurat_obj = seu
    ), file = paste0("results/seurat_results_", basename(data_path)))

    # Notifica completamento
    send_telegram_message(paste("Completata analisi Seurat per:", data_path))

    return(hotspots)

  }, error = function(e) {
    send_telegram_message(paste("Errore nell'analisi Seurat per:", data_path, "-", e$message))
    stop(e)
  })
}

# Configurazione parallelizzazione
plan(multisession, workers = min(3, availableCores()))

# Path dei dati simulati
data_paths <- list.files("data", pattern = "simulated_.*_correlation.rds", full.names = TRUE)

# Esegui analisi in parallelo
results <- future_lapply(data_paths, analyze_seurat, future.seed = TRUE)

# Cleanup
plan(sequential)
