#!/usr/bin/env Rscript

# Script per la simulazione di dati di trascrittomica spaziale (aggiornato)
# Genera dati simulati con caratteristiche più realistiche

library(tidyverse)
library(MASS)
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

# Parametri della simulazione
n_cells <- 100  # Numero di cellule (ridotto per test iniziale)
n_genes <- 10   # Numero di geni (ridotto per test iniziale)
grid_size <- 100 # Dimensione della griglia spaziale

# Funzione per generare coordinate spaziali
generate_spatial_coordinates <- function(n_cells, grid_size) {
  data.frame(
    x = runif(n_cells, 0, grid_size),
    y = runif(n_cells, 0, grid_size)
  )
}

# Funzione per generare matrice di correlazione spaziale
generate_spatial_correlation <- function(coords, correlation_level) {
  dist_matrix <- as.matrix(dist(coords))

  # Parametri per diversi livelli di correlazione
  range_param <- switch(correlation_level,
                        "high" = 20,
                        "medium" = 50,
                        "low" = 80
  )

  # Matrice di correlazione spaziale usando funzione esponenziale
  cor_matrix <- exp(-dist_matrix / range_param)
  diag(cor_matrix) <- 1
  return(cor_matrix)
}

# Funzione per generare espressione genica con dipendenza spaziale
generate_expression_data <- function(n_cells, n_genes, spatial_cor, correlation_level) {
  # Parametri per diversi livelli di correlazione dell'espressione
  gene_cor <- switch(correlation_level,
                     "high" = 0.8,
                     "medium" = 0.5,
                     "low" = 0.2
  )

  # Matrice di varianza-covarianza per i geni
  gene_sigma <- matrix(gene_cor, n_genes, n_genes)
  diag(gene_sigma) <- 1

  # Generazione dati con mvrnorm
  expression_data <- mvrnorm(n_cells, mu = rep(0, n_genes), Sigma = gene_sigma)

  # Aggiunta della componente spaziale
  for (i in 1:n_genes) {
    spatial_effect <- mvrnorm(1, mu = rep(0, n_cells), Sigma = spatial_cor)
    expression_data[, i] <- expression_data[, i] + spatial_effect
  }

  # Trasformazione in conteggi usando distribuzione negativa binomiale
  expression_data <- exp(expression_data)
  dispersion <- 0.5  # Dispersione per distribuzione NB
  expression_data <- matrix(rnbinom(length(expression_data), mu = expression_data, size = dispersion),
                            nrow = n_cells)

  # Introduzione di zeri sparsi (dropout rate)
  dropout_rate <- 0.4
  expression_data[sample(length(expression_data), size = length(expression_data) * dropout_rate)] <- 0

  return(expression_data)
}

# Funzione principale per la simulazione
simulate_spatial_transcriptomics <- function(correlation_level) {
  tryCatch({
    send_telegram_message(paste("Iniziata simulazione per livello di correlazione:",
                                correlation_level))

    # Genera coordinate
    coords <- generate_spatial_coordinates(n_cells, grid_size)

    # Genera correlazione spaziale
    spatial_cor <- generate_spatial_correlation(coords, correlation_level)

    # Genera dati di espressione
    expression_data <- generate_expression_data(n_cells, n_genes,
                                                spatial_cor, correlation_level)

    # Creazione dei veri hotspot per validazione
    true_hotspots <- kmeans(coords, centers = 3)$cluster

    # Preparazione output
    result <- list(
      coordinates = coords,
      expression = expression_data,
      spatial_correlation = spatial_cor,
      true_hotspots = true_hotspots,
      parameters = list(
        n_cells = n_cells,
        n_genes = n_genes,
        correlation_level = correlation_level
      )
    )

    # Salva i risultati
    saveRDS(result,
            file = paste0("data/simulated_", correlation_level, "_correlation.rds"))

    send_telegram_message(paste("Completata simulazione per livello di correlazione:",
                                correlation_level))

    return(result)

  }, error = function(e) {
    send_telegram_message(paste("Errore nella simulazione:", e$message))
    stop(e)
  })
}

# Configurazione parallelizzazione
plan(multisession, workers = min(3, availableCores()))

# Esegui simulazioni per i tre livelli di correlazione
results <- future_lapply(c("high", "medium", "low"), simulate_spatial_transcriptomics,
                         future.seed = TRUE)

# Cleanup
plan(sequential)
