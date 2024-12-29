#!/usr/bin/env Rscript

# Script per l'analisi geospaziale con GAM (metodo geospaziale con errori gaussiani)

library(tidyverse)
library(mgcv)  # Per GAM
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

# Funzione per analizzare un dataset con GAM
analyze_gam <- function(data_path) {
  tryCatch({
    # Notifica inizio analisi
    send_telegram_message(paste("Iniziata analisi GAM geospaziale per:", data_path))

    # Carica il dataset
    data <- readRDS(data_path)

    # Prepara i dati per il modello GAM
    df <- data.frame(
      expression = as.vector(data$expression),
      x = rep(data$coordinates$x, ncol(data$expression)),
      y = rep(data$coordinates$y, ncol(data$expression)),
      gene = factor(rep(1:ncol(data$expression), each = nrow(data$expression)))
    )

    # Fit del modello GAM
    gam_model <- gam(expression ~ s(x, y, bs = "tp") + gene,
                     data = df,
                     family = gaussian(),
                     method = "REML")

    # Calcola i valori predetti
    predicted_values <- predict(gam_model, type = "response")

    # Identifica hotspot basati sui valori predetti
    hotspots <- kmeans(matrix(predicted_values, ncol = ncol(data$expression)), centers = 3)$cluster

    # Salva i risultati
    saveRDS(list(
      hotspots = hotspots,
      model = gam_model,
      scores = predicted_values
    ), file = paste0("results/gam_geospatial_results_", basename(data_path)))

    # Notifica completamento
    send_telegram_message(paste("Completata analisi GAM geospaziale per:", data_path))

    return(hotspots)

  }, error = function(e) {
    send_telegram_message(paste("Errore nell'analisi GAM geospaziale per:", data_path, "-", e$message))
    stop(e)
  })
}

# Configurazione parallelizzazione
plan(multisession, workers = min(3, availableCores()))

# Path dei dati simulati
data_paths <- list.files("data", pattern = "simulated_.*_correlation.rds", full.names = TRUE)

# Esegui analisi in parallelo
results <- future_lapply(data_paths, analyze_gam, future.seed = TRUE)

# Cleanup
plan(sequential)
