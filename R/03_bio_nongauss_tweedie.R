#!/usr/bin/env Rscript

# Script per l'analisi con il modello Tweedie

library(tidyverse)
library(mgcv)  # Per GAM con famiglia Tweedie
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

# Funzione per analizzare un dataset con il modello Tweedie
analyze_tweedie <- function(data_path) {
  tryCatch({
    # Notifica inizio analisi
    send_telegram_message(paste("Iniziata analisi Tweedie per:", data_path))

    # Carica il dataset
    data <- readRDS(data_path)

    # Prepara i dati per il modello Tweedie
    df <- data.frame(
      expression = as.vector(data$expression),
      x = rep(data$coordinates$x, ncol(data$expression)),
      y = rep(data$coordinates$y, ncol(data$expression)),
      gene = factor(rep(1:ncol(data$expression), each = nrow(data$expression)))
    )

    # Fit del modello GAM con famiglia Tweedie
    tweedie_model <- gam(expression ~ s(x, y, bs = "tp") + gene,
                         data = df,
                         family = tw(link = "log"),
                         method = "REML")

    # Calcola i valori predetti
    predicted_values <- predict(tweedie_model, type = "response")

    # Identifica hotspot basati sui valori predetti
    hotspots <- kmeans(matrix(predicted_values, ncol = ncol(data$expression)), centers = 3)$cluster

    # Salva i risultati
    saveRDS(list(
      hotspots = hotspots,
      model = tweedie_model,
      scores = predicted_values
    ), file = paste0("results/tweedie_results_", basename(data_path)))

    # Notifica completamento
    send_telegram_message(paste("Completata analisi Tweedie per:", data_path))

    return(hotspots)

  }, error = function(e) {
    send_telegram_message(paste("Errore nell'analisi Tweedie per:", data_path, "-", e$message))
    stop(e)
  })
}

# Configurazione parallelizzazione
plan(multisession, workers = min(3, availableCores()))

# Path dei dati simulati
data_paths <- list.files("data", pattern = "simulated_.*_correlation.rds", full.names = TRUE)

# Esegui analisi in parallelo
results <- future_lapply(data_paths, analyze_tweedie, future.seed = TRUE)

# Cleanup
plan(sequential)
