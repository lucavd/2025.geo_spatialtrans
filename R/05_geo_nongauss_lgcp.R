#!/usr/bin/env Rscript

# Script per l'analisi con Log-Gaussian Cox Process (LGCP)

library(spatstat)
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

# Funzione per analizzare un dataset con LGCP
analyze_lgcp <- function(data_path) {
  tryCatch({
    # Notifica inizio analisi
    send_telegram_message(paste("Iniziata analisi LGCP per:", data_path))

    # Carica il dataset
    data <- readRDS(data_path)

    # Estrai coordinate e dati di espressione
    coords <- data.frame(x = data$coordinates$x, y = data$coordinates$y)
    expression <- data$expression

    # Crea un oggetto ppp (point pattern) per spatstat
    ppp_data <- ppp(coords$x, coords$y, window = owin(range(coords$x), range(coords$y)))

    # Aggiungi covariata di espressione genica
    covariate <- as.im(matrix(expression, nrow = length(unique(coords$y)), ncol = length(unique(coords$x))),
                       W = owin(range(coords$x), range(coords$y)))

    # Stima del modello LGCP
    lgcp_model <- kppm(ppp_data ~ covariate, clusters = "LGCP")

    # Predizione dell'intensità
    intensity_pred <- predict(lgcp_model, type = "intensity")

    # Identificazione degli hotspot
    threshold <- quantile(intensity_pred, 0.95)
    hotspots <- intensity_pred > threshold

    # Salva i risultati
    saveRDS(list(
      hotspots = hotspots,
      model = lgcp_model,
      intensity = intensity_pred
    ), file = paste0("results/lgcp_results_", basename(data_path)))

    # Notifica completamento
    send_telegram_message(paste("Completata analisi LGCP per:", data_path))

    return(hotspots)

  }, error = function(e) {
    send_telegram_message(paste("Errore nell'analisi LGCP per:", data_path, "-", e$message))
    stop(e)
  })
}

# Configurazione parallelizzazione
plan(multisession, workers = min(3, availableCores()))

# Path dei dati simulati
data_paths <- list.files("data", pattern = "simulated_.*_correlation.rds", full.names = TRUE)

# Esegui analisi in parallelo
results <- future_lapply(data_paths, analyze_lgcp, future.seed = TRUE)

# Cleanup
plan(sequential)
