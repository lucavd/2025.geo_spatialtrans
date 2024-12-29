#!/usr/bin/env Rscript

# Script per calcolare e confrontare le metriche tra i metodi su tutti i livelli di correlazione

library(tidyverse)
library(telegram.bot)

# Configura Telegram Bot
bot <- Bot(token = Sys.getenv("R_TELEGRAM_BOT_lucavdRbot"))
chat_id <- Sys.getenv("R_TELEGRAM_BOT_lucavdRbot_chatID")

# Livelli di correlazione
correlation_levels <- c("high", "medium", "low")

# Inizializza dataframe per raccogliere tutte le metriche
all_metrics <- tibble()

# Analizza ciascun livello di correlazione
for (correlation_level in correlation_levels) {

  # Carica i risultati per il livello di correlazione corrente
  seurat_results <- readRDS(paste0("results/seurat_results_simulated_", correlation_level, "_correlation.rds"))
  tweedie_results <- readRDS(paste0("results/tweedie_results_simulated_", correlation_level, "_correlation.rds"))
  gam_results <- readRDS(paste0("results/gam_geospatial_results_simulated_", correlation_level, "_correlation.rds"))
  rf_results <- readRDS(paste0("results/rf_spatial_results_simulated_", correlation_level, "_correlation.rds"))

  # True hotspots (da simulazione)
  data <- readRDS(paste0("data/simulated_", correlation_level, "_correlation.rds"))
  true_hotspots <- data$true_hotspots

  # Funzione per calcolare FP e FN
  calculate_fp_fn <- function(predicted, true) {
    fp <- sum(predicted != true)
    fn <- sum(!predicted %in% true)
    return(list(fp = fp, fn = fn))
  }

  # Calcolo delle metriche per Seurat
  seurat_predicted <- seurat_results$hotspots
  seurat_metrics <- calculate_fp_fn(seurat_predicted, true_hotspots)
  tp_seurat <- sum(seurat_predicted == true_hotspots)
  precision_seurat <- tp_seurat / (tp_seurat + seurat_metrics$fp)
  recall_seurat <- tp_seurat / (tp_seurat + seurat_metrics$fn)
  f1_seurat <- 2 * (precision_seurat * recall_seurat) / (precision_seurat + recall_seurat)

  # Calcolo delle metriche per Tweedie
  tweedie_predicted <- tweedie_results$hotspots
  tweedie_metrics <- calculate_fp_fn(tweedie_predicted, true_hotspots)
  tp_tweedie <- sum(tweedie_predicted == true_hotspots)
  precision_tweedie <- tp_tweedie / (tp_tweedie + tweedie_metrics$fp)
  recall_tweedie <- tp_tweedie / (tp_tweedie + tweedie_metrics$fn)
  f1_tweedie <- 2 * (precision_tweedie * recall_tweedie) / (precision_tweedie + recall_tweedie)

  # Calcolo delle metriche per GAM
  gam_predicted <- gam_results$hotspots
  gam_metrics <- calculate_fp_fn(gam_predicted, true_hotspots)
  tp_gam <- sum(gam_predicted == true_hotspots)
  precision_gam <- tp_gam / (tp_gam + gam_metrics$fp)
  recall_gam <- tp_gam / (tp_gam + gam_metrics$fn)
  f1_gam <- 2 * (precision_gam * recall_gam) / (precision_gam + recall_gam)

  # Calcolo delle metriche per RF
  rf_predicted <- rf_results$hotspots
  rf_metrics <- calculate_fp_fn(rf_predicted, true_hotspots)
  tp_rf <- sum(rf_predicted == true_hotspots)
  precision_rf <- tp_rf / (tp_rf + rf_metrics$fp)
  recall_rf <- tp_rf / (tp_rf + rf_metrics$fn)
  f1_rf <- 2 * (precision_rf * recall_rf) / (precision_rf + recall_rf)

  # Creazione di un dataframe di confronto per il livello corrente
  metrics <- tibble(
    Correlation_Level = correlation_level,
    Method = c("Seurat", "Tweedie", "GAM", "RF"),
    True_Positives = c(tp_seurat, tp_tweedie, tp_gam, tp_rf),
    False_Positives = c(seurat_metrics$fp, tweedie_metrics$fp, gam_metrics$fp, rf_metrics$fp),
    False_Negatives = c(seurat_metrics$fn, tweedie_metrics$fn, gam_metrics$fn, rf_metrics$fn),
    Precision = c(precision_seurat, precision_tweedie, precision_gam, precision_rf),
    Recall = c(recall_seurat, recall_tweedie, recall_gam, recall_rf),
    F1_Score = c(f1_seurat, f1_tweedie, f1_gam, f1_rf)
  )

  # Aggiungi i risultati al dataframe complessivo
  all_metrics <- bind_rows(all_metrics, metrics)
}

# Salva tutte le metriche
write_csv(all_metrics, "results/all_correlation_levels_metrics.csv")

# Invia la tabella delle metriche su Telegram
tryCatch({
  bot$send_document(chat_id = chat_id, document = "results/all_correlation_levels_metrics.csv")
}, error = function(e) {
  warning("Errore nell'invio della tabella su Telegram: ", e$message)
})

# Grafico delle metriche per tutti i livelli di correlazione
metrics_plot <- all_metrics %>%
  pivot_longer(cols = c(Precision, Recall, F1_Score), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Method, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Correlation_Level) +
  theme_minimal() +
  labs(title = "Confronto delle metriche tra i metodi per diversi livelli di correlazione", y = "Value")

ggsave("results/all_correlation_levels_metrics_plot.png", plot = metrics_plot)

# Invia il grafico delle metriche su Telegram
tryCatch({
  bot$send_photo(chat_id = chat_id, photo = "results/all_correlation_levels_metrics_plot.png")
}, error = function(e) {
  warning("Errore nell'invio del grafico delle metriche su Telegram: ", e$message)
})

# Grafico per FP e FN
errors_plot <- all_metrics %>%
  pivot_longer(cols = c(False_Positives, False_Negatives), names_to = "Error_Type", values_to = "Count") %>%
  ggplot(aes(x = Method, y = Count, fill = Error_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Correlation_Level) +
  theme_minimal() +
  labs(title = "Confronto dei falsi positivi e negativi tra i metodi per diversi livelli di correlazione", y = "Count")

ggsave("results/all_correlation_levels_fp_fn_errors.png", plot = errors_plot)

# Invia il grafico degli errori su Telegram
tryCatch({
  bot$send_photo(chat_id = chat_id, photo = "results/all_correlation_levels_fp_fn_errors.png")
}, error = function(e) {
  warning("Errore nell'invio del grafico degli errori su Telegram: ", e$message)
})
