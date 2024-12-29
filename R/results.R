#!/usr/bin/env Rscript
# Script per calcolare e confrontare le metriche tra i metodi su tutti i livelli di correlazione

library(tidyverse)
library(telegram.bot)

# Configura Telegram Bot
bot <- Bot(token = Sys.getenv("R_TELEGRAM_BOT_lucavdRbot"))
chat_id <- Sys.getenv("R_TELEGRAM_BOT_lucavdRbot_chatID")

# Funzione per inviare messaggi e file su Telegram
send_telegram <- function(content, type = "message") {
  tryCatch({
    if (type == "message") {
      bot$send_message(chat_id = chat_id, text = content)
    } else if (type == "document") {
      bot$send_document(chat_id = chat_id, document = content)
    } else if (type == "photo") {
      bot$send_photo(chat_id = chat_id, photo = content)
    }
  }, error = function(e) {
    warning("Errore nell'invio su Telegram: ", e$message)
  })
}

# Funzione per calcolare FP e FN
calculate_fp_fn <- function(predicted, true) {
  fp <- sum(predicted != true)
  fn <- sum(!predicted %in% true)
  return(list(fp = fp, fn = fn))
}

# Funzione per il tema personalizzato
custom_theme <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(color = "black", face = "bold")
    )
}

# Inizializza dataframe per le metriche
all_metrics <- tibble()

# Analizza ciascun livello di correlazione
for(correlation_level in c("low", "medium", "high")) {
  # Carica i dati
  data <- readRDS(paste0("data/simulated_", correlation_level, "_correlation.rds"))
  true_hotspots <- data$true_hotspots

  # Carica i risultati per ogni metodo
  seurat_results <- readRDS(paste0("results/seurat_results_simulated_", correlation_level, "_correlation.rds"))
  tweedie_results <- readRDS(paste0("results/tweedie_results_simulated_", correlation_level, "_correlation.rds"))
  gam_results <- readRDS(paste0("results/gam_geospatial_results_simulated_", correlation_level, "_correlation.rds"))
  rf_results <- readRDS(paste0("results/ranger_results_simulated_", correlation_level, "_correlation.rds"))
  gnn_results <- readRDS(paste0("results/gnn_results_simulated_", correlation_level, "_correlation.rds"))

  # Lista dei metodi e risultati
  methods_list <- list(
    Seurat = seurat_results$hotspots,
    Tweedie = tweedie_results$hotspots,
    GAM = gam_results$hotspots,
    RF = rf_results$hotspots,
    GNN = gnn_results$hotspots
  )

  # Calcola metriche per ogni metodo
  metrics_list <- lapply(names(methods_list), function(method) {
    predicted <- methods_list[[method]]
    metrics <- calculate_fp_fn(predicted, true_hotspots)
    tp <- sum(predicted == true_hotspots)
    precision <- tp / (tp + metrics$fp)
    recall <- tp / (tp + metrics$fn)
    f1 <- 2 * (precision * recall) / (precision + recall)

    tibble(
      Correlation_Level = correlation_level,
      Method = method,
      True_Positives = tp,
      False_Positives = metrics$fp,
      False_Negatives = metrics$fn,
      Precision = precision,
      Recall = recall,
      F1_Score = f1
    )
  })

  # Combina i risultati
  all_metrics <- bind_rows(all_metrics, bind_rows(metrics_list))
}

# Colori personalizzati
custom_colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F")

# Grafico delle metriche principali
metrics_plot <- all_metrics %>%
  pivot_longer(cols = c(Precision, Recall, F1_Score), names_to = "Metric", values_to = "Value") %>%
  mutate(Correlation_Level = factor(Correlation_Level, levels = c("low", "medium", "high"))) %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Metric ~ Correlation_Level) +
  scale_fill_manual(values = custom_colors) +
  custom_theme() +
  labs(
    title = "Confronto delle metriche tra i metodi",
    y = "Valore",
    x = "Metodo"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Salva il grafico delle metriche
ggsave(
  "results/all_correlation_levels_metrics_plot.png",
  plot = metrics_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Grafico per errori (FP e FN)
errors_plot <- all_metrics %>%
  pivot_longer(
    cols = c(False_Positives, False_Negatives),
    names_to = "Error_Type",
    values_to = "Count"
  ) %>%
  mutate(Correlation_Level = factor(Correlation_Level, levels = c("low", "medium", "high"))) %>%
  ggplot(aes(x = Method, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Error_Type ~ Correlation_Level) +
  scale_fill_manual(values = custom_colors) +
  custom_theme() +
  labs(
    title = "Confronto dei falsi positivi e negativi tra i metodi",
    y = "Conteggio",
    x = "Metodo"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Salva il grafico degli errori
ggsave(
  "results/all_correlation_levels_fp_fn_errors.png",
  plot = errors_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Grafico di confronto F1-Score
f1_plot <- all_metrics %>%
  mutate(Correlation_Level = factor(Correlation_Level, levels = c("low", "medium", "high"))) %>%
  ggplot(aes(x = Correlation_Level, y = F1_Score, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  custom_theme() +
  labs(
    title = "Confronto F1-Score per livello di correlazione",
    y = "F1-Score",
    x = "Livello di Correlazione"
  )

# Salva il grafico F1-Score
ggsave(
  "results/f1_score_comparison.png",
  plot = f1_plot,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Salva le metriche in CSV
write_csv(all_metrics, "results/all_correlation_levels_metrics.csv")

# Invia risultati su Telegram
send_telegram("Analisi completata! Invio risultati...")
send_telegram("results/all_correlation_levels_metrics.csv", "document")
send_telegram("results/all_correlation_levels_metrics_plot.png", "photo")
send_telegram("results/all_correlation_levels_fp_fn_errors.png", "photo")
send_telegram("results/f1_score_comparison.png", "photo")

# Return list of results for targets
list(
  metrics = all_metrics,
  plots = list(
    metrics = "results/all_correlation_levels_metrics_plot.png",
    errors = "results/all_correlation_levels_fp_fn_errors.png",
    f1_comparison = "results/f1_score_comparison.png"
  )
)
