##############################
#  LIBRERIE E CONFIGURAZIONI
##############################
library(data.table)     # Per gestire grandi tabelle più velocemente
library(tidyverse)      # Per funzioni ausiliarie (se le preferisci)
library(spaMM)          # Per il modello Negative Binomial + Matern
library(MASS)           # Per la generazione di dati NB
library(furrr)          # Per la parallelizzazione via purrr
library(telegram.bot)   # Per notifiche Telegram
library(yardstick)

# Per sicurezza:
conflicted::conflicts_prefer(dplyr::filter)

# Imposta parallelizzazione
avail_thr <- parallel::detectCores(logical=FALSE) - 20L
plan(multisession, workers = max(avail_thr, 1L))

# Imposta il bot Telegram
bot <- Bot(token = "R_TELEGRAM_BOT_lucavdRbot")
chat_id <- "118334609"


##############################
#  NOTIFICA DI INIZIO
##############################
bot$sendMessage(
  chat_id = chat_id,
  text = "Il bot è connesso correttamente! Analisi partita!"
)


##############################
#  SIMULAZIONE DEI DATI
##############################
set.seed(3011)

# Parametri ridotti per velocizzare
n_cells  <- 2000    # ridotto da 10000
n_genes  <- 200     # ridotto da 1000

hotspot_genes_prop <- 0.05   # % di geni attivi negli hotspot
hotspot_intensity  <- 3      # fattore moltiplicativo nell’hotspot
p_zero_inflation   <- 0.2    # % di zeri extra
nb_size            <- 2      # parametro Negative Binomial

# 1. Generazione delle coordinate delle 2000 cellule
df_cells <- tibble(
  cell_id = 1:n_cells,
  x       = runif(n_cells, 0, 100),
  y       = runif(n_cells, 0, 100),
  # 3 hotspot (esempio)
  hotspot = (x < 20 & y < 20) | (x > 80 & y < 20) | (x < 20 & y > 80)
)

# 2. Simulazione delle espressioni per 200 geni
df_data <- map_dfr(1:n_genes, function(g) {

  # Indica se questo gene è "hotspot" o no
  is_hotspot_gene <- (runif(1) < hotspot_genes_prop)

  # Intensità media (baseline)
  baseline_mu <- rexp(1, rate = 1)  # intensità casuale
  multiplier  <- if (is_hotspot_gene) hotspot_intensity else 1

  # Creiamo un data.frame/tibble per questo gene
  df_tmp <- df_cells %>%
    mutate(
      mu          = baseline_mu * if_else(hotspot, multiplier, 1.0),
      counts_nb   = rnbinom(n(), size = nb_size, mu = mu),
      counts_zinf = if_else(runif(n()) < p_zero_inflation, 0L, counts_nb),
      gene_id     = g
    ) %>%
    dplyr::select(cell_id, gene_id, x, y, hotspot, counts_zinf)

  df_tmp
})

# Convertiamo in data.table per velocizzare le operazioni
DT_data <- as.data.table(df_data)
setkey(DT_data, gene_id)  # chiave su gene_id per subset più rapido


##############################
#  FUNZIONI DI VALUTAZIONE
##############################
evaluate_methods <- function(df_sub) {
  # Modello standard (LM su log(count+1))
  fit_lm <- lm(log(counts_zinf + 1) ~ x + y, data = df_sub)
  pred_lm <- predict(fit_lm, newdata = df_sub)

  # Modello spaziale NB + Matern
  fit_spa <- fitme(
    counts_zinf ~ x + y + Matern(1 | x + y),
    data   = df_sub,
    family = negbin2()
  )
  pred_spa <- predict(fit_spa, newdata = df_sub)

  # -----------------------------
  # 1) Calcolo TPR, FDR, F1 usando la SOGGETTIVA soglia mediana
  # -----------------------------
  thr_lm  <- median(pred_lm)
  thr_spa <- median(pred_spa)

  df_class <- df_sub %>%
    mutate(
      pred_hotspot_lm  = (pred_lm  > thr_lm),
      pred_hotspot_spa = (pred_spa > thr_spa)
    )

  # Calcolo TPR (recall) e FDR -> precision = 1 - FDR
  tpr_lm  <- mean(df_class$pred_hotspot_lm  & df_class$hotspot) / mean(df_class$hotspot)
  fdr_lm  <- mean(df_class$pred_hotspot_lm  & (!df_class$hotspot)) / mean(df_class$pred_hotspot_lm)
  # F1 = 2 * (precision * recall) / (precision + recall)
  # precision_lm = 1 - fdr_lm
  precision_lm <- 1 - fdr_lm
  f1_lm <- 2 * (precision_lm * tpr_lm) / (precision_lm + tpr_lm)

  tpr_spa <- mean(df_class$pred_hotspot_spa & df_class$hotspot) / mean(df_class$hotspot)
  fdr_spa <- mean(df_class$pred_hotspot_spa & (!df_class$hotspot)) / mean(df_class$pred_hotspot_spa)
  precision_spa <- 1 - fdr_spa
  f1_spa <- 2 * (precision_spa * tpr_spa) / (precision_spa + tpr_spa)

  # -----------------------------
  # 2) Calcolo AUC (ROC e PR) con yardstick
  #    Per farlo correttamente occorrono i dati "completi":
  #    - etichetta binaria (hotspot) come fattore
  #    - predizione continua
  # -----------------------------
  df_preds <- df_sub %>%
    mutate(
      hotspot_factor = factor(hotspot, levels = c(FALSE, TRUE)),  # verità
      .pred_lm       = pred_lm,
      .pred_spa      = pred_spa
    )

  # ROC-AUC e PR-AUC per LM
  # Nota: yardstick vuole la sintassi `roc_auc(data, truth, .pred)`.
  roc_lm <- roc_auc(df_preds, truth = hotspot_factor, .pred_lm)
  pr_lm  <- pr_auc(df_preds, truth = hotspot_factor, .pred_lm)

  # ROC-AUC e PR-AUC per SPA
  roc_spa <- roc_auc(df_preds, truth = hotspot_factor, .pred_spa)
  pr_spa  <- pr_auc(df_preds, truth = hotspot_factor, .pred_spa)

  # Convertiamo i risultati yardstick in vettori
  auc_roc_lm <- roc_lm$.estimate
  auc_pr_lm  <- pr_lm$.estimate
  auc_roc_spa <- roc_spa$.estimate
  auc_pr_spa  <- pr_spa$.estimate

  # -----------------------------
  # 3) Restituiamo un tibble con tutte le metriche
  # -----------------------------
  tibble(
    TPR_lm  = tpr_lm,
    FDR_lm  = fdr_lm,
    F1_lm   = f1_lm,
    AUC_ROC_lm = auc_roc_lm,
    AUC_PR_lm  = auc_pr_lm,

    TPR_spa  = tpr_spa,
    FDR_spa  = fdr_spa,
    F1_spa   = f1_spa,
    AUC_ROC_spa = auc_roc_spa,
    AUC_PR_spa  = auc_pr_spa
  )
}



##############################
#  ANALISI A BLOCCHI (CHUNK)
##############################
# Supponiamo di voler dividere i 200 geni in chunk di 50
genes_chunks <- split(unique(DT_data$gene_id),
                      ceiling(seq_along(unique(DT_data$gene_id)) / 50))

# Creiamo un oggetto per accumulare i risultati
df_eval_all <- NULL

# Cicliamo su ciascun chunk (per evitare di caricare troppi dati in RAM alla volta)
for (ich in seq_along(genes_chunks)) {
  chunk_genes <- genes_chunks[[ich]]

  # Estrai dal data.table solo i geni di questo chunk
  DT_chunk <- DT_data[J(chunk_genes), nomatch = 0]  # subset veloce grazie a setkey

  # Suddividiamo la tabella in liste di subset (uno per gene)
  list_subsets <- split(DT_chunk, by = "gene_id")

  # Possiamo parallelizzare la mappatura su list_subsets:
  chunk_res <- future_map_dfr(list_subsets, function(df_sub) {
    tryCatch(
      evaluate_methods(df_sub),
      error = function(e) tibble(error = as.character(e))
    )
  })

  # Aggiungiamo colonna con il gene_id (se manca)
  # In evaluate_methods non lo aggiungiamo, quindi lo aggiungiamo adesso
  chunk_res <- chunk_res %>%
    mutate(gene_id = as.numeric(names(list_subsets)))

  # Accumula i risultati
  df_eval_all <- bind_rows(df_eval_all, chunk_res)

  # Messaggio di avanzamento su console (o su Telegram se vuoi)
  cat(sprintf("Chunk %d di %d completato. Geni: %s\n",
              ich, length(genes_chunks),
              paste(chunk_genes, collapse = ", ")))
}


##############################
#  SALVATAGGIO RISULTATI
##############################
tryCatch(
  {
    saveRDS(df_eval_all, "df_eval_all.rds")

    # Notifica di successo
    bot$sendMessage(
      chat_id = chat_id,
      text = "Il file df_eval_all.rds è stato salvato con successo!"
    )
  },
  error = function(e) {
    # Notifica di errore
    bot$sendMessage(
      chat_id = chat_id,
      text = paste("Errore durante il processo:", e$message)
    )
  }
)

# Fine script
