library(tidyverse)
library(spaMM)
library(MASS)

conflicted::conflicts_prefer(dplyr::filter)
avail_thr <- parallel::detectCores(logical=FALSE) - 20L
bot <- telegram.bot::Bot(token = "1089304740:AAHA_buNHuJ2XMTHvCzrPgFaiZDApafHCWo")
chat_id <- "118334609"

bot$sendMessage(chat_id = chat_id, text = "Il bot è connesso correttamente! Analisi partita!")


set.seed(123)

n_cells <- 10000
n_genes <- 1000
hotspot_genes_prop <- 0.05   # % di geni attivi negli hotspot
hotspot_intensity <- 3       # fattore moltiplicativo nell’hotspot
p_zero_inflation  <- 0.2     # % di zeri extra
nb_size           <- 2       # parametro Negative Binomial

# 1. Generazione delle coordinate delle 10000 cellule
df_cells <- tibble(
  cell_id = 1:n_cells,
  x       = runif(n_cells, 0, 100),
  y       = runif(n_cells, 0, 100),
  # 3 hotspot (esempio)
  hotspot = (x < 20 & y < 20) | (x > 80 & y < 20) | (x < 20 & y > 80)
)

# 2. Simulazione delle espressioni per 1000 geni
#    Ciascun gene ha una probabilità "hotspot_genes_prop" di essere up-regulated nei hotspot
df_data <- map_dfr(1:n_genes, function(g) {

  # Indica se questo gene è "hotspot" o no
  is_hotspot_gene <- (runif(1) < hotspot_genes_prop)

  # Intensità media
  baseline_mu <- rexp(1, rate = 1)  # intensità casuale
  multiplier  <- if (is_hotspot_gene) hotspot_intensity else 1

  df_tmp <- df_cells %>%
    mutate(
      # A livello di log, nelle zone hotspot aumenta l'espressione
      mu = baseline_mu * if_else(hotspot, multiplier, 1.0),
      # Simulazione Negative Binomial
      counts_nb = rnbinom(n(), size = nb_size, mu = mu),
      # Aggiunta zero-inflation
      counts_zinf = if_else(runif(n()) < p_zero_inflation, 0L, counts_nb),
      gene_id = g
    ) %>%
    dplyr::select(cell_id, gene_id, x, y, hotspot, counts_zinf)

  df_tmp
})

# Il dataset finale (df_data) contiene ~10 milioni di righe (10000 cellule x 1000 geni).
# Per scopi illustrativi, ci si può limitare a un subset nelle prove preliminari.

# 3. Analisi: confronto tra (a) metodo standard vs. (b) metodo avanzato.
#    (a) Metodo standard: log(count+1) ~ coordinate, ignorando la natura di conteggio
#    (b) Metodo avanzato: Negative Binomial + correlazione spaziale (Matern).

# Esempio: focus su un singolo gene per illustrare il flusso base:
df_gene_ex <- df_data %>% filter(gene_id == 1)

# (a) Modello “standard”
fit_lm <- lm(log(counts_zinf + 1) ~ x + y, data = df_gene_ex)

# (b) Modello avanzato (Negative Binomial + Matern)
fit_spa <- fitme(
  counts_zinf ~ x + y + Matern(1 | x + y),
  data   = df_gene_ex,
  family = negbin2(),
  control.HLfit = list(NbThreads=max(avail_thr, 1L))
)

summary(fit_lm)
summary(fit_spa)

# 4. Metriche di successo (identificazione hotspot)
#    Per una valutazione globale, si iterano i modelli su tutti i geni o su un campione di geni:
#    - Si stima un punteggio di "hotspot detection" (ad es. espressione predetta nella zona)
#    - Si etichettano i geni/cellule come "hotspot" o "non-hotspot" secondo una soglia
#    - Calcolo TPR e FDR confrontando con la verità (df_data$hotspot e la designazione del gene)
#    Qui mostriamo schematicamente come potremmo procedere per un subset di geni.

evaluate_methods <- function(df_sub) {
  # Modello standard
  fit_lm <- lm(log(counts_zinf + 1) ~ x + y, data = df_sub)
  pred_lm <- predict(fit_lm, newdata = df_sub)

  # Modello spaziale
  fit_spa <- fitme(
    counts_zinf ~ x + y + Matern(1 | x + y),
    data   = df_sub,
    family = nbinom2()
  )
  pred_spa <- predict(fit_spa, newdata = df_sub)

  # Esempio: classifica come "hotspot predetto" se la predizione > mediana
  thr_lm  <- median(pred_lm)
  thr_spa <- median(pred_spa)

  df_sub <- df_sub %>%
    mutate(
      pred_hotspot_lm  = (pred_lm  > thr_lm),
      pred_hotspot_spa = (pred_spa > thr_spa)
    )

  # Calcolo TPR e FDR per i due metodi
  tpr_lm  <- mean(df_sub$pred_hotspot_lm  & df_sub$hotspot) / mean(df_sub$hotspot)
  fdr_lm  <- mean(df_sub$pred_hotspot_lm  & (!df_sub$hotspot)) / mean(df_sub$pred_hotspot_lm)

  tpr_spa <- mean(df_sub$pred_hotspot_spa & df_sub$hotspot) / mean(df_sub$hotspot)
  fdr_spa <- mean(df_sub$pred_hotspot_spa & (!df_sub$hotspot)) / mean(df_sub$pred_hotspot_spa)

  tibble(
    TPR_lm  = tpr_lm,
    FDR_lm  = fdr_lm,
    TPR_spa = tpr_spa,
    FDR_spa = fdr_spa
  )
}

# Esecuzione con gestione degli errori
tryCatch(
  {
    # Esempio di valutazione su un campione di 10 geni
    df_eval <- df_data %>%
      group_by(gene_id) %>%
      sample_n(1) %>% # per sicurezza: estrai un rigo per gene solo a scopo di illustrazione
      filter(gene_id %in% sample(1:n_genes, 10)) %>%
      group_by(gene_id) %>%
      nest() %>%
      mutate(metrics = map(data, ~ tryCatch(
        evaluate_methods(.x),
        error = function(e) tibble(error = as.character(e))
      ))) %>%
      unnest(cols = c(metrics))

    # Salvataggio di df_eval in un file .rds
    saveRDS(df_eval, "df_eval.rds")

    # Notifica di successo
    bot$sendMessage(
      chat_id = chat_id,
      text = "Il file df_eval.rds è stato salvato con successo!"
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



