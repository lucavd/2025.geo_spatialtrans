library(dplyr)

df_eval_all %>%
  # In caso ci fossero righe con errore, le escludiamo
  filter(!("error" %in% names(df_eval_all)) | is.na(error)) %>%
  summarise(
    TPR_lm_media  = mean(TPR_lm, na.rm = TRUE),
    FDR_lm_media  = mean(FDR_lm, na.rm = TRUE),
    TPR_spa_media = mean(TPR_spa, na.rm = TRUE),
    FDR_spa_media = mean(FDR_spa, na.rm = TRUE)
  )

  # In genere, più alta TPR e più bassa FDR indicano prestazioni migliori.

# Comparazione per Gene
# contare quanti geni risultano “migliori” per l’uno o l’altro metodo

df_eval_compare <- df_eval_all %>%
  filter(!("error" %in% names(df_eval_all)) | is.na(error)) %>%
  mutate(
    better_TPR = case_when(
      TPR_spa > TPR_lm ~ "spa",
      TPR_spa < TPR_lm ~ "lm",
      TRUE             ~ "tie"
    ),
    better_FDR = case_when(
      FDR_spa < FDR_lm ~ "spa",
      FDR_spa > FDR_lm ~ "lm",
      TRUE             ~ "tie"
    )
  )

head(df_eval_compare)

df_eval_compare %>%
  count(better_TPR)

df_eval_compare %>%
  count(better_FDR)


library(ggplot2)

# Scatter plot TPR
ggplot(df_eval_compare, aes(x = TPR_lm, y = TPR_spa)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Confronto TPR: LM vs. SPA",
       x = "TPR metodo LM",
       y = "TPR metodo SPA")

# Scatter plot FDR
ggplot(df_eval_compare, aes(x = FDR_lm, y = FDR_spa)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Confronto FDR: LM vs. SPA",
       x = "FDR metodo LM",
       y = "FDR metodo SPA")

# Dai due grafici si possono trarre alcune considerazioni generali:
#
#   1. **TPR (True Positive Rate)**
#   - Il grafico “Confronto TPR: LM vs. SPA” mostra per ogni gene (un punto nel plot) la coppia \((\text{TPR}_\text{LM}, \text{TPR}_\text{SPA})\).
# - Se un punto sta **sopra** la linea rossa diagonale (cioè \(\text{TPR}_\text{SPA} > \text{TPR}_\text{LM}\)), significa che il metodo SPA (modello spaziale) ha una TPR più alta di LM per quel gene.
# - Nel tuo grafico, la maggior parte dei punti appare **sopra** la diagonale, indicando che spesso SPA raggiunge una TPR maggiore rispetto a LM.
#
# 2. **FDR (False Discovery Rate)**
#   - Il secondo grafico riporta \((\text{FDR}_\text{LM}, \text{FDR}_\text{SPA})\).
# - Se un punto sta **sopra** la diagonale (\(\text{FDR}_\text{SPA} > \text{FDR}_\text{LM}\)), il metodo SPA ha un FDR più alto, quindi introduce più “falsi positivi” in proporzione.
# - Nel grafico, molti punti sono effettivamente **sopra** la diagonale, quindi in media SPA tende ad avere FDR superiore a LM.
#
# In altre parole:
#   - **SPA** sembrerebbe identificare più veri positivi (TPR più alta) ma, in parallelo, potrebbe marcare anche più falsi positivi (FDR più alto).
# - **LM** cattura meno veri positivi, ma tende ad avere un FDR più contenuto.
#
# Questo riflette il classico **trade-off**: se un metodo è più sensibile (TPR), spesso paga lo scotto di avere più falsi positivi (FDR). La “bontà” di un metodo dipende dunque da cosa è più importante per la tua analisi:
#   - **Se vuoi massimizzare la sensibilità** (meglio non perdere geni veri-positivi), potresti preferire SPA (magari cercando di regolare la soglia di predizione).
# - **Se vuoi tenere bassa la FDR** (pochi falsi positivi), allora LM potrebbe essere più sicuro, a discapito della sensibilità.
#
# Puoi anche calcolare e confrontare altre metriche (ad es. F1-score, precision, AUC, ecc.) oppure variare la soglia di classificazione (ad esempio anziché usare la **mediana** come threshold, valutare soglie diverse) per vedere se è possibile trovare un compromesso migliore tra TPR e FDR.
