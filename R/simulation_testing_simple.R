#!/usr/bin/env Rscript

# ===========================
#  simulation_testing.R
# ===========================
# Confronta un dataset simulato di trascrittomica spaziale con un dataset reale Visium (stxBrain).
# Estende i test precedenti aggiungendo una verifica della struttura di correlazione genica.

# --------------------------
# Librerie
# --------------------------
library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)  # per combinare plot

# --------------------------
# 1. Funzione di caricamento dati reali Visium
# --------------------------
load_reference_data <- function() {
  InstallData("stxBrain")   # Scarica i dati se non già presenti
  ref_obj <- LoadData("stxBrain", type = "anterior1")

  # Preprocessing semplice
  ref_obj <- SCTransform(ref_obj, assay = "Spatial",
                         variable.features.n = 2000,
                         method = "glmGamPoi",
                         verbose = FALSE)
  return(ref_obj)
}

# --------------------------
# 2. Funzione di caricamento dati simulati
# --------------------------
load_simulated_data <- function(sim_path) {
  sim_list <- readRDS(sim_path)
  # sim_list contiene coordinates, expression, ecc.
  # Creiamo un SeuratObject (usando la matrice di conteggio in formato gene x cell)
  expr_mat <- t(sim_list$expression)
  colnames(expr_mat) <- paste0("cell_", seq_len(ncol(expr_mat)))
  rownames(expr_mat) <- paste0("gene_", seq_len(nrow(expr_mat)))

  sim_obj <- CreateSeuratObject(counts = expr_mat)
  return(sim_obj)
}

# --------------------------
# 3. Confronto distribuzioni di base (già esistente)
# --------------------------
compare_expression_distributions <- function(sim_obj, ref_obj, n_genes = 40) {
  sim_counts <- GetAssayData(sim_obj, slot = "counts")
  ref_counts <- GetAssayData(ref_obj, slot = "counts")

  set.seed(123)
  sim_genes <- rownames(sim_counts)
  ref_genes <- rownames(ref_counts)

  # Selezioniamo geni da confrontare
  common_genes <- intersect(sim_genes, ref_genes)
  if(length(common_genes) >= n_genes) {
    chosen_genes <- sample(common_genes, n_genes)
  } else {
    # Non ci sono geni "veramente" in comune (naming differente). Campioniamo separatamente
    chosen_genes <- sample(sim_genes, min(n_genes, length(sim_genes)))
    ref_genes_sub <- sample(ref_genes, min(n_genes, length(ref_genes)))
  }

  calc_stats <- function(mat, gene_subset) {
    mat_sub <- mat[gene_subset, , drop = FALSE]
    mean_vals <- rowMeans(mat_sub)
    var_vals  <- apply(mat_sub, 1, var)
    cv_vals   <- sqrt(var_vals)/mean_vals
    zero_vals <- rowSums(mat_sub == 0)/ncol(mat_sub)

    tibble(
      gene      = gene_subset,
      mean_expr = mean_vals,
      var_expr  = var_vals,
      cv_expr   = cv_vals,
      zero_frac = zero_vals
    )
  }

  if (exists("ref_genes_sub")) {
    # Caso "no overlapping gene names"
    sim_stats <- calc_stats(sim_counts, chosen_genes)
    ref_stats <- calc_stats(ref_counts, ref_genes_sub)
    sim_stats$dataset <- "Simulato"
    ref_stats$dataset <- "Reale"
    # Rinominiamo i geni per evitare collisione
    sim_stats$gene <- paste0(sim_stats$gene, "_SIM")
    ref_stats$gene <- paste0(ref_stats$gene, "_REF")
    combined_stats <- bind_rows(sim_stats, ref_stats)
  } else {
    # Abbiamo chosen_genes comuni
    sim_stats <- calc_stats(sim_counts, chosen_genes) %>%
      mutate(dataset = "Simulato")
    ref_stats <- calc_stats(ref_counts, chosen_genes) %>%
      mutate(dataset = "Reale")
    combined_stats <- bind_rows(sim_stats, ref_stats)
  }

  # Plot: Media vs Varianza (log scale)
  p_mv <- ggplot(combined_stats, aes(x = mean_expr, y = var_expr, color = dataset)) +
    geom_point(alpha = 0.7) +
    scale_x_log10() + scale_y_log10() +
    theme_minimal() +
    labs(x = "Media espressione (log10)", y = "Varianza (log10)",
         title = "Confronto media-varianza")

  # Plot: Media vs Zero fraction
  p_zero <- ggplot(combined_stats, aes(x = mean_expr, y = zero_frac, color = dataset)) +
    geom_point(alpha = 0.7) +
    scale_x_log10() +
    theme_minimal() +
    labs(x = "Media espressione (log10)", y = "Frazione di zeri",
         title = "Confronto dropout (zeri) in funzione della media")

  out <- list(
    stats = combined_stats,
    plot_meanvar = p_mv,
    plot_zerofrac = p_zero
  )
  return(out)
}

# --------------------------
# 4. NUOVO: Struttura di correlazione genica
# --------------------------
compare_gene_correlations <- function(sim_obj, ref_obj, top_n = 50) {
  # Per essere sicuri di avere pattern di correlazione confrontabili,
  # scegliamo i geni più variabili in ciascun dataset (o in comune, se i nomi coincidono).

  # Identifica le feature variabili
  sim_obj <- FindVariableFeatures(sim_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  ref_obj <- FindVariableFeatures(ref_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  sim_top_genes <- head(VariableFeatures(sim_obj), top_n)
  ref_top_genes <- head(VariableFeatures(ref_obj), top_n)

  # Se i nomi dei geni simulati non hanno corrispondenza con i reali,
  # potremmo ugualmente confrontare le strutture di correlazione "internamente"
  # (cioè sim vs. ref su set di geni *diversi*), purché la lunghezza coincida.

  sim_counts <- GetAssayData(sim_obj, slot = "data")  # slot "data" = log-normalized
  ref_counts <- GetAssayData(ref_obj, slot = "data")

  # Per semplicità, se i nomi coincidono usiamo quelli in comune
  common_genes <- intersect(sim_top_genes, ref_top_genes)
  if (length(common_genes) >= 5) {
    # Calcoliamo la correlazione su geni davvero corrispondenti
    sim_mat <- as.matrix(sim_counts[common_genes, , drop = FALSE])
    ref_mat <- as.matrix(ref_counts[common_genes, , drop = FALSE])
  } else {
    # Altrimenti, prendiamo i primi top_n dal simulato e i primi top_n dal reale
    sim_mat <- as.matrix(sim_counts[sim_top_genes, , drop = FALSE])
    ref_mat <- as.matrix(ref_counts[ref_top_genes, , drop = FALSE])
    # Notare che stiamo confrontando struttura di correlazione di due insiemi genici diversi.
    # Ha senso se vogliamo solo vedere se "statisticamente" emergono pattern simili.
  }

  # Matrici di correlazione gene-gene
  cor_sim <- cor(t(sim_mat), method = "pearson", use = "pairwise.complete.obs")
  cor_ref <- cor(t(ref_mat), method = "pearson", use = "pairwise.complete.obs")

  # Creiamo un indice di somiglianza: correlazione tra i "flatten" delle due matrici
  # (solo parte triangolare superiore per non duplicare).
  upper_ind <- upper.tri(cor_sim, diag = FALSE)
  sim_flat <- cor_sim[upper_ind]
  ref_flat <- cor_ref[upper_ind]  # Assumendo abbiano dimensioni identiche

  # Se le dimensioni differiscono, bisogna gestirlo. In caso di geni diversi,
  # useremo la dimensione minima comune (già gestito con “common_genes” o con top_genes).
  # In situazioni più complesse, servirebbe un matching.
  if (!all(dim(cor_sim) == dim(cor_ref))) {
    warning("Le dimensioni delle matrici di correlazione differiscono, non si può confrontare direttamente la parte triangolare. Uso un approccio ridotto.")

    # Riduciamo le dimensioni a min(n_sim_genes, n_ref_genes)
    g_sim <- nrow(cor_sim)
    g_ref <- nrow(cor_ref)
    g_min <- min(g_sim, g_ref)
    cor_sim <- cor_sim[1:g_min, 1:g_min]
    cor_ref <- cor_ref[1:g_min, 1:g_min]
    upper_ind <- upper.tri(cor_sim, diag = FALSE)
    sim_flat <- cor_sim[upper_ind]
    ref_flat <- cor_ref[upper_ind]
  }

  # Calcoliamo la correlazione tra i vettori "flatten" di correlazioni
  cor_of_cor <- cor(sim_flat, ref_flat, method = "spearman", use = "complete.obs")

  # Creiamo due heatmap semplificate
  # Per la visualizzazione, usiamo ggplot su un melt delle matrici di correlazione
  melt_cor <- function(cor_mat, dataset_name) {
    df <- as.data.frame(as.table(cor_mat))
    colnames(df) <- c("Gene1", "Gene2", "Corr")
    df$Dataset <- dataset_name
    return(df)
  }

  sim_cor_df <- melt_cor(cor_sim, "Simulato")
  ref_cor_df <- melt_cor(cor_ref, "Reale")

  # Per evitare di disegnare un sacco di quadratini, ci limitiamo a un subset se è molto grande
  # Oppure saltiamo la maschera upper-tri per intero
  if (nrow(sim_cor_df) > 20000) {
    sim_cor_df <- sim_cor_df[1:20000, ]
  }
  if (nrow(ref_cor_df) > 20000) {
    ref_cor_df <- ref_cor_df[1:20000, ]
  }

  p_sim <- ggplot(sim_cor_df, aes(x = Gene1, y = Gene2, fill = Corr)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", limits = c(-1,1)) +
    theme_minimal(base_size = 9) +
    theme(axis.text = element_blank()) +
    labs(title = "Matrice di correlazione (Simulato)")

  p_ref <- ggplot(ref_cor_df, aes(x = Gene1, y = Gene2, fill = Corr)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", limits = c(-1,1)) +
    theme_minimal(base_size = 9) +
    theme(axis.text = element_blank()) +
    labs(title = "Matrice di correlazione (Reale)")

  # Output
  out_plot <- p_sim + p_ref
  return(list(
    correlation_of_correlations = cor_of_cor,
    cor_matrix_sim = cor_sim,
    cor_matrix_ref = cor_ref,
    plot = out_plot
  ))
}

# --------------------------
# 5. Funzione principale (estesa)
# --------------------------
run_simulation_testing <- function(sim_path = "data/simulated_image_correlation.rds",
                                   output_dir = "results/sim_tests/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Caricamento dati simulati...")
  sim_obj <- load_simulated_data(sim_path)

  message("Caricamento dati reali (Visium stxBrain)...")
  ref_obj <- load_reference_data()

  # 5a) Confronto distribuzioni di base
  message("Confronto distribuzioni di espressione...")
  comparison <- compare_expression_distributions(sim_obj, ref_obj, n_genes = 40)
  ggsave(
    filename = file.path(output_dir, "mean_var_comparison.png"),
    plot = comparison$plot_meanvar,
    width = 6, height = 5
  )
  ggsave(
    filename = file.path(output_dir, "zero_fraction_comparison.png"),
    plot = comparison$plot_zerofrac,
    width = 6, height = 5
  )
  readr::write_csv(comparison$stats, file.path(output_dir, "expression_stats.csv"))

  # 5b) NUOVO: Struttura di correlazione genica
  message("Confronto della struttura di correlazione genica...")
  cor_test <- compare_gene_correlations(sim_obj, ref_obj, top_n = 50)

  # Salviamo la correlazione tra matrici di correlazione
  cat("Correlazione tra le due matrici di correlazione (spearman): ",
      cor_test$correlation_of_correlations, "\n")

  # Salva la heatmap comparativa
  ggsave(
    filename = file.path(output_dir, "gene_correlation_comparison.png"),
    plot = cor_test$plot,
    width = 10, height = 5
  )

  # Salva il valore di "cor_of_cor" in un file di testo, per praticità
  writeLines(
    paste("Correlation of correlation matrices (spearman):", cor_test$correlation_of_correlations),
    con = file.path(output_dir, "correlation_of_correlations.txt")
  )

  message("Test completato! Risultati salvati in: ", output_dir)

  # Possiamo restituire in R la lista di risultati, se l’utente lancia la funzione da console
  return(list(
    distribution_comparison = comparison,
    correlation_test = cor_test
  ))
}

# --------------------------
# Esempio di utilizzo:
# --------------------------
# run_simulation_testing("data/simulated_easy_correlation.rds", "results/sim_tests/")
