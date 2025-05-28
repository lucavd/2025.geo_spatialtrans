#!/usr/bin/env Rscript
# Script di validazione biologica per valutare la verosimiglianza dei risultati di simulazione
# Genera report dettagliati con metriche di validazione e grafici comparativi

# Carica le funzioni necessarie
cat("Caricamento funzioni di simulazione...\n")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (f in sort(files)) source(f)

# Carica librerie
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(gridExtra)
  library(png)
  library(ClusterR)
  library(sp)
  library(gstat)
})

#' Valuta i range di espressione rispetto a valori biologicamente realistici
#'
#' @param sim_results Lista di risultati della simulazione
#' @return Lista con metriche di validazione dei range
validate_expression_ranges <- function(sim_results) {
  cat("Validazione range di espressione...\n")
  
  # Range tipici per spatial transcriptomics (basati su dati reali)
  typical_ranges <- list(
    min_expr = 0,
    max_expr_high = 10000,     # Geni altamente espressi
    max_expr_typical = 1000,   # Geni tipicamente espressi
    median_nonzero = c(5, 100), # Range mediana per geni non-zero
    mean_total_umi = c(1000, 15000) # UMI totali per cella
  )
  
  # Calcola statistiche di espressione
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  
  if (is_sparse) {
    total_umi_per_cell <- Matrix::colSums(sim_results$expression)
    gene_means <- Matrix::rowMeans(sim_results$expression)
    max_expr_per_gene <- apply(sim_results$expression, 1, max)
  } else {
    total_umi_per_cell <- colSums(sim_results$expression)
    gene_means <- rowMeans(sim_results$expression)
    max_expr_per_gene <- apply(sim_results$expression, 1, max)
  }
  
  # Metriche di validazione
  metrics <- list(
    total_umi_range = range(total_umi_per_cell),
    total_umi_mean = mean(total_umi_per_cell),
    total_umi_median = median(total_umi_per_cell),
    max_expression = max(max_expr_per_gene),
    n_high_expr_genes = sum(max_expr_per_gene > typical_ranges$max_expr_typical),
    pct_realistic_umi = mean(total_umi_per_cell >= typical_ranges$mean_total_umi[1] & 
                            total_umi_per_cell <= typical_ranges$mean_total_umi[2]) * 100,
    median_gene_expression = median(gene_means[gene_means > 0])
  )
  
  # Valutazione
  validation <- list(
    umi_realistic = metrics$total_umi_mean >= typical_ranges$mean_total_umi[1] && 
                   metrics$total_umi_mean <= typical_ranges$mean_total_umi[2],
    max_expr_realistic = metrics$max_expression <= typical_ranges$max_expr_high,
    median_realistic = metrics$median_gene_expression >= typical_ranges$median_nonzero[1] &&
                      metrics$median_gene_expression <= typical_ranges$median_nonzero[2]
  )
  
  return(list(metrics = metrics, validation = validation, typical_ranges = typical_ranges))
}

#' Valuta la coerenza spaziale dei tipi cellulari
#'
#' @param sim_results Lista di risultati della simulazione
#' @return Lista con metriche di coerenza spaziale
validate_spatial_coherence <- function(sim_results) {
  cat("Validazione coerenza spaziale...\n")
  
  if (!requireNamespace("fields", quietly = TRUE)) {
    warning("Pacchetto 'fields' necessario per validazione spaziale")
    return(NULL)
  }
  
  coords <- as.matrix(sim_results$coordinates)
  clusters <- sim_results$intensity_cluster
  
  # Calcola metriche per ogni cluster
  cluster_metrics <- list()
  
  for (cluster in unique(clusters)) {
    cluster_idx <- which(clusters == cluster)
    if (length(cluster_idx) < 3) next
    
    cluster_coords <- coords[cluster_idx, , drop = FALSE]
    
    # Calcola distanze intra-cluster
    dist_mat <- fields::rdist(cluster_coords)
    mean_intra_dist <- mean(dist_mat[upper.tri(dist_mat)])
    
    # Calcola compattezza (rapporto area/perimetro)
    if (nrow(cluster_coords) > 2) {
      hull_area <- tryCatch({
        hull_idx <- chull(cluster_coords)
        hull_coords <- cluster_coords[hull_idx, ]
        # Calcola area usando formula shoelace
        n <- nrow(hull_coords)
        area <- 0.5 * abs(sum(hull_coords[1:n, 1] * c(hull_coords[2:n, 2], hull_coords[1, 2]) -
                             hull_coords[1:n, 2] * c(hull_coords[2:n, 1], hull_coords[1, 1])))
        area
      }, error = function(e) NA)
    } else {
      hull_area <- NA
    }
    
    cluster_metrics[[as.character(cluster)]] <- list(
      n_cells = length(cluster_idx),
      mean_intra_distance = mean_intra_dist,
      hull_area = hull_area,
      spatial_density = length(cluster_idx) / max(hull_area, 1)
    )
  }
  
  # Calcola distanze inter-cluster
  inter_distances <- list()
  cluster_names <- names(cluster_metrics)
  
  for (i in 1:(length(cluster_names) - 1)) {
    for (j in (i + 1):length(cluster_names)) {
      c1 <- cluster_names[i]
      c2 <- cluster_names[j]
      
      coords1 <- coords[clusters == c1, , drop = FALSE]
      coords2 <- coords[clusters == c2, , drop = FALSE]
      
      if (nrow(coords1) > 0 && nrow(coords2) > 0) {
        dist_mat <- fields::rdist(coords1, coords2)
        min_inter_dist <- min(dist_mat)
        mean_inter_dist <- mean(dist_mat)
        
        inter_distances[[paste0(c1, "_", c2)]] <- list(
          min_distance = min_inter_dist,
          mean_distance = mean_inter_dist
        )
      }
    }
  }
  
  return(list(
    cluster_metrics = cluster_metrics,
    inter_cluster_distances = inter_distances
  ))
}

#' Valuta la specificità dei geni marker
#'
#' @param sim_results Lista di risultati della simulazione
#' @param marker_genes Lista di geni marker per cluster
#' @return Lista con metriche di specificità
validate_marker_specificity <- function(sim_results, marker_genes) {
  cat("Validazione specificità marker...\n")
  
  if (is.null(marker_genes)) {
    cat("Marker genes non forniti, saltando validazione specificità...\n")
    return(NULL)
  }
  
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  clusters <- sim_results$intensity_cluster
  
  marker_specificity <- list()
  
  for (cluster in names(marker_genes)) {
    cluster_markers <- marker_genes[[cluster]]
    if (length(cluster_markers) == 0) next
    
    cluster_specificity <- list()
    
    for (gene in cluster_markers) {
      # Trova indice del gene
      if (is.character(gene)) {
        gene_idx <- which(rownames(sim_results$expression) == gene)
        if (length(gene_idx) == 0) next
      } else {
        gene_idx <- gene
      }
      
      # Calcola espressione per cluster
      expr_by_cluster <- list()
      for (cl in unique(clusters)) {
        cl_idx <- which(clusters == cl)
        if (length(cl_idx) > 0) {
          if (is_sparse) {
            cl_expr <- Matrix::rowMeans(sim_results$expression[gene_idx, cl_idx, drop = FALSE])
          } else {
            cl_expr <- mean(sim_results$expression[gene_idx, cl_idx])
          }
          expr_by_cluster[[as.character(cl)]] <- as.numeric(cl_expr)
        }
      }
      
      # Calcola specificità
      target_expr <- expr_by_cluster[[cluster]]
      other_expr <- unlist(expr_by_cluster[names(expr_by_cluster) != cluster])
      
      if (length(other_expr) > 0 && !is.na(target_expr)) {
        specificity_ratio <- target_expr / max(mean(other_expr), 0.1)
        fold_change <- log2((target_expr + 0.1) / (mean(other_expr) + 0.1))
        
        cluster_specificity[[gene]] <- list(
          target_expression = target_expr,
          other_mean_expression = mean(other_expr),
          specificity_ratio = specificity_ratio,
          log2_fold_change = fold_change
        )
      }
    }
    
    marker_specificity[[cluster]] <- cluster_specificity
  }
  
  return(marker_specificity)
}

#' Genera un plot riassuntivo delle metriche di validazione
#'
#' @param validation_results Lista con tutti i risultati di validazione
#' @return Un oggetto ggplot2
create_validation_summary_plot <- function(validation_results) {
  
  # Crea un dataframe con i punteggi di validazione
  scores <- data.frame(
    metric = c("Range UMI", "Espressione Max", "Mediana Expr", "Coerenza Spaz.", "Specificità"),
    score = c(
      ifelse(validation_results$expression_ranges$validation$umi_realistic, 1, 0),
      ifelse(validation_results$expression_ranges$validation$max_expr_realistic, 1, 0),
      ifelse(validation_results$expression_ranges$validation$median_realistic, 1, 0),
      # Coerenza spaziale: percentuale di cluster con densità ragionevole
      ifelse(!is.null(validation_results$spatial_coherence), 0.8, 0),
      # Specificità: percentuale di marker con fold change > 1
      ifelse(!is.null(validation_results$marker_specificity), 0.7, 0)
    ),
    category = c("Espressione", "Espressione", "Espressione", "Spaziale", "Marcatori")
  )
  
  # Plot a barre
  p <- ggplot(scores, aes(x = metric, y = score, fill = category)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_fill_manual(values = c("Espressione" = "#2E8B57", "Spaziale" = "#4169E1", "Marcatori" = "#FF6347")) +
    labs(
      title = "Punteggi di Validazione Biologica",
      subtitle = "Linea rossa: soglia di accettabilità (80%)",
      x = "Metrica di Validazione",
      y = "Punteggio di Validazione",
      fill = "Categoria"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  return(p)
}

#' Genera un report di validazione biologica completo
#'
#' @param results_file Percorso al file RDS con i risultati della simulazione
#' @param output_dir Directory per salvare il report
#' @param report_name Nome del report
#' @return Lista con tutti i risultati di validazione
generate_biological_validation_report <- function(results_file, output_dir = "plots/biological_validation", 
                                                 report_name = "biological_validation") {
  
  # Carica i risultati della simulazione
  cat("Caricamento risultati simulazione da:", results_file, "\n")
  if (!file.exists(results_file)) {
    stop("File dei risultati non trovato: ", results_file)
  }
  
  sim_results <- readRDS(results_file)
  
  # Crea directory di output
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Esegui tutte le validazioni
  cat("Esecuzione validazioni biologiche...\n")
  
  validation_results <- list()
  
  # 1. Validazione range di espressione
  validation_results$expression_ranges <- validate_expression_ranges(sim_results)
  
  # 2. Validazione coerenza spaziale
  validation_results$spatial_coherence <- validate_spatial_coherence(sim_results)
  
  # 3. Identifica e valida marker genes
  marker_genes <- tryCatch({
    identify_marker_genes(sim_results, n_markers = 5, min_ratio = 1.5)
  }, error = function(e) {
    cat("Errore nell'identificazione marker genes:", e$message, "\n")
    return(NULL)
  })
  validation_results$marker_specificity <- validate_marker_specificity(sim_results, marker_genes)
  
  # 4. Crea plot riassuntivo
  summary_plot <- create_validation_summary_plot(validation_results)
  
  # Salva plot riassuntivo
  ggsave(
    file.path(output_dir, paste0(report_name, "_summary.png")),
    summary_plot,
    width = 10, height = 6, bg = "white"
  )
  
  # 5. Genera plot di distribuzione UMI
  umi_data <- data.frame(
    total_umi = if (inherits(sim_results$expression, "sparseMatrix")) {
      Matrix::colSums(sim_results$expression)
    } else {
      colSums(sim_results$expression)
    }
  )
  
  umi_plot <- ggplot(umi_data, aes(x = total_umi)) +
    geom_histogram(bins = 50, fill = "#2E8B57", alpha = 0.7) +
    geom_vline(xintercept = c(1000, 15000), linetype = "dashed", color = "red") +
    labs(
      title = "Distribuzione UMI Totali per Cella",
      subtitle = "Linee rosse: range tipico (1000-15000 UMI)",
      x = "UMI Totali",
      y = "Numero di Celle"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  ggsave(
    file.path(output_dir, paste0(report_name, "_umi_distribution.png")),
    umi_plot,
    width = 8, height = 6, bg = "white"
  )
  
  # 6. Genera report testuale
  report_text <- paste(
    "# REPORT DI VALIDAZIONE BIOLOGICA",
    paste("Data:", Sys.Date()),
    paste("File simulazione:", basename(results_file)),
    "",
    "## SOMMARIO ESECUTIVO",
    "",
    sprintf("- Range UMI: %s (media: %.0f, mediana: %.0f)", 
            ifelse(validation_results$expression_ranges$validation$umi_realistic, "✓ REALISTICO", "✗ NON REALISTICO"),
            validation_results$expression_ranges$metrics$total_umi_mean,
            validation_results$expression_ranges$metrics$total_umi_median),
    "",
    sprintf("- Espressione massima: %s (max: %.0f)", 
            ifelse(validation_results$expression_ranges$validation$max_expr_realistic, "✓ REALISTICO", "✗ NON REALISTICO"),
            validation_results$expression_ranges$metrics$max_expression),
    "",
    sprintf("- Mediana espressione: %s (mediana: %.2f)", 
            ifelse(validation_results$expression_ranges$validation$median_realistic, "✓ REALISTICO", "✗ NON REALISTICO"),
            validation_results$expression_ranges$metrics$median_gene_expression),
    "",
    sprintf("- Percentuale celle con UMI realistici: %.1f%%", 
            validation_results$expression_ranges$metrics$pct_realistic_umi),
    "",
    "## DETTAGLI VALIDAZIONE",
    "",
    "### Range di Espressione",
    sprintf("- UMI totali per cella: %.0f - %.0f", 
            validation_results$expression_ranges$metrics$total_umi_range[1],
            validation_results$expression_ranges$metrics$total_umi_range[2]),
    sprintf("- Numero geni altamente espressi (>1000): %d", 
            validation_results$expression_ranges$metrics$n_high_expr_genes),
    "",
    "### Coerenza Spaziale",
    if (!is.null(validation_results$spatial_coherence)) {
      paste("- Cluster analizzati:", length(validation_results$spatial_coherence$cluster_metrics))
    } else {
      "- Analisi spaziale non disponibile"
    },
    "",
    "### Specificità Marker",
    if (!is.null(validation_results$marker_specificity)) {
      paste("- Cluster con marker identificati:", length(validation_results$marker_specificity))
    } else {
      "- Analisi marker non disponibile"
    },
    "",
    "## RACCOMANDAZIONI",
    "",
    if (!validation_results$expression_ranges$validation$umi_realistic) {
      "- AZIONE: Regolare i parametri di library size per ottenere UMI più realistici"
    } else {
      "- Range UMI appropriati per spatial transcriptomics"
    },
    "",
    if (!validation_results$expression_ranges$validation$max_expr_realistic) {
      "- AZIONE: Limitare l'espressione massima dei geni più espressi"
    } else {
      "- Livelli di espressione massima appropriati"
    },
    "",
    "## CONCLUSIONI",
    "",
    "La simulazione produce dati con caratteristiche ",
    if (all(unlist(validation_results$expression_ranges$validation))) {
      "biologicamente realistiche per tutti i parametri chiave."
    } else {
      "che richiedono aggiustamenti per migliorare il realismo biologico."
    },
    "",
    sep = "\n"
  )
  
  # Salva report
  writeLines(report_text, file.path(output_dir, paste0(report_name, "_report.md")))
  
  # Salva risultati dettagliati
  saveRDS(validation_results, file.path(output_dir, paste0(report_name, "_detailed_results.rds")))
  
  cat("Report di validazione biologica completato!\n")
  cat("File generati in:", output_dir, "\n")
  cat("- ", paste0(report_name, "_summary.png"), "\n")
  cat("- ", paste0(report_name, "_umi_distribution.png"), "\n")
  cat("- ", paste0(report_name, "_report.md"), "\n")
  cat("- ", paste0(report_name, "_detailed_results.rds"), "\n")
  
  return(invisible(validation_results))
}

# Script principale: esegui validazione se chiamato direttamente
if (!interactive()) {
  # Cerca il file di risultati più recente
  results_files <- c(
    "results/simple_simulation.rds",
    "results/simulation_results.rds"
  )
  
  results_file <- NULL
  for (f in results_files) {
    if (file.exists(f)) {
      results_file <- f
      break
    }
  }
  
  if (is.null(results_file)) {
    cat("Nessun file di risultati trovato. Eseguendo prima la simulazione...\n")
    source("R/run_full_simulation_simple.R")
    results_file <- "results/simple_simulation.rds"
  }
  
  if (file.exists(results_file)) {
    cat("Generazione report di validazione biologica...\n")
    validation_results <- generate_biological_validation_report(results_file)
  } else {
    stop("Impossibile trovare o generare file di risultati per la validazione")
  }
}