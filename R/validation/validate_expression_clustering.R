#' Validazione clustering basato su espressione
#'
#' Valida la qualità del clustering basato su profili di espressione
#' analizzando separazione, coerenza spaziale e qualità biologica
#'
#' @param clustering_result Risultato di expression_based_clustering()
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param output_dir Directory per salvare i risultati
#' @param prefix Prefisso per i file di output
#' @return Lista con metriche di validazione
#' @export
validate_expression_clustering <- function(
  clustering_result,
  cell_df = NULL,
  output_dir = "R/validation",
  prefix = "expression_clustering"
) {
  
  # Setup
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (is.null(cell_df)) {
    cell_df <- clustering_result$cell_df
  }
  
  pca_coords <- clustering_result$pca_coords
  clusters <- clustering_result$clusters
  expr_matrix <- clustering_result$expression_matrix
  
  cat("=== VALIDAZIONE CLUSTERING BASATO SU ESPRESSIONE ===\n")
  
  # 1. Metriche di clustering generale
  n_clusters <- length(unique(clusters))
  n_cells <- nrow(cell_df)
  cluster_sizes <- table(clusters)
  
  cat("Numero di cluster:", n_clusters, "\n")
  cat("Numero di celle:", n_cells, "\n")
  cat("Dimensioni cluster (min/max):", min(cluster_sizes), "/", max(cluster_sizes), "\n")
  
  # 2. Silhouette score su spazio PCA
  cat("Calcolo silhouette score...\n")
  silhouette_score <- calculate_silhouette_score(pca_coords, clusters)
  
  # 3. Coerenza spaziale
  cat("Calcolo coerenza spaziale...\n")
  spatial_coherence <- calculate_spatial_coherence(cell_df, clusters)
  
  # 4. Separazione nell'espressione
  cat("Calcolo separazione espressione...\n")
  expression_separation <- calculate_expression_separation(expr_matrix, clusters)
  
  # 5. Marker gene enrichment
  cat("Calcolo marker gene enrichment...\n")
  marker_enrichment <- calculate_marker_enrichment(expr_matrix, clusters)
  
  # 6. Visualizzazioni
  cat("Generazione visualizzazioni...\n")
  
  # Plot PCA con cluster
  plot_pca_clusters(pca_coords, clusters, 
                   file.path(output_dir, paste0(prefix, "_pca_clusters.png")))
  
  # Plot spaziale
  plot_spatial_clusters(cell_df, 
                       file.path(output_dir, paste0(prefix, "_spatial_clusters.png")))
  
  # Plot silhouette
  plot_silhouette_analysis(pca_coords, clusters,
                          file.path(output_dir, paste0(prefix, "_silhouette.png")))
  
  # Plot marker heatmap (temporaneamente disabilitato per debug)
  tryCatch({
    plot_marker_heatmap(expr_matrix, clusters, marker_enrichment$top_markers,
                       file.path(output_dir, paste0(prefix, "_marker_heatmap.png")))
  }, error = function(e) {
    cat("Warning: Marker heatmap skipped due to error:", e$message, "\n")
  })
  
  # 7. Compila risultati
  validation_results <- list(
    n_clusters = n_clusters,
    n_cells = n_cells,
    cluster_sizes = cluster_sizes,
    silhouette_score = silhouette_score,
    spatial_coherence = spatial_coherence,
    expression_separation = expression_separation,
    marker_enrichment = marker_enrichment,
    timestamp = Sys.time()
  )
  
  # Salva risultati
  saveRDS(validation_results, file.path(output_dir, paste0(prefix, "_validation_results.rds")))
  
  # Report riassuntivo
  generate_validation_report(validation_results, output_dir, prefix)
  
  cat("Validazione completata. Risultati salvati in:", output_dir, "\n")
  
  return(validation_results)
}

#' Calcola silhouette score
calculate_silhouette_score <- function(pca_coords, clusters) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    warning("Pacchetto 'cluster' non disponibile")
    return(NA)
  }
  
  if (length(unique(clusters)) < 2) {
    return(NA)
  }
  
  # Calcola distanze
  dist_matrix <- dist(pca_coords)
  
  # Silhouette analysis
  sil <- cluster::silhouette(clusters, dist_matrix)
  
  list(
    average = mean(sil[, "sil_width"]),
    per_cluster = tapply(sil[, "sil_width"], clusters, mean),
    silhouette_matrix = sil
  )
}

#' Calcola coerenza spaziale
calculate_spatial_coherence <- function(cell_df, clusters) {
  # Per ogni cluster, calcola quanto sono spazialmente compatti
  coherence_scores <- numeric(length(unique(clusters)))
  names(coherence_scores) <- unique(clusters)
  
  for (cl in unique(clusters)) {
    cl_cells <- cell_df[clusters == cl, ]
    
    if (nrow(cl_cells) < 3) {
      coherence_scores[as.character(cl)] <- NA
      next
    }
    
    # Calcola centroide
    centroid <- c(mean(cl_cells$x), mean(cl_cells$y))
    
    # Distanza media dal centroide
    distances <- sqrt((cl_cells$x - centroid[1])^2 + (cl_cells$y - centroid[2])^2)
    mean_dist <- mean(distances)
    
    # Normalizza per area occupata
    area <- (max(cl_cells$x) - min(cl_cells$x)) * (max(cl_cells$y) - min(cl_cells$y))
    coherence_scores[as.character(cl)] <- mean_dist / sqrt(area)
  }
  
  list(
    per_cluster = coherence_scores,
    average = mean(coherence_scores, na.rm = TRUE)
  )
}

#' Calcola separazione nell'espressione
calculate_expression_separation <- function(expr_matrix, clusters) {
  # Assicurati che la matrice sia nel formato corretto (geni x celle)
  if (ncol(expr_matrix) < nrow(expr_matrix)) {
    expr_matrix <- t(expr_matrix)
  }
  
  # Between-cluster sum of squares vs within-cluster sum of squares
  unique_clusters <- unique(clusters)
  cluster_means <- matrix(0, nrow = nrow(expr_matrix), ncol = length(unique_clusters))
  colnames(cluster_means) <- as.character(unique_clusters)
  
  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    cl_indices <- which(clusters == cl)
    if (length(cl_indices) == 1) {
      cluster_means[, i] <- expr_matrix[, cl_indices]
    } else {
      cluster_means[, i] <- rowMeans(expr_matrix[, cl_indices, drop = FALSE])
    }
  }
  
  global_mean <- rowMeans(expr_matrix)
  
  # Between-cluster SS
  between_ss <- 0
  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    n_cl <- sum(clusters == cl)
    cl_mean <- cluster_means[, i]
    between_ss <- between_ss + n_cl * sum((cl_mean - global_mean)^2)
  }
  
  # Within-cluster SS
  within_ss <- 0
  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    cl_indices <- which(clusters == cl)
    cl_mean <- cluster_means[, i]
    for (j in cl_indices) {
      within_ss <- within_ss + sum((expr_matrix[, j] - cl_mean)^2)
    }
  }
  
  # F-statistic like measure
  f_stat <- (between_ss / (length(unique(clusters)) - 1)) / (within_ss / (ncol(expr_matrix) - length(unique(clusters))))
  
  list(
    between_ss = between_ss,
    within_ss = within_ss,
    f_statistic = f_stat,
    separation_ratio = between_ss / (between_ss + within_ss)
  )
}

#' Calcola arricchimento marker genes
calculate_marker_enrichment <- function(expr_matrix, clusters) {
  # Assicurati che la matrice sia nel formato corretto
  if (ncol(expr_matrix) < nrow(expr_matrix)) {
    expr_matrix <- t(expr_matrix)
  }
  
  n_genes <- nrow(expr_matrix)
  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene_", 1:n_genes)
  }
  
  top_markers <- list()
  
  for (cl in unique(clusters)) {
    cl_indices <- which(clusters == cl)
    other_indices <- which(clusters != cl)
    
    if (length(cl_indices) < 3 || length(other_indices) < 3) {
      top_markers[[as.character(cl)]] <- character(0)
      next
    }
    
    # Calcola fold change
    cl_means <- rowMeans(expr_matrix[, cl_indices, drop = FALSE])
    other_means <- rowMeans(expr_matrix[, other_indices, drop = FALSE])
    
    # Evita divisione per zero
    other_means[other_means == 0] <- 1e-6
    fold_changes <- cl_means / other_means
    
    # Seleziona top marker (alto fold change)
    top_indices <- order(fold_changes, decreasing = TRUE)[1:min(10, length(fold_changes))]
    top_markers[[as.character(cl)]] <- gene_names[top_indices]
  }
  
  list(
    top_markers = top_markers,
    n_markers_per_cluster = sapply(top_markers, length)
  )
}

#' Plot PCA con cluster
plot_pca_clusters <- function(pca_coords, clusters, output_file) {
  library(ggplot2)
  
  df <- data.frame(
    PC1 = pca_coords[, 1],
    PC2 = pca_coords[, 2],
    Cluster = factor(clusters)
  )
  
  p <- ggplot(df, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 1, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Clustering su spazio PCA",
         x = "Prima componente principale",
         y = "Seconda componente principale") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(output_file, p, width = 8, height = 6, bg = "white")
}

#' Plot cluster spaziali
plot_spatial_clusters <- function(cell_df, output_file) {
  library(ggplot2)
  
  p <- ggplot(cell_df, aes(x = x, y = y, color = intensity_cluster)) +
    geom_point(size = 1, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Cluster spaziali",
         x = "Coordinata X", 
         y = "Coordinata Y",
         color = "Cluster") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(output_file, p, width = 8, height = 6, bg = "white")
}

#' Plot silhouette analysis
plot_silhouette_analysis <- function(pca_coords, clusters, output_file) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    return(NULL)
  }
  
  library(ggplot2)
  
  dist_matrix <- dist(pca_coords)
  sil <- cluster::silhouette(clusters, dist_matrix)
  
  # Converti in dataframe per ggplot
  sil_df <- data.frame(
    cluster = factor(sil[, "cluster"]),
    neighbor = sil[, "neighbor"],
    sil_width = sil[, "sil_width"]
  )
  
  sil_df$index <- 1:nrow(sil_df)
  
  p <- ggplot(sil_df, aes(x = index, y = sil_width, fill = cluster)) +
    geom_col() +
    facet_wrap(~cluster, scales = "free_x") +
    theme_minimal() +
    labs(title = paste("Silhouette Analysis (avg =", round(mean(sil[, "sil_width"]), 3), ")"),
         x = "Indice cella", 
         y = "Silhouette width") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(output_file, p, width = 10, height = 6, bg = "white")
}

#' Plot heatmap marker genes
plot_marker_heatmap <- function(expr_matrix, clusters, top_markers, output_file) {
  library(ggplot2)
  library(dplyr)
  
  # Assicurati formato corretto
  if (ncol(expr_matrix) < nrow(expr_matrix)) {
    expr_matrix <- t(expr_matrix)
  }
  
  # Seleziona subset di marker genes
  all_markers <- unique(unlist(top_markers))
  if (length(all_markers) == 0) {
    return(NULL)
  }
  
  # Limita a massimo 50 marker per visualizzazione
  if (length(all_markers) > 50) {
    all_markers <- all_markers[1:50]
  }
  
  # Subset matrice
  marker_indices <- which(rownames(expr_matrix) %in% all_markers)
  if (length(marker_indices) == 0) {
    marker_indices <- 1:min(50, nrow(expr_matrix))
  }
  
  marker_matrix <- expr_matrix[marker_indices, , drop = FALSE]
  
  # Calcola medie per cluster
  unique_clusters <- unique(clusters)
  cluster_means <- matrix(0, nrow = nrow(marker_matrix), ncol = length(unique_clusters))
  colnames(cluster_means) <- as.character(unique_clusters)
  
  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    cl_indices <- which(clusters == cl)
    if (length(cl_indices) == 1) {
      cluster_means[, i] <- marker_matrix[, cl_indices]
    } else {
      cluster_means[, i] <- rowMeans(marker_matrix[, cl_indices, drop = FALSE])
    }
  }
  
  # Converti per ggplot  
  if (length(unique_clusters) > 0 && nrow(marker_matrix) > 0) {
    heatmap_df <- expand.grid(
      Gene = rownames(marker_matrix),
      Cluster = colnames(cluster_means),
      stringsAsFactors = FALSE
    )
    heatmap_df$Expression <- as.vector(cluster_means)
  } else {
    # Se non ci sono marker, crea un dataframe vuoto
    heatmap_df <- data.frame(
      Gene = character(0),
      Cluster = character(0),
      Expression = numeric(0)
    )
  }
  
  if (nrow(heatmap_df) > 0) {
    p <- ggplot(heatmap_df, aes(x = Cluster, y = Gene, fill = Expression)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      labs(title = "Heatmap Marker Genes",
           x = "Cluster", 
           y = "Gene") +
      theme(axis.text.y = element_text(size = 8),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA))
  } else {
    # Plot vuoto se non ci sono dati
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No marker genes found", size = 5) +
      theme_minimal() +
      labs(title = "Heatmap Marker Genes - No Data",
           x = "Cluster", 
           y = "Gene") +
      theme(panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA))
  }
  
  ggsave(output_file, p, width = 8, height = 10, bg = "white")
}

#' Genera report di validazione
generate_validation_report <- function(validation_results, output_dir, prefix) {
  report_file <- file.path(output_dir, paste0(prefix, "_validation_report.md"))
  
  report_lines <- c(
    "# Validation Report - Expression-based Clustering",
    "",
    paste("**Timestamp:**", validation_results$timestamp),
    "",
    "## Summary Statistics",
    "",
    paste("- **Number of clusters:**", validation_results$n_clusters),
    paste("- **Number of cells:**", validation_results$n_cells),
    paste("- **Average silhouette score:**", round(validation_results$silhouette_score$average, 3)),
    paste("- **Average spatial coherence:**", round(validation_results$spatial_coherence$average, 3)),
    paste("- **Expression separation ratio:**", round(validation_results$expression_separation$separation_ratio, 3)),
    "",
    "## Cluster Sizes",
    "",
    paste(names(validation_results$cluster_sizes), ":", validation_results$cluster_sizes, collapse = "\n"),
    "",
    "## Quality Metrics",
    "",
    "### Silhouette Score by Cluster",
    "",
    paste(names(validation_results$silhouette_score$per_cluster), ":", 
          round(validation_results$silhouette_score$per_cluster, 3), collapse = "\n"),
    "",
    "### Spatial Coherence by Cluster",
    "",
    paste(names(validation_results$spatial_coherence$per_cluster), ":", 
          round(validation_results$spatial_coherence$per_cluster, 3), collapse = "\n"),
    "",
    "## Files Generated",
    "",
    paste0("- `", prefix, "_pca_clusters.png` - PCA visualization with clusters"),
    paste0("- `", prefix, "_spatial_clusters.png` - Spatial distribution of clusters"),
    paste0("- `", prefix, "_silhouette.png` - Silhouette analysis"),
    paste0("- `", prefix, "_marker_heatmap.png` - Marker gene heatmap"),
    paste0("- `", prefix, "_validation_results.rds` - Complete validation results"),
    ""
  )
  
  writeLines(report_lines, report_file)
}