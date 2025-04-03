#!/usr/bin/env Rscript

# Script: analyze_and_compare_clusters.R
# Descrizione:
#  - Carica un dataset RDS di trascrittomica single-cell simulata,
#    con coordinate spaziali e "intensity_cluster" come etichette vere.
#  - Esegue un'analisi con Seurat (normalizzazione, PCA, clustering, UMAP)
#  - Esegue clustering HDBSCAN su UMAP
#  - Confronta i cluster predetti rispetto a quelli reali (intensity_cluster)
#    usando l'algoritmo Hungarian (solve_LSAP del package 'clue') per ottenere
#    una corrispondenza ottimale e calcolare TP/FP/FN/Precision in modo coerente.

# -----------------------------------------------------------------------------
# Caricamento librerie
# -----------------------------------------------------------------------------
suppressMessages({
  library(tidyverse)    # Per dplyr, ggplot2, etc.
  library(Seurat)       # Per analisi scRNA-seq
  library(future)       # Parallelizzazione
  library(future.apply) # Parallelizzazione
  library(patchwork)    # Per combinare i plot
  library(dbscan)       # Per HDBSCAN
  library(cluster)      # silhouette, etc.
  library(factoextra)   # visualizzazioni di clustering
  library(mclust)       # ARI, AMI
  library(Matrix)       # Gestione matrici sparse
  library(clue)         # Hungarian algorithm (solve_LSAP)
})

# -----------------------------------------------------------------------------
# Funzione di calcolo TP, FP, FN, Precision
# -----------------------------------------------------------------------------
calculate_metrics <- function(true_labels, predicted_labels, method_name) {
  # Convertiamo a fattori per garantire che i livelli siano espliciti
  true_factors      <- factor(true_labels)
  predicted_factors <- factor(predicted_labels)

  all_classes <- union(levels(true_factors), levels(predicted_factors))
  true_factors      <- factor(true_factors,      levels = all_classes)
  predicted_factors <- factor(predicted_factors, levels = all_classes)

  # Creiamo la tabella di contingenza
  contingency <- table(true_factors, predicted_factors)

  # TP per classe (diagonale)
  TP <- diag(contingency)
  # FP per classe (somma colonna - diagonale)
  FP <- colSums(contingency) - TP
  # FN per classe (somma riga - diagonale)
  FN <- rowSums(contingency) - TP
  # Precision per classe
  Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)

  # Tabella per singola classe
  metrics_per_class <- tibble(
    Class           = rownames(contingency),
    True_Positives  = TP,
    False_Positives = FP,
    False_Negatives = FN,
    Precision       = Precision
  )

  # Metriche totali (sommando su tutte le classi)
  total_TP <- sum(TP)
  total_FP <- sum(FP)
  total_FN <- sum(FN)
  total_Precision <- ifelse((total_TP + total_FP) > 0,
                            total_TP / (total_TP + total_FP),
                            NA)

  metrics_total <- tibble(
    Class           = "Total",
    True_Positives  = total_TP,
    False_Positives = total_FP,
    False_Negatives = total_FN,
    Precision       = total_Precision
  )

  # Combiniamo e aggiungiamo la colonna 'Method'
  metrics <- bind_rows(metrics_per_class, metrics_total) %>%
    mutate(Method = method_name) %>%
    dplyr::select(Method, everything())

  return(metrics)
}

# -----------------------------------------------------------------------------
# Funzione per la mappatura Hungarian (righe = cluster predetti, colonne = cluster veri)
# -----------------------------------------------------------------------------
assign_clusters_hungarian <- function(true_labels, predicted_labels) {
  # 1) Costruiamo la matrice di contingenza INVERTITA:
  #    righe = cluster predetti, colonne = cluster veri
  cont_mat <- table(predicted_labels, true_labels)

  # 2) Matrice di costi (massimo - overlap)
  max_val  <- max(cont_mat)
  cost_mat <- max_val - cont_mat

  # 3) Risoluzione dell'assegnamento ottimo
  #    Ora le righe (predicted) <= colonne (true) se abbiamo meno cluster predetti che veri
  assignment <- solve_LSAP(cost_mat)

  # row_names = cluster predetti, col_names = cluster veri
  row_names <- rownames(cont_mat)  # cluster predetti
  col_names <- colnames(cont_mat)  # cluster veri

  # 4) Costruiamo un vettore di mapping: "cluster_predetto" -> "cluster_vero"
  #    assignment[i] = j significa: la riga i (predetto i) è assegnata alla colonna j (vero j)
  cluster_map <- character(length(row_names))
  names(cluster_map) <- row_names

  for (i in seq_along(row_names)) {
    col_assigned <- assignment[i]
    cluster_map[row_names[i]] <- col_names[col_assigned]
  }

  # 5) Convertiamo i cluster predetti in character
  predicted_labels_char <- as.character(predicted_labels)
  predicted_mapped      <- character(length(predicted_labels_char))

  # 6) Assegniamo l'etichetta vera corrispondente
  for (pred_cl in row_names) {
    new_label <- cluster_map[pred_cl]
    predicted_mapped[predicted_labels_char == pred_cl] <- new_label
  }

  # Volendo, si può rimappare in factor con i livelli = col_names
  # predicted_mapped <- factor(predicted_mapped, levels = col_names)

  return(predicted_mapped)
}

# -----------------------------------------------------------------------------
# Funzione principale
# -----------------------------------------------------------------------------
analyze_and_compare_clusters <- function(
    rds_path          = NULL,
    output_path       = NULL,
    k_cell_types      = NULL,
    seurat_resolution = NULL,
    hdbscan_minPts    = NULL
) {
  # 1) Parallelizzazione
  plan(multisession, workers = min(3, availableCores()))

  # 2) Verifica file RDS
  if (!file.exists(rds_path)) {
    stop(paste("Errore: il file", rds_path, "non esiste."))
  }
  cat("Analizzando il file:", rds_path, "\n")

  # 3) Caricamento dataset
  data <- readRDS(rds_path)

  # 4) Check 'intensity_cluster'
  if (!"intensity_cluster" %in% names(data)) {
    stop("Errore: 'intensity_cluster' non trovato nei dati simulati.")
  }

  # k_cell_types
  if (is.null(k_cell_types)) {
    if ("parameters" %in% names(data) && "k_cell_types" %in% names(data$parameters)) {
      k_cell_types <- data$parameters$k_cell_types
      cat("Numero di tipi cellulari letto dal file RDS:", k_cell_types, "\n")
    } else {
      k_cell_types <- length(unique(data$intensity_cluster))
      cat("Numero di tipi cellulari dedotto dai dati:", k_cell_types, "\n")
    }
  } else {
    cat("Numero di tipi cellulari specificato dall'utente:", k_cell_types, "\n")
  }

  # 5) IDs delle cellule
  cell_ids <- paste("Cell", seq_len(nrow(data$expression)), sep="_")

  # 6) Matrice di espressione
  cat("Creazione della matrice di espressione...\n")
  expression_matrix <- Matrix::Matrix(t(data$expression), sparse = TRUE)

  # 7) Crea oggetto Seurat
  cat("Creazione dell'oggetto Seurat...\n")
  seu <- CreateSeuratObject(
    counts   = expression_matrix,
    meta.data = data.frame(
      cell              = cell_ids,
      x                 = data$coordinates$x,
      y                 = data$coordinates$y,
      intensity_cluster = data$intensity_cluster
    )
  )
  rownames(seu@meta.data) <- cell_ids

  # 8) Plot diagnostici pre-filtering
  cat("Generazione plot diagnostici pre-filtering...\n")
  plot1 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
  plot2 <- VlnPlot(seu, features = "nCount_RNA",   pt.size = 0.1) + NoLegend()
  print(plot1 + plot2)

  # 9) Filtro geni a bassa espressione
  cat("Filtraggio geni a bassa espressione...\n")
  counts_matrix <- GetAssayData(seu, assay = "RNA", slot = "counts")
  genes_to_keep <- rownames(seu)[Matrix::rowSums(counts_matrix) > 0]
  seu <- subset(seu, features = genes_to_keep)

  cat(sprintf("Numero di geni dopo il filtraggio: %d\n", nrow(seu)))
  cat(sprintf("Numero di cellule dopo il filtraggio: %d\n", ncol(seu)))

  # 10) Plot diagnostici post-filtering
  cat("Generazione plot diagnostici post-filtering...\n")
  plot3 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
  plot4 <- VlnPlot(seu, features = "nCount_RNA",   pt.size = 0.1) + NoLegend()
  print(plot3 + plot4)

  # 11) Normalizzazione con SCTransform
  cat("Normalizzazione dei dati con SCTransform...\n")
  seu <- SCTransform(seu, verbose = FALSE)
  DefaultAssay(seu) <- "SCT"

  plot5 <- VlnPlot(seu, features = "nFeature_SCT", pt.size = 0.1) + NoLegend()
  plot6 <- VlnPlot(seu, features = "nCount_SCT",   pt.size = 0.1) + NoLegend()
  print(plot5 + plot6)

  # 12) PCA
  cat("Esecuzione PCA...\n")
  seu <- RunPCA(seu, assay = "SCT", npcs = 30, approximate = FALSE)
  cat("PCA completata.\n")
  print(seu[["pca"]], dims = 1:5, nfeatures = 5)

  elbow_plot <- ElbowPlot(seu, ndims = 30) + scale_x_continuous(breaks = seq(0, 30, 5))
  print(elbow_plot)

  # 13) Clustering con Seurat
  cat("Esecuzione clustering con Seurat...\n")
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = seurat_resolution, verbose = FALSE)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("Errore: 'seurat_clusters' non è stato creato da FindClusters.")
  } else {
    num_seurat_clusters <- length(unique(seu@meta.data$seurat_clusters))
    cat("Clustering Seurat completato. # cluster =", num_seurat_clusters, "\n")
    if (num_seurat_clusters < 2) {
      warning("Avviso: meno di 2 cluster identificati. Aumentare la risoluzione.")
    }
  }

  # 14) UMAP
  cat("Esecuzione UMAP...\n")
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30,
                 n.neighbors = 5, min.dist = 0.25)
  DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") +
    ggtitle("Seurat Clusters")

  # 15) Silhouette Score (su UMAP)
  cat("Calcolo Silhouette Score...\n")
  distance_matrix <- dist(Embeddings(seu[["umap"]])[, 1:2])
  clusters       <- seu@meta.data$seurat_clusters
  silhouette_vals <- silhouette(as.numeric(clusters), distance_matrix)
  seu@meta.data$silhouette_score <- silhouette_vals[, 3]
  fviz_silhouette(silhouette_vals, label = FALSE, print.summary = TRUE)

  # 16) Plot spaziale (x,y)
  cat("Plot spaziale cluster Seurat...\n")
  ggplot(seu@meta.data, aes(x = x, y = y, color = as.factor(seurat_clusters))) +
    geom_point(alpha = 0.8) +
    labs(title = "Spatial Plot - Seurat Clusters",
         x = "Spatial X", y = "Spatial Y") +
    theme_minimal()

  # 17) Plot proporzioni cluster
  cat("Plot proporzioni Seurat...\n")
  total_cells <- ncol(seu)
  cell_counts <- table(seu@meta.data$seurat_clusters) %>%
    as.data.frame() %>%
    rename(Cluster = Var1, Count = Freq) %>%
    mutate(Percentage = (Count / total_cells) * 100)

  ggplot(cell_counts, aes(x = Cluster, y = Percentage, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack",
             colour = "black", size = 0.1, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)),
              position = position_stack(vjust = 0.5), size = 3) +
    labs(x = "Cluster", y = "Percentage of Cells") +
    theme_minimal() +
    scale_y_continuous(breaks = seq(0, 100, by = 10))

  # 18) HDBSCAN su UMAP
  cat("HDBSCAN su UMAP...\n")
  umap_embeddings <- Embeddings(seu, "umap")
  cluster_results <- hdbscan(umap_embeddings, minPts = hdbscan_minPts)
  cat("Risultati HDBSCAN:\n")
  print(table(cluster_results$cluster))

  seu <- AddMetaData(seu, metadata = cluster_results$cluster, col.name = "HDBSCAN")

  # 19) Assegnamento ottimo Hungarian: Seurat
  cat("Assegnamento ottimo (Hungarian) per Seurat...\n")
  true_labels        <- seu@meta.data$intensity_cluster
  pred_seurat        <- seu@meta.data$seurat_clusters
  pred_seurat_mapped <- assign_clusters_hungarian(true_labels, pred_seurat)

  seu@meta.data$SeuratMapped <- pred_seurat_mapped

  # 20) Assegnamento ottimo (Hungarian) per HDBSCAN
  cat("Assegnamento ottimo (Hungarian) per HDBSCAN...\n")
  pred_hdbscan           <- seu@meta.data$HDBSCAN
  pred_hdbscan_mapped    <- assign_clusters_hungarian(true_labels, pred_hdbscan)
  seu@meta.data$HDBSCANMapped <- pred_hdbscan_mapped

  # 21) Tabelle di contingenza post-mappatura
  cat("Tabella di contingenza - Seurat (mapped) vs True:\n")
  cont_seurat_mapped <- table(True = true_labels, Predicted = pred_seurat_mapped)
  print(cont_seurat_mapped)

  cat("Tabella di contingenza - HDBSCAN (mapped) vs True:\n")
  cont_hdbscan_mapped <- table(True = true_labels, Predicted = pred_hdbscan_mapped)
  print(cont_hdbscan_mapped)

  # 22) Calcolo delle metriche di performance
  cat("\nCalcolo metriche di performance...\n")
  metrics_seurat    <- calculate_metrics(true_labels, pred_seurat_mapped,  "Seurat_Hungarian")
  metrics_hdbscan   <- calculate_metrics(true_labels, pred_hdbscan_mapped, "HDBSCAN_Hungarian")
  metrics_combined  <- bind_rows(metrics_seurat, metrics_hdbscan)

  cat("\nMetriche di Performance Combinate:\n")
  print(metrics_combined)

  # 23) Salviamo le metriche su CSV
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  write_csv(metrics_combined, output_path)
  cat("Metriche salvate in:", output_path, "\n")

  # 24) Salvataggio oggetto Seurat
  output_rds <- sub("\\.csv$", ".rds", output_path)
  output_rds <- sub("metrics_comparison", "seurat_hdbscan_combined", output_rds)
  dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(seu, file = output_rds)
  cat("Oggetto Seurat salvato in:", output_rds, "\n")

  # Ripristina la parallelizzazione
  plan(sequential)

  return(metrics_combined)
}

# ============================
# Esempio di utilizzo
# ============================
 analyze_and_compare_clusters(
   rds_path          = "results/hard/simulated_granuloma_hard.rds",
   output_path       = "results/metrics_comparison_hard_hungarian.csv",
   seurat_resolution = 0.4,
   hdbscan_minPts    = 7
 )
