#!/usr/bin/env Rscript

# Script: analyze_and_compare_clusters.R
# Descrizione: Analizza i dati simulati con Seurat e HDBSCAN, rinomina i cluster come "cells_a", "cells_b", "cells_c",
#              e calcola True Positives (TP), False Positives (FP), False Negatives (FN) e Precision confrontando con i cluster reali.

# ========================================
# Caricamento librerie
# ========================================
suppressMessages({
  library(tidyverse)    # Per funzioni d'aiuto (dplyr, ggplot2, ecc.)
  library(Seurat)        # Per analisi di trascrittomica a singola cellula
  library(future)        # Per parallelizzazione
  library(future.apply)  # Per parallelizzazione
  library(patchwork)     # Per combinare plot
  library(dbscan)        # Per HDBSCAN
  library(cluster)       # Per silhouette/clustering o altri metodi
  library(factoextra)    # Per visualizzare risultati di clustering
  library(mclust)        # Per calcolare ARI e AMI
  library(Matrix)        # Per gestire matrici sparse
})

# ========================================
# 1) Impostazione della Parallelizzazione
# ========================================
plan(multisession, workers = min(3, availableCores()))  # Adatta il numero di workers secondo il tuo sistema

# ========================================
# 2) Path dei dati simulati
# ========================================
data_path <- "data/simulated_image_correlation.rds"

# ========================================
# 3) Verifica che il file esista
# ========================================
if (!file.exists(data_path)) {
  stop(paste("Errore: il file", data_path, "non esiste."))
}

cat("Analizzando il file:", data_path, "\n")

# ========================================
# 4) Caricamento del dataset
# ========================================
data <- readRDS(data_path)

# ========================================
# 5) Verifica che 'intensity_cluster' sia presente
# ========================================
if (!"intensity_cluster" %in% names(data)) {
  stop("Errore: 'intensity_cluster' non trovato nei dati simulati.")
}

# ========================================
# 6) Genera gli identificatori delle cellule
# ========================================
cell_ids <- paste("Cell", seq_len(nrow(data$expression)), sep="_")

# ========================================
# 7) Crea la matrice di espressione
# ========================================
cat("Creazione della matrice di espressione...\n")
expression_matrix <- Matrix::Matrix(t(data$expression), sparse = TRUE)  # Seurato si aspetta righe = geni, colonne = cellule

# ========================================
# 8) Crea l'oggetto Seurat includendo 'intensity_cluster' nel meta.data
# ========================================
cat("Creazione dell'oggetto Seurat...\n")
seu <- CreateSeuratObject(
  counts = expression_matrix,
  meta.data = data.frame(
    cell = cell_ids,
    x = data$coordinates$x,
    y = data$coordinates$y,
    intensity_cluster = data$intensity_cluster  # Aggiungi questa riga
  )
)

# Imposta i nomi delle righe
rownames(seu@meta.data) <- cell_ids

# ========================================
# 9) Plot diagnostici per metriche pre-filtering
# ========================================
cat("Generazione dei plot diagnostici pre-filtering...\n")
plot1 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
plot2 <- VlnPlot(seu, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
print(plot1 + plot2)

# ========================================
# 10) Filtro geni che hanno espressione pari a 0 in tutte le cellule
# ========================================
cat("Filtraggio dei geni a bassa espressione...\n")
counts_matrix <- GetAssayData(seu, assay = "RNA", slot = "counts")
genes_to_keep <- rownames(seu)[Matrix::rowSums(counts_matrix) > 0]
seu <- subset(seu, features = genes_to_keep)

# Verifica numero di geni e celle rimasti
cat(sprintf("Numero di geni dopo il filtraggio: %d\n", nrow(seu)))
cat(sprintf("Numero di cellule dopo il filtraggio: %d\n", ncol(seu)))

# ========================================
# 11) Plot diagnostici per metriche post-filtering
# ========================================
cat("Generazione dei plot diagnostici post-filtering...\n")
plot3 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
plot4 <- VlnPlot(seu, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
print(plot3 + plot4)

# ========================================
# 12) Normalizzazione con SCTransform
# ========================================
cat("Normalizzazione dei dati con SCTransform...\n")
seu <- SCTransform(seu, verbose=FALSE)

# Imposta l'assay di default
DefaultAssay(seu) <- "SCT"

# Ispezione delle metriche dopo la normalizzazione
plot5 <- VlnPlot(seu, features = "nFeature_SCT", pt.size = 0.1) + NoLegend()
plot6 <- VlnPlot(seu, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
print(plot5 + plot6)

# ========================================
# 13) PCA
# ========================================
cat("Esecuzione della PCA...\n")
seu <- RunPCA(
  object = seu,
  assay = "SCT",
  npcs = 30,
  approximate = FALSE
)

# Verifica PCA
cat("PCA completata.\n")
print(seu[["pca"]], dims = 1:5, nfeatures = 5)

# Selezione delle PC da utilizzare
elbow_plot <- ElbowPlot(seu, ndims = 30) +
  scale_x_continuous(breaks = seq(0, 30, 5))
print(elbow_plot)

# ========================================
# 14) Rinomina dei cluster reali basati sulla dimensione
# ========================================
cat("Rinomina dei cluster reali basati sulla dimensione...\n")
true_cluster_sizes <- table(seu@meta.data$intensity_cluster)
true_clusters_sorted <- sort(true_cluster_sizes, decreasing = TRUE)
true_sorted_names <- names(true_clusters_sorted)

# Assegna etichette "cells_a", "cells_b", "cells_c" basate sulla dimensione
num_true_clusters <- min(length(true_sorted_names), 5)
true_labels <- c("cells_a", "cells_b", "cells_c", "cells_d", "cells_e")[1:num_true_clusters]

# Crea un vettore di mappatura per i cluster reali
true_cluster_renamed <- setNames(true_labels, true_sorted_names[1:num_true_clusters])

# Assegna "noise" per eventuali cluster rimanenti
if (length(true_sorted_names) > num_true_clusters) {
  true_cluster_renamed <- c(true_cluster_renamed, setNames(rep("noise", length(true_sorted_names) - num_true_clusters),
                                                           true_sorted_names[(num_true_clusters + 1):length(true_sorted_names)]))
}

# Crea una nuova colonna per i cluster reali rinominati
seu@meta.data$RenamedTrue <- true_cluster_renamed[as.character(seu@meta.data$intensity_cluster)]

# ========================================
# 15) Clustering con Seurat
# ========================================
cat("Esecuzione del clustering con Seurat...\n")
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seurat_resolution <- 0.25  # Aumentata da 0.25 a 0.5
seu <- FindClusters(seu, resolution = seurat_resolution, verbose = FALSE)

# Verifica presenza dei cluster
if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("Errore: 'seurat_clusters' non è stato creato durante FindClusters.")
} else {
  num_seurat_clusters <- length(unique(seu@meta.data$seurat_clusters))
  cat("Clustering con Seurat completato. Numero di cluster identificati:", num_seurat_clusters, "\n")
  if (num_seurat_clusters < 2) {
    warning("Avviso: meno di 2 cluster identificati da Seurat. Considera di aumentare la risoluzione.")
  }
}

# Rinomina dei cluster Seurat basati sulla dimensione
cat("Rinomina dei cluster Seurat basati sulla dimensione...\n")
seurat_cluster_sizes <- table(seu@meta.data$seurat_clusters)
seurat_clusters_sorted <- sort(seurat_cluster_sizes, decreasing = TRUE)
seurat_sorted_names <- names(seurat_clusters_sorted)

# Assegna etichette "cells_a", "cells_b", "cells_c" basate sulla dimensione
num_seurat_labels <- min(length(seurat_sorted_names), 5)
seurat_labels <- c("cells_a", "cells_b", "cells_c", "cells_d", "cells_e")[1:num_seurat_labels]

# Crea un vettore di mappatura per i cluster Seurat
seurat_cluster_renamed <- setNames(seurat_labels, seurat_sorted_names[1:num_seurat_labels])

# Assegna "noise" per eventuali cluster rimanenti
if (length(seurat_sorted_names) > num_seurat_labels) {
  seurat_cluster_renamed <- c(seurat_cluster_renamed, setNames(rep("noise", length(seurat_sorted_names) - num_seurat_labels),
                                                               seurat_sorted_names[(num_seurat_labels + 1):length(seurat_sorted_names)]))
}

# Crea una nuova colonna per i cluster Seurat rinominati
seu@meta.data$RenamedSeurat <- seurat_cluster_renamed[as.character(seu@meta.data$seurat_clusters)]

# Verifica che tutte le etichette siano state assegnate correttamente
if (any(is.na(seu@meta.data$RenamedSeurat))) {
  warning("Alcuni cluster Seurat non sono stati rinominati correttamente.")
}

# ========================================
# 16) Run UMAP con Seurat
# ========================================
cat("Esecuzione di UMAP con Seurat...\n")
seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, n.neighbors = 5, min.dist = 0.25)
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") + ggtitle("Seurat Clusters")

# ========================================
# 17) Validazione qualità del clustering con Silhouette Score
# ========================================
cat("Calcolo del Silhouette Score...\n")
distance_matrix <- dist(Embeddings(seu[['umap']])[, 1:2])
clusters <- seu@meta.data$seurat_clusters
silhouette_vals <- silhouette(as.numeric(clusters), distance_matrix)
seu@meta.data$silhouette_score <- silhouette_vals[, 3]

# Plot Silhouette
fviz_silhouette(silhouette_vals, label = FALSE, print.summary = TRUE)

# ========================================
# 18) Plot clusters nello spazio xy
# ========================================
cat("Creazione del plot spaziale dei cluster Seurat...\n")
ggplot(seu@meta.data, aes(x = x, y = y, color = as.factor(seurat_clusters))) +
  geom_point(alpha = 0.8) +
  labs(title = "Spatial Plot of Simulated Data - Seurat Clusters", x = "Spatial X", y = "Spatial Y", color = "Cluster") +
  theme_minimal()

# ========================================
# 19) Plot proporzioni clusters
# ========================================
cat("Creazione del plot delle proporzioni dei cluster Seurat...\n")
total_cells <- ncol(seu)
cell_counts <- table(seu@meta.data$seurat_clusters) %>%
  as.data.frame() %>%
  rename(Cluster = Var1, Count = Freq) %>%
  mutate(Percentage = (Count / total_cells) * 100)

ggplot(cell_counts, aes(x = Cluster, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", size = 0.1, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), position = position_stack(vjust = 0.5), size = 3) +
  labs(x = "Cluster", y = "Percentage of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))

# ========================================
# 20) Applicazione di HDBSCAN sui dati UMAP
# ========================================
cat("Applicazione di HDBSCAN sui dati UMAP...\n")
umap_embeddings <- Embeddings(seu, "umap")

# Parametri di HDBSCAN (iniziali)
hdbscan_minPts <- 7  # Questo parametro può essere modificato nella fase di tuning

# Esegui HDBSCAN
cluster_results <- hdbscan(umap_embeddings, minPts = hdbscan_minPts)

# Verifica risultato HDBSCAN
cat("Risultati HDBSCAN:\n")
print(table(cluster_results$cluster))

# Aggiungi i cluster HDBSCAN all'oggetto Seurat
seu <- AddMetaData(seu, metadata = cluster_results$cluster, col.name = "HDBSCAN")

# Rinomina dei cluster HDBSCAN basati sulla dimensione
cat("Rinomina dei cluster HDBSCAN basati sulla dimensione...\n")
hdbscan_cluster_sizes <- table(cluster_results$cluster)
hdbscan_clusters_sorted <- sort(hdbscan_cluster_sizes, decreasing = TRUE)
hdbscan_sorted_names <- names(hdbscan_clusters_sorted)

# Assegna etichette "cells_a", "cells_b", "cells_c" basate sulla dimensione
num_hdbscan_labels <- min(length(hdbscan_sorted_names), 5)
hdbscan_labels <- c("cells_a", "cells_b", "cells_c", "cells_d", "cells_e")[1:num_hdbscan_labels]

# Crea un vettore di mappatura per i cluster HDBSCAN
hdbscan_cluster_renamed <- setNames(hdbscan_labels, hdbscan_sorted_names[1:num_hdbscan_labels])

# Assegna "noise" per eventuali cluster rimanenti
if (length(hdbscan_sorted_names) > num_hdbscan_labels) {
  hdbscan_cluster_renamed <- c(hdbscan_cluster_renamed, setNames(rep("noise", length(hdbscan_sorted_names) - num_hdbscan_labels),
                                                                 hdbscan_sorted_names[(num_hdbscan_labels + 1):length(hdbscan_sorted_names)]))
}

# Crea una nuova colonna per i cluster HDBSCAN rinominati
seu@meta.data$RenamedHDBSCAN <- hdbscan_cluster_renamed[as.character(seu@meta.data$HDBSCAN)]

# Verifica che tutte le etichette siano state assegnate correttamente
if (any(is.na(seu@meta.data$RenamedHDBSCAN))) {
  warning("Alcuni cluster HDBSCAN non sono stati rinominati correttamente.")
}

# ========================================
# 21) Creazione delle tabelle di contingenza
# ========================================
cat("Creazione delle tabelle di contingenza...\n")
# Creazione della tabella di contingenza con i cluster HDBSCAN rinominati
contingency_table_renamed_hdbscan <- table(seu@meta.data$intensity_cluster, seu@meta.data$RenamedHDBSCAN)
cat("Tabella di contingenza (intensity_cluster vs Renamed HDBSCAN):\n")
print(contingency_table_renamed_hdbscan)

# Creazione della tabella di contingenza con i cluster Seurat rinominati
contingency_table_renamed_seurat <- table(seu@meta.data$intensity_cluster, seu@meta.data$RenamedSeurat)
cat("Tabella di contingenza (intensity_cluster vs Renamed Seurat):\n")
print(contingency_table_renamed_seurat)

# ========================================
# 22) Salva l'oggetto Seurat con entrambi i cluster rinominati
# ========================================
cat("Salvataggio dell'oggetto Seurat con i cluster rinominati...\n")
saveRDS(seu, file = "results/seurat_hdbscan_combined_simulated_image_correlation.rds")
cat("Analisi completata e oggetto Seurat salvato in 'results/seurat_hdbscan_combined_simulated_image_correlation.rds'.\n")

# ========================================
# 23) Calcolo di TP, FP, FN e Precision
# ========================================
cat("Calcolo delle metriche di performance per i cluster Seurat e HDBSCAN...\n")

# Estrai le etichette vere e predette
true_labels <- as.character(seu@meta.data$RenamedTrue)
predicted_seurat <- as.character(seu@meta.data$RenamedSeurat)
predicted_hdbscan <- as.character(seu@meta.data$RenamedHDBSCAN)

# Funzione per calcolare TP, FP, FN e Precision per ogni metodo
calculate_metrics <- function(true_labels, predicted_labels, method_name) {
  # Converti le etichette in fattori per gestire tutte le classi
  true_factors <- factor(true_labels)
  predicted_factors <- factor(predicted_labels)

  # Trova tutte le classi presenti in entrambe le etichette
  all_classes <- union(levels(true_factors), levels(predicted_factors))

  # Ricrea le etichette come fattori con tutti i livelli
  true_factors <- factor(true_labels, levels = all_classes)
  predicted_factors <- factor(predicted_labels, levels = all_classes)

  # Crea la tabella di contingenza
  contingency <- table(true_factors, predicted_factors)

  # Stampa la tabella di contingenza per debug
  cat(sprintf("\nTabella di contingenza per %s:\n", method_name))
  print(contingency)

  # Calcola True Positives (TP) per ogni classe
  TP <- diag(contingency)

  # Calcola False Positives (FP) per ogni classe
  FP <- colSums(contingency) - TP

  # Calcola False Negatives (FN) per ogni classe
  FN <- rowSums(contingency) - TP

  # Calcola Precisione per ogni classe
  Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)

  # Crea un tibble con i risultati per ogni classe
  metrics_per_class <- tibble(
    Class = rownames(contingency),
    True_Positives = TP,
    False_Positives = FP,
    False_Negatives = FN,
    Precision = Precision
  )

  # Calcola le metriche totali
  total_TP <- sum(TP)
  total_FP <- sum(FP)
  total_FN <- sum(FN)
  total_Precision <- ifelse((total_TP + total_FP) > 0, total_TP / (total_TP + total_FP), NA)

  metrics_total <- tibble(
    Class = "Total",
    True_Positives = total_TP,
    False_Positives = total_FP,
    False_Negatives = total_FN,
    Precision = total_Precision
  )

  # Combina le metriche per classe e totali
  metrics <- bind_rows(metrics_per_class, metrics_total) %>%
    mutate(Method = method_name) %>%
    dplyr::select(Method, everything())  # Specifica dplyr::select

  return(metrics)
}

# Calcola le metriche per Seurat
cat("\nCalcolo di True Positives (TP), False Positives (FP), False Negatives (FN) e Precision per Seurat:\n")
metrics_seurat <- calculate_metrics(true_labels, predicted_seurat, "Seurat")

# Calcola le metriche per HDBSCAN
cat("\nCalcolo di True Positives (TP), False Positives (FP) e False Negatives (FN) e Precision per HDBSCAN:\n")
metrics_hdbscan <- calculate_metrics(true_labels, predicted_hdbscan, "HDBSCAN")

# Combina i risultati
metrics_combined <- bind_rows(metrics_seurat, metrics_hdbscan)

# Stampa i risultati
cat("\nMetriche di Performance Combinate:\n")
print(metrics_combined)

# Salva i risultati in un file CSV
write_csv(metrics_combined, "results/metrics_comparison.csv")
cat("Confronto completato e metriche salvate in 'results/metrics_comparison_colon_40k100gene5clust.csv'.\n")

# ========================================
# 24) Ripristina la parallelizzazione
# ========================================
plan(sequential)
