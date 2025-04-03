#!/usr/bin/env Rscript

# Script: analyze_and_compare_clusters_label_free.R
# Descrizione:
#  - Carica un dataset RDS di trascrittomica single-cell simulata,
#    con coordinate spaziali e "intensity_cluster" come etichette vere.
#  - Esegue un'analisi con Seurat (normalizzazione, PCA, clustering, UMAP)
#  - Esegue clustering HDBSCAN su UMAP
#  - Confronta i cluster predetti con quelli reali usando l'ARI (da mclust)
#    e l'AMI (da aricode), metriche "label-free" che non richiedono
#    di rinominare i cluster.

# ---------------------------
# Caricamento librerie
# ---------------------------
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(future)
  library(future.apply)
  library(patchwork)
  library(dbscan)     # per HDBSCAN
  library(cluster)    # silhouette se ti interessa
  library(factoextra) # per eventuali visualizzazioni di clustering
  library(mclust)     # Adjusted Rand Index (adjustedRandIndex)
  library(Matrix)

  # Per AMI, usiamo il pacchetto aricode
  # Se non ce l'hai installato, devi prima fare:
  # install.packages("aricode")
  library(aricode)    # Per AMI(true, pred)
})

analyze_and_compare_clusters_label_free <- function(
    rds_path          = NULL,
    output_path       = NULL,
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

  # 4) Verifica 'intensity_cluster'
  if (!"intensity_cluster" %in% names(data)) {
    stop("Errore: 'intensity_cluster' non trovato nei dati simulati.")
  }
  true_labels <- data$intensity_cluster

  # 5) Crea IDs
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
      intensity_cluster = true_labels
    )
  )
  rownames(seu@meta.data) <- cell_ids

  # (Plot diagnostici pre-filtering opzionali)
  cat("Generazione plot diagnostici pre-filtering...\n")
  plot1 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
  plot2 <- VlnPlot(seu, features = "nCount_RNA",   pt.size = 0.1) + NoLegend()
  print(plot1 + plot2)

  # 8) Filtro geni a bassa espressione
  cat("Filtraggio geni a bassa espressione...\n")
  counts_matrix <- GetAssayData(seu, assay = "RNA", slot = "counts")
  genes_to_keep <- rownames(seu)[Matrix::rowSums(counts_matrix) > 0]
  seu <- subset(seu, features = genes_to_keep)

  cat(sprintf("Numero di geni dopo il filtraggio: %d\n", nrow(seu)))
  cat(sprintf("Numero di cellule dopo il filtraggio: %d\n", ncol(seu)))

  # 9) Plot diagnostici post-filtering
  cat("Generazione plot diagnostici post-filtering...\n")
  plot3 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
  plot4 <- VlnPlot(seu, features = "nCount_RNA",   pt.size = 0.1) + NoLegend()
  print(plot3 + plot4)

  # 10) Normalizzazione con SCTransform
  cat("Normalizzazione dei dati con SCTransform...\n")
  seu <- SCTransform(seu, verbose = FALSE)
  DefaultAssay(seu) <- "SCT"

  plot5 <- VlnPlot(seu, features = "nFeature_SCT", pt.size = 0.1) + NoLegend()
  plot6 <- VlnPlot(seu, features = "nCount_SCT",   pt.size = 0.1) + NoLegend()
  print(plot5 + plot6)

  # 11) PCA
  cat("Esecuzione PCA...\n")
  seu <- RunPCA(seu, assay = "SCT", npcs = 30, approximate = FALSE)
  cat("PCA completata.\n")
  print(seu[["pca"]], dims = 1:5, nfeatures = 5)

  elbow_plot <- ElbowPlot(seu, ndims = 30) + scale_x_continuous(breaks = seq(0, 30, 5))
  print(elbow_plot)

  # 12) Clustering Seurat
  cat("Clustering con Seurat (resolution =", seurat_resolution, ")...\n")
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = seurat_resolution, verbose = FALSE)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("Errore: 'seurat_clusters' non è stato creato.")
  } else {
    num_seurat_clusters <- length(unique(seu@meta.data$seurat_clusters))
    cat("Clustering Seurat completato. # cluster =", num_seurat_clusters, "\n")
  }

  # 13) UMAP
  cat("Esecuzione UMAP...\n")
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, n.neighbors = 5, min.dist = 0.25)
  DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") +
    ggtitle("Seurat Clusters")

  # 14) Silhouette (opzionale)
  cat("Calcolo Silhouette Score...\n")
  distance_matrix <- dist(Embeddings(seu[["umap"]])[, 1:2])
  clusters       <- seu@meta.data$seurat_clusters
  silhouette_vals <- silhouette(as.numeric(clusters), distance_matrix)
  seu@meta.data$silhouette_score <- silhouette_vals[, 3]
  fviz_silhouette(silhouette_vals, label = FALSE, print.summary = TRUE)

  # 15) HDBSCAN su UMAP
  cat("HDBSCAN su UMAP...\n")
  umap_embeddings <- Embeddings(seu, "umap")
  cluster_results <- hdbscan(umap_embeddings, minPts = hdbscan_minPts)
  cat("Risultati HDBSCAN:\n")
  print(table(cluster_results$cluster))

  seu <- AddMetaData(seu, metadata = cluster_results$cluster, col.name = "HDBSCAN")

  # 16) Calcolo ARI e AMI (label-free)
  cat("Calcolo ARI e AMI...\n")
  seurat_clusters  <- seu@meta.data$seurat_clusters
  hdbscan_clusters <- seu@meta.data$HDBSCAN

  # ARI da mclust
  ari_seurat   <- adjustedRandIndex(true_labels, seurat_clusters)
  ari_hdbscan  <- adjustedRandIndex(true_labels, hdbscan_clusters)

  # AMI da aricode
  ami_seurat   <- AMI(true_labels, seurat_clusters)
  ami_hdbscan  <- AMI(true_labels, hdbscan_clusters)

  metrics_df <- tibble(
    Method = c("Seurat", "HDBSCAN"),
    ARI    = c(ari_seurat, ari_hdbscan),
    AMI    = c(ami_seurat, ami_hdbscan)
  )
  cat("Risultati ARI/AMI:\n")
  print(metrics_df)

  # 17) Salviamo su CSV
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  write_csv(metrics_df, output_path)
  cat("Metriche ARI/AMI salvate in:", output_path, "\n")

  # 18) Salvataggio Seurat
  output_rds <- sub("\\.csv$", ".rds", output_path)
  output_rds <- sub("metrics_ari_ami", "seurat_hdbscan_ari_ami", output_rds)
  dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(seu, file = output_rds)
  cat("Oggetto Seurat salvato in:", output_rds, "\n")

  # Fine e ritorno in sequenziale
  plan(sequential)
  return(metrics_df)
}

# ============================
# Esempio di utilizzo
# ============================
 analyze_and_compare_clusters_label_free(
   rds_path          = "results/easy/simulated_granuloma_easy.rds",
   output_path       = "results/metrics_ari_ami_easy_0_1.csv",
   seurat_resolution = 0.25,
   hdbscan_minPts    = 7
 )
