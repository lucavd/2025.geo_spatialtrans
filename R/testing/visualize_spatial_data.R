#!/usr/bin/env Rscript
# Script per visualizzazioni esplorative dei dati spatial transcriptomics
# IMPORTANTE: Eseguire DOPO full_test_visHD.R per visualizzare dati già validati biologicamente
# Usage: Rscript visualize_spatial_data.R

library(ggplot2)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(Matrix)  # Per gestire matrici sparse

# Configurazione
output_dir <- "results/visualizations"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Carica i dati salvati da full_test_visHD.R
cat("Caricamento dati salvati...\n")

# Verifica che i dati siano stati validati biologicamente
if (!file.exists("R/testing/full_test_result.rds")) {
  stop("ERRORE: full_test_result.rds non trovato. Eseguire prima full_test_visHD.R")
}
if (!file.exists("R/testing/full_simulation_data.rds")) {
  stop("ERRORE: full_simulation_data.rds non trovato. Eseguire prima full_test_visHD.R")
}

final_data <- readRDS("R/testing/full_simulation_data.rds")
full_result <- readRDS("R/testing/full_test_result.rds")

# Verifica che i dati siano stati validati con successo
if (full_result$status != "PASS") {
  stop("ERRORE: I dati non hanno superato la validazione biologica. Status: ", full_result$status)
}

cat("✓ Dati biologicamente validati caricati con successo\n")
cat("✓ Validazione timestamp:", as.character(full_result$timestamp), "\n")

# Estrai i componenti necessari
simulated_data <- final_data$expression  # Matrice sparsa dgCMatrix
cell_df <- data.frame(
  X = final_data$coordinates$x,
  Y = final_data$coordinates$y,
  cluster = final_data$clusters
)
params <- list(
  cfg = final_data$config,
  diff_cfg = full_result$difficulty_config
)

# Verifica che simulated_data sia una matrice
if (!is.matrix(simulated_data) && !inherits(simulated_data, "Matrix")) {
  stop("simulated_data deve essere una matrice o Matrix sparsa")
}

# Estrai informazioni base
n_cells <- ncol(simulated_data)
n_genes <- nrow(simulated_data)
spatial_coords <- as.data.frame(cell_df[, c("X", "Y")])
cell_types <- as.factor(cell_df$cluster)

cat(sprintf("Dataset: %d celle, %d geni, %d tipi cellulari\n", 
            n_cells, n_genes, length(unique(cell_types))))

# 1. DISTRIBUZIONE SPAZIALE DEI TIPI CELLULARI
cat("\n1. Generazione mappa tipi cellulari...\n")
p1 <- ggplot(cbind(spatial_coords, cluster = cell_types), aes(x = X, y = Y, color = cluster)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_color_brewer(palette = "Set3", name = "Cell Type") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "white"),
        legend.position = "right") +
  labs(title = "Distribuzione Spaziale dei Tipi Cellulari",
       subtitle = sprintf("%d celle, %d tipi", n_cells, length(unique(cell_types))))

ggsave(file.path(output_dir, "01_cell_types_spatial.png"), p1, 
       width = 10, height = 8, dpi = 300, bg = "white")

# 2. DENSITÀ CELLULARE PER REGIONE
cat("2. Generazione heatmap densità cellulare...\n")
# Crea griglia esagonale per densità
library(hexbin)
hb <- hexbin(spatial_coords$X, spatial_coords$Y, xbins = 50)
hex_df <- data.frame(hcell2xy(hb), count = hb@count)

p2 <- ggplot(hex_df, aes(x = x, y = y, fill = count)) +
  geom_hex(stat = "identity") +
  scale_fill_viridis(name = "Cell\nDensity", option = "plasma") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey95")) +
  labs(title = "Densità Cellulare Spaziale",
       subtitle = "Aggregazione esagonale")

ggsave(file.path(output_dir, "02_cell_density.png"), p2, 
       width = 8, height = 8, dpi = 300)

# 3. DISTRIBUZIONE UMI TOTALI
cat("3. Generazione distribuzione UMI...\n")
umi_per_cell <- colSums(simulated_data)
spatial_umi <- cbind(spatial_coords, umi = umi_per_cell)

p3a <- ggplot(spatial_umi, aes(x = X, y = Y, color = log10(umi + 1))) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_viridis(name = "log10(UMI+1)", option = "magma") +
  coord_fixed() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "black")) +
  labs(title = "Distribuzione Spaziale UMI Totali")

# Istogramma UMI per tipo cellulare
umi_df <- data.frame(umi = umi_per_cell, cell_type = cell_types)
p3b <- ggplot(umi_df, aes(x = umi, fill = cell_type)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_x_log10() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Distribuzione UMI per Tipo Cellulare",
       x = "UMI totali (log scale)", y = "Numero di celle")

p3 <- grid.arrange(p3a, p3b, ncol = 2)
ggsave(file.path(output_dir, "03_umi_distribution.png"), p3, 
       width = 16, height = 8, dpi = 300)

# 4. GENI MARKER BIOLOGICAMENTE DEFINITI
cat("4. Visualizzazione geni marker secondo configurazione biologica...\n")

# Estrai configurazione marker dalla validazione biologica
marker_config <- params$diff_cfg$marker_params
n_marker_per_type <- marker_config$marker_genes_per_type  # 25
fold_expected <- marker_config$marker_expression_fold      # 2.5
overlap_expected <- marker_config$marker_overlap_fold     # 0.05

cat(sprintf("Configurazione biologica: %d marker/tipo, fold=%.1fx, overlap=%.2f\n", 
            n_marker_per_type, fold_expected, overlap_expected))

# USA L'ASSEGNAZIONE SEQUENZIALE BIOLOGICA (come in 06b_expression_baseline.R)
biological_markers <- list()
n_types <- length(unique(cell_types))
for (k in 1:min(n_types, 8)) {  # Max 8 tipi
  start_idx <- (k - 1) * n_marker_per_type + 1
  end_idx <- min(k * n_marker_per_type, n_genes)
  if (start_idx <= end_idx) {
    biological_markers[[k]] <- start_idx:end_idx
    cat(sprintf("Tipo %d: Geni %d-%d (%d marker)\n", k, start_idx, end_idx, length(start_idx:end_idx)))
  }
}

# Calcola medie per validazione
mean_expr_by_type <- matrix(0, nrow = n_genes, ncol = n_types)
for (i in 1:n_types) {
  cells_type_i <- which(cell_types == i)
  if (length(cells_type_i) > 0) {
    mean_expr_by_type[, i] <- rowMeans(simulated_data[, cells_type_i, drop = FALSE])
  }
}

# VALIDAZIONE: Verifica che i marker biologici siano effettivamente più espressi
marker_validation <- data.frame()
for (k in 1:length(biological_markers)) {
  marker_genes <- biological_markers[[k]]
  if (length(marker_genes) > 0) {
    # Media marker vs media altri geni in questo tipo
    marker_mean <- mean(mean_expr_by_type[marker_genes, k])
    non_marker_mean <- mean(mean_expr_by_type[-marker_genes, k])
    fold_observed <- marker_mean / (non_marker_mean + 1e-6)
    
    marker_validation <- rbind(marker_validation, data.frame(
      tipo = k,
      fold_observed = fold_observed,
      fold_expected = fold_expected,
      marker_mean = marker_mean,
      non_marker_mean = non_marker_mean
    ))
  }
}

print(marker_validation)
top_markers <- biological_markers  # Usa marker biologici invece di quelli calcolati

# Visualizza espressione spaziale dei marker
marker_plots <- list()
for (ct in seq_len(min(4, length(unique(cell_types))))) {  # Primi 4 tipi cellulari
  gene_idx <- top_markers[[ct]][1]  # Top marker
  expr_values <- simulated_data[gene_idx, ]
  
  p <- ggplot(cbind(spatial_coords, expr = expr_values), 
              aes(x = X, y = Y, color = expr)) +
    geom_point(size = 0.2, alpha = 0.7) +
    scale_color_viridis(name = "Expression", option = "viridis") +
    coord_fixed() +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "black")) +
    labs(title = sprintf("Marker Cell Type %d (Gene %d)", ct, gene_idx))
  
  marker_plots[[ct]] <- p
}

p4 <- do.call(grid.arrange, c(marker_plots, ncol = 2))
ggsave(file.path(output_dir, "04_marker_genes_spatial.png"), p4, 
       width = 12, height = 12, dpi = 300)

# 5. HEATMAP ESPRESSIONE GENI MARKER
cat("5. Generazione heatmap geni marker...\n")
# Seleziona tutti i marker identificati
all_markers <- unique(unlist(top_markers))
marker_expr <- simulated_data[all_markers, ]

# Ordina celle per tipo
cell_order <- order(cell_types)
marker_expr_ordered <- marker_expr[, cell_order]

# Subsample per visualizzazione (max 2000 celle per dataset grandi)
if (ncol(marker_expr_ordered) > 2000) {
  # Campiona in modo bilanciato per tipo cellulare
  sample_idx <- c()
  cells_per_type <- 250  # Max celle per tipo
  unique_types <- unique(cell_types[cell_order])
  
  for (ct in unique_types) {
    ct_indices <- which(cell_types[cell_order] == ct)
    n_sample <- min(length(ct_indices), cells_per_type)
    if (n_sample > 0) {
      sample_idx <- c(sample_idx, sample(ct_indices, n_sample))
    }
  }
  
  # Verifica validità degli indici
  sample_idx <- sample_idx[sample_idx <= ncol(marker_expr_ordered)]
  sample_idx <- sample_idx[sample_idx > 0]
  
  if (length(sample_idx) > 0) {
    marker_expr_plot <- marker_expr_ordered[, sample_idx, drop = FALSE]
    cell_types_plot <- cell_types[cell_order][sample_idx]
    umi_plot <- umi_per_cell[cell_order][sample_idx]
  } else {
    # Fallback: usa prime 2000 celle
    marker_expr_plot <- marker_expr_ordered[, 1:min(2000, ncol(marker_expr_ordered)), drop = FALSE]
    cell_types_plot <- cell_types[cell_order][1:min(2000, length(cell_types[cell_order]))]
    umi_plot <- umi_per_cell[cell_order][1:min(2000, length(umi_per_cell[cell_order]))]
  }
} else {
  marker_expr_plot <- marker_expr_ordered
  cell_types_plot <- cell_types[cell_order]
  umi_plot <- umi_per_cell[cell_order]
}

# Crea annotazioni per il subset
# Genera nomi colonne se non presenti
if (is.null(colnames(marker_expr_plot))) {
  colnames(marker_expr_plot) <- paste0("Cell_", seq_len(ncol(marker_expr_plot)))
}

col_annotation <- data.frame(
  CellType = as.factor(cell_types_plot),
  row.names = colnames(marker_expr_plot)
)

# Heatmap con pheatmap
library(pheatmap)
library(RColorBrewer)

# Normalizza per riga (z-score)
marker_expr_scaled <- t(scale(t(as.matrix(marker_expr_plot))))
marker_expr_scaled[is.na(marker_expr_scaled)] <- 0

# Colori per annotazioni - assicurati che corrispondano ai livelli effettivi
actual_cell_types <- sort(unique(as.character(col_annotation$CellType)))
n_colors_needed <- length(actual_cell_types)
color_palette <- if (n_colors_needed <= 8) {
  brewer.pal(max(3, n_colors_needed), "Set3")[1:n_colors_needed]
} else {
  colorRampPalette(brewer.pal(8, "Set3"))(n_colors_needed)
}

ann_colors <- list(
  CellType = setNames(color_palette, actual_cell_types)
)

tryCatch({
  # Crea heatmap con gestione corretta del device
  heatmap_file <- file.path(output_dir, "05_marker_heatmap.png")
  
  # Usa pheatmap con filename per evitare problemi di device
  pheatmap(marker_expr_scaled,
           show_colnames = FALSE,
           show_rownames = TRUE,
           annotation_col = col_annotation,
           annotation_colors = ann_colors,
           clustering_cols = FALSE,
           clustering_rows = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Espressione Geni Marker Biologici per Tipo Cellulare",
           filename = heatmap_file,
           width = 12, height = 8)
  
  cat("✓ Heatmap salvata:", heatmap_file, "\n")
}, error = function(e) {
  cat("Errore in heatmap:", e$message, "\n")
  # Fallback: prova con ggsave
  tryCatch({
    library(pheatmap)
    p_heat <- pheatmap(marker_expr_scaled,
                       show_colnames = FALSE,
                       show_rownames = TRUE,
                       annotation_col = col_annotation,
                       annotation_colors = ann_colors,
                       clustering_cols = FALSE,
                       clustering_rows = TRUE,
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       main = "Espressione Geni Marker Biologici per Tipo Cellulare",
                       silent = TRUE)
    ggsave(file.path(output_dir, "05_marker_heatmap.png"), p_heat, 
           width = 12, height = 8, dpi = 300)
  }, error = function(e2) {
    cat("Fallback heatmap fallito:", e2$message, "\n")
  })
})

# 6. STATISTICHE DI QUALITÀ
cat("7. Generazione statistiche di qualità...\n")
# Calcola metriche
genes_per_cell <- colSums(simulated_data > 0)
zero_inflation <- sum(simulated_data == 0) / (n_cells * n_genes)
detection_rate <- rowMeans(simulated_data > 0)

# Panel statistiche
p7a <- ggplot(data.frame(genes = genes_per_cell, cell_type = cell_types), 
              aes(x = cell_type, y = genes, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Geni Rilevati per Tipo Cellulare", y = "Numero di geni")

p7b <- ggplot(data.frame(detection = detection_rate), aes(x = detection)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribuzione Detection Rate dei Geni",
       x = "Frazione di celle con espressione > 0", y = "Numero di geni")

# Scatter UMI vs geni
p7c <- ggplot(data.frame(umi = umi_per_cell, genes = genes_per_cell, 
                         cell_type = cell_types), 
              aes(x = umi, y = genes, color = cell_type)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_log10() +
  scale_color_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "UMI vs Geni Rilevati", x = "UMI totali (log)", y = "Geni rilevati")

p7 <- grid.arrange(p7a, p7b, p7c, ncol = 2)
ggsave(file.path(output_dir, "07_quality_metrics.png"), p7, 
       width = 12, height = 10, dpi = 300)

# 8. VALIDAZIONE LIBRARY SIZE E BIOLOGICA
cat("8. Validazione parametri biologici target...\n")

# Estrai parametri target dalla configurazione
library_target <- params$diff_cfg$cell_specific_params$library_size_params$mean_library_size
dropout_range <- params$diff_cfg$dropout_params$dropout_range
sparsity_target <- mean(dropout_range) * 100

cat(sprintf("Target biologici: Library=%.0f UMI, Dropout=%.0f-%.0f%%, Sparsity~%.0f%%\n", 
            library_target, dropout_range[1]*100, dropout_range[2]*100, sparsity_target))

# Calcola metriche osservate
observed_library <- mean(umi_per_cell)
observed_sparsity <- sum(simulated_data == 0) / (n_cells * n_genes) * 100
observed_cv <- sd(umi_per_cell) / mean(umi_per_cell)

# 9. VALIDAZIONE FOLD-CHANGE MARKER SPECIFICI
cat("9. Validazione fold-change marker per ogni tipo cellulare...\n")

# Crea plot per validation fold-change
fold_validation_data <- data.frame()
for (k in 1:length(biological_markers)) {
  if (length(biological_markers[[k]]) > 0) {
    marker_genes <- biological_markers[[k]][1:min(5, length(biological_markers[[k]]))]  # Max 5 marker per tipo
    
    for (g in marker_genes) {
      # Calcola espressione media di questo gene in questo tipo vs altri tipi
      expr_in_type <- mean(simulated_data[g, cell_types == k])
      expr_in_others <- mean(simulated_data[g, cell_types != k])
      fold_obs <- expr_in_type / (expr_in_others + 1e-6)
      
      fold_validation_data <- rbind(fold_validation_data, data.frame(
        gene = paste0("Gene_", g),
        cell_type = factor(paste0("Tipo_", k)),
        fold_observed = fold_obs,
        fold_expected = fold_expected,
        expr_in_type = expr_in_type,
        expr_in_others = expr_in_others,
        marker_for_type = k
      ))
    }
  }
}

p8a <- ggplot(fold_validation_data, aes(x = cell_type, y = fold_observed)) +
  geom_boxplot(alpha = 0.7, fill = "lightblue") +
  geom_hline(yintercept = fold_expected, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 1, y = fold_expected + 0.2, label = paste("Target:", fold_expected, "x"), 
           color = "red", hjust = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Validazione Fold-Change Geni Marker",
       subtitle = sprintf("Ogni boxplot mostra fold-change dei marker per quel tipo (%d marker/tipo)", n_marker_per_type),
       x = "Tipo Cellulare", y = "Fold-Change Osservato")

# 10. CONFRONTO ESPRESSIONE MARKER: SPECIFICO vs ALTRI TIPI  
cat("10. Confronto espressione marker: specifico vs altri tipi...\n")

# Crea un plot che mostra chiaramente la specificità dei marker
specificity_data <- data.frame()
for (k in 1:min(4, length(biological_markers))) {  # Primi 4 tipi
  if (length(biological_markers[[k]]) > 0) {
    marker_genes <- biological_markers[[k]][1:2]  # Primi 2 marker per tipo
    
    for (g in marker_genes) {
      # Aggiungi dati per il tipo specifico
      cells_this_type <- which(cell_types == k)
      if (length(cells_this_type) > 0) {
        specificity_data <- rbind(specificity_data, data.frame(
          expr = simulated_data[g, cells_this_type],
          gene = sprintf("Marker_Tipo%d_Gene%d", k, g),
          condition = paste0("Tipo_", k, " (specifico)"),
          is_specific = TRUE,
          marker_for_type = k
        ))
      }
      
      # Aggiungi dati per tutti gli altri tipi
      cells_other_types <- which(cell_types != k)
      if (length(cells_other_types) > 0) {
        specificity_data <- rbind(specificity_data, data.frame(
          expr = simulated_data[g, sample(cells_other_types, min(1000, length(cells_other_types)))],
          gene = sprintf("Marker_Tipo%d_Gene%d", k, g),
          condition = "Altri tipi",
          is_specific = FALSE,
          marker_for_type = k
        ))
      }
    }
  }
}

p8b <- ggplot(specificity_data, aes(x = condition, y = log2(expr + 1), fill = is_specific)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.8) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "red"), 
                    name = "Specificità", labels = c("Altri tipi", "Tipo specifico")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Specificità Geni Marker: Espressione in Tipo Specifico vs Altri",
       subtitle = "I marker dovrebbero essere altamente espressi nel loro tipo specifico (rosso)",
       x = "Condizione", y = "log2(Expression + 1)")

# Combina i plot
p8 <- grid.arrange(p8a, p8b, nrow = 2, heights = c(1, 2))
ggsave(file.path(output_dir, "08_marker_validation_specificity.png"), p8, 
       width = 16, height = 12, dpi = 300, bg = "white")

# RIEPILOGO CON VALIDAZIONE BIOLOGICA
cat("\n=== VISUALIZZAZIONI GENERATE CON VALIDAZIONE BIOLOGICA ===\n")
cat("Output directory:", output_dir, "\n")
cat("1. 01_cell_types_spatial.png - Mappa spaziale dei tipi cellulari\n")
cat("2. 02_cell_density.png - Heatmap densità cellulare\n")
cat("3. 03_umi_distribution.png - Distribuzione UMI spaziale e per tipo\n")
cat("4. 04_marker_genes_spatial.png - Espressione spaziale marker biologici\n")
cat("5. 05_marker_heatmap.png - Heatmap marker con assegnazione biologica\n")
cat("6. 07_quality_metrics.png - Metriche di qualità dataset\n")
cat("7. 08_marker_validation_specificity.png - Validazione fold-change e specificità marker\n")

# REPORT VALIDAZIONE FINALE
cat("\n=== REPORT VALIDAZIONE BIOLOGICA ===\n")
if (exists("marker_validation") && nrow(marker_validation) > 0) {
  avg_fold_observed <- mean(marker_validation$fold_observed, na.rm = TRUE)
  cat(sprintf("✓ Fold-change marker medio: %.2fx (target: %.1fx)\n", avg_fold_observed, fold_expected))
  
  fold_ok <- abs(avg_fold_observed - fold_expected) < 1.0
  cat(sprintf("✓ Validazione fold-change: %s\n", ifelse(fold_ok, "PASS", "REVIEW")))
}

validation_ok <- abs(observed_library - library_target) < (library_target * 0.3)
cat(sprintf("✓ Library size: %.0f UMI (target: %.0f) - %s\n", 
            observed_library, library_target, ifelse(validation_ok, "PASS", "REVIEW")))

sparsity_ok <- abs(observed_sparsity - sparsity_target) < 20
cat(sprintf("✓ Sparsity: %.1f%% (target: ~%.0f%%) - %s\n", 
            observed_sparsity, sparsity_target, ifelse(sparsity_ok, "PASS", "REVIEW")))

# Salva anche statistiche riassuntive
stats_summary <- list(
  n_cells = n_cells,
  n_genes = n_genes,
  n_cell_types = length(unique(cell_types)),
  zero_inflation = zero_inflation,
  mean_umi_per_cell = mean(umi_per_cell),
  mean_genes_per_cell = mean(genes_per_cell),
  validation_info = list(
    validated_by = "full_test_visHD.R",
    validation_timestamp = full_result$timestamp,
    validation_status = full_result$status
  )
)

saveRDS(stats_summary, file.path(output_dir, "stats_summary.rds"))
cat("\nStatistiche salvate in stats_summary.rds\n")
cat("\nNOTA: Validazione biologica già completata da full_test_visHD.R\n")
