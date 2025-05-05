#' Plot Spaziale dell'Espressione dei Geni Marker
#'
#' Visualizza la distribuzione spaziale dell'espressione di un gene specifico.
#' 
#' @param sim_results Lista di risultati della simulazione
#' @param gene_id ID o nome del gene marker da visualizzare
#' @param log_transform Booleano, se TRUE applica log(x+1) ai valori di espressione
#' @param title Titolo opzionale del grafico
#' @param color_scale Scala di colori da utilizzare (default: viridis)
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_minimal labs
#' @export
plot_marker_gene_spatial <- function(sim_results, gene_id, log_transform = TRUE, 
                                    title = NULL, color_scale = "viridis") {
  # Estrazione delle coordinate spaziali
  coords <- sim_results$coordinates
  
  # Verifica che la matrice di espressione sia nel formato corretto (geni x celle)
  if (!is.matrix(sim_results$expression)) {
    stop("L'oggetto expression deve essere una matrice")
  }
  
  if (ncol(sim_results$expression) != nrow(coords)) {
    warning("La matrice di espressione non ha il formato corretto (geni x celle). Tentativo di trasposizione...")
    # Transpone la matrice se ha la dimensione sbagliata
    if (nrow(sim_results$expression) == nrow(coords)) {
      sim_results$expression <- t(sim_results$expression)
    } else {
      stop(paste("Dimensioni incompatibili: matrice expression", 
                 dim(sim_results$expression)[1], "x", dim(sim_results$expression)[2], 
                 ", numero di celle", nrow(coords)))
    }
  }
  
  # Estrazione dell'espressione del gene
  if (is.character(gene_id)) {
    gene_idx <- which(rownames(sim_results$expression) == gene_id)
    if (length(gene_idx) == 0) {
      stop("Gene ID not found in the expression matrix")
    }
  } else if (is.numeric(gene_id)) {
    gene_idx <- gene_id
    if (gene_idx > nrow(sim_results$expression)) {
      stop("Gene index out of bounds")
    }
  } else {
    stop("gene_id must be either a character string or a numeric index")
  }
  
  # Estrazione dell'espressione
  gene_expr <- sim_results$expression[gene_idx, ]
  
  # Log-trasformazione se richiesta
  if (log_transform) {
    gene_expr <- log1p(gene_expr)
  }
  
  # Verifica che intensity_cluster sia un factor
  if (!is.factor(sim_results$intensity_cluster)) {
    sim_results$intensity_cluster <- factor(sim_results$intensity_cluster)
  }
  
  # Creazione di un dataframe per il plot
  plot_df <- data.frame(
    x = coords$x,
    y = coords$y,
    expression = gene_expr,
    cluster = sim_results$intensity_cluster
  )
  
  # Creazione del titolo se non fornito
  if (is.null(title)) {
    if (is.character(gene_id)) {
      title <- paste("Espressione spaziale del gene", gene_id)
    } else {
      title <- paste("Espressione spaziale del gene", rownames(sim_results$expression)[gene_idx])
    }
  }
  
  # Creazione del plot
  p <- ggplot(plot_df, aes(x = x, y = y, color = expression)) +
    geom_point(size = 1) +
    scale_color_viridis_c(option = color_scale) +
    theme_minimal() +
    labs(title = title, 
         x = "Coordinata X", 
         y = "Coordinata Y", 
         color = "Espressione") +
    theme(aspect.ratio = 1)
  
  return(p)
}

#' Plot della Relazione Media-Varianza dei Geni
#'
#' Visualizza la relazione tra media e varianza dell'espressione genica, 
#' evidenziando la sovradispersione tipica dei dati trascrittomici.
#' 
#' @param sim_results Lista di risultati della simulazione
#' @param log_scale Booleano, se TRUE usa scala logaritmica per entrambi gli assi
#' @param highlight_markers Booleano, se TRUE evidenzia i geni marker
#' @param title Titolo opzionale del grafico
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_point geom_abline scale_x_log10 scale_y_log10 
#'                     theme_minimal labs
#' @importFrom stats var
#' @export
plot_mean_variance_relationship <- function(sim_results, log_scale = TRUE, 
                                          highlight_markers = TRUE, title = NULL) {
  # Calcolo della media e varianza per ogni gene
  gene_means <- rowMeans(sim_results$expression)
  gene_vars <- apply(sim_results$expression, 1, var)
  
  # Identificazione dei geni marker, se richiesto
  markers <- NULL
  if (highlight_markers && !is.null(sim_results$parameters$marker_params$marker_genes_per_type)) {
    # Estrazione dei geni marker dai parametri
    markers <- sim_results$parameters$marker_params$marker_genes_per_type
    markers <- unlist(markers)
  }
  
  # Creazione di un dataframe per il plot
  plot_df <- data.frame(
    gene = rownames(sim_results$expression),
    mean = gene_means,
    variance = gene_vars,
    is_marker = rownames(sim_results$expression) %in% markers
  )
  
  # Creazione del titolo se non fornito
  if (is.null(title)) {
    title <- "Relazione Media-Varianza dell'Espressione Genica"
  }
  
  # Creazione del plot base
  p <- ggplot(plot_df, aes(x = mean, y = variance)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey") +
    theme_minimal() +
    labs(title = title, 
         x = "Media dell'Espressione", 
         y = "Varianza dell'Espressione")
  
  # Aggiunta dei punti, eventualmente colorati per marker
  if (highlight_markers && !is.null(markers)) {
    p <- p + geom_point(aes(color = is_marker), alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"),
                         labels = c("FALSE" = "Non marker", "TRUE" = "Marker"))
  } else {
    p <- p + geom_point(alpha = 0.7)
  }
  
  # Applicazione della scala logaritmica se richiesta
  if (log_scale) {
    p <- p + scale_x_log10() + scale_y_log10()
  }
  
  return(p)
}

#' Plot della Relazione Dropout-Espressione Media
#'
#' Visualizza la relazione tra il tasso di dropout (frazione di zeri) e
#' l'espressione media dei geni, un aspetto critico delle tecnologie trascrittomiche.
#'
#' @param sim_results Lista di risultati della simulazione
#' @param log_mean Booleano, se TRUE applica log(x+1) alle medie
#' @param highlight_markers Booleano, se TRUE evidenzia i geni marker
#' @param title Titolo opzionale del grafico
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_x_log10
#'                     theme_minimal labs
#' @export
plot_dropout_vs_expression <- function(sim_results, log_mean = TRUE, 
                                     highlight_markers = TRUE, title = NULL) {
  # Calcolo della frazione di zeri per gene
  zero_fraction <- rowMeans(sim_results$expression == 0)
  
  # Calcolo dell'espressione media per gene
  gene_means <- rowMeans(sim_results$expression)
  if (log_mean) {
    gene_means <- log1p(gene_means)
  }
  
  # Identificazione dei geni marker, se richiesto
  markers <- NULL
  if (highlight_markers && !is.null(sim_results$parameters$marker_params$marker_genes_per_type)) {
    # Estrazione dei geni marker dai parametri
    markers <- sim_results$parameters$marker_params$marker_genes_per_type
    markers <- unlist(markers)
  }
  
  # Creazione di un dataframe per il plot
  plot_df <- data.frame(
    gene = rownames(sim_results$expression),
    mean_expression = gene_means,
    zero_fraction = zero_fraction,
    is_marker = rownames(sim_results$expression) %in% markers
  )
  
  # Creazione del titolo se non fornito
  if (is.null(title)) {
    title <- "Relazione Dropout-Espressione Media"
  }
  
  # Creazione del plot base
  p <- ggplot(plot_df, aes(x = mean_expression, y = zero_fraction)) +
    geom_smooth(method = "loess", se = TRUE, color = "blue", alpha = 0.2) +
    theme_minimal() +
    labs(title = title, 
         x = ifelse(log_mean, "Log(Media dell'Espressione + 1)", "Media dell'Espressione"), 
         y = "Frazione di Zeri (Dropout Rate)")
  
  # Aggiunta dei punti, eventualmente colorati per marker
  if (highlight_markers && !is.null(markers)) {
    p <- p + geom_point(aes(color = is_marker), alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"),
                         labels = c("FALSE" = "Non marker", "TRUE" = "Marker"))
  } else {
    p <- p + geom_point(alpha = 0.7)
  }
  
  return(p)
}

#' Plot UMAP/t-SNE Colorato per Tipo Cellulare Simulato
#'
#' Crea una visualizzazione ridotta dimensionalmente (UMAP o t-SNE) dei profili 
#' di espressione, colorata per i tipi cellulari simulati.
#'
#' @param sim_results Lista di risultati della simulazione
#' @param method Metodo di riduzione dimensionale ("umap" o "tsne")
#' @param n_genes_hvg Numero di geni altamente variabili da utilizzare
#' @param perplexity Parametro di perplexity per t-SNE
#' @param title Titolo opzionale del grafico
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs
#' @export
plot_dimensionality_reduction <- function(sim_results, method = "umap", 
                                        n_genes_hvg = 2000, perplexity = 30, 
                                        title = NULL) {
  # Verifica della disponibilità dei pacchetti necessari
  if (method == "umap") {
    if (!requireNamespace("umap", quietly = TRUE)) {
      stop("Il pacchetto 'umap' è necessario per questo metodo. Installalo con install.packages('umap')")
    }
  } else if (method == "tsne") {
    if (!requireNamespace("Rtsne", quietly = TRUE)) {
      stop("Il pacchetto 'Rtsne' è necessario per questo metodo. Installalo con install.packages('Rtsne')")
    }
  } else {
    stop("Metodo non riconosciuto. Scegli tra 'umap' o 'tsne'")
  }
  
  # Selezione dei geni più variabili
  genes_var <- apply(sim_results$expression, 1, var)
  hvg_idx <- order(genes_var, decreasing = TRUE)[1:min(n_genes_hvg, length(genes_var))]
  
  # Preparazione della matrice di espressione (trasposizione per avere celle x geni)
  expr_matrix <- t(sim_results$expression[hvg_idx, ])
  
  # Applicazione della riduzione dimensionale
  reduced_dims <- NULL
  if (method == "umap") {
    umap_result <- umap::umap(expr_matrix)
    reduced_dims <- umap_result$layout
    colnames(reduced_dims) <- c("UMAP1", "UMAP2")
  } else if (method == "tsne") {
    tsne_result <- Rtsne::Rtsne(expr_matrix, perplexity = perplexity, check_duplicates = FALSE)
    reduced_dims <- tsne_result$Y
    colnames(reduced_dims) <- c("tSNE1", "tSNE2")
  }
  
  # Creazione di un dataframe per il plot
  plot_df <- data.frame(
    reduced_dims,
    cluster = sim_results$intensity_cluster
  )
  
  # Creazione del titolo se non fornito
  if (is.null(title)) {
    if (method == "umap") {
      title <- "Visualizzazione UMAP dei Profili di Espressione"
    } else {
      title <- "Visualizzazione t-SNE dei Profili di Espressione"
    }
  }
  
  # Creazione del plot
  dim_names <- colnames(reduced_dims)
  p <- ggplot(plot_df, aes_string(x = dim_names[1], y = dim_names[2], color = "cluster")) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(title = title, 
         x = dim_names[1], 
         y = dim_names[2], 
         color = "Tipo Cellulare")
  
  return(p)
}

#' Identifica Automaticamente i Geni Marker per ogni Tipo Cellulare
#'
#' @param sim_results Lista di risultati della simulazione
#' @param n_markers Numero di marker per tipo cellulare da identificare
#' @return Lista con i geni marker per ogni tipo cellulare
#' @importFrom stats aggregate
#' @export
identify_marker_genes <- function(sim_results, n_markers = 5) {
  # Estrazione dei tipi cellulari unici
  cell_types <- levels(sim_results$intensity_cluster)
  
  # Creazione della matrice di espressione cellula x gene
  expr_t <- t(sim_results$expression)
  expr_df <- data.frame(expr_t, cluster = sim_results$intensity_cluster)
  
  # Calcolo dell'espressione media per ogni cluster e gene
  mean_expr_by_cluster <- list()
  for (cluster in cell_types) {
    # Subset delle celle appartenenti al cluster
    cluster_cells <- expr_df$cluster == cluster
    # Calcolo della media per ogni gene
    cluster_means <- colMeans(expr_df[cluster_cells, -ncol(expr_df)])
    mean_expr_by_cluster[[cluster]] <- cluster_means
  }
  
  # Identificazione dei geni marker per ogni cluster
  marker_genes <- list()
  for (i in seq_along(cell_types)) {
    cluster <- cell_types[i]
    # Calcolo del rapporto di espressione rispetto agli altri cluster
    ratios <- list()
    for (j in seq_along(cell_types)) {
      if (i != j) {
        other_cluster <- cell_types[j]
        # Rapporto di espressione tra il cluster corrente e un altro cluster
        ratio <- mean_expr_by_cluster[[cluster]] / (mean_expr_by_cluster[[other_cluster]] + 1e-10)
        ratios[[other_cluster]] <- ratio
      }
    }
    
    # Media dei rapporti rispetto a tutti gli altri cluster
    avg_ratio <- Reduce(`+`, ratios) / length(ratios)
    
    # Selezione dei top n_markers geni con il rapporto più alto
    marker_idx <- order(avg_ratio, decreasing = TRUE)[1:min(n_markers, length(avg_ratio))]
    marker_genes[[cluster]] <- names(avg_ratio)[marker_idx]
  }
  
  return(marker_genes)
}

#' Identifica i Geni Altamente Variabili
#'
#' @param expression_matrix Matrice di espressione (geni x celle)
#' @param n_hvg Numero di geni da selezionare
#' @return Vettore con i nomi/indici dei geni selezionati
#' @importFrom stats var
#' @export
find_highly_variable_genes <- function(expression_matrix, n_hvg = 2000) {
  # Calcolo della varianza per ogni gene
  gene_var <- apply(expression_matrix, 1, var)
  
  # Ordinamento dei geni per varianza
  hvg_idx <- order(gene_var, decreasing = TRUE)[1:min(n_hvg, length(gene_var))]
  
  return(rownames(expression_matrix)[hvg_idx])
}

#' Plot della distribuzione spaziale dell'espressione per tutti i cluster
#'
#' Crea un pannello di plot che mostrano l'espressione spaziale dei geni marker
#' per ogni cluster, con un gene marker rappresentativo per ogni cluster.
#'
#' @param sim_results Lista di risultati della simulazione
#' @param marker_list Lista di geni marker per ogni cluster
#' @param title Titolo del grafico
#' @param log_transform Applicare trasformazione log all'espressione
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_viridis_c theme_minimal labs
#' @export
plot_spatial_expression_panel <- function(sim_results, marker_list = NULL, 
                                          title = "Distribuzione spaziale dell'espressione genica",
                                          subtitle = "Un gene marker rappresentativo per ogni cluster",
                                          log_transform = TRUE) {
  if (is.null(marker_list)) {
    # Identificazione automatica dei geni marker se non specificati
    marker_list <- identify_marker_genes(sim_results, n_markers = 1)
  }
  
  # Estrazione delle coordinate spaziali
  coords <- sim_results$coordinates
  
  # Verifica che la matrice di espressione sia nel formato corretto (geni x celle)
  if (!is.matrix(sim_results$expression)) {
    stop("L'oggetto expression deve essere una matrice")
  }
  
  if (ncol(sim_results$expression) != nrow(coords)) {
    warning("La matrice di espressione non ha il formato corretto (geni x celle). Tentativo di trasposizione...")
    # Transpone la matrice se ha la dimensione sbagliata
    if (nrow(sim_results$expression) == nrow(coords)) {
      sim_results$expression <- t(sim_results$expression)
    } else {
      stop(paste("Dimensioni incompatibili: matrice expression", 
                 dim(sim_results$expression)[1], "x", dim(sim_results$expression)[2], 
                 ", numero di celle", nrow(coords)))
    }
  }
  
  # Verifica che intensity_cluster sia un factor
  if (!is.factor(sim_results$intensity_cluster)) {
    sim_results$intensity_cluster <- factor(sim_results$intensity_cluster)
  }
  
  # Crea un dataframe per il plot combinato
  plot_df <- data.frame()
  
  # Per ogni cluster, aggiungi il gene marker rappresentativo
  for (cluster in names(marker_list)) {
    if (length(marker_list[[cluster]]) > 0) {
      gene_id <- marker_list[[cluster]][1]  # Prendi il primo gene marker
      
      # Trova l'indice del gene
      if (is.character(gene_id)) {
        gene_idx <- which(rownames(sim_results$expression) == gene_id)
        if (length(gene_idx) == 0) {
          warning(paste("Gene marker", gene_id, "per cluster", cluster, "non trovato nella matrice"))
          next
        }
      } else {
        gene_idx <- gene_id
        if (gene_idx > nrow(sim_results$expression)) {
          warning(paste("Indice gene", gene_idx, "fuori dai limiti della matrice"))
          next
        }
      }
      
      # Estrai l'espressione
      gene_expr <- sim_results$expression[gene_idx, ]
      
      # Log-trasformazione se richiesta
      if (log_transform) {
        gene_expr <- log1p(gene_expr)
      }
      
      # Nome del plot
      marker_name <- ifelse(is.character(gene_id), gene_id, 
                           paste0("gene_", gene_idx))
      
      # Aggiungi al dataframe
      cluster_df <- data.frame(
        x = coords$x,
        y = coords$y,
        expression = gene_expr,
        cluster = sim_results$intensity_cluster,
        marker = paste0("marker_cluster_", cluster)
      )
      
      plot_df <- rbind(plot_df, cluster_df)
    }
  }
  
  # Crea il plot
  if (nrow(plot_df) > 0) {
    p <- ggplot(plot_df, aes(x = x, y = y, color = expression)) +
      geom_point(size = 1) +
      facet_wrap(~ marker) +
      scale_color_viridis_c(option = "viridis") +
      theme_minimal() +
      labs(title = title, 
           subtitle = subtitle,
           x = "Coordinata X", 
           y = "Coordinata Y", 
           color = "Espressione")
    
    return(p)
  } else {
    warning("Nessun dato per il plot")
    return(NULL)
  }
}

#' Plot delle distanze intra-cluster
#'
#' Visualizza le distanze medie all'interno di ogni cluster
#' colorando i punti in base alla distanza media dagli altri punti del cluster.
#'
#' @param sim_results Lista di risultati della simulazione
#' @param n_points Numero massimo di punti da visualizzare per cluster (per efficienza)
#' @param title Titolo del grafico
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_viridis_c theme_minimal labs
#' @export
plot_distances_by_cluster <- function(sim_results, n_points = 5000,
                                     title = "Distanza media intra-cluster",
                                     subtitle = "Distribuzione spaziale per cluster") {
  # Verifica requisiti
  if (!requireNamespace("fields", quietly = TRUE)) {
    warning("Il pacchetto 'fields' è necessario per il calcolo delle distanze")
    return(NULL)
  }
  
  if (!is.factor(sim_results$intensity_cluster)) {
    sim_results$intensity_cluster <- factor(sim_results$intensity_cluster)
  }
  
  # Preparazione del dataframe per il plot
  plot_df <- data.frame()
  
  # Per ogni cluster
  for (clust in levels(sim_results$intensity_cluster)) {
    # Seleziona tutti i punti di questo cluster
    idx <- which(sim_results$intensity_cluster == clust)
    
    # Campiona per efficienza
    if (length(idx) > n_points) {
      set.seed(42)  # Per riproducibilità
      idx <- sample(idx, n_points)
    }
    
    # Estrai coordinate per questo subset
    coords <- as.matrix(sim_results$coordinates[idx, ])
    
    # Calcola matrice di distanza
    dist_mat <- fields::rdist(coords)
    
    # Calcola distanza media per ogni punto
    mean_dist <- rowMeans(dist_mat)
    
    # Crea dataframe per questo cluster
    cluster_df <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      mean_distance = mean_dist,
      cluster_label = paste0("cells_", tolower(clust))
    )
    
    # Aggiungi al dataframe principale
    plot_df <- rbind(plot_df, cluster_df)
  }
  
  # Crea il plot
  if (nrow(plot_df) > 0) {
    p <- ggplot(plot_df, aes(x = x, y = y, color = mean_distance)) +
      geom_point(size = 1.2) +
      facet_wrap(~ cluster_label) +
      scale_color_viridis_c(option = "inferno", direction = -1) +
      theme_minimal() +
      labs(title = title, 
           subtitle = subtitle,
           x = "x", 
           y = "y", 
           color = "Distanza media")
    
    return(p)
  } else {
    warning("Nessun dato per il plot")
    return(NULL)
  }
}

#' Plot della distribuzione dell'espressione genica
#'
#' Visualizza gli istogrammi di distribuzione dei conteggi di espressione,
#' separando geni stabili e variabili o marker e non marker.
#'
#' @param sim_results Lista di risultati della simulazione
#' @param marker_genes Vettore di ID dei geni marker 
#' @param title Titolo del grafico
#' @param bins Numero di bin per gli istogrammi
#' @param max_count Limite superiore dei conteggi da mostrare
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap theme_minimal labs scale_fill_manual
#' @export
plot_expression_distribution <- function(sim_results, marker_genes = NULL,
                                        title = "Distribuzione dell'espressione genica",
                                        subtitle = "Confronto tra geni stabili (sub-Poisson) e variabili (Negative Binomial)",
                                        bins = 50, max_count = 200) {
  # Verifica che la matrice di espressione sia nel formato corretto
  if (!is.matrix(sim_results$expression)) {
    stop("L'oggetto expression deve essere una matrice")
  }
  
  # Se marker_genes è NULL, cerca di identificarli
  if (is.null(marker_genes) && !is.null(sim_results$parameters$marker_params$marker_genes_per_type)) {
    marker_list <- sim_results$parameters$marker_params$marker_genes_per_type
    marker_genes <- unlist(marker_list)
  }
  
  # Calcola statistiche per ogni gene
  gene_means <- rowMeans(sim_results$expression)
  gene_vars <- apply(sim_results$expression, 1, var)
  
  # Calcola la dispersione (varianza/media) per distinguere i tipi di geni
  dispersion <- gene_vars / pmax(gene_means, 1e-5)
  
  # Dividi geni in stabili (sotto-Poisson) e variabili (sopra-Poisson)
  stable_genes <- dispersion <= 1.2  # Soglia arbitraria vicina a 1
  variable_genes <- dispersion > 1.2
  
  # Crea categorie
  if (!is.null(marker_genes)) {
    is_marker <- rownames(sim_results$expression) %in% marker_genes
    
    # Dividi in 4 categorie
    gene_categories <- rep("Variable (Non-marker)", nrow(sim_results$expression))
    gene_categories[stable_genes & !is_marker] <- "Stable (Non-marker)"
    gene_categories[stable_genes & is_marker] <- "Stable (Marker)"
    gene_categories[variable_genes & is_marker] <- "Variable (Marker)"
    
    # Ordine personalizzato
    category_order <- c("Stable (Non-marker)", "Stable (Marker)", 
                        "Variable (Non-marker)", "Variable (Marker)")
  } else {
    # Dividi in 2 categorie
    gene_categories <- ifelse(stable_genes, "Stable (sub-Poisson)", "Variable (Negative Binomial)")
    category_order <- c("Stable (sub-Poisson)", "Variable (Negative Binomial)")
  }
  
  # Trasponi per avere conteggi × categorie
  expr_t <- t(sim_results$expression)
  expression_df <- data.frame()
  
  # Per ogni categoria
  for (category in unique(gene_categories)) {
    cat_idx <- which(gene_categories == category)
    if (length(cat_idx) > 0) {
      # Estrai conteggi per questa categoria
      counts <- as.vector(expr_t[, cat_idx])
      
      # Limita ai conteggi non-zero e sotto max_count
      counts <- counts[counts > 0 & counts <= max_count]
      
      # Aggiungi al dataframe
      if (length(counts) > 0) {
        df <- data.frame(
          counts = counts,
          type = category
        )
        expression_df <- rbind(expression_df, df)
      }
    }
  }
  
  # Crea il plot
  if (nrow(expression_df) > 0) {
    # Converti a fattore con ordine corretto
    expression_df$type <- factor(expression_df$type, levels = category_order)
    
    # Crea il plot
    p <- ggplot(expression_df, aes(x = counts, fill = type)) +
      geom_histogram(bins = bins) +
      facet_wrap(~ type, scales = "free_y", ncol = 1) +
      theme_minimal() +
      labs(title = title, 
           subtitle = subtitle,
           x = "Counts", 
           y = "Frequenza",
           fill = "Tipo di Gene") +
      scale_fill_manual(values = c("Stable (sub-Poisson)" = "#8FBC8F", 
                                  "Variable (Negative Binomial)" = "#87CEFA",
                                  "Stable (Non-marker)" = "#8FBC8F", 
                                  "Stable (Marker)" = "#32CD32",
                                  "Variable (Non-marker)" = "#87CEFA", 
                                  "Variable (Marker)" = "#1E90FF"))
    
    return(p)
  } else {
    warning("Nessun dato per il plot")
    return(NULL)
  }
}

#' Genera e Salva Tutti i Plot di Validazione per un Dataset Simulato
#'
#' @param sim_results Lista di risultati della simulazione
#' @param marker_genes Vettore di ID dei geni marker da visualizzare (default: NULL, seleziona automaticamente)
#' @param n_markers Numero di marker per tipo cellulare se marker_genes è NULL
#' @param output_dir Directory in cui salvare i plot
#' @param file_prefix Prefisso per i nomi dei file
#' @param file_format Formato dei file ("png", "pdf", "svg", etc.)
#' @param width Larghezza dei plot in pollici
#' @param height Altezza dei plot in pollici
#' @param skip_dim_reduction Se TRUE, salta la creazione del plot di riduzione dimensionale
#' @return Invisibilmente, una lista con tutti gli oggetti ggplot2
#' @importFrom ggplot2 ggsave
#' @export
generate_validation_plots <- function(sim_results, marker_genes = NULL, n_markers = 5,
                                    output_dir = "plots", file_prefix = "validation", 
                                    file_format = "png", width = 8, height = 6,
                                    skip_dim_reduction = FALSE) {
  # Verifica e crea la directory se necessario
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Identificazione automatica dei geni marker se non specificati
  marker_list <- NULL
  if (is.null(marker_genes)) {
    marker_list <- identify_marker_genes(sim_results, n_markers = n_markers)
    marker_genes <- unique(unlist(marker_list))
  }
  
  # Lista per contenere tutti i plot
  all_plots <- list()
  
  # 1. Plot della relazione media-varianza
  mean_var_plot <- plot_mean_variance_relationship(sim_results)
  all_plots[["mean_variance"]] <- mean_var_plot
  filename <- file.path(output_dir, paste0(file_prefix, "_mean_variance.", file_format))
  ggsave(filename, mean_var_plot, width = width, height = height)
  
  # 2. Plot della relazione dropout-espressione
  dropout_plot <- plot_dropout_vs_expression(sim_results)
  all_plots[["dropout"]] <- dropout_plot
  filename <- file.path(output_dir, paste0(file_prefix, "_dropout.", file_format))
  ggsave(filename, dropout_plot, width = width, height = height)
  
  # 3. Plot dimensionalità ridotta (UMAP) - solo se richiesto
  if (!skip_dim_reduction) {
    dim_reduction_plot <- tryCatch({
      plot_dimensionality_reduction(sim_results)
    }, error = function(e) {
      warning("Could not create dimensionality reduction plot: ", e$message)
      return(NULL)
    })
    
    if (!is.null(dim_reduction_plot)) {
      all_plots[["dim_reduction"]] <- dim_reduction_plot
      filename <- file.path(output_dir, paste0(file_prefix, "_dim_reduction.", file_format))
      ggsave(filename, dim_reduction_plot, width = width, height = height)
    }
  }
  
  # 4. Plot spaziali per i geni marker individuali
  marker_plots <- list()
  for (gene_id in marker_genes) {
    marker_plot <- tryCatch({
      plot_marker_gene_spatial(sim_results, gene_id)
    }, error = function(e) {
      warning("Could not create marker plot for gene ", gene_id, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(marker_plot)) {
      marker_plots[[gene_id]] <- marker_plot
      # Crea un nome file sicuro (sostituzione di caratteri non validi)
      safe_gene_id <- gsub("[^a-zA-Z0-9]", "_", gene_id)
      filename <- file.path(output_dir, paste0(file_prefix, "_marker_", safe_gene_id, ".", file_format))
      ggsave(filename, marker_plot, width = width, height = height)
    }
  }
  all_plots[["marker_plots"]] <- marker_plots
  
  # 5. Plot di espressione spaziale per tutti i cluster
  spatial_expr_plot <- tryCatch({
    plot_spatial_expression_panel(sim_results, marker_list = marker_list)
  }, error = function(e) {
    warning("Could not create spatial expression panel: ", e$message)
    return(NULL)
  })
  
  if (!is.null(spatial_expr_plot)) {
    all_plots[["spatial_expression"]] <- spatial_expr_plot
    filename <- file.path(output_dir, paste0(file_prefix, "_spatial_expression.", file_format))
    ggsave(filename, spatial_expr_plot, width = width * 1.5, height = height * 1.2)
  }
  
  # 6. Plot delle distanze intra-cluster
  distances_plot <- tryCatch({
    plot_distances_by_cluster(sim_results)
  }, error = function(e) {
    warning("Could not create distances by cluster plot: ", e$message)
    return(NULL)
  })
  
  if (!is.null(distances_plot)) {
    all_plots[["distances_by_cluster"]] <- distances_plot
    filename <- file.path(output_dir, paste0(file_prefix, "_distances_by_cluster.", file_format))
    ggsave(filename, distances_plot, width = width * 1.5, height = height * 1.2)
  }
  
  # 7. Plot delle distribuzioni di espressione
  expr_dist_plot <- tryCatch({
    plot_expression_distribution(sim_results, marker_genes = marker_genes)
  }, error = function(e) {
    warning("Could not create expression distribution plot: ", e$message)
    return(NULL)
  })
  
  if (!is.null(expr_dist_plot)) {
    all_plots[["expression_distribution"]] <- expr_dist_plot
    filename <- file.path(output_dir, paste0(file_prefix, "_expression_distribution.", file_format))
    ggsave(filename, expr_dist_plot, width = width, height = height * 1.2)
  }
  
  return(invisible(all_plots))
}