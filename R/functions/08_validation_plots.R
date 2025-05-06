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
  
  # Verifica se la matrice è sparsa o densa
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  
  # Per matrici sparse, non serve controllare se è una matrice standard
  if (!is_sparse && !is.matrix(sim_results$expression)) {
    stop("L'oggetto expression deve essere una matrice o una matrice sparsa")
  }
  
  # Controllo dimensioni (funziona sia per matrici sparse che dense)
  if (ncol(sim_results$expression) != nrow(coords)) {
    warning("La matrice di espressione non ha il formato corretto (geni x celle). Tentativo di trasposizione...")
    # Trasponi la matrice se ha la dimensione sbagliata
    if (nrow(sim_results$expression) == nrow(coords)) {
      sim_results$expression <- t(sim_results$expression)
    } else {
      stop(paste("Dimensioni incompatibili: matrice expression", 
                 nrow(sim_results$expression), "x", ncol(sim_results$expression), 
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
  
  # Estrazione dell'espressione (funziona sia per matrici sparse che dense)
  if (is_sparse) {
    # Converti solo questa riga specifica da sparsa a vettore numerico standard
    gene_expr <- as.numeric(as.matrix(sim_results$expression[gene_idx, , drop = FALSE]))
  } else {
    gene_expr <- sim_results$expression[gene_idx, ]
  }
  
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
      gene_name <- if (!is.null(rownames(sim_results$expression)) && gene_idx <= length(rownames(sim_results$expression))) {
        rownames(sim_results$expression)[gene_idx]
      } else {
        paste0("gene_", gene_idx)
      }
      title <- paste("Espressione spaziale del gene", gene_name)
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
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
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
                                          highlight_markers = TRUE, title = NULL,
                                          max_points = 5000) {
  # Verifica se la matrice è sparsa
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  
  # Calcolo della media e varianza per ogni gene, ottimizzato per matrici sparse
  if (is_sparse) {
    # Per matrici sparse, usa le funzioni di Matrix più efficienti
    require(Matrix)
    gene_means <- Matrix::rowMeans(sim_results$expression)
    
    # Per la varianza, possiamo usare una formula equivalente:
    # Var(X) = E(X^2) - (E(X))^2
    # Calcola E(X^2) in modo efficiente per matrici sparse
    expr2 <- sim_results$expression^2
    
    # Usa rowMeans su X^2
    mean_squared <- Matrix::rowMeans(expr2)
    
    # Calcola la varianza usando la formula
    gene_vars <- mean_squared - gene_means^2
  } else {
    # Per matrici dense, usa le funzioni standard di R
    gene_means <- rowMeans(sim_results$expression)
    gene_vars <- apply(sim_results$expression, 1, var)
  }
  
  # Identificazione dei geni marker, se richiesto
  markers <- NULL
  if (highlight_markers && !is.null(sim_results$parameters$marker_params$marker_genes_per_type)) {
    # Estrazione dei geni marker dai parametri
    markers <- sim_results$parameters$marker_params$marker_genes_per_type
    markers <- unlist(markers)
  }
  
  # Campionamento per dataset grandi
  if (length(gene_means) > max_points) {
    set.seed(42)  # Per riproducibilità
    if (!is.null(markers) && length(markers) > 0) {
      # Assicurati di mantenere tutti i marker
      marker_indices <- which(rownames(sim_results$expression) %in% markers)
      non_marker_indices <- setdiff(1:length(gene_means), marker_indices)
      
      # Campiona dai non-marker
      n_to_sample <- min(max_points - length(marker_indices), length(non_marker_indices))
      if (n_to_sample > 0) {
        sampled_non_markers <- sample(non_marker_indices, n_to_sample)
        indices_to_keep <- c(marker_indices, sampled_non_markers)
      } else {
        indices_to_keep <- marker_indices[1:min(max_points, length(marker_indices))]
      }
    } else {
      # Campionamento casuale
      indices_to_keep <- sample(1:length(gene_means), max_points)
    }
    
    gene_means <- gene_means[indices_to_keep]
    gene_vars <- gene_vars[indices_to_keep]
    
    if (!is.null(rownames(sim_results$expression))) {
      gene_names <- rownames(sim_results$expression)[indices_to_keep]
    } else {
      gene_names <- paste0("gene_", indices_to_keep)
    }
  } else {
    gene_names <- rownames(sim_results$expression)
    if (is.null(gene_names)) {
      gene_names <- paste0("gene_", seq_along(gene_means))
    }
  }
  
  # Creazione di un dataframe per il plot
  plot_df <- data.frame(
    gene = gene_names,
    mean = as.numeric(gene_means),  # Assicura conversione a vettore numerico
    variance = as.numeric(gene_vars),  # Assicura conversione a vettore numerico
    is_marker = gene_names %in% markers
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
         y = "Varianza dell'Espressione") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
  # Aggiunta dei punti, eventualmente colorati per marker
  if (highlight_markers && !is.null(markers) && any(plot_df$is_marker)) {
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
#' @param max_points Numero massimo di punti da visualizzare (per efficienza)
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_x_log10
#'                     theme_minimal labs
#' @export
plot_dropout_vs_expression <- function(sim_results, log_mean = TRUE, 
                                     highlight_markers = TRUE, title = NULL,
                                     max_points = 5000) {
  # Verifica se la matrice è sparsa
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  
  # Calcolo della frazione di zeri e medie, ottimizzato per matrici sparse
  if (is_sparse) {
    require(Matrix)
    
    # Per matrici sparse, il calcolo degli zeri è efficiente usando caratteristiche della matrice sparsa
    # La frazione di zeri è 1 - (numero di elementi non-zero / numero totale di elementi)
    nnz_per_row <- Matrix::rowSums(sim_results$expression != 0)
    total_cols <- ncol(sim_results$expression)
    zero_fraction <- 1 - (nnz_per_row / total_cols)
    
    # Calcolo media per riga
    gene_means <- Matrix::rowMeans(sim_results$expression)
  } else {
    # Per matrici dense, usa le funzioni standard
    zero_fraction <- rowMeans(sim_results$expression == 0)
    gene_means <- rowMeans(sim_results$expression)
  }
  
  # Log-trasformazione delle medie, se richiesto
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
  
  # Campionamento per dataset grandi
  if (length(gene_means) > max_points) {
    set.seed(42)  # Per riproducibilità
    if (!is.null(markers) && length(markers) > 0) {
      # Assicurati di mantenere tutti i marker
      marker_indices <- which(rownames(sim_results$expression) %in% markers)
      non_marker_indices <- setdiff(1:length(gene_means), marker_indices)
      
      # Campiona dai non-marker
      n_to_sample <- min(max_points - length(marker_indices), length(non_marker_indices))
      if (n_to_sample > 0) {
        sampled_non_markers <- sample(non_marker_indices, n_to_sample)
        indices_to_keep <- c(marker_indices, sampled_non_markers)
      } else {
        indices_to_keep <- marker_indices[1:min(max_points, length(marker_indices))]
      }
    } else {
      # Campionamento casuale
      indices_to_keep <- sample(1:length(gene_means), max_points)
    }
    
    gene_means <- gene_means[indices_to_keep]
    zero_fraction <- zero_fraction[indices_to_keep]
    
    if (!is.null(rownames(sim_results$expression))) {
      gene_names <- rownames(sim_results$expression)[indices_to_keep]
    } else {
      gene_names <- paste0("gene_", indices_to_keep)
    }
  } else {
    gene_names <- rownames(sim_results$expression)
    if (is.null(gene_names)) {
      gene_names <- paste0("gene_", seq_along(gene_means))
    }
  }
  
  # Creazione di un dataframe per il plot
  plot_df <- data.frame(
    gene = gene_names,
    mean_expression = as.numeric(gene_means),  # Assicura conversione a vettore numerico
    zero_fraction = as.numeric(zero_fraction),  # Assicura conversione a vettore numerico
    is_marker = gene_names %in% markers
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
         y = "Frazione di Zeri (Dropout Rate)") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
  # Aggiunta dei punti, eventualmente colorati per marker
  if (highlight_markers && !is.null(markers) && any(plot_df$is_marker)) {
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
         color = "Tipo Cellulare") +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
  return(p)
}

#' Identifica Automaticamente i Geni Marker per ogni Tipo Cellulare
#'
#' @param sim_results Lista di risultati della simulazione
#' @param n_markers Numero di marker per tipo cellulare da identificare
#' @param min_ratio Rapporto minimo di espressione tra cluster target e altri cluster
#' @return Lista con i geni marker per ogni tipo cellulare
#' @importFrom stats aggregate
#' @export
identify_marker_genes <- function(sim_results, n_markers = 5, min_ratio = 1.5) {
  # Estrazione dei tipi cellulari unici
  cell_types <- levels(sim_results$intensity_cluster)
  
  # Verifica se la matrice è sparsa
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  
  # Calcolo dell'espressione media per ogni cluster e gene in modo efficiente
  # usando un metodo diverso per matrici sparse
  if (is_sparse) {
    require(Matrix)
    
    # Preparazione dell'espressione media per cluster
    mean_expr_by_cluster <- list()
    
    # Ottieni le dimensioni senza trasporre l'intera matrice
    n_genes <- nrow(sim_results$expression)
    n_cells <- ncol(sim_results$expression)
    
    for (cluster in cell_types) {
      # Indici delle celle appartenenti a questo cluster
      idx <- which(sim_results$intensity_cluster == cluster)
      
      if (length(idx) > 0) {
        # Selezione solo delle colonne relative a questo cluster
        # e calcolo della media per gene
        cluster_expr <- sim_results$expression[, idx, drop = FALSE]
        mean_expr_by_cluster[[cluster]] <- Matrix::rowMeans(cluster_expr)
      } else {
        mean_expr_by_cluster[[cluster]] <- rep(0, n_genes)
      }
    }
  } else {
    # Metodo originale per matrici dense
    expr_t <- t(sim_results$expression)
    expr_df <- data.frame(expr_t, cluster = sim_results$intensity_cluster)
    
    mean_expr_by_cluster <- list()
    for (cluster in cell_types) {
      cluster_cells <- expr_df$cluster == cluster
      if (sum(cluster_cells) > 0) {
        cluster_means <- colMeans(expr_df[cluster_cells, -ncol(expr_df), drop = FALSE])
        mean_expr_by_cluster[[cluster]] <- cluster_means
      } else {
        mean_expr_by_cluster[[cluster]] <- rep(0, ncol(expr_df) - 1)
      }
    }
  }
  
  # Identificazione dei geni marker per ogni cluster con criterio più restrittivo
  marker_genes <- list()
  for (i in seq_along(cell_types)) {
    cluster <- cell_types[i]
    # Calcolo del rapporto di espressione rispetto agli altri cluster
    ratios <- list()
    for (j in seq_along(cell_types)) {
      if (i != j) {
        other_cluster <- cell_types[j]
        # Rapporto di espressione tra il cluster corrente e un altro cluster
        # Aggiunta di pseudocount per evitare divisione per zero
        ratio <- (mean_expr_by_cluster[[cluster]] + 0.1) / (mean_expr_by_cluster[[other_cluster]] + 0.1)
        ratios[[other_cluster]] <- ratio
      }
    }
    
    # Media dei rapporti rispetto a tutti gli altri cluster
    if (length(ratios) > 0) {
      avg_ratio <- Reduce(`+`, ratios) / length(ratios)
      
      # Assicurati che avg_ratio sia un vettore numerico standard
      avg_ratio <- as.numeric(avg_ratio)  
      
      # Assegna nomi se necessario
      if (is.null(names(avg_ratio))) {
        if (!is.null(rownames(sim_results$expression))) {
          names(avg_ratio) <- rownames(sim_results$expression)
        } else {
          names(avg_ratio) <- paste0("gene_", seq_along(avg_ratio))
        }
      }
      
      # Seleziona solo geni con un rapporto minimo di espressione
      high_ratio_genes <- which(avg_ratio >= min_ratio)
      
      # Ordina per rapporto e seleziona i top n_markers
      if (length(high_ratio_genes) > 0) {
        ordered_genes <- high_ratio_genes[order(avg_ratio[high_ratio_genes], decreasing = TRUE)]
        top_genes <- ordered_genes[1:min(n_markers, length(ordered_genes))]
        marker_genes[[cluster]] <- names(avg_ratio)[top_genes]
      } else {
        # Se non ci sono geni con rapporto sufficientemente alto, usa i top n comunque
        top_genes <- order(avg_ratio, decreasing = TRUE)[1:min(n_markers, length(avg_ratio))]
        marker_genes[[cluster]] <- names(avg_ratio)[top_genes]
      }
    } else {
      marker_genes[[cluster]] <- character(0)
    }
  }
  
  return(marker_genes)
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
    marker_list <- identify_marker_genes(sim_results, n_markers = 1, min_ratio = 1.5)
  }
  
  # Estrazione delle coordinate spaziali
  coords <- sim_results$coordinates
  
  # Verifica se la matrice è sparsa o densa
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  
  # Per matrici sparse, non serve controllare se è una matrice standard
  if (!is_sparse && !is.matrix(sim_results$expression)) {
    stop("L'oggetto expression deve essere una matrice o una matrice sparsa")
  }
  
  # Controllo dimensioni (funziona sia per matrici sparse che dense)
  if (ncol(sim_results$expression) != nrow(coords)) {
    warning("La matrice di espressione non ha il formato corretto (geni x celle). Tentativo di trasposizione...")
    # Trasponi la matrice se ha la dimensione sbagliata
    if (nrow(sim_results$expression) == nrow(coords)) {
      sim_results$expression <- t(sim_results$expression)
    } else {
      stop(paste("Dimensioni incompatibili: matrice expression", 
                 nrow(sim_results$expression), "x", ncol(sim_results$expression), 
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
      
      # Estrazione dell'espressione (funziona sia per matrici sparse che dense)
      if (is_sparse) {
        # Converti solo questa riga specifica da sparsa a vettore numerico standard
        gene_expr <- as.numeric(as.matrix(sim_results$expression[gene_idx, , drop = FALSE]))
      } else {
        gene_expr <- sim_results$expression[gene_idx, ]
      }
      
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
           color = "Espressione") +
      theme(panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA))
    
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
plot_distances_by_cluster <- function(sim_results, n_points = 2000,
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
           color = "Distanza media") +
      theme(panel.background = element_rect(fill = "white", colour = NA),
            plot.background = element_rect(fill = "white", colour = NA))
    
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
#' @param max_sample Numero massimo di punti da campionare (per efficienza)
#' @return Un oggetto ggplot2
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap theme_minimal labs scale_fill_manual
#' @export
plot_expression_distribution <- function(sim_results, marker_genes = NULL,
                                        title = "Distribuzione dell'espressione genica",
                                        subtitle = "Confronto tra geni a bassa e alta dispersione",
                                        bins = 30, max_count = 50, max_sample = 20000) {
  # Verifica se la matrice è sparsa o densa
  is_sparse <- inherits(sim_results$expression, "sparseMatrix")
  
  # Per matrici sparse, non serve controllare se è una matrice standard
  if (!is_sparse && !is.matrix(sim_results$expression)) {
    stop("L'oggetto expression deve essere una matrice o una matrice sparsa")
  }
  
  # Calcolo della media e varianza per ogni gene, ottimizzato per matrici sparse
  if (is_sparse) {
    require(Matrix)
    # Calcolo media
    gene_means <- Matrix::rowMeans(sim_results$expression)
    
    # Per la varianza, utilizziamo la formula: Var(X) = E(X^2) - (E(X))^2
    expr2 <- sim_results$expression^2
    mean_squared <- Matrix::rowMeans(expr2)
    gene_vars <- as.numeric(mean_squared - gene_means^2)
    gene_means <- as.numeric(gene_means)
  } else {
    # Per matrici dense, usa le funzioni standard
    gene_means <- rowMeans(sim_results$expression)
    gene_vars <- apply(sim_results$expression, 1, var)
  }
  
  # Calcola la dispersione e classifica i geni
  dispersion <- gene_vars / pmax(gene_means, 1e-5)
  dispersion_median <- median(dispersion)
  
  # Categorie basate sulla dispersione
  gene_categories <- ifelse(dispersion <= dispersion_median, 
                          "Bassa dispersione", "Alta dispersione")
  category_order <- c("Bassa dispersione", "Alta dispersione")
  
  # Campionamento per dataset grandi
  if (nrow(sim_results$expression) > 1000 || ncol(sim_results$expression) > 10000) {
    # Campiona un sottoinsieme di celle e geni
    set.seed(42)  # Per riproducibilità
    
    # Seleziona fino a 200 geni casuali
    n_genes_sample <- min(200, nrow(sim_results$expression))
    gene_sample <- sample(1:nrow(sim_results$expression), n_genes_sample)
    
    # Seleziona un campione limitato di celle
    n_cells_sample <- min(5000, ncol(sim_results$expression))
    cell_sample <- sample(1:ncol(sim_results$expression), n_cells_sample)
    
    # Crea un dataframe per i valori di espressione campionati
    expr_values <- c()
    categories <- c()
    
    # Per ogni categoria, estrai e campiona valori
    for (category in category_order) {
      # Genes in this category
      cat_genes <- gene_sample[gene_categories[gene_sample] == category]
      
      if (length(cat_genes) > 0) {
        # Extract expression values for these genes
        for (gene_idx in cat_genes) {
          # Estrai valori in modo efficiente per matrice sparsa o densa
          if (is_sparse) {
            values <- as.numeric(as.matrix(sim_results$expression[gene_idx, cell_sample, drop = FALSE]))
          } else {
            values <- sim_results$expression[gene_idx, cell_sample]
          }
          
          values <- values[values > 0 & values <= max_count]
          
          # Limita ulteriormente se necessario
          if (length(values) > 5000) {
            values <- sample(values, 5000)
          }
          
          if (length(values) > 0) {
            expr_values <- c(expr_values, values)
            categories <- c(categories, rep(category, length(values)))
          }
        }
      }
    }
    
    # Costruisci il dataframe per il plot
    plot_data <- data.frame(
      counts = expr_values,
      type = categories
    )
  } else {
    # Se il dataset è piccolo, estrai tutti i valori
    # Ma gestisci in modo diverso per matrici sparse
    if (is_sparse) {
      plot_data <- data.frame()
      
      # Estrazione per categoria in modo efficiente
      for (category in category_order) {
        cat_idx <- which(gene_categories == category)
        
        if (length(cat_idx) > 0) {
          # Seleziona solo queste righe della matrice sparsa
          submatrix <- sim_results$expression[cat_idx, , drop = FALSE]
          
          # Converti a matrice densa (solo per queste righe)
          dense_submatrix <- as.matrix(submatrix)
          
          # Estrai i valori non-zero che soddisfano il criterio
          counts <- as.vector(dense_submatrix)
          counts <- counts[counts > 0 & counts <= max_count]
          
          if (length(counts) > 0) {
            df <- data.frame(
              counts = counts,
              type = category
            )
            plot_data <- rbind(plot_data, df)
          }
        }
      }
    } else {
      # Per matrici dense, usa l'approccio originale
      expr_t <- t(sim_results$expression)
      plot_data <- data.frame()
      
      for (category in category_order) {
        cat_idx <- which(gene_categories == category)
        if (length(cat_idx) > 0) {
          counts <- as.vector(expr_t[, cat_idx])
          counts <- counts[counts > 0 & counts <= max_count]
          
          if (length(counts) > 0) {
            df <- data.frame(
              counts = counts,
              type = category
            )
            plot_data <- rbind(plot_data, df)
          }
        }
      }
    }
  }
  
  # Limita la dimensione totale del dataframe
  if (nrow(plot_data) > max_sample) {
    set.seed(42)
    plot_data <- plot_data[sample(1:nrow(plot_data), max_sample), ]
  }
  
  # Verifica se ci sono dati per il plot
  if (nrow(plot_data) == 0) {
    warning("Nessun dato disponibile per il plot di distribuzione dell'espressione")
    return(NULL)
  }
  
  # Converti a fattore con ordine corretto
  plot_data$type <- factor(plot_data$type, levels = category_order)
  
  # Crea il plot
  p <- ggplot(plot_data, aes(x = counts, fill = type)) +
    geom_histogram(bins = bins, position = "identity", alpha = 0.7) +
    facet_wrap(~ type, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(title = title, 
         subtitle = subtitle,
         x = "Conteggi di Espressione", 
         y = "Frequenza",
         fill = "Tipo di Gene") +
    scale_fill_manual(values = c("Bassa dispersione" = "#4DAF4A", 
                                "Alta dispersione" = "#377EB8")) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
  
  return(p)
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
    marker_list <- identify_marker_genes(sim_results, n_markers = n_markers, min_ratio = 1.5)
    marker_genes <- unique(unlist(marker_list))
  }
  
  # Lista per contenere tutti i plot
  all_plots <- list()
  
  # 1. Plot della relazione media-varianza
  mean_var_plot <- plot_mean_variance_relationship(sim_results, max_points = 5000)
  all_plots[["mean_variance"]] <- mean_var_plot
  filename <- file.path(output_dir, paste0(file_prefix, "_mean_variance.", file_format))
  ggsave(filename, mean_var_plot, width = width, height = height, bg = "white")
  
  # 2. Plot della relazione dropout-espressione
  dropout_plot <- plot_dropout_vs_expression(sim_results, max_points = 5000)
  all_plots[["dropout"]] <- dropout_plot
  filename <- file.path(output_dir, paste0(file_prefix, "_dropout.", file_format))
  ggsave(filename, dropout_plot, width = width, height = height, bg = "white")
  
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
      ggsave(filename, dim_reduction_plot, width = width, height = height, bg = "white")
    }
  }
  
  # 4. Plot spaziali per i geni marker individuali (limitati a massimo 3 per evitare sovraccarico)
  marker_plots <- list()
  if (length(marker_genes) > 3) {
    marker_genes <- marker_genes[1:3]
  }
  
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
      filename <- file.path(output_dir, paste0(file_prefix, "_marker_gene_", safe_gene_id, ".", file_format))
      ggsave(filename, marker_plot, width = width, height = height, bg = "white")
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
    ggsave(filename, spatial_expr_plot, width = width * 1.5, height = height * 1.2, bg = "white")
  }
  
  # 6. Plot delle distanze intra-cluster
  distances_plot <- tryCatch({
    plot_distances_by_cluster(sim_results, n_points = 2000)
  }, error = function(e) {
    warning("Could not create distances by cluster plot: ", e$message)
    return(NULL)
  })
  
  if (!is.null(distances_plot)) {
    all_plots[["distances_by_cluster"]] <- distances_plot
    filename <- file.path(output_dir, paste0(file_prefix, "_distances_by_cluster.", file_format))
    ggsave(filename, distances_plot, width = width * 1.5, height = height * 1.2, bg = "white")
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
    ggsave(filename, expr_dist_plot, width = width, height = height * 1.2, bg = "white")
  }
  
  # Crea un file README.md nella directory
  readme_text <- paste(
    "# Plot di Validazione",
    "",
    "Questa directory contiene i plot di validazione generati per valutare la qualità della simulazione spaziale trascrittomica. Di seguito una breve descrizione di ciascun tipo di plot:",
    "",
    "## Plot Principali",
    "",
    "1. spatial_clusters.png: Visualizzazione spaziale dei cluster cellulari identificati",
    "2. mean_variance.png: Relazione tra media e varianza dell'espressione genica",
    "3. dropout.png: Relazione tra tasso di dropout e livello di espressione genica",
    "4. expression_distribution.png: Distribuzione dell'espressione genica per geni a bassa e alta dispersione",
    "",
    "## Plot Specifici",
    "",
    "1. marker_gene_*.png: Visualizzazione spaziale dell'espressione di geni marker specifici",
    "2. *_distances_by_cluster.png: Distanza media tra le celle all'interno di ciascun cluster",
    "3. *_spatial_expression.png: Pannello che mostra l'espressione spaziale di un gene marker per ogni cluster",
    "",
    sep = "\n"
  )
  
  writeLines(readme_text, file.path(output_dir, "README.md"))
  
  return(invisible(all_plots))
}