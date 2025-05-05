test_that("plot_marker_gene_spatial funziona correttamente", {
  # Creazione di un dataset di simulazione semplificato
  n_cells <- 100
  n_genes <- 50
  k_cell_types <- 3
  
  # Coordinate spaziali casuali
  set.seed(123)
  coordinates <- data.frame(
    x = runif(n_cells, 0, 100),
    y = runif(n_cells, 0, 100)
  )
  
  # Cluster casuali
  intensity_cluster <- factor(sample(1:k_cell_types, n_cells, replace = TRUE))
  levels(intensity_cluster) <- paste0("Type_", 1:k_cell_types)
  
  # Matrice di espressione casuale
  expression <- matrix(rpois(n_cells * n_genes, lambda = 5), nrow = n_genes, ncol = n_cells)
  rownames(expression) <- paste0("gene_", 1:n_genes)
  
  # Risultati simulati
  sim_results <- list(
    coordinates = coordinates,
    intensity_cluster = intensity_cluster,
    expression = expression,
    parameters = list(
      marker_params = list(
        marker_genes_per_type = list(
          Type_1 = c("gene_1", "gene_2"),
          Type_2 = c("gene_3", "gene_4"),
          Type_3 = c("gene_5", "gene_6")
        )
      )
    )
  )
  
  # Test con nome gene
  p1 <- plot_marker_gene_spatial(sim_results, "gene_1")
  expect_s3_class(p1, "ggplot")
  
  # Test con indice gene
  p2 <- plot_marker_gene_spatial(sim_results, 1)
  expect_s3_class(p2, "ggplot")
  
  # Test con log_transform = FALSE
  p3 <- plot_marker_gene_spatial(sim_results, "gene_1", log_transform = FALSE)
  expect_s3_class(p3, "ggplot")
  
  # Test con titolo personalizzato
  p4 <- plot_marker_gene_spatial(sim_results, "gene_1", title = "Test Plot")
  expect_s3_class(p4, "ggplot")
  expect_equal(p4$labels$title, "Test Plot")
  
  # Test con scala di colori personalizzata
  p5 <- plot_marker_gene_spatial(sim_results, "gene_1", color_scale = "magma")
  expect_s3_class(p5, "ggplot")
  
  # Test errore per gene inesistente
  expect_error(plot_marker_gene_spatial(sim_results, "gene_999"))
  expect_error(plot_marker_gene_spatial(sim_results, 999))
})

test_that("plot_mean_variance_relationship funziona correttamente", {
  # Usa lo stesso dataset di simulazione
  n_cells <- 100
  n_genes <- 50
  k_cell_types <- 3
  
  # Coordinate spaziali casuali
  set.seed(123)
  coordinates <- data.frame(
    x = runif(n_cells, 0, 100),
    y = runif(n_cells, 0, 100)
  )
  
  # Cluster casuali
  intensity_cluster <- factor(sample(1:k_cell_types, n_cells, replace = TRUE))
  levels(intensity_cluster) <- paste0("Type_", 1:k_cell_types)
  
  # Matrice di espressione casuale con Negative Binomial per simulare sovradispersione
  expression <- matrix(rnbinom(n_cells * n_genes, mu = 5, size = 2), nrow = n_genes, ncol = n_cells)
  rownames(expression) <- paste0("gene_", 1:n_genes)
  
  # Risultati simulati
  sim_results <- list(
    coordinates = coordinates,
    intensity_cluster = intensity_cluster,
    expression = expression,
    parameters = list(
      marker_params = list(
        marker_genes_per_type = list(
          Type_1 = c("gene_1", "gene_2"),
          Type_2 = c("gene_3", "gene_4"),
          Type_3 = c("gene_5", "gene_6")
        )
      )
    )
  )
  
  # Test base
  p1 <- plot_mean_variance_relationship(sim_results)
  expect_s3_class(p1, "ggplot")
  
  # Test senza scala logaritmica
  p2 <- plot_mean_variance_relationship(sim_results, log_scale = FALSE)
  expect_s3_class(p2, "ggplot")
  
  # Test senza evidenziazione dei marker
  p3 <- plot_mean_variance_relationship(sim_results, highlight_markers = FALSE)
  expect_s3_class(p3, "ggplot")
  
  # Test con titolo personalizzato
  p4 <- plot_mean_variance_relationship(sim_results, title = "Test Plot")
  expect_s3_class(p4, "ggplot")
  expect_equal(p4$labels$title, "Test Plot")
})

test_that("plot_dropout_vs_expression funziona correttamente", {
  # Usa lo stesso dataset di simulazione ma con dropout
  n_cells <- 100
  n_genes <- 50
  k_cell_types <- 3
  
  # Coordinate spaziali casuali
  set.seed(123)
  coordinates <- data.frame(
    x = runif(n_cells, 0, 100),
    y = runif(n_cells, 0, 100)
  )
  
  # Cluster casuali
  intensity_cluster <- factor(sample(1:k_cell_types, n_cells, replace = TRUE))
  levels(intensity_cluster) <- paste0("Type_", 1:k_cell_types)
  
  # Matrice di espressione con dropout
  expression <- matrix(rnbinom(n_cells * n_genes, mu = 5, size = 1), nrow = n_genes, ncol = n_cells)
  # Introduzione di dropout in base all'espressione media
  gene_means <- rowMeans(expression)
  for (i in 1:n_genes) {
    dropout_rate <- exp(-0.5 * gene_means[i])
    dropout_mask <- runif(n_cells) < dropout_rate
    expression[i, dropout_mask] <- 0
  }
  rownames(expression) <- paste0("gene_", 1:n_genes)
  
  # Risultati simulati
  sim_results <- list(
    coordinates = coordinates,
    intensity_cluster = intensity_cluster,
    expression = expression,
    parameters = list(
      marker_params = list(
        marker_genes_per_type = list(
          Type_1 = c("gene_1", "gene_2"),
          Type_2 = c("gene_3", "gene_4"),
          Type_3 = c("gene_5", "gene_6")
        )
      )
    )
  )
  
  # Test base
  p1 <- plot_dropout_vs_expression(sim_results)
  expect_s3_class(p1, "ggplot")
  
  # Test senza log delle medie
  p2 <- plot_dropout_vs_expression(sim_results, log_mean = FALSE)
  expect_s3_class(p2, "ggplot")
  
  # Test senza evidenziazione dei marker
  p3 <- plot_dropout_vs_expression(sim_results, highlight_markers = FALSE)
  expect_s3_class(p3, "ggplot")
  
  # Test con titolo personalizzato
  p4 <- plot_dropout_vs_expression(sim_results, title = "Test Plot")
  expect_s3_class(p4, "ggplot")
  expect_equal(p4$labels$title, "Test Plot")
})

test_that("identify_marker_genes funziona correttamente", {
  # Crea un dataset di simulazione con espressione differenziale chiara
  n_cells <- 99  # Usando 99 cellule per avere 33 per gruppo
  n_genes <- 50
  k_cell_types <- 3
  
  # Coordinate spaziali casuali
  set.seed(123)
  coordinates <- data.frame(
    x = runif(n_cells, 0, 100),
    y = runif(n_cells, 0, 100)
  )
  
  # Cluster
  cell_per_type <- n_cells / k_cell_types
  intensity_cluster <- factor(rep(1:k_cell_types, each = cell_per_type))
  levels(intensity_cluster) <- paste0("Type_", 1:k_cell_types)
  
  # Matrice di espressione con marker chiari
  expression <- matrix(rpois(n_cells * n_genes, lambda = 2), nrow = n_genes, ncol = n_cells)
  # Imposta alcuni geni come marker per ciascun tipo
  for (type in 1:k_cell_types) {
    cells_idx <- which(intensity_cluster == paste0("Type_", type))
    marker_genes <- ((type-1)*5 + 1):(type*5)  # 5 marker per tipo
    # Aumenta l'espressione dei marker nel tipo specifico
    expression[marker_genes, cells_idx] <- expression[marker_genes, cells_idx] * 10
  }
  rownames(expression) <- paste0("gene_", 1:n_genes)
  
  # Risultati simulati
  sim_results <- list(
    coordinates = coordinates,
    intensity_cluster = intensity_cluster,
    expression = expression
  )
  
  # Test identificazione marker
  markers <- identify_marker_genes(sim_results, n_markers = 5)
  
  # Verifica che ci sia una lista per ogni tipo cellulare
  expect_equal(length(markers), k_cell_types)
  
  # Verifica che ogni tipo abbia il numero corretto di marker
  for (type in 1:k_cell_types) {
    expect_equal(length(markers[[paste0("Type_", type)]]), 5)
  }
  
  # Verifica che almeno alcuni dei marker identificati siano quelli attesi
  for (type in 1:k_cell_types) {
    expected_markers <- paste0("gene_", ((type-1)*5 + 1):(type*5))
    type_markers <- markers[[paste0("Type_", type)]]
    # Dovrebbe esserci una sovrapposizione significativa tra i marker attesi e quelli trovati
    overlap <- sum(type_markers %in% expected_markers)
    expect_gte(overlap, 2)  # almeno 2 dei 5 marker dovrebbero essere identificati
  }
})

test_that("find_highly_variable_genes funziona correttamente", {
  # Crea una matrice di espressione con variabilità controllata
  n_cells <- 100
  n_genes <- 200
  
  # Matrice di base
  set.seed(123)
  expression <- matrix(rpois(n_cells * n_genes, lambda = 5), nrow = n_genes, ncol = n_cells)
  
  # Rendi alcuni geni molto più variabili
  hvg_idx <- 1:20
  for (i in hvg_idx) {
    # Aumenta la variabilità di questi geni
    expression[i, ] <- rnbinom(n_cells, mu = 5, size = 0.5)
  }
  
  rownames(expression) <- paste0("gene_", 1:n_genes)
  
  # Test identificazione HVG
  hvg <- find_highly_variable_genes(expression, n_hvg = 20)
  
  # Verifica che il numero di geni restituiti sia corretto
  expect_equal(length(hvg), 20)
  
  # Verifica che almeno alcuni dei geni altamente variabili siano stati identificati
  hvg_names <- paste0("gene_", hvg_idx)
  overlap <- sum(hvg %in% hvg_names)
  expect_gte(overlap, 10)  # almeno la metà dei geni altamente variabili dovrebbe essere identificata
})

test_that("generate_validation_plots funziona correttamente", {
  # Modifica la funzione per testare senza dipendenze UMAP/t-SNE
  # Test semplice per verificare che la funzione non generi errori
  # Crea un dataset di simulazione minimo
  n_cells <- 50
  n_genes <- 30
  k_cell_types <- 2
  
  # Coordinate spaziali casuali
  set.seed(123)
  coordinates <- data.frame(
    x = runif(n_cells, 0, 100),
    y = runif(n_cells, 0, 100)
  )
  
  # Cluster casuali
  intensity_cluster <- factor(sample(1:k_cell_types, n_cells, replace = TRUE))
  levels(intensity_cluster) <- paste0("Type_", 1:k_cell_types)
  
  # Matrice di espressione
  expression <- matrix(rpois(n_cells * n_genes, lambda = 3), nrow = n_genes, ncol = n_cells)
  rownames(expression) <- paste0("gene_", 1:n_genes)
  
  # Risultati simulati
  sim_results <- list(
    coordinates = coordinates,
    intensity_cluster = intensity_cluster,
    expression = expression,
    parameters = list(
      marker_params = list(
        marker_genes_per_type = list(
          Type_1 = c("gene_1", "gene_2"),
          Type_2 = c("gene_3", "gene_4")
        )
      )
    )
  )
  
  # Crea una directory temporanea per i test
  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_plots")
  
  # Test della funzione wrapper con skip_dim_reduction = TRUE per evitare la dipendenza da umap/Rtsne
  plots <- tryCatch({
    generate_validation_plots(
      sim_results,
      marker_genes = c("gene_1", "gene_3"),
      output_dir = test_dir,
      file_prefix = "test",
      file_format = "png",
      width = 6,
      height = 4,
      skip_dim_reduction = TRUE
    )
  }, error = function(e) {
    fail(paste("generate_validation_plots ha generato un errore:", e$message))
    return(NULL)
  })
  
  # Verifica che la funzione abbia restituito una lista di plot
  expect_type(plots, "list")
  expect_true("mean_variance" %in% names(plots))
  expect_true("dropout" %in% names(plots))
  
  # Verifica che i file siano stati creati
  expect_true(file.exists(file.path(test_dir, "test_mean_variance.png")))
  expect_true(file.exists(file.path(test_dir, "test_dropout.png")))
  expect_true(file.exists(file.path(test_dir, "test_marker_gene_1.png")))
  expect_true(file.exists(file.path(test_dir, "test_marker_gene_3.png")))
  
  # Pulizia
  unlink(test_dir, recursive = TRUE)
})