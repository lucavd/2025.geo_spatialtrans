# Test per moduli espressione (06*)

test_that("initialize_expression_params valida correttamente i parametri", {
  # Test con parametri standard completi per evitare avvisi
  params <- initialize_expression_params(
    marker_params = list(
      marker_genes_per_type = 10,
      marker_expression_fold = 1.5,
      marker_overlap_fold = 0.2
    ),
    dropout_params = list(
      dropout_range = c(0.2, 0.5),
      dispersion_range = c(2.0, 1.0)
    ),
    random_seed = 123
  )
  
  # Verifica la struttura del risultato
  expect_type(params, "list")
  expect_true(all(c("marker_params", "spatial_params", "dropout_params", 
                    "library_size_params", "cell_specific_params", 
                    "hybrid_params", "random_seed") %in% names(params)))
  
  # Verifica che i parametri abbiano i valori corretti
  expect_equal(params$marker_params$marker_genes_per_type, 10)
  expect_equal(params$dropout_params$dropout_range, c(0.2, 0.5))
  expect_equal(params$random_seed, 123)
  
  # Test con parametri non validi - in questo test VOGLIAMO le avvertenze
  # Per segnalare a testthat che queste avvertenze sono attese, usiamo suppressWarnings
  suppressWarnings(
    expect_warning(
      bad_params <- initialize_expression_params(
        marker_params = list(marker_genes_per_type = -5),
        random_seed = 123
      )
    )
  )
  
  # Verifica che i parametri siano stati corretti
  expect_equal(bad_params$marker_params$marker_genes_per_type, 10)
})

test_that("generate_baseline_expression crea profili di base corretti", {
  # Test di base
  mean_expr <- generate_baseline_expression(
    n_genes = 20,
    k_cell_types = 3,
    marker_params = list(
      marker_genes_per_type = 5,
      marker_expression_fold = 1.5,
      marker_overlap_fold = 0.2
    ),
    random_seed = 123
  )
  
  # Verifica la struttura del risultato
  expect_type(mean_expr, "list")
  expect_equal(length(mean_expr), 3)
  expect_equal(length(mean_expr[[1]]), 20)
  
  # Verifica che i marker siano potenziati
  # Primo tipo cellulare: primi 5 geni hanno espressione elevata
  first_type_markers <- mean_expr[[1]][1:5]
  expect_true(all(first_type_markers > 3))  # 2 (baseline) + 1.5 (fold)
  
  # Secondo tipo: geni 6-10 hanno espressione elevata
  second_type_markers <- mean_expr[[2]][6:10]
  expect_true(all(second_type_markers > 3))
  
  # Verifica overlapping nel secondo tipo cellulare
  # I primi 5 geni (marker del primo tipo) hanno espressione parziale
  overlap_expr <- mean_expr[[2]][1:5]
  expect_true(all(overlap_expr > 2 & overlap_expr < 3))  # 2 + 0.2
})

test_that("calculate_spatial_distances calcola correttamente le distanze", {
  # Crea un dataset semplice di test
  test_df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5),
    intensity_cluster = factor(c(rep(1, 12), rep(2, 13)))
  )
  
  # Calcola le distanze
  dist_result <- calculate_spatial_distances(test_df, random_seed = 123)
  
  # Verifica la struttura del risultato
  expect_type(dist_result, "list")
  expect_true(all(c("dist_mat", "mean_dist", "local_density") %in% names(dist_result)))
  
  # Verifica dimensioni
  expect_equal(dim(dist_result$dist_mat), c(25, 25))
  expect_equal(length(dist_result$mean_dist), 25)
  expect_equal(length(dist_result$local_density), 25)
  
  # Le distanze diagonali dovrebbero essere zero
  expect_equal(as.numeric(diag(dist_result$dist_mat)), rep(0, 25))
  
  # Le densità locali dovrebbero essere tra 0 e 1
  expect_true(all(dist_result$local_density >= 0 & dist_result$local_density <= 1))
})

test_that("generate_library_sizes produce dimensioni libreria realistiche", {
  # Crea un dataset semplice di test
  test_df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5),
    intensity_cluster = factor(c(rep(1, 12), rep(2, 13)))
  )
  
  # Genera dimensioni libreria
  lib_size <- generate_library_sizes(
    test_df,
    library_size_params = list(
      mean_library_size = 1000,
      library_size_cv = 0.2,
      spatial_effect_on_library = 0.5,
      cell_type_effect = TRUE
    ),
    spatial_params = list(spatial_range = 10),
    k_cell_types = 2,
    random_seed = 123
  )
  
  # Verifica dimensioni output
  expect_equal(length(lib_size), 25)
  
  # Verifica che i valori siano tutti positivi
  expect_true(all(lib_size > 0))
  
  # Verifica che la media sia approssimativamente quella desiderata
  expect_true(abs(mean(lib_size) - 1000) < 300)  # Permettiamo una certa variabilità
})

test_that("calculate_dropout_probabilities genera dropout corretto", {
  # Crea un dataset semplice di test con boundary_dist
  test_df <- data.frame(
    x = rep(1:5, each = 5),
    y = rep(1:5, times = 5),
    intensity_cluster = factor(c(rep(1, 12), rep(2, 13))),
    boundary_dist = runif(25)
  )
  
  # Media distanza (simulata)
  mean_dist <- runif(25)
  
  # Test con gradient_regions = FALSE
  dropout_prob_nogradient <- calculate_dropout_probabilities(
    test_df,
    mean_dist,
    dropout_params = list(dropout_range = c(0.1, 0.5)),
    spatial_params = list(gradient_regions = FALSE)
  )
  
  # Verifica risultati
  expect_equal(length(dropout_prob_nogradient), 25)
  expect_true(all(dropout_prob_nogradient >= 0.1 & dropout_prob_nogradient <= 0.5))
  
  # Test con gradient_regions = TRUE
  dropout_prob_gradient <- calculate_dropout_probabilities(
    test_df,
    mean_dist,
    dropout_params = list(dropout_range = c(0.1, 0.5)),
    spatial_params = list(gradient_regions = TRUE)
  )
  
  # Verifica risultati
  expect_equal(length(dropout_prob_gradient), 25)
  expect_true(all(dropout_prob_gradient >= 0.1 & dropout_prob_gradient <= 0.5))
  
  # I due metodi dovrebbero produrre risultati diversi
  expect_false(identical(dropout_prob_nogradient, dropout_prob_gradient))
})

test_that("generate_gene_modules crea moduli di geni correlati", {
  # Parametri di test
  n_genes <- 20
  n_cells <- 30
  
  # Genera moduli genici - versione base
  modules_result <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = TRUE,
      n_gene_modules = 4,
      module_correlation = 0.7
    ),
    random_seed = 123
  )
  
  # Verifica la struttura del risultato
  expect_type(modules_result, "list")
  expect_true(all(c("gene_modules", "module_noise") %in% names(modules_result)))
  
  # Verifica i moduli genici
  expect_equal(length(modules_result$gene_modules), 4)
  
  # Verifica il rumore
  expect_equal(dim(modules_result$module_noise), c(30, 20))
  
  # Test con moduli disabilitati
  modules_disabled <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = FALSE,
      n_gene_modules = 4,
      module_correlation = 0.7
    ),
    random_seed = 123
  )
  
  expect_null(modules_disabled$gene_modules)
  expect_null(modules_disabled$module_noise)
})

test_that("generate_gene_modules crea moduli avanzati con fattori latenti", {
  # Parametri di test
  n_genes <- 50
  n_cells <- 30
  
  # Genera moduli genici con funzionalità avanzate e tutti i parametri esplicitamente definiti
  modules_result <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = TRUE,
      n_gene_modules = 4,
      module_correlation = 0.7,
      module_hierarchical = TRUE,
      module_overlap = 0.1,
      module_size_distribution = "exponential",
      n_latent_factors = 3,
      module_network_density = 0.2,
      latent_factor_strength = 0.8
    ),
    random_seed = 123
  )
  
  # Verifica la struttura del risultato estesa
  expect_type(modules_result, "list")
  expect_true(all(c("gene_modules", "module_noise", "latent_factors", "module_network") %in% 
                  names(modules_result)))
  
  # Verifica i fattori latenti
  expect_true(!is.null(modules_result$latent_factors))
  expect_equal(dim(modules_result$latent_factors), c(30, 3))
  
  # Verifica la rete di moduli
  expect_true(!is.null(modules_result$module_network))
  expect_true(is.matrix(modules_result$module_network))
  
  # Verifica che con module_hierarchical=TRUE ci siano più moduli del numero iniziale
  expect_true(length(modules_result$gene_modules) >= 4)
  
  # Verifica distribuzione esponenziale - dovremmo avere moduli di dimensioni diverse
  module_sizes <- sapply(modules_result$gene_modules, length)
  expect_true(length(unique(module_sizes)) > 1)
})

test_that("generate_gene_modules gestisce correttamente la sovrapposizione tra moduli", {
  n_genes <- 50
  n_cells <- 30
  
  # Test con sovrapposizione dei moduli alta
  modules_overlap_high <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = TRUE,
      n_gene_modules = 4,
      module_correlation = 0.7,
      module_hierarchical = FALSE,
      module_overlap = 0.3,  # Sovrapposizione elevata
      n_latent_factors = 2,
      module_network_density = 0.2,
      latent_factor_strength = 0.8,
      module_size_distribution = "uniform"
    ),
    random_seed = 123
  )
  
  # Test con sovrapposizione dei moduli bassa
  modules_overlap_low <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = TRUE,
      n_gene_modules = 4,
      module_correlation = 0.7,
      module_hierarchical = FALSE,
      module_overlap = 0.0,  # Nessuna sovrapposizione
      n_latent_factors = 2,
      module_network_density = 0.2,
      latent_factor_strength = 0.8,
      module_size_distribution = "uniform"
    ),
    random_seed = 123
  )
  
  # Calcola i totali dei geni nei moduli
  overlap_high_total <- sum(sapply(modules_overlap_high$gene_modules, length))
  overlap_low_total <- sum(sapply(modules_overlap_low$gene_modules, length))
  
  # Conta i geni unici (non duplicati) nei moduli
  overlap_high_unique <- length(unique(unlist(modules_overlap_high$gene_modules)))
  overlap_low_unique <- length(unique(unlist(modules_overlap_low$gene_modules)))
  
  # Con overlap elevato, il rapporto tra geni totali e unici dovrebbe essere > 1
  expect_gt(overlap_high_total / overlap_high_unique, 1.0)
  
  # Con overlap basso o nullo, il rapporto dovrebbe essere vicino a 1
  expect_lte(overlap_low_total / overlap_low_unique, 1.1)
  
  # Dovrebbe esserci più sovrapposizione nel set con overlap alto
  expect_gt(overlap_high_total / overlap_high_unique, 
            overlap_low_total / overlap_low_unique)
})