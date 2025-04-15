# Carica manualmente i file in ordine
files <- list.files("R/functions", full.names = TRUE, pattern = ".*\\.R$")
sorted_files <- sort(files)
for (file in sorted_files) { 
  source(file) 
}

cat("\nTest con dimensioni e parametri del modulo variabili\n")
# Test con dimensioni variabili
test_sizes <- function() {
  n_genes <- 50
  n_cells <- 30
  
  # Con distribuzione uniforme
  uniform_modules <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = TRUE,
      n_gene_modules = 5,
      module_correlation = 0.7,
      module_size_distribution = "uniform"
    ),
    random_seed = 123
  )
  
  # Con distribuzione esponenziale
  exp_modules <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = TRUE,
      n_gene_modules = 5,
      module_correlation = 0.7,
      module_size_distribution = "exponential"
    ),
    random_seed = 123
  )
  
  # Verifica la differenza nelle dimensioni
  uniform_sizes <- sapply(uniform_modules$gene_modules, length)
  exp_sizes <- sapply(exp_modules$gene_modules, length)
  
  cat("Dimensione moduli con distribuzione uniforme:", uniform_sizes, "\n")
  cat("Deviazione standard dimensioni uniformi:", sd(uniform_sizes), "\n")
  cat("Dimensione moduli con distribuzione esponenziale:", exp_sizes, "\n")
  cat("Deviazione standard dimensioni esponenziali:", sd(exp_sizes), "\n")
}

test_sizes()

cat("\nTest con fattori latenti e network di moduli\n")
test_latent_factors <- function() {
  n_genes <- 50
  n_cells <- 30
  
  # Con fattori latenti e networking attivo
  advanced_modules <- generate_gene_modules(
    n_genes,
    n_cells,
    cell_specific_params = list(
      use_gene_modules = TRUE,
      n_gene_modules = 5,
      module_correlation = 0.7,
      module_hierarchical = TRUE,
      n_latent_factors = 3,
      module_network_density = 0.3
    ),
    random_seed = 123
  )
  
  # Verifica numero di moduli dopo la gerarchizzazione
  cat("Numero di moduli dopo gerarchizzazione:", length(advanced_modules$gene_modules), "\n")
  cat("Dovrebbe essere >= 5\n")
  
  # Verifica la rete di moduli
  network_density <- sum(advanced_modules$module_network > 0) / 
                     (length(advanced_modules$gene_modules)^2 - length(advanced_modules$gene_modules))
  cat("Densità effettiva della rete:", network_density, "\n")
  
  # Verifica i fattori latenti
  cat("Dimensioni matrice fattori latenti:", dim(advanced_modules$latent_factors), "\n")
}

test_latent_factors()
