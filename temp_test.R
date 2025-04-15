# Carica manualmente i file in ordine
files <- list.files("R/functions", full.names = TRUE, pattern = ".*\\.R$")
sorted_files <- sort(files)
for (file in sorted_files) { 
  cat("Caricamento file:", file, "\n")
  source(file) 
}

cat("\nVerifica interfaccia generate_gene_modules\n")
# Esegui un test semplice
n_genes <- 20
n_cells <- 30
  
# Genera moduli genici - versione base
res <- generate_gene_modules(
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

cat("Risultato generate_gene_modules:\n")
str(res)
