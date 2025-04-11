#' Genera moduli di geni co-espressi
#'
#' Crea gruppi di geni co-espressi con correlazione controllata.
#'
#' @param n_genes Numero totale di geni
#' @param n_cells Numero di celle
#' @param cell_specific_params Parametri per i moduli genici
#' @param random_seed Seed per riproducibilità
#' @return Lista con i moduli genici e il rumore correlato
#' @export
generate_gene_modules <- function(
  n_genes,
  n_cells,
  cell_specific_params = list(
    use_gene_modules = TRUE,
    n_gene_modules = 5,
    module_correlation = 0.7
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Inizializza strutture di output
  gene_modules <- NULL
  module_noise <- NULL
  
  # Crea moduli di geni se richiesto
  if (cell_specific_params$use_gene_modules) {
    # Calcola il numero di geni per modulo
    genes_per_module <- ceiling(n_genes / cell_specific_params$n_gene_modules)
    
    # Assegna geni ai moduli
    gene_modules <- list()
    for (m in 1:cell_specific_params$n_gene_modules) {
      start_idx <- (m-1) * genes_per_module + 1
      end_idx <- min(m * genes_per_module, n_genes)
      gene_modules[[m]] <- start_idx:end_idx
    }
    
    # Crea rumore correlato per ogni modulo
    module_noise <- matrix(0, nrow = n_cells, ncol = n_genes)
    
    # Genera rumore base per ogni modulo
    base_module_noise <- matrix(
      rnorm(cell_specific_params$n_gene_modules * n_cells),
      nrow = n_cells, ncol = cell_specific_params$n_gene_modules
    )
    
    # Applica il rumore del modulo a ciascun gene appartenente al modulo
    for (m in 1:length(gene_modules)) {
      module_genes <- gene_modules[[m]]
      # Assegna lo stesso rumore base a tutti i geni del modulo, scalato per la correlazione
      module_noise[, module_genes] <- base_module_noise[, m] * cell_specific_params$module_correlation
    }
  }
  
  return(list(
    gene_modules = gene_modules,
    module_noise = module_noise
  ))
}