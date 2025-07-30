#' Genera moduli di geni co-espressi
#'
#' Crea gruppi di geni co-espressi con correlazione controllata e struttura di rete.
#' Implementa un modello avanzato di co-espressione che include:
#' - Moduli gerarchici con sotto-moduli
#' - Interazioni tra moduli diversi (cross-talk)
#' - Distribuzione realistica delle dimensioni dei moduli
#' - Possibilità di assegnare geni a più moduli
#' - Fattori latenti che controllano l'attivazione dei moduli
#'
#' @param n_genes Numero totale di geni
#' @param n_cells Numero di celle
#' @param cell_specific_params Parametri per i moduli genici
#' @param random_seed Seed per riproducibilità
#' @return Lista con i moduli genici, il rumore correlato e i fattori latenti
#' @importFrom stats rnorm rbinom rpois runif
#' @export
generate_gene_modules <- function(
  n_genes,
  n_cells,
  cell_specific_params = list(
    use_gene_modules = TRUE,
    n_gene_modules = 5,
    module_correlation = 0.7,
    module_hierarchical = FALSE,
    module_overlap = 0.1,
    module_size_distribution = "exponential",
    n_latent_factors = 3,
    module_network_density = 0.2,
    latent_factor_strength = 0.8
  ),
  random_seed = 123
) {
  # Imposta il seed per riproducibilità
  set.seed(random_seed)
  
  # Inizializza strutture di output
  gene_modules <- NULL
  module_noise <- NULL
  latent_factors <- NULL
  module_network <- NULL
  
  # Crea moduli di geni se richiesto
  use_gene_modules <- ifelse(is.null(cell_specific_params$use_gene_modules), TRUE, cell_specific_params$use_gene_modules)
  if (use_gene_modules) {
    # Determinazione del numero di moduli
    n_modules <- ifelse(is.null(cell_specific_params$n_gene_modules), 5, cell_specific_params$n_gene_modules)
    
    # Creazione dei moduli con distribuzione di dimensioni realistica
    gene_modules <- list()
    
    module_size_distribution <- ifelse(is.null(cell_specific_params$module_size_distribution), "exponential", cell_specific_params$module_size_distribution)
    if (module_size_distribution == "exponential") {
      # Distribuzione esponenziale: molti moduli piccoli, pochi moduli grandi
      # come nei pathway biologici reali
      module_sizes <- rpois(n_modules, lambda = n_genes / (n_modules * 2))
      module_sizes <- pmax(module_sizes, 5)  # Dimensione minima di un modulo
      module_sizes <- round(module_sizes / sum(module_sizes) * n_genes * 0.8)  # Lascia il 20% dei geni non assegnati
    } else if (module_size_distribution == "uniform") {
      # Distribuzione uniforme (come nell'implementazione originale)
      genes_per_module <- ceiling(n_genes / n_modules)
      module_sizes <- rep(genes_per_module, n_modules)
      module_sizes[n_modules] <- n_genes - sum(module_sizes[1:(n_modules-1)])
    } else {
      # Default a poisson con media personalizzabile
      module_sizes <- rpois(n_modules, lambda = n_genes / n_modules)
      module_sizes <- pmax(module_sizes, 3)  # Dimensione minima di un modulo
      # Riscala per assicurarsi che non superino n_genes
      module_sizes <- round(module_sizes / sum(module_sizes) * n_genes * 0.8)
    }
    
    # Assegna geni ai moduli (versione base)
    all_genes <- 1:n_genes
    available_genes <- all_genes
    
    for (m in 1:n_modules) {
      # Se non ci sono più geni disponibili, esci dal ciclo
      if (length(available_genes) == 0) break
      
      # Limita la dimensione del modulo ai geni disponibili
      size <- min(module_sizes[m], length(available_genes))
      
      # Seleziona i geni per questo modulo
      if (size > 0) {
        module_genes <- sample(available_genes, size)
        gene_modules[[m]] <- module_genes
        
        # Rimuovi i geni assegnati (se non desideriamo sovrapposizioni)
        module_overlap <- ifelse(is.null(cell_specific_params$module_overlap), 0.1, cell_specific_params$module_overlap)
        
        if (module_overlap < 0.01) {
          available_genes <- setdiff(available_genes, module_genes)
        } else {
          # Con sovrapposizione, rimuoviamo solo una parte dei geni
          n_remove <- round(size * (1 - module_overlap))
          if (n_remove > 0 && n_remove <= length(module_genes)) {
            genes_to_remove <- sample(module_genes, n_remove)
            available_genes <- setdiff(available_genes, genes_to_remove)
          }
        }
      }
    }
    
    # Se abilitata la gerarchia, crea moduli gerarchici
    module_hierarchical <- ifelse(is.null(cell_specific_params$module_hierarchical), FALSE, cell_specific_params$module_hierarchical)
    if (module_hierarchical && n_modules >= 3) {
      # Seleziona alcuni moduli per creare sotto-moduli
      n_parent_modules <- max(1, round(n_modules / 3))
      parent_modules <- sample(1:n_modules, n_parent_modules)
      
      for (parent in parent_modules) {
        # Crea 2-3 sotto-moduli dal modulo genitore
        n_children <- sample(2:3, 1)
        parent_genes <- gene_modules[[parent]]
        
        if (length(parent_genes) >= n_children * 3) {  # Solo se ci sono abbastanza geni
          # Dividi i geni del genitore nei sotto-moduli
          child_genes <- split(parent_genes, cut(seq_along(parent_genes), n_children))
          
          # Aggiungi i sotto-moduli alla lista
          for (i in 1:length(child_genes)) {
            gene_modules[[length(gene_modules) + 1]] <- unlist(child_genes[i])
          }
        }
      }
    }
    
    # Creazione della rete di interazione tra moduli
    n_total_modules <- length(gene_modules)
    module_network <- matrix(0, nrow = n_total_modules, ncol = n_total_modules)
    
    # Generiamo connessioni tra moduli con densità specificata
    module_network_density <- ifelse(is.null(cell_specific_params$module_network_density), 0.2, cell_specific_params$module_network_density)
    if (module_network_density > 0) {
      for (i in 1:(n_total_modules-1)) {
        for (j in (i+1):n_total_modules) {
          if (runif(1) < module_network_density) {
            # Assegna un peso di interazione casuale tra 0.1 e 0.5
            interaction_weight <- runif(1, 0.1, 0.5)
            module_network[i, j] <- interaction_weight
            module_network[j, i] <- interaction_weight  # Matrice simmetrica
          }
        }
      }
    }
    
    # Genera fattori latenti che influenzano i moduli
    n_latent <- ifelse(is.null(cell_specific_params$n_latent_factors), 3, cell_specific_params$n_latent_factors)
    latent_factor_strength <- ifelse(is.null(cell_specific_params$latent_factor_strength), 0.8, cell_specific_params$latent_factor_strength)
    if (n_latent > 0) {
      # Matrice di coefficienti che mappa fattori latenti ai moduli
      latent_to_module <- matrix(0, nrow = n_total_modules, ncol = n_latent)
      
      # Ogni modulo è influenzato da 1-2 fattori latenti
      for (m in 1:n_total_modules) {
        n_factors <- sample(1:min(2, n_latent), 1)
        factor_indices <- sample(1:n_latent, n_factors)
        latent_to_module[m, factor_indices] <- runif(n_factors, 0.5, 1.0)
      }
      
      # Genera valori dei fattori latenti per ogni cellula
      latent_factors <- matrix(rnorm(n_cells * n_latent), nrow = n_cells, ncol = n_latent)
      
      # Inizializza matrice di rumore correlato
      module_noise <- matrix(0, nrow = n_cells, ncol = n_genes)
      
      # Calcola il contributo dei fattori latenti per ogni modulo
      module_activities <- latent_factors %*% t(latent_to_module)
      
      # Applica il network tra moduli per propagare gli effetti
      if (sum(module_network) > 0) {
        # Normalizza la matrice di rete
        network_norm <- sweep(module_network, 1, rowSums(module_network) + 1e-10, "/")
        
        # Propaga l'attivazione attraverso la rete (2 passi)
        orig_activities <- module_activities
        for (step in 1:2) {
          module_activities <- 0.7 * orig_activities + 0.3 * (module_activities %*% network_norm)
        }
      }
      
      # Applica l'attività dei moduli ai geni
      for (m in 1:n_total_modules) {
        module_genes <- gene_modules[[m]]
        if (length(module_genes) > 0) {
          # L'intensità dell'effetto può variare per diversi geni nello stesso modulo
          module_correlation <- ifelse(is.null(cell_specific_params$module_correlation), 0.7, cell_specific_params$module_correlation)
          gene_weights <- runif(length(module_genes), 
                              module_correlation * 0.5,
                              module_correlation * 1.5)
          
          # Applica l'attività del modulo a ciascun gene con intensità variabile
          for (i in 1:length(module_genes)) {
            g <- module_genes[i]
            # Aggiungi l'attività del modulo scalata per il peso del gene
            module_noise[, g] <- module_noise[, g] + 
                               module_activities[, m] * gene_weights[i] * 
                               latent_factor_strength
          }
        }
      }
      
      # Aggiungi rumore individuale specifico del gene
      gene_specific_noise <- matrix(
        rnorm(n_cells * n_genes, 0, 0.2),
        nrow = n_cells, ncol = n_genes
      )
      
      # Combina il rumore correlato con il rumore specifico del gene
      module_noise <- module_noise + gene_specific_noise
    } else {
      # Versione originale se non si usano fattori latenti
      module_noise <- matrix(0, nrow = n_cells, ncol = n_genes)
      base_module_noise <- matrix(
        rnorm(n_total_modules * n_cells),
        nrow = n_cells, ncol = n_total_modules
      )
      
      for (m in 1:n_total_modules) {
        module_genes <- gene_modules[[m]]
        if (length(module_genes) > 0) {
          module_noise[, module_genes] <- module_noise[, module_genes] + 
                                        base_module_noise[, m] * cell_specific_params$module_correlation
        }
      }
    }
  }
  
  return(list(
    gene_modules = gene_modules,
    module_noise = module_noise,
    latent_factors = latent_factors,
    module_network = module_network
  ))
}