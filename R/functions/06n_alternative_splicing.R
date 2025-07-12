#' Adattamento spaziale dello splicing alternativo
#'
#' Giustificazione: Le analisi di Ding et al. (Genome Biology 2020) e Stewart et al. 
#' (Nature Communications 2021) dimostrano pattern spaziali di splicing alternativo 
#' nei tessuti, con implicazioni funzionali rilevanti.
#'
#' Questo modulo implementa:
#' - Modelli varianti di splicing per lo stesso gene con pattern spaziali distinti
#' - Incorpori regolazione coordinata dello splicing in specifici domini tissutali
#' - Implementi correlazioni tra splicing alternativo e microambienti locali

#' Genera varianti di splicing per i geni
#'
#' Crea varianti di splicing alternativo per un sottoinsieme di geni,
#' definendo i loro pattern spaziali basati sulla posizione nel tessuto.
#'
#' @param cell_df Dataframe delle cellule con coordinate e cluster
#' @param gene_info Informazioni sui geni (opzionale)
#' @param n_genes Numero totale di geni
#' @param splicing_params Parametri per lo splicing alternativo
#' @param random_seed Seed per riproducibilità
#' @return Lista con definizione delle varianti di splicing
#' @importFrom stats runif rbinom rnorm
#' @importFrom sp coordinates
#' @importFrom gstat gstat predict vgm
#' @export
generate_splicing_variants <- function(
  cell_df,
  gene_info = NULL,
  n_genes,
  splicing_params = list(
    fraction_genes_with_variants = 0.15,  # Frazione di geni con splicing alternativo
    n_variants_per_gene = c(2, 3),        # Range di varianti per gene (min, max)
    variant_effect_strength = c(0.3, 1.5),  # Range di intensità dell'effetto (min, max)
    spatial_regulation = 0.7,             # Intensità della regolazione spaziale (0-1)
    cell_type_regulation = 0.6,           # Intensità della regolazione specifica per tipo cellulare
    coordinated_splicing_groups = 3,      # Numero di gruppi di splicing coordinato
    n_splicing_regulators = 6,            # Numero di "regolatori" di splicing
    splicing_noise = 0.2                  # Rumore nel pattern di splicing
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  fraction_genes_with_variants <- splicing_params$fraction_genes_with_variants
  n_variants_per_gene <- splicing_params$n_variants_per_gene
  variant_effect_strength <- splicing_params$variant_effect_strength
  spatial_regulation <- splicing_params$spatial_regulation
  cell_type_regulation <- splicing_params$cell_type_regulation
  coordinated_splicing_groups <- splicing_params$coordinated_splicing_groups
  n_splicing_regulators <- splicing_params$n_splicing_regulators
  splicing_noise <- splicing_params$splicing_noise
  
  # Valori predefiniti se non specificati
  if (is.null(fraction_genes_with_variants)) fraction_genes_with_variants <- 0.15
  if (is.null(n_variants_per_gene)) n_variants_per_gene <- c(2, 3)
  if (is.null(variant_effect_strength)) variant_effect_strength <- c(0.3, 1.5)
  if (is.null(spatial_regulation)) spatial_regulation <- 0.7
  if (is.null(cell_type_regulation)) cell_type_regulation <- 0.6
  if (is.null(coordinated_splicing_groups)) coordinated_splicing_groups <- 3
  if (is.null(n_splicing_regulators)) n_splicing_regulators <- 6
  if (is.null(splicing_noise)) splicing_noise <- 0.2
  
  # Dimensioni
  n_cells <- nrow(cell_df)
  
  # Seleziona geni con varianti di splicing
  n_spliced_genes <- round(n_genes * fraction_genes_with_variants)
  spliced_genes <- sample(1:n_genes, n_spliced_genes)
  
  # Genera campi dei regolatori di splicing (pattern spaziali che influenzano lo splicing)
  splicing_regulators <- matrix(0, nrow = n_cells, ncol = n_splicing_regulators)
  
  # Converti cell_df in oggetto spatial per generare campi spaziali
  sp_df <- cell_df
  sp::coordinates(sp_df) <- ~ x + y
  
  for (r in 1:n_splicing_regulators) {
    # Definisci un modello di correlazione spaziale per questo regolatore
    range_param <- mean(c(diff(range(cell_df$x)), diff(range(cell_df$y)))) * runif(1, 0.1, 0.3)
    
    # Crea un campo gaussiano correlato spazialmente come "regolatore"
    regulator_gp <- gstat::gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                          beta = 0, model = gstat::vgm(psill = 1.0, 
                                                 range = range_param, 
                                                 model = "Exp"), 
                          nmax = 30)
    
    # Genera e normalizza il campo
    reg_field <- gstat::predict(regulator_gp, newdata = sp_df, nsim = 1)$sim1
    splicing_regulators[, r] <- scale(reg_field)
  }
  
  # Organizza i geni in gruppi di splicing coordinato
  splicing_groups <- sample(1:coordinated_splicing_groups, n_spliced_genes, replace = TRUE)
  
  # Crea lista di informazioni su varianti di splicing
  splicing_info <- list()
  
  # Per ogni gene con splicing alternativo, genera informazioni sulle varianti
  for (i in 1:n_spliced_genes) {
    gene_id <- spliced_genes[i]
    group_id <- splicing_groups[i]
    
    # Determina il numero di varianti per questo gene
    n_variants <- sample(n_variants_per_gene[1]:n_variants_per_gene[2], 1)
    
    # Crea liste per dettagli delle varianti
    variant_ids <- paste0(gene_id, ".", 1:n_variants)
    variant_names <- paste0("Var", 1:n_variants)
    
    # Genera "pesi" di regolatori per questo gene (quanto ogni regolatore influenza ogni variante)
    # Questi determinano come i vari regolatori spaziali influenzano le diverse varianti
    regulator_weights <- matrix(runif(n_variants * n_splicing_regulators, -1, 1), 
                              nrow = n_variants, ncol = n_splicing_regulators)
    
    # Per geni nello stesso gruppo, simili pattern di regolazione (splicing coordinato)
    if (i > 1 && any(splicing_groups[1:(i-1)] == group_id)) {
      # Trova un gene già processato nello stesso gruppo
      template_idx <- sample(which(splicing_groups[1:(i-1)] == group_id), 1)
      template_gene <- spliced_genes[template_idx]
      
      # Copia i pesi di regolazione dal gene template, aggiungendo variazione minore
      template_weights <- splicing_info[[paste0("gene_", template_gene)]]$regulator_weights
      base_weights <- template_weights[sample(1:nrow(template_weights), n_variants, replace = TRUE), ]
      regulator_weights <- base_weights + matrix(rnorm(n_variants * n_splicing_regulators, 0, 0.3), 
                                              nrow = n_variants, ncol = n_splicing_regulators)
    }
    
    # Effetto del tipo cellulare sullo splicing (if cell type information is available)
    cell_type_effects <- NULL
    if ("intensity_cluster" %in% colnames(cell_df)) {
      cell_types <- sort(unique(cell_df$intensity_cluster))
      n_cell_types <- length(cell_types)
      
      # Genera pesi per ogni tipo cellulare e variante
      cell_type_effects <- matrix(runif(n_variants * n_cell_types, -1, 1), 
                                nrow = n_variants, ncol = n_cell_types)
      
      # Per geni nello stesso gruppo, pattern simili di regolazione per tipo cellulare
      if (i > 1 && any(splicing_groups[1:(i-1)] == group_id)) {
        template_idx <- sample(which(splicing_groups[1:(i-1)] == group_id), 1)
        template_gene <- spliced_genes[template_idx]
        template_effects <- splicing_info[[paste0("gene_", template_gene)]]$cell_type_effects
        
        if (!is.null(template_effects)) {
          # Usa pesi simili con piccola variazione
          base_effects <- template_effects[sample(1:nrow(template_effects), n_variants, replace = TRUE), ]
          cell_type_effects <- base_effects + matrix(rnorm(n_variants * n_cell_types, 0, 0.3), 
                                                  nrow = n_variants, ncol = n_cell_types)
        }
      }
    }
    
    # Genera informazioni sull'effetto delle varianti
    variant_effects <- runif(n_variants, 
                           variant_effect_strength[1], 
                           variant_effect_strength[2])
    
    # Bilanciamento: vettore di probabilità di default tra le varianti (somma a 1)
    base_probs <- rep(1/n_variants, n_variants)
    
    # Memorizza tutte le informazioni per questo gene
    splicing_info[[paste0("gene_", gene_id)]] <- list(
      gene_id = gene_id,
      group_id = group_id,
      n_variants = n_variants,
      variant_ids = variant_ids,
      variant_names = variant_names,
      variant_effects = variant_effects,
      regulator_weights = regulator_weights,
      cell_type_effects = cell_type_effects,
      base_probs = base_probs
    )
  }
  
  # Crea una tabella di metadati di splicing
  splicing_metadata <- data.frame(
    gene_id = integer(),
    has_variants = logical(),
    n_variants = integer(),
    splicing_group = integer(),
    stringsAsFactors = FALSE
  )
  
  for (g in 1:n_genes) {
    info_key <- paste0("gene_", g)
    if (info_key %in% names(splicing_info)) {
      info <- splicing_info[[info_key]]
      splicing_metadata <- rbind(splicing_metadata, data.frame(
        gene_id = g,
        has_variants = TRUE,
        n_variants = info$n_variants,
        splicing_group = info$group_id,
        stringsAsFactors = FALSE
      ))
    } else {
      splicing_metadata <- rbind(splicing_metadata, data.frame(
        gene_id = g,
        has_variants = FALSE,
        n_variants = 0,
        splicing_group = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(list(
    splicing_info = splicing_info,
    splicing_regulators = splicing_regulators,
    splicing_metadata = splicing_metadata,
    params = splicing_params
  ))
}

#' Calcola propensioni di splicing per ogni cellula
#'
#' Calcola la propensione di ciascuna cellula a generare ciascuna variante
#' di splicing, basata su regolatori spaziali e specifici del tipo cellulare.
#'
#' @param cell_df Dataframe delle cellule con coordinate e cluster
#' @param splicing_variants Risultato della funzione generate_splicing_variants
#' @param splicing_params Parametri per il calcolo delle propensioni
#' @param random_seed Seed per riproducibilità
#' @return Matrice di propensioni di splicing per ogni cellula/gene
#' @importFrom stats plogis
#' @export
calculate_splicing_propensity <- function(
  cell_df,
  splicing_variants,
  splicing_params = list(
    spatial_regulation = 0.7,             # Intensità della regolazione spaziale (0-1)
    cell_type_regulation = 0.6,           # Intensità della regolazione specifica per tipo cellulare
    splicing_noise = 0.2                  # Rumore nel pattern di splicing
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Estrai parametri
  spatial_regulation <- splicing_params$spatial_regulation
  cell_type_regulation <- splicing_params$cell_type_regulation
  splicing_noise <- splicing_params$splicing_noise
  
  # Usa i parametri dalla lista di varianti se non specificati
  if (is.null(spatial_regulation))
    spatial_regulation <- splicing_variants$params$spatial_regulation
  if (is.null(cell_type_regulation)) 
    cell_type_regulation <- splicing_variants$params$cell_type_regulation
  if (is.null(splicing_noise))
    splicing_noise <- splicing_variants$params$splicing_noise
  
  # Dimensioni
  n_cells <- nrow(cell_df)
  
  # Estrai elenco di geni con splicing alternativo
  splicing_info <- splicing_variants$splicing_info
  splicing_regulators <- splicing_variants$splicing_regulators
  
  # Crea una lista per memorizzare le propensioni di splicing per ogni gene
  splicing_propensity <- list()
  
  # Per ogni gene con splicing alternativo
  for (gene_key in names(splicing_info)) {
    gene_info <- splicing_info[[gene_key]]
    gene_id <- gene_info$gene_id
    n_variants <- gene_info$n_variants
    
    # Pesi di regolazione per questo gene
    regulator_weights <- gene_info$regulator_weights
    
    # Calcola l'effetto dei regolatori spaziali per ogni variante in ogni cellula
    # Dimensione: n_cells x n_variants
    spatial_effects <- splicing_regulators %*% t(regulator_weights)
    
    # Inizializza la matrice di propensione di splicing
    propensity <- matrix(0, nrow = n_cells, ncol = n_variants)
    
    # Aggiungi l'effetto spaziale
    propensity <- propensity + spatial_regulation * spatial_effects
    
    # Aggiungi l'effetto del tipo cellulare, se disponibile
    if (!is.null(gene_info$cell_type_effects) && "intensity_cluster" %in% colnames(cell_df)) {
      cell_types <- sort(unique(cell_df$intensity_cluster))
      cell_type_effects <- gene_info$cell_type_effects
      
      # Per ogni tipo cellulare
      for (ct in 1:length(cell_types)) {
        cell_type <- cell_types[ct]
        cells_of_type <- which(cell_df$intensity_cluster == cell_type)
        
        # Applica l'effetto specifico del tipo cellulare alle celle di quel tipo
        if (length(cells_of_type) > 0) {
          for (v in 1:n_variants) {
            propensity[cells_of_type, v] <- propensity[cells_of_type, v] + 
              cell_type_regulation * cell_type_effects[v, ct]
          }
        }
      }
    }
    
    # Aggiungi rumore casuale
    if (splicing_noise > 0) {
      noise <- matrix(rnorm(n_cells * n_variants, 0, splicing_noise), nrow = n_cells, ncol = n_variants)
      propensity <- propensity + noise
    }
    
    # Trasforma la propensione in probabilità usando una funzione softmax
    # Per garantire che le probabilità sommino a 1 per ogni cellula
    propensity_exp <- exp(propensity)
    probabilities <- propensity_exp / rowSums(propensity_exp)
    
    # Memorizza le probabilità per questo gene
    splicing_propensity[[gene_key]] <- probabilities
  }
  
  return(splicing_propensity)
}

#' Genera espressione di varianti di splicing
#'
#' Genera matrici di espressione per le varianti di splicing
#' basate sulla propensione di splicing e sull'espressione di base.
#'
#' @param expression_matrix Matrice di espressione originale
#' @param splicing_variants Risultato della funzione generate_splicing_variants
#' @param splicing_propensity Risultato della funzione calculate_splicing_propensity
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici di espressione modificate per varianti di splicing
#' @importFrom stats rmultinom
#' @export
generate_splicing_expression <- function(
  expression_matrix,
  splicing_variants,
  splicing_propensity,
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Dimensioni
  n_cells <- nrow(expression_matrix)
  
  # Estrai informazioni sulle varianti e propensioni
  splicing_info <- splicing_variants$splicing_info
  
  # Inizializza una matrice per l'espressione totale modificata (summing up all variants)
  modified_expression <- expression_matrix
  
  # Crea una lista per memorizzare l'espressione di ogni variante di splicing
  variant_expression <- list()
  
  # Per ogni gene con splicing alternativo
  for (gene_key in names(splicing_info)) {
    gene_info <- splicing_info[[gene_key]]
    gene_id <- gene_info$gene_id
    n_variants <- gene_info$n_variants
    variant_effects <- gene_info$variant_effects
    
    # Ignora geni che non esistono nella matrice di espressione
    if (gene_id > ncol(expression_matrix)) {
      next
    }
    
    # Propensione a ciascuna variante per ogni cellula
    gene_propensity <- splicing_propensity[[gene_key]]
    
    # Estrai espressione originale per questo gene
    gene_expression <- expression_matrix[, gene_id]
    
    # Crea matrice per l'espressione delle varianti di questo gene
    gene_variants_expr <- matrix(0, nrow = n_cells, ncol = n_variants)
    
    # Per ogni cellula
    for (cell in 1:n_cells) {
      expr_level <- gene_expression[cell]
      
      if (expr_level > 0) {
        # Proporzioni di splicing per questa cellula
        cell_proportions <- gene_propensity[cell, ]
        
        # Se abbiamo conteggi interi (UMI-based), usiamo il multinomiale per assegnare alle varianti
        if (expr_level >= 5 && floor(expr_level) == expr_level) {
          # Genera conteggi multinomiali per le varianti
          variant_counts <- rmultinom(n = 1, size = expr_level, prob = cell_proportions)
          gene_variants_expr[cell, ] <- variant_counts
        } else {
          # Per espressione non intera o molto bassa, ripartisci proporzionalmente
          gene_variants_expr[cell, ] <- expr_level * cell_proportions
        }
        
        # Applica gli effetti specifici delle varianti (aumento/riduzione dell'espressione)
        for (v in 1:n_variants) {
          gene_variants_expr[cell, v] <- gene_variants_expr[cell, v] * variant_effects[v]
        }
      }
    }
    
    # Memorizza l'espressione delle varianti
    for (v in 1:n_variants) {
      variant_id <- gene_info$variant_ids[v]
      variant_expression[[variant_id]] <- gene_variants_expr[, v]
    }
    
    # Modifica l'espressione totale con la somma delle varianti
    modified_expression[, gene_id] <- rowSums(gene_variants_expr)
  }
  
  # Converti la lista di espressione varianti in una matrice
  # Prepara un elenco di tutti i variant_id e loro gene_id
  all_variants <- data.frame(
    variant_id = character(),
    gene_id = integer(),
    variant_index = integer(),
    stringsAsFactors = FALSE
  )
  
  # Costruisci all_variants solo per varianti effettivamente presenti, senza rbind massivo
  all_variants_list <- list()
  idx <- 1
  for (gene_key in names(splicing_info)) {
    gene_info <- splicing_info[[gene_key]]
    gene_id <- gene_info$gene_id
    for (v in 1:gene_info$n_variants) {
      variant_id <- gene_info$variant_ids[v]
      all_variants_list[[idx]] <- data.frame(
        variant_id = variant_id,
        gene_id = gene_id,
        variant_index = v,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  all_variants <- do.call(rbind, all_variants_list)
  
  # Crea una matrice per l'espressione di tutte le varianti
  # Crea variant_matrix come matrice sparsa
  library(Matrix)
  variant_matrix <- Matrix(0, nrow = nrow(modified_expression), ncol = nrow(all_variants), sparse = TRUE)
  colnames(variant_matrix) <- all_variants$variant_id
  
  # Riempie la matrice con l'espressione delle varianti
  for (i in 1:nrow(all_variants)) {
    variant_id <- all_variants$variant_id[i]
    if (variant_id %in% names(variant_expression)) {
      variant_matrix[, i] <- variant_expression[[variant_id]]
    }
  }
  
  return(list(
    modified_expression = modified_expression,  # Espressione totale modificata
    variant_matrix = variant_matrix,            # Matrice con espressione di ogni variante
    variant_metadata = all_variants             # Metadati delle varianti
  ))
}

#' Genera modello completo di splicing alternativo
#'
#' Funzione wrapper che esegue l'intero processo di generazione
#' di varianti di splicing con pattern spaziali.
#'
#' @param cell_df Dataframe delle cellule con coordinate e cluster
#' @param expression_matrix Matrice di espressione originale
#' @param gene_info Informazioni sui geni (opzionale)
#' @param splicing_params Parametri per lo splicing alternativo
#' @param random_seed Seed per riproducibilità
#' @return Lista con matrici di espressione e informazioni di splicing
#' @export
generate_alternative_splicing <- function(
  cell_df,
  expression_matrix,
  gene_info = NULL,
  splicing_params = list(
    fraction_genes_with_variants = 0.15,  # Frazione di geni con splicing alternativo
    n_variants_per_gene = c(2, 3),        # Range di varianti per gene (min, max)
    variant_effect_strength = c(0.3, 1.5),  # Range di intensità dell'effetto (min, max)
    spatial_regulation = 0.7,             # Intensità della regolazione spaziale (0-1)
    cell_type_regulation = 0.6,           # Intensità della regolazione specifica per tipo cellulare
    coordinated_splicing_groups = 3,      # Numero di gruppi di splicing coordinato
    n_splicing_regulators = 6,            # Numero di "regolatori" di splicing
    splicing_noise = 0.2                  # Rumore nel pattern di splicing
  ),
  random_seed = 123
) {
  # Imposta seed per riproducibilità
  set.seed(random_seed)
  
  # Dimensioni
  n_genes <- ncol(expression_matrix)
  
  # 1. Genera varianti di splicing
  splicing_variants <- generate_splicing_variants(
    cell_df = cell_df,
    gene_info = gene_info,
    n_genes = n_genes,
    splicing_params = splicing_params,
    random_seed = random_seed
  )
  
  # 2. Calcola propensioni di splicing per ogni cellula
  splicing_propensity <- calculate_splicing_propensity(
    cell_df = cell_df,
    splicing_variants = splicing_variants,
    splicing_params = splicing_params,
    random_seed = random_seed
  )
  
  # 3. Genera espressione di varianti di splicing
  splicing_expression <- generate_splicing_expression(
    expression_matrix = expression_matrix,
    splicing_variants = splicing_variants,
    splicing_propensity = splicing_propensity,
    random_seed = random_seed
  )
  
  # Restituisci risultati
  return(list(
    modified_expression = splicing_expression$modified_expression,  # Espressione totale modificata
    variant_matrix = splicing_expression$variant_matrix,            # Matrice con espressione di ogni variante
    variant_metadata = splicing_expression$variant_metadata,        # Metadati delle varianti
    splicing_variants = splicing_variants,                          # Informazioni sulle varianti
    splicing_propensity = splicing_propensity,                      # Propensioni di splicing
    splicing_regulators = splicing_variants$splicing_regulators,    # Regolatori di splicing
    params = splicing_params                                        # Parametri utilizzati
  ))
}