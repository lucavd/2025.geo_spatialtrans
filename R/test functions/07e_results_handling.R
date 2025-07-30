#' Salva i risultati della simulazione
#'
#' Prepara e salva i risultati della simulazione di dati di trascrittomica spaziale.
#'
#' @param result Risultati della simulazione
#' @param output_path Path dove salvare i risultati
#' @param save_module_data Se TRUE, salva anche i dati dei moduli biologici aggiuntivi
#' @return Path dove i risultati sono stati salvati
#' @export
save_simulation_results <- function(
  result,
  output_path,
  save_module_data = TRUE
) {
  # Crea directory se non esiste
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  
  # Verifica se ci sono dati dei moduli biologici aggiuntivi
  has_module_data <- !is.null(result$module_data) && length(result$module_data) > 0
  
  # Se richiesto, salva anche i dati dei moduli in file separati
  if (save_module_data && has_module_data) {
    module_dir <- file.path(dirname(output_path), "module_data")
    dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
    
    module_paths <- list()
    
    # Salva i dati di ogni modulo in un file separato
    for (module_name in names(result$module_data)) {
      module_path <- file.path(module_dir, paste0(module_name, ".rds"))
      saveRDS(result$module_data[[module_name]], file = module_path)
      module_paths[[module_name]] <- module_path
    }
    
    # Aggiunge i percorsi dei file dei moduli al risultato
    result$module_paths <- module_paths
  }
  
  # Salva il risultato principale
  saveRDS(result, file = output_path)
  
  return(output_path)
}

#' Carica risultati di simulazione salvati
#'
#' Carica risultati precedentemente salvati di una simulazione.
#'
#' @param file_path Path del file salvato
#' @param load_module_data Se TRUE, carica anche i dati dei moduli biologici aggiuntivi
#' @return Risultati della simulazione
#' @export
load_simulation_results <- function(
  file_path,
  load_module_data = TRUE
) {
  if (!file.exists(file_path)) {
    stop("File non trovato: ", file_path)
  }
  
  result <- readRDS(file_path)
  
  # Carica i dati dei moduli se richiesto e se esistono i percorsi
  if (load_module_data && !is.null(result$module_paths)) {
    for (module_name in names(result$module_paths)) {
      module_path <- result$module_paths[[module_name]]
      if (file.exists(module_path)) {
        result$module_data[[module_name]] <- readRDS(module_path)
      } else {
        warning("File del modulo non trovato: ", module_path)
      }
    }
  }
  
  return(result)
}

#' Visualizza i risultati della simulazione con rappresentazioni specifiche per moduli
#'
#' Crea visualizzazioni dei risultati di simulazione, inclusi gli effetti dei moduli biologici aggiuntivi.
#'
#' @param result Risultati della simulazione
#' @param output_dir Directory dove salvare le visualizzazioni
#' @param prefix Prefisso per i nomi dei file
#' @param plot_modules Quali moduli visualizzare (NULL = tutti)
#' @return Lista con percorsi dei plot generati
#' @export
visualize_simulation_results <- function(
  result,
  output_dir = "results",
  prefix = "sim",
  plot_modules = NULL
) {
  # Crea directory se non esiste
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Lista per salvare i percorsi dei plot
  plot_paths <- list()
  
  # Plot base dell'espressione genica
  expression_plot <- visualize_spatial_expression(
    result$coordinates,
    result$expression,
    result$intensity_cluster
  )
  
  # Salva il plot base
  expression_path <- file.path(output_dir, paste0(prefix, "_expression.png"))
  ggplot2::ggsave(expression_path, expression_plot, width = 10, height = 8)
  plot_paths$expression <- expression_path
  
  # Determina quali moduli visualizzare
  available_modules <- names(result$module_data)
  if (is.null(plot_modules)) {
    plot_modules <- available_modules
  } else {
    plot_modules <- intersect(plot_modules, available_modules)
  }
  
  # Visualizzazioni specifiche per ogni modulo
  for (module_name in plot_modules) {
    if (module_name == "lr_interactions" && !is.null(result$module_data$lr_interactions)) {
      # Visualizzazione delle interazioni ligando-recettore
      lr_plot <- visualize_lr_interactions(
        result$coordinates,
        result$module_data$lr_interactions$signaling_effects
      )
      lr_plot_path <- file.path(output_dir, paste0(prefix, "_lr_interactions.png"))
      ggplot2::ggsave(lr_plot_path, lr_plot, width = 10, height = 8)
      plot_paths$lr_interactions <- lr_plot_path
    }
    
    if (module_name == "temporal_dynamics" && !is.null(result$module_data$temporal_dynamics)) {
      # Visualizzazione della pseudotime e della velocità RNA
      temporal_plot <- visualize_pseudotime(
        result$coordinates,
        result$module_data$temporal_dynamics$pseudotime
      )
      temporal_plot_path <- file.path(output_dir, paste0(prefix, "_pseudotime.png"))
      ggplot2::ggsave(temporal_plot_path, temporal_plot, width = 10, height = 8)
      plot_paths$pseudotime <- temporal_plot_path
      
      if (!is.null(result$module_data$temporal_dynamics$velocity)) {
        velocity_plot <- visualize_rna_velocity(
          result$coordinates,
          result$module_data$temporal_dynamics$velocity
        )
        velocity_plot_path <- file.path(output_dir, paste0(prefix, "_rna_velocity.png"))
        ggplot2::ggsave(velocity_plot_path, velocity_plot, width = 10, height = 8)
        plot_paths$rna_velocity <- velocity_plot_path
      }
    }
    
    if (module_name == "alternative_splicing" && !is.null(result$module_data$alternative_splicing)) {
      # Visualizzazione dello splicing alternativo
      if (length(result$module_data$alternative_splicing$genes_with_variants) > 0) {
        splicing_plot <- visualize_alternative_splicing(
          result$coordinates,
          result$module_data$alternative_splicing$variant_matrices
        )
        splicing_plot_path <- file.path(output_dir, paste0(prefix, "_splicing.png"))
        ggplot2::ggsave(splicing_plot_path, splicing_plot, width = 10, height = 8)
        plot_paths$alternative_splicing <- splicing_plot_path
      }
    }
    
    if (module_name == "anisotropic_patterns" && !is.null(result$module_data$anisotropic_patterns)) {
      # Visualizzazione delle strutture anisotropiche
      aniso_plot <- visualize_anisotropic_structures(
        result$coordinates,
        result$module_data$anisotropic_patterns$structure_mask
      )
      aniso_plot_path <- file.path(output_dir, paste0(prefix, "_anisotropic.png"))
      ggplot2::ggsave(aniso_plot_path, aniso_plot, width = 10, height = 8)
      plot_paths$anisotropic_patterns <- aniso_plot_path
    }
    
    if (module_name == "microenvironment_3d" && !is.null(result$module_data$microenvironment_3d)) {
      # Visualizzazione del microambiente 3D
      micro3d_plot <- visualize_3d_layers(
        result$coordinates,
        result$module_data$microenvironment_3d$z_positions
      )
      micro3d_plot_path <- file.path(output_dir, paste0(prefix, "_3d_layers.png"))
      ggplot2::ggsave(micro3d_plot_path, micro3d_plot, width = 10, height = 8)
      plot_paths$microenvironment_3d <- micro3d_plot_path
    }
  }
  
  return(plot_paths)
}

# Funzioni di utilità per visualizzazione

#' Visualizza l'espressione genica spaziale
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_minimal labs
visualize_spatial_expression <- function(coordinates, expression, clusters) {
  # Crea un dataframe per la visualizzazione
  df <- data.frame(
    x = coordinates$x,
    y = coordinates$y,
    cluster = clusters
  )
  
  # Visualizza la distribuzione spaziale dei cluster
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster)) +
    ggplot2::geom_point() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Spatial Expression Patterns",
                  x = "X Coordinate", 
                  y = "Y Coordinate",
                  color = "Cell Type")
  
  return(p)
}

#' Visualizza interazioni ligando-recettore
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_minimal labs
visualize_lr_interactions <- function(coordinates, signaling_effects) {
  # Crea un dataframe per la visualizzazione
  df <- data.frame(
    x = coordinates$x,
    y = coordinates$y,
    signal_strength = rowMeans(signaling_effects, na.rm = TRUE)
  )
  
  # Visualizza la forza del segnale delle interazioni L-R
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = signal_strength)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Ligand-Receptor Signaling Effects",
                  x = "X Coordinate", 
                  y = "Y Coordinate",
                  color = "Signal Strength")
  
  return(p)
}

#' Visualizza pseudotime
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_minimal labs
visualize_pseudotime <- function(coordinates, pseudotime) {
  # Crea un dataframe per la visualizzazione
  df <- data.frame(
    x = coordinates$x,
    y = coordinates$y,
    pseudotime = pseudotime
  )
  
  # Visualizza il pseudotime
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = pseudotime)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "viridis") +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Pseudotime Trajectory",
                  x = "X Coordinate", 
                  y = "Y Coordinate",
                  color = "Pseudotime")
  
  return(p)
}

#' Visualizza RNA velocity
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_y_reverse theme_minimal labs
visualize_rna_velocity <- function(coordinates, velocity) {
  # Crea un dataframe per la visualizzazione della velocità
  df <- data.frame(
    x = coordinates$x,
    y = coordinates$y,
    vx = velocity$vx,
    vy = velocity$vy,
    velocity_magnitude = sqrt(velocity$vx^2 + velocity$vy^2)
  )
  
  # Normalizza i vettori di velocità per una migliore visualizzazione
  scale_factor <- 5 / max(df$velocity_magnitude, na.rm = TRUE)
  df$vx_scaled <- df$vx * scale_factor
  df$vy_scaled <- df$vy * scale_factor
  
  # Visualizza i vettori di velocità RNA
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_segment(
      ggplot2::aes(x = x, y = y, xend = x + vx_scaled, yend = y - vy_scaled, color = velocity_magnitude),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"))
    ) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y), alpha = 0.3) +
    ggplot2::scale_color_viridis_c(option = "magma") +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "RNA Velocity Vectors",
                  x = "X Coordinate", 
                  y = "Y Coordinate",
                  color = "Velocity Magnitude")
  
  return(p)
}

#' Visualizza splicing alternativo
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c facet_wrap theme_minimal labs
visualize_alternative_splicing <- function(coordinates, variant_matrices) {
  # Seleziona il primo gene con varianti per visualizzazione
  if (length(variant_matrices) == 0) {
    return(NULL)
  }
  
  # Prendi il primo gene con varianti
  gene_name <- names(variant_matrices)[1]
  variant_matrix <- variant_matrices[[gene_name]]
  
  # Crea un dataframe per visualizzare il rapporto tra varianti
  if (ncol(variant_matrix) >= 2) {
    df <- data.frame(
      x = coordinates$x,
      y = coordinates$y,
      variant_ratio = variant_matrix[, 1] / (variant_matrix[, 1] + variant_matrix[, 2])
    )
    
    # Visualizza il rapporto tra le varianti
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = variant_ratio)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_viridis_c(option = "cividis") +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Alternative Splicing for Gene:", gene_name),
                    x = "X Coordinate", 
                    y = "Y Coordinate",
                    color = "Variant Ratio")
    
    return(p)
  } else {
    return(NULL)
  }
}

#' Visualizza strutture anisotropiche
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_minimal labs
visualize_anisotropic_structures <- function(coordinates, structure_mask) {
  # Crea un dataframe per la visualizzazione
  df <- data.frame(
    x = coordinates$x,
    y = coordinates$y
  )
  
  # Se structure_mask è una lista, prendi la prima struttura
  if (is.list(structure_mask) && length(structure_mask) > 0) {
    df$structure_value <- structure_mask[[1]]
  } else {
    df$structure_value <- structure_mask
  }
  
  # Visualizza le strutture anisotropiche
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = structure_value)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "inferno") +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Anisotropic Structures",
                  x = "X Coordinate", 
                  y = "Y Coordinate",
                  color = "Structure Value")
  
  return(p)
}

#' Visualizza gli strati 3D
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_minimal labs
visualize_3d_layers <- function(coordinates, z_positions) {
  # Crea un dataframe per la visualizzazione
  df <- data.frame(
    x = coordinates$x,
    y = coordinates$y,
    z = z_positions
  )
  
  # Visualizza le posizioni Z
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = z)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "turbo") +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "3D Microenvironment Layers",
                  x = "X Coordinate", 
                  y = "Y Coordinate",
                  color = "Z Position")
  
  return(p)
}