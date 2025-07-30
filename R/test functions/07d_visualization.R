#' Crea visualizzazione dei risultati della simulazione
#'
#' Genera plot dei risultati della simulazione di dati di trascrittomica spaziale.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param config Configurazione di base della simulazione
#' @param difficulty_config Configurazione di difficoltà
#' @param expression_results Risultati dell'espressione genica
#' @return Oggetto ggplot con la visualizzazione
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_y_reverse coord_fixed theme_minimal labs
#' @export
create_simulation_plots <- function(
  cell_df,
  config,
  difficulty_config,
  expression_results = NULL
) {
  # Crea il plot in base alla modalità
  if (config$grid_mode) {
    # Per griglia, usiamo geom_tile
    p <- ggplot(cell_df, aes(x = x, y = y, fill = intensity_cluster)) +
      geom_tile(width = config$grid_resolution, height = config$grid_resolution) +
      scale_y_reverse() +
      coord_fixed() +
      theme_minimal() +
      labs(title = sprintf("Visium HD (%dum) - Livello difficoltà: %s", 
                           config$grid_resolution, 
                           difficulty_config$difficulty_level),
           subtitle = sprintf("Griglia %dum, %d bin, %d geni (marker/tipo: %d, fold: %.1f)",
                             config$grid_resolution, 
                             nrow(cell_df), 
                             config$n_genes,
                             difficulty_config$marker_params$marker_genes_per_type,
                             difficulty_config$marker_params$marker_expression_fold),
           fill = "Cell Type")
  } else {
    # Per sampling casuale, usiamo punti
    p <- ggplot(cell_df, aes(x = x, y = y, color = intensity_cluster)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_y_reverse() +
      coord_fixed() +
      theme_minimal() +
      labs(title = sprintf("Distribuzione spaziale (livello difficoltà: %s)", 
                           difficulty_config$difficulty_level),
           subtitle = sprintf("Threshold: %.2f, Geni marker: %d per tipo, Fold-change: %.1f",
                             config$threshold_value,
                             difficulty_config$marker_params$marker_genes_per_type,
                             difficulty_config$marker_params$marker_expression_fold),
           color = "Cell Type")
  }
  
  return(p)
}

#' Genera e salva plot della simulazione
#'
#' Crea e opzionalmente salva i plot dei risultati della simulazione.
#'
#' @param cell_df Dataframe delle celle con coordinate e cluster
#' @param config Configurazione di base della simulazione
#' @param difficulty_config Configurazione di difficoltà
#' @param output_plot Path dove salvare il plot (NULL per non salvare)
#' @param expression_results Risultati dell'espressione genica (opzionale)
#' @return Oggetto ggplot con la visualizzazione
#' @importFrom ggplot2 ggsave
#' @export
generate_and_save_plots <- function(
  cell_df,
  config,
  difficulty_config,
  output_plot = NULL,
  expression_results = NULL
) {
  # Crea il plot principale
  p <- create_simulation_plots(
    cell_df, config, difficulty_config, expression_results
  )
  
  # Mostra il plot
  print(p)
  
  # Salva il plot se richiesto
  if (!is.null(output_plot)) {
    # Crea directory se necessario
    output_dir <- dirname(output_plot)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    # Calcola l'aspect ratio della griglia spaziale
    x_range <- range(cell_df$x, na.rm = TRUE)
    y_range <- range(cell_df$y, na.rm = TRUE)
    aspect_ratio <- diff(y_range) / diff(x_range)
    # Dimensione base (inches) per il lato maggiore
    base_size <- 8
    if (aspect_ratio >= 1) {
      height_in <- base_size
      width_in  <- base_size / aspect_ratio
    } else {
      width_in  <- base_size
      height_in <- base_size * aspect_ratio
    }
    # Salva il plot con dimensioni esplicite per evitare ritagli
    ggplot2::ggsave(
      filename = output_plot,
      plot = p,
      device = "png",
      dpi = 300,
      width = width_in,
      height = height_in,
      units = "in",
      bg = "white"
    )
  }
  
  return(p)
}