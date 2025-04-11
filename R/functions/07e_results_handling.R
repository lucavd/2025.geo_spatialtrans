#' Salva i risultati della simulazione
#'
#' Prepara e salva i risultati della simulazione di dati di trascrittomica spaziale.
#'
#' @param result Risultati della simulazione
#' @param output_path Path dove salvare i risultati
#' @return Path dove i risultati sono stati salvati
#' @export
save_simulation_results <- function(
  result,
  output_path
) {
  # Crea directory se non esiste
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  
  # Salva il risultato
  saveRDS(result, file = output_path)
  
  return(output_path)
}

#' Carica risultati di simulazione salvati
#'
#' Carica risultati precedentemente salvati di una simulazione.
#'
#' @param file_path Path del file salvato
#' @return Risultati della simulazione
#' @export
load_simulation_results <- function(
  file_path
) {
  if (!file.exists(file_path)) {
    stop("File non trovato: ", file_path)
  }
  
  result <- readRDS(file_path)
  return(result)
}