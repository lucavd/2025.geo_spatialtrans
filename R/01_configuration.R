#' Configurazione per ottimizzare la performance delle simulazioni
#'
#' Imposta i parametri ottimali per la parallelizzazione e
#' gestione della memoria durante le simulazioni.
#'
#' @param max_memory_gb Dimensione massima allocabile in GB (default: 50)
#' @param workers Numero di worker in parallelo (default: 16)
#' @return Invisibly restituisce il numero di worker configurati
#' @importFrom future plan multisession nbrOfWorkers availableCores
#' @export
setup_performance <- function(max_memory_gb = 50, workers = 16) {
  options(future.globals.maxSize = max_memory_gb * 1024^2)
  options(future.fork.enable = FALSE)
  options(future.gc = TRUE)
  plan(multisession, workers = workers)
  
  invisible(future::nbrOfWorkers())
}