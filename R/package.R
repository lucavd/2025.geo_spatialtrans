#' geo.spatialtrans: Un Pacchetto per la Simulazione di Dati di Trascrittomica Spaziale
#'
#' Framework avanzato per la simulazione di dati di trascrittomica spaziale
#' che incorpora complessi pattern biologici, correlazione spaziale
#' e artefatti tecnici realistici.
#'
#' @section Funzioni principali:
#' \itemize{
#'   \item \code{\link{simulate_spatial_transcriptomics}}: Funzione principale per la simulazione completa
#'   \item \code{\link{generate_expression_profiles}}: Genera profili di espressione per celle/spot
#'   \item \code{\link{prepare_image}}: Prepara un'immagine per la simulazione
#'   \item \code{\link{cluster_image}}: Esegue clustering sull'immagine
#'   \item \code{\link{create_sampling_grid}}: Crea griglia di campionamento
#'   \item \code{\link{generate_synthetic_tissue}}: Genera immagini sintetiche di tessuto
#'   \item \code{\link{analyze_and_compare_clusters}}: Analizza e confronta metodi di clustering
#' }
#'
#' @section Componenti modulari:
#' \itemize{
#'   \item Configurazione (01)
#'   \item Helper functions (02)
#'   \item Elaborazione immagini (03)
#'   \item Clustering (04)
#'   \item Grid sampling (05a-05c)
#'   \item Generazione profili espressione (06a-06k)
#'   \item Pipeline di simulazione (07a-07f)
#' }
#'
#' @docType package
#' @name geo.spatialtrans
NULL