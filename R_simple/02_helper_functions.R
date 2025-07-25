#' Normalizza un vettore tra 0 e 1
#'
#' @param x Vettore numerico da normalizzare
#' @return Vettore normalizzato tra 0 e 1
#' @examples
#' scale01(c(1, 2, 3, 4, 5))
#' @export
scale01 <- function(x) {
  if (max(x) == min(x)) return(rep(0.5, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

#' Normalizza un vettore tra 0 e 1 (versione vettorizzata)
#'
#' @param x Vettore numerico da normalizzare
#' @return Vettore normalizzato tra 0 e 1
#' @noRd
scale01_vec <- function(x) {
  if (all(x == x[1])) return(rep(0.5, length(x)))
  (x - min(x)) / (max(x) - min(x))
}