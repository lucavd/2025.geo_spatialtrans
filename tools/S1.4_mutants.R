# tools/S1.4_mutants.R — S1.4: ambiente di derive_cell_geometry() e mutanti per C-mut (pre-registrazione, commit c72c01a).
#  M1 clamp rimosso; M2 beta ignorato (= 1); M3 sigma ignorato (= 0); M4 POI -> centroide del territorio;
#  M5 area non corretta per il poligono (fattore k = 1).
dcg_env <- function(file = "R/04b4_derive_cell_geometry.R") { e <- new.env(); sys.source(file, envir = e); e }
dcg_mutate <- function(e, m = "none") {
  if (m == "none") return(e)
  if (m == "M1") {
    e$.dcg_clamp <- function(a_eq, ...) rep(1, length(a_eq))
  } else if (m == "M2") {
    f <- e$.dcg_target_area
    e$.dcg_target_area <- function(a_t, ct, ti, z) { ct$nucleus_area_slope <- 1; f(a_t, ct, ti, z) }
  } else if (m == "M3") {
    f <- e$.dcg_target_area
    e$.dcg_target_area <- function(a_t, ct, ti, z) f(a_t, ct, ti, 0 * z)
  } else if (m == "M4") {
    e$.dcg_center <- function(terr, cen, mode, prec) {
      cc <- sf::st_coordinates(suppressWarnings(sf::st_centroid(terr))); list(x = cc[, 1], y = cc[, 2]) }
  } else if (m == "M5") {
    e$.dcg_poly_factor <- function(nv) 1
  } else stop("mutante sconosciuto: ", m)
  e
}
