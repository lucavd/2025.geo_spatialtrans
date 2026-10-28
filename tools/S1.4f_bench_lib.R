# tools/S1.4f_bench_lib.R - S1.4f: nucleo orientato in parte lungo il territorio (leggi V0, V1, V2), cataloghi dei bracci, centri del banco A
# senza i doppioni al bordo (BL-174). Definizioni: docs/steps/S1.4f.md. Richiede: R/04b*.R, tools/R3_voronoi_metrics.R, tools/S1.4e_bench_lib.R.

#' Leggi della concentrazione: rho_i = 1 - exp(-lambda w_i) (derive_cell_geometry(), colonne nucleus_orient_conc e nucleus_orient_law).
S14F_LAWS <- c(V0 = "constant", V1 = "elongation", V2 = "compression")
#' Griglia di lambda della taratura (piu' ET = lambda infinito, preso da S1.4e).
S14F_GRID <- c(0, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10)

#' Catalogo di un braccio orientato in parte: l'ellisse di ET (allungamento vero del tessuto, nucleus_orient = "territory") con lambda e legge.
law_catalog <- function(ct, A, law, lambda) {
  stopifnot(law %in% names(S14F_LAWS), length(lambda) == 1, is.finite(lambda), lambda >= 0)
  ct <- arm_catalog(ct, "ET", A); ct$nucleus_orient_conc <- lambda; ct$nucleus_orient_law <- S14F_LAWS[[law]]; ct
}

#' BL-174: toglie dai centri del banco A (output di centers_A()) i centri del segmentatore entro r um da un punto di Luca.
s14f_dedup <- function(cc, r = 1.0) {
  co <- which(cc$src == "CO"); man <- which(cc$src == "manual")
  if (!length(co) || !length(man)) return(cc)
  dm <- rep(Inf, length(co))
  for (k in man) dm <- pmin(dm, (cc$x[co] - cc$x[k])^2 + (cc$y[co] - cc$y[k])^2)
  drop <- co[dm < r^2]
  if (length(drop)) cc[-drop, , drop = FALSE] else cc
}
