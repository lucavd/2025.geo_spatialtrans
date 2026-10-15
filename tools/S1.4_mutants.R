# tools/S1.4_mutants.R — S1.4: ambiente di derive_cell_geometry() e mutanti per C-mut (pre-registrazione, commit c72c01a).
#  M1 clamp rimosso; M2 beta ignorato (= 1); M3 sigma ignorato (= 0); M4 POI -> centroide del territorio;
#  M5 area non corretta per il poligono (fattore k = 1). X1–X9: mutanti aggiunti dalla revisione avversariale (CODE-19).
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
  } else if (m %in% sprintf("X%d", 1:9)) {
    # mutanti aggiuntivi della revisione avversariale (revisore «code», CODE-19), testo sostituito nel sorgente
    mt <- function(fname, from, to) { src <- deparse(get(fname, e), width.cutoff = 500L)
      new <- gsub(from, to, src, fixed = TRUE); if (identical(src, new)) stop("sostituzione non applicata: ", fname, " / ", from)
      f <- eval(parse(text = new)); environment(f) <- e; assign(fname, f, envir = e) }
    switch(m,
      X1 = mt(".dcg_clamp", "for (it in seq_len(n_iter))", "for (it in seq_len(0L))"),           # niente bisezione
      X2 = mt(".dcg_clamp", "sc/cf", "sc"),                                                        # ellisse senza margine 1/cf
      X3 = { f <- e$.dcg_territory_axis; e$.dcg_territory_axis <- function(g) -f(g) },            # segno dell'asse
      X4 = mt("derive_cell_geometry", "z <- stats::rnorm(n)", "u <- stats::runif(n); z <- stats::rnorm(n); u <- u"),  # ordine del flusso
      X5 = mt(".dcg_territory_axis", "ifelse(hole, -1, 1) * ", ""),                                # buchi non sottratti
      X6 = mt(".dcg_clamp", "pmin(1, cf * d/", "pmin(1, d/"),                                      # clamp del cerchio a 1.0 d
      X7 = mt("derive_cell_geometry", "a_n <- s^2 * a_tgt", "a_n <- a_tgt"),                       # area non scalata
      X8 = mt(".dcg_orientation", "th <- pi * u - pi/2", "th <- pi * u"),                          # angolo in [0, pi)
      X9 = { e$.dcg_center <- function(terr, cen, mode, prec) { p <- polylabelr::poi(terr, precision = 1)
               list(x = vapply(p, function(q) q$x, 0), y = vapply(p, function(q) q$y, 0)) } })     # POI grossolano
  } else stop("mutante sconosciuto: ", m)
  e
}
