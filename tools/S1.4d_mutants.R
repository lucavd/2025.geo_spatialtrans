# tools/S1.4d_mutants.R — S1.4d: mutanti del modello "clip" di derive_cell_geometry() per il cancello G-2
# (pre-registrazione docs/steps/S1.4d.md, commit 5118bce). Ogni mutante sostituisce un frammento del sorgente;
# se la sostituzione non si applica e' un errore (il mutante non esisterebbe).
#  M1 tetto ignorato; M2 ritaglio contro il territorio del vicino; M3 senza rientro; M4 rientro con modello di
#  precisione e senza guardia sui vuoti; M5 nessun ingrandimento (solo ritaglio); M6 sf senza WKB; M7 flusso invertito.
dcg_env_d <- function(file = "R/04b4_derive_cell_geometry.R") { e <- new.env(); sys.source(file, envir = e); e }
dcg_mutate_d <- function(e, m = "none") {
  if (m == "none") return(e)
  mt <- function(fname, from, to) {
    src <- deparse(get(fname, e), width.cutoff = 500L)
    new <- gsub(from, to, src, fixed = TRUE)
    if (identical(src, new)) stop("sostituzione non applicata: ", fname, " / ", from)
    f <- eval(parse(text = new)); environment(f) <- e; assign(fname, f, envir = e)
  }
  switch(m,
    M1 = mt(".dcg_run_clip", "a_star <- pmin(a_tgt, phi * pmin(a_t, geos::geos_area(tin)))", "a_star <- a_tgt"),
    M2 = mt(".dcg_run_clip", "kp, nv, tin, a_star", "kp, nv, tin[c(seq_along(tin)[-1], 1L)], a_star"),
    M3 = mt(".dcg_inset", "-inset", "-0 * inset"),
    M4 = { mt(".dcg_inset", "geos::geos_buffer(geos::geos_set_precision(tg, 0), -inset", "geos::geos_buffer(tg, -inset")
           mt(".dcg_inset", "if (any(emp))", "if (FALSE)") },
    M5 = mt(".dcg_clip_solve", "need <- which(area < a_star * (1 - tol))", "need <- integer(0)"),
    M6 = mt(".dcg_to_sfc", "sf::st_as_sfc(geos::geos_write_wkb(g))", "sf::st_as_sfc(g)"),
    M7 = mt("derive_cell_geometry", "z <- stats::rnorm(n)", "u <- stats::runif(n); z <- stats::rnorm(n); u <- u"),
    stop("mutante sconosciuto: ", m))
  e
}
