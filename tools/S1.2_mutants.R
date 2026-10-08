# tools/S1.2_mutants.R — mutanti di seed_centroids() per misurare il potere dei check C (BL-048)
# Ogni mutante sostituisce UN helper in una copia isolata dell'ambiente di R/04b2_seed_centroids.R.
# Uso: E <- sc_env(); E <- sc_mutate(E, "M3"); E$seed_centroids(...)

sc_env <- function(file = "R/04b2_seed_centroids.R") {
  E <- new.env(parent = globalenv())
  sys.source(file, envir = E)
  E
}

sc_mutants <- list(
  # M1: campionamento sul riquadro, senza test di appartenenza al poligono
  M1 = list(.sc_in_polygon = function(poly, x, y) rep(TRUE, length(x))),
  # M2: ricerca dei vicini solo nella propria cella della griglia
  M2 = list(.sc_neighbour_offsets = function(nxp) 0L),
  # M3: arrotondamento deterministico del numero di cellule per regione
  M3 = list(.sc_round_regions = function(lambda) as.integer(round(lambda))),
  # M3b: largest remainder deterministico per tipo
  M3b = list(.sc_allocate_types = function(n, fraction) {
    e <- n * fraction; fl <- floor(e + 1e-9); r <- n - sum(fl)
    if (r > 0) { o <- order(-(e - fl), seq_along(e))[seq_len(r)]; fl[o] <- fl[o] + 1 }
    as.integer(fl)
  }),
  # M5: set.seed() interno senza ripristino dello stato globale
  M5 = list(.sc_with_seed = function(seed, expr) { set.seed(seed); expr }),
  # M6: densita' di miscela aritmetica (design v1.1) invece di armonica (D3)
  M6 = list(.sc_mix_density = function(fraction, density) sum(fraction * density)),
  # M4: Bridson (2007, doi:10.1145/1278780.1278807) troncato a n, un tipo, regione rettangolare
  M4 = list(.sc_place_all = function(tasks, polys, rdf, type_df, max_attempts_factor) {
    X <- numeric(0); Y <- numeric(0); TI <- integer(0); RI <- integer(0)
    failed <- integer(nrow(tasks)); atts <- integer(nrow(tasks))
    for (t in seq_len(nrow(tasks))) {
      bb <- sf::st_bbox(polys[tasks$region_idx[t]]); n <- tasks$n[t]
      d <- type_df$min_dist_um[tasks$type_idx[t]]; cs <- d / sqrt(2)
      nx <- ceiling((bb[["xmax"]] - bb[["xmin"]]) / cs) + 1; ny <- ceiling((bb[["ymax"]] - bb[["ymin"]]) / cs) + 1
      G <- matrix(0L, nx, ny)
      gi <- function(x) floor((x - bb[["xmin"]]) / cs) + 1; gj <- function(y) floor((y - bb[["ymin"]]) / cs) + 1
      px <- stats::runif(1, bb[["xmin"]], bb[["xmax"]]); py <- stats::runif(1, bb[["ymin"]], bb[["ymax"]])
      G[gi(px), gj(py)] <- 1L; active <- 1L; att <- 1L
      while (length(active) && length(px) < n) {
        a <- active[as.integer(stats::runif(1) * length(active)) + 1L]
        found <- FALSE
        for (k in 1:30) {
          att <- att + 1L
          rr <- d * sqrt(stats::runif(1, 1, 4)); th <- stats::runif(1, 0, 2 * pi)
          x <- px[a] + rr * cos(th); y <- py[a] + rr * sin(th)
          if (x < bb[["xmin"]] || x > bb[["xmax"]] || y < bb[["ymin"]] || y > bb[["ymax"]]) next
          i <- gi(x); j <- gj(y)
          ii <- max(1, i - 2):min(nx, i + 2); jj <- max(1, j - 2):min(ny, j + 2)
          ids <- G[ii, jj]; ids <- ids[ids > 0L]
          if (length(ids) && any((px[ids] - x)^2 + (py[ids] - y)^2 < d^2)) next
          px <- c(px, x); py <- c(py, y); G[i, j] <- length(px); active <- c(active, length(px)); found <- TRUE
          break
        }
        if (!found) active <- active[active != a]
      }
      X <- c(X, px); Y <- c(Y, py); TI <- c(TI, rep(tasks$type_idx[t], length(px)))
      RI <- c(RI, rep(tasks$region_id[t], length(px))); failed[t] <- n - length(px); atts[t] <- att
    }
    list(x = X, y = Y, type_idx = TI, region_id = RI, task_failed = failed, task_attempts = atts,
         grid_cell_um = NA_real_, hardcore = TRUE)
  })
)

sc_mutate <- function(E, id) {
  E2 <- sc_env()
  for (nm in names(sc_mutants[[id]])) {
    f <- sc_mutants[[id]][[nm]]
    environment(f) <- E2
    assign(nm, f, envir = E2)
  }
  E2
}
