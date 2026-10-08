# tools/S1.3_variants.R — motore alternativo per il confronto dei motori di S1.3 (addendum C-10).
# Il pacchetto ha UN solo motore (.tv_engine()); qui l'altro sostituisce .tv_tiles() in una copia dell'ambiente,
# come i mutanti. Uso: E <- tv_variant(tv_env(), "deldir"); E$tessellate_voronoi(...)
tv_variants <- list(
  deldir = list(
    .tv_engine = function() "deldir",
    .tv_tiles = function(x, y, rw) {
      n <- length(x)
      if (n == 1) return(sf::st_sfc(.tv_rect(rw)))
      dd <- deldir::deldir(x, y, rw = rw, round = FALSE)
      tl <- deldir::tile.list(dd)
      pt <- vapply(tl, function(t) t$ptNum, 1L)
      if (length(tl) != n || any(sort(pt) != seq_len(n))) stop("tessellate_voronoi(): deldir ha scartato o duplicato generatori.", call. = FALSE)
      tl <- tl[order(pt)]
      sf::st_sfc(lapply(tl, function(t) sf::st_polygon(list(cbind(c(t$x, t$x[1]), c(t$y, t$y[1]))))))
    }),
  geos = list()
)
tv_variant <- function(E, id) {
  if (!id %in% names(tv_variants)) stop("variante sconosciuta: ", id)
  for (nm in names(tv_variants[[id]])) { f <- tv_variants[[id]][[nm]]; environment(f) <- E; assign(nm, f, envir = E) }
  E
}
