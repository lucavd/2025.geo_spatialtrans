# tools/S1.3_mutants.R — mutanti di tessellate_voronoi() per misurare il potere dei check C (BL-048)
# Ogni mutante sostituisce UN helper in una copia isolata dell'ambiente di R/04b3_tessellate_voronoi.R.
# Uso: E <- tv_env(); E <- tv_mutate(E, "M3"); E$tessellate_voronoi(...)

tv_env <- function(file = "R/04b3_tessellate_voronoi.R") {
  E <- new.env(parent = globalenv())
  sys.source(file, envir = E)
  E
}

tv_mutants <- list(
  # M1: tile permutati di una posizione (associazione tile <-> generatore sbagliata)
  M1 = list(.tv_tiles = function(x, y, rw) {
    tl <- TV_ORIG$.tv_tiles(x, y, rw); n <- length(tl)
    if (n > 1) tl[c(2:n, 1)] else tl
  }),
  # M2: nessun ritaglio sulla regione
  M2 = list(.tv_clip = function(tl, reg) list(geoms = unclass(tl), clipped = rep(FALSE, length(tl)))),
  # M3: testo del design §5.3: Voronoi globale, ritaglio sulla regione del proprio centroide
  M3 = list(.tv_point_sets = function(ridx, nreg) list(list(idx = seq_along(ridx)))),
  # M4: Chaikin senza intersezione con il territorio non smussato
  M4 = list(.tv_smooth = function(geoms, x, y, n_iter) {
    out <- lapply(geoms, function(g) {
      if (inherits(g, "POLYGON")) sf::st_polygon(lapply(unclass(g), TV_ORIG$.tv_chaikin_ring, n_iter = n_iter)) else
        sf::st_multipolygon(lapply(unclass(g), function(p) lapply(p, TV_ORIG$.tv_chaikin_ring, n_iter = n_iter)))
    })
    list(geoms = out, n_dropped = 0L)
  }),
  # M5: frammenti orfani scartati (area persa)
  M5 = list(.tv_fragments = function(geoms, x, y, ridx) {
    n <- length(geoms); lost <- integer(n)
    for (i in seq_len(n)) {
      g <- geoms[[i]]
      if (inherits(g, "MULTIPOLYGON") && length(unclass(g)) > 1) {
        pcs <- TV_ORIG$.tv_pieces(g)
        own <- which(lengths(sf::st_intersects(sf::st_sfc(pcs), sf::st_sfc(sf::st_point(c(x[i], y[i]))))) > 0)[1]
        if (is.na(own)) own <- 1L
        geoms[[i]] <- pcs[[own]]; lost[i] <- length(pcs) - 1L
      }
    }
    list(geoms = geoms, lost = lost, gained = integer(n), n_fragments = sum(lost), n_reassigned = 0L, n_snapped = 0L,
         n_snap_ok = 0L, n_snap_failed = 0L, n_isolated_cells = 0L)
  }),
  # M7 (dopo la revisione, RA-code-12): orfano al confine condiviso PIU' CORTO invece che al piu' lungo
  M7 = list(.tv_pick = function(len, ids, lmin = 1e-6) { ok <- len > lmin; if (!any(ok)) return(NA_integer_); ids[ok & len == min(len[ok])][1] }),
  # M6: region_id sfasato di una posizione (indice della regione sbagliato)
  M6 = list(.tv_region_index = function(region_id, ids) { m <- match(region_id, ids); (m %% length(ids)) + 1L })
)

tv_mutate <- function(E, id) {
  if (identical(id, "none")) return(E)
  if (!id %in% names(tv_mutants)) stop("mutante sconosciuto: ", id)
  orig <- new.env(); for (nm in ls(E, all.names = TRUE)) assign(nm, get(nm, envir = E), envir = orig)
  assign("TV_ORIG", orig, envir = globalenv())
  for (nm in names(tv_mutants[[id]])) {
    f <- tv_mutants[[id]][[nm]]; environment(f) <- E; assign(nm, f, envir = E)
  }
  E
}
