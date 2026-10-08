#' Tassella le regioni in territori cellulari di Voronoi (Step 1, S1.3)
#'
#' Per ogni regione di `extract_regions()`: diagramma di Voronoi dei soli
#' centroidi della regione (decisione D-S1.3.1, 2026-10-08), ritagliato sulla
#' regione. Se il ritaglio spezza una cella, il pezzo che contiene il generatore
#' resta alla cellula e ogni pezzo orfano passa al territorio della stessa
#' regione con cui condivide il confine più lungo (D-S1.3.2). Con
#' `corner_smoothing > 0` ogni territorio è smussato con Chaikin e intersecato
#' con il territorio non smussato (D-S1.3.3): le lacune sono spazio
#' extracellulare, riportate in `region_check$gap_area`.
#' Design: `docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md` §5.3,
#' con le deviazioni dichiarate in `results/S1.3/S1.3_preregistration.md`.
#'
#' Convenzione spaziale: µm, y verso il basso, come `extract_regions()` e
#' `seed_centroids()` (punti in coordinate continue, regioni sui lati dei pixel).
#'
#' @param centroids data.frame con `cell_id` (unico), `region_id`, `x`, `y`
#'   (µm), ad esempio `seed_centroids()$centroids`. Ogni centroide deve stare
#'   nella propria regione; nessun duplicato di posizione nella stessa regione.
#' @param regions output di `extract_regions()`: lista con `region_df`
#'   (`region_id`) e `region_polygons` (sfc POLYGON o MULTIPOLYGON, µm), nello
#'   stesso ordine.
#' @param corner_smoothing 0 = Voronoi puro; in (0, 1] → Chaikin con
#'   `max(1, round(3 * corner_smoothing))` iterazioni.
#' @param keep_tiles se TRUE restituisce anche i tile non ritagliati (diagnostica,
#'   regola d'interna di R3).
#' @param verbose stampa un riepilogo.
#' @return lista con `cell_territories` (sfc, ordine di `cell_id` crescente),
#'   `territory_df` (`cell_id`, `region_id`, `territory_area`, `tile_area`,
#'   `clipped`, `n_pieces_lost`, `n_pieces_gained`, `smooth_loss`),
#'   `region_check` (`region_id`, `area_um2`, `sum_territory`, `gap_area`,
#'   `n_cells`), `tiles` (sfc o NULL), `info`.
#' @export
tessellate_voronoi <- function(centroids,
                               regions,
                               corner_smoothing = 0,
                               keep_tiles = FALSE,
                               verbose = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  a <- .tv_check_args(centroids, regions, corner_smoothing, keep_tiles)
  cen <- a$centroids
  rdf <- regions$region_df
  polys <- regions$region_polygons
  n <- nrow(cen); nreg <- nrow(rdf)
  n_iter <- .tv_n_iter(corner_smoothing)
  ridx <- .tv_region_index(cen$region_id, rdf$region_id)

  geoms <- vector("list", n)
  tile_sfg <- if (keep_tiles) vector("list", n) else NULL
  tile_area <- rep(NA_real_, n); clipped <- rep(FALSE, n)
  for (s in .tv_point_sets(ridx, nreg)) {
    k <- s$idx
    rw <- .tv_window(polys[unique(ridx[k])])
    tl <- .tv_tiles(cen$x[k], cen$y[k], rw)
    tile_area[k] <- as.numeric(sf::st_area(tl))
    if (keep_tiles) tile_sfg[k] <- unclass(tl)
    for (r in unique(ridx[k])) {
      kk <- which(ridx[k] == r)
      cl <- .tv_clip(tl[kk], polys[r])
      geoms[k[kk]] <- cl$geoms
      clipped[k[kk]] <- cl$clipped
    }
  }

  fr <- .tv_fragments(geoms, cen$x, cen$y, ridx)
  geoms <- fr$geoms

  smooth_loss <- numeric(n); n_smooth_dropped <- 0L
  if (n_iter > 0) {
    sm <- .tv_smooth(geoms, cen$x, cen$y, n_iter)
    smooth_loss <- vapply(seq_len(n), function(i) .tv_area(geoms[[i]]) - .tv_area(sm$geoms[[i]]), 0)
    geoms <- sm$geoms; n_smooth_dropped <- sm$n_dropped
  }

  fin <- .tv_finalise(geoms)
  terr <- fin$sfc
  t_area <- as.numeric(sf::st_area(terr))
  tdf <- data.frame(cell_id = cen$cell_id, region_id = cen$region_id, territory_area = t_area,
                    tile_area = tile_area, clipped = clipped, n_pieces_lost = fr$lost,
                    n_pieces_gained = fr$gained, smooth_loss = smooth_loss)
  r_area <- if (nreg) as.numeric(sf::st_area(polys)) else numeric(0)
  s_terr <- vapply(seq_len(nreg), function(r) sum(t_area[ridx == r]), 0)
  n_cells <- tabulate(ridx, nreg)
  rcheck <- data.frame(region_id = rdf$region_id, area_um2 = r_area, sum_territory = s_terr,
                       gap_area = ifelse(n_cells > 0, r_area - s_terr, NA_real_), n_cells = n_cells)
  info <- list(n_cells = n, n_regions = nreg, n_regions_without_cells = sum(n_cells == 0),
               engine = .tv_engine(), corner_smoothing = corner_smoothing, n_iter = n_iter,
               n_fragments = fr$n_fragments, n_fragments_reassigned = fr$n_reassigned, n_snapped = fr$n_snapped,
               n_multipart = fin$n_multipart, n_repaired = fin$n_repaired,
               n_smooth_pieces_dropped = n_smooth_dropped,
               elapsed_s = proc.time()[["elapsed"]] - t0)
  if (verbose) {
    message(sprintf("tessellate_voronoi(): %d territori in %d regioni (%s), %d frammenti riassegnati, %d multiparte, %d riparati, %.1f s",
                    n, nreg, info$engine, info$n_fragments_reassigned, info$n_multipart, info$n_repaired, info$elapsed_s))
  }
  list(cell_territories = terr, territory_df = tdf, region_check = rcheck,
       tiles = if (keep_tiles) sf::st_sfc(tile_sfg) else NULL, info = info)
}

# ---- helper interni ---------------------------------------------------------

#' Validazione degli argomenti; restituisce i centroidi ordinati per cell_id.
#' @keywords internal
.tv_check_args <- function(centroids, regions, corner_smoothing, keep_tiles) {
  if (!is.data.frame(centroids)) stop("tessellate_voronoi(): `centroids` deve essere un data.frame.", call. = FALSE)
  need <- c("cell_id", "region_id", "x", "y")
  miss <- setdiff(need, names(centroids))
  if (length(miss)) stop("tessellate_voronoi(): centroids senza colonne ", paste(miss, collapse = ", "), call. = FALSE)
  if (!is.list(regions) || is.null(regions$region_df) || is.null(regions$region_polygons)) {
    stop("tessellate_voronoi(): `regions` deve essere l'output di extract_regions() (region_df + region_polygons).", call. = FALSE)
  }
  rdf <- regions$region_df
  if (!"region_id" %in% names(rdf)) stop("tessellate_voronoi(): region_df senza region_id.", call. = FALSE)
  if (length(regions$region_polygons) != nrow(rdf)) stop("tessellate_voronoi(): region_polygons e region_df hanno lunghezze diverse.", call. = FALSE)
  if (is.factor(rdf$region_id) || is.factor(centroids$region_id)) stop("tessellate_voronoi(): region_id e' un factor: passare interi o caratteri.", call. = FALSE)
  if (anyDuplicated(rdf$region_id)) stop("tessellate_voronoi(): region_id duplicati in region_df.", call. = FALSE)
  if (!is.numeric(corner_smoothing) || length(corner_smoothing) != 1 || !is.finite(corner_smoothing) ||
      corner_smoothing < 0 || corner_smoothing > 1) {
    stop("tessellate_voronoi(): corner_smoothing deve essere un numero in [0, 1].", call. = FALSE)
  }
  if (!is.logical(keep_tiles) || length(keep_tiles) != 1 || is.na(keep_tiles)) stop("tessellate_voronoi(): keep_tiles deve essere TRUE/FALSE.", call. = FALSE)
  if (!is.numeric(centroids$x) || !is.numeric(centroids$y) || any(!is.finite(centroids$x)) || any(!is.finite(centroids$y))) {
    stop("tessellate_voronoi(): x e y devono essere numerici finiti.", call. = FALSE)
  }
  if (anyDuplicated(centroids$cell_id)) stop("tessellate_voronoi(): cell_id duplicati.", call. = FALSE)
  bad_r <- !(centroids$region_id %in% rdf$region_id)
  if (any(bad_r)) stop(sprintf("tessellate_voronoi(): %d centroidi con region_id assente da region_df (es. %s).",
                               sum(bad_r), as.character(centroids$region_id[which(bad_r)[1]])), call. = FALSE)
  polys <- regions$region_polygons
  if (nrow(rdf) && any(sf::st_is_empty(polys))) stop("tessellate_voronoi(): region_polygons contiene geometrie vuote.", call. = FALSE)
  gt <- as.character(sf::st_geometry_type(polys))
  if (any(!gt %in% c("POLYGON", "MULTIPOLYGON"))) stop("tessellate_voronoi(): region_polygons deve contenere POLYGON o MULTIPOLYGON.", call. = FALSE)
  cen <- centroids[order(centroids$cell_id, method = "radix"), c("cell_id", "region_id", "x", "y")]
  rownames(cen) <- NULL
  ridx <- match(cen$region_id, rdf$region_id)
  dup <- duplicated(data.frame(ridx, cen$x, cen$y))
  if (any(dup)) stop(sprintf("tessellate_voronoi(): %d centroidi duplicati (stessa posizione nella stessa regione).", sum(dup)), call. = FALSE)
  for (r in unique(ridx)) {
    k <- which(ridx == r)
    pts <- sf::st_as_sf(data.frame(x = cen$x[k], y = cen$y[k]), coords = c("x", "y"))
    inside <- lengths(sf::st_intersects(sf::st_geometry(pts), polys[r])) > 0
    if (any(!inside)) stop(sprintf("tessellate_voronoi(): %d centroidi fuori dalla propria regione (region_id %s).",
                                   sum(!inside), as.character(rdf$region_id[r])), call. = FALSE)
  }
  list(centroids = cen)
}

#' Numero di iterazioni di Chaikin (D-S1.3.3).
#' @keywords internal
.tv_n_iter <- function(cs) if (cs == 0) 0L else as.integer(max(1, round(3 * cs)))

#' Indice della regione (riga di region_df) per ogni centroide.
#' @keywords internal
.tv_region_index <- function(region_id, ids) match(region_id, ids)

#' Insiemi di punti tassellati insieme: uno per regione (D-S1.3.1).
#' @keywords internal
.tv_point_sets <- function(ridx, nreg) {
  lapply(Filter(function(r) any(ridx == r), seq_len(nreg)), function(r) list(idx = which(ridx == r)))
}

#' Riquadro di calcolo: bbox delle regioni allargato di max(1 µm, 10 % del lato).
#' @keywords internal
.tv_window <- function(p) {
  bb <- sf::st_bbox(p)
  pad <- max(1, 0.1 * max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]]))
  c(bb[["xmin"]] - pad, bb[["xmax"]] + pad, bb[["ymin"]] - pad, bb[["ymax"]] + pad)
}

#' Poligono del riquadro rw = c(xmin, xmax, ymin, ymax).
#' @keywords internal
.tv_rect <- function(rw) sf::st_polygon(list(rbind(c(rw[1], rw[3]), c(rw[2], rw[3]), c(rw[2], rw[4]), c(rw[1], rw[4]), c(rw[1], rw[3]))))

#' Motore del diagramma (uno solo nel pacchetto; confronto in tools/S1.3_variants.R).
#' @keywords internal
.tv_engine <- function() "geos"

#' Tile di Voronoi nel riquadro, nell'ordine dei punti (sfc POLYGON): GEOS
#' (`sf::st_voronoi`, `point_order = TRUE`), celle riportate al riquadro.
#' @keywords internal
.tv_tiles <- function(x, y, rw) {
  n <- length(x)
  if (n == 1) return(sf::st_sfc(.tv_rect(rw)))
  env <- .tv_rect(rw)
  v <- sf::st_voronoi(sf::st_multipoint(cbind(x, y)), envelope = env, point_order = TRUE)
  cells <- sf::st_collection_extract(sf::st_sfc(v), "POLYGON")
  if (length(cells) != n) stop("tessellate_voronoi(): GEOS ha restituito ", length(cells), " celle per ", n, " generatori.", call. = FALSE)
  # le celle GEOS si estendono oltre l'envelope: le si riporta al riquadro, come deldir
  envs <- sf::st_sfc(env)
  inside <- lengths(sf::st_within(cells, envs)) > 0
  if (any(!inside)) {
    out <- sf::st_intersection(cells[!inside], envs)
    cells[!inside] <- sf::st_cast(out, "POLYGON")
  }
  cells
}

#' Ritaglio dei tile sulla regione: solo i tile non contenuti vengono intersecati.
#' @keywords internal
.tv_clip <- function(tl, reg) {
  within <- lengths(sf::st_within(tl, reg)) > 0
  geoms <- unclass(tl)
  nb <- which(!within)
  if (length(nb)) {
    rg <- reg[[1]]
    for (i in nb) {
      tb <- sf::st_bbox(tl[i])
      loc <- suppressWarnings(sf::st_crop(sf::st_sfc(rg), tb))
      g <- if (length(loc) == 0 || all(sf::st_is_empty(loc))) sf::st_polygon() else
        suppressWarnings(sf::st_intersection(tl[i], loc))
      geoms[[i]] <- .tv_polys_only(g)
    }
  }
  list(geoms = geoms, clipped = !within)
}

#' Solo la parte poligonale di una geometria (POLYGON, MULTIPOLYGON o POLYGON vuoto).
#' @keywords internal
.tv_polys_only <- function(g) {
  if (inherits(g, "sfc")) g <- if (length(g)) g[[1]] else sf::st_polygon()
  if (inherits(g, c("POLYGON", "MULTIPOLYGON"))) return(g)
  if (inherits(g, "GEOMETRYCOLLECTION")) {
    p <- Filter(function(z) inherits(z, c("POLYGON", "MULTIPOLYGON")), unclass(g))
    if (!length(p)) return(sf::st_polygon())
    u <- sf::st_union(sf::st_sfc(p))
    return(u[[1]])
  }
  sf::st_polygon()
}

#' Area di una geometria sfg.
#' @keywords internal
.tv_area <- function(g) if (sf::st_is_empty(g)) 0 else as.numeric(sf::st_area(sf::st_sfc(g)))

#' Pezzi POLYGON di una geometria.
#' @keywords internal
.tv_pieces <- function(g) {
  if (inherits(g, "POLYGON")) return(list(g))
  lapply(unclass(g), function(p) sf::st_polygon(p))
}

#' Frammenti (D-S1.3.2): il pezzo con il generatore resta; gli orfani vanno al
#' territorio della stessa regione con il confine condiviso più lungo.
#' @keywords internal
.tv_fragments <- function(geoms, x, y, ridx) {
  n <- length(geoms); lost <- integer(n); gained <- integer(n)
  multi <- which(vapply(geoms, function(g) inherits(g, "MULTIPOLYGON") && length(unclass(g)) > 1, TRUE))
  orphans <- list()
  for (i in multi) {
    pcs <- .tv_pieces(geoms[[i]])
    pt <- sf::st_sfc(sf::st_point(c(x[i], y[i])))
    own <- which(lengths(sf::st_intersects(sf::st_sfc(pcs), pt)) > 0)
    if (!length(own)) own <- which.max(vapply(pcs, .tv_area, 0))
    own <- own[1]
    geoms[[i]] <- pcs[[own]]
    for (j in setdiff(seq_along(pcs), own)) orphans[[length(orphans) + 1]] <- list(g = pcs[[j]], from = i)
    lost[i] <- length(pcs) - 1L
  }
  n_frag <- length(orphans); n_re <- 0L; n_snap <- 0L
  if (n_frag) {
    pending <- seq_len(n_frag)
    for (pass in 1:10) {
      progress <- FALSE
      for (o in pending) {
        og <- orphans[[o]]$g; r <- ridx[orphans[[o]]$from]
        cand <- which(ridx == r)
        bbs <- sf::st_sfc(lapply(geoms[cand], function(g) if (sf::st_is_empty(g)) sf::st_polygon() else sf::st_as_sfc(sf::st_bbox(g))[[1]]))
        hit <- cand[lengths(sf::st_intersects(bbs, sf::st_sfc(sf::st_as_sfc(sf::st_bbox(og))[[1]]))) > 0]
        if (!length(hit)) next
        ob <- sf::st_boundary(sf::st_sfc(og))
        len <- vapply(hit, function(j) {
          # tolleranza 1e-7 µm: i vertici condivisi di due tile deldir differiscono all'ultima cifra (S1.3, prova di fumo)
          s <- suppressWarnings(sf::st_intersection(ob, sf::st_buffer(sf::st_boundary(sf::st_sfc(geoms[[j]])), 1e-7)))
          if (!length(s) || all(sf::st_is_empty(s))) return(0)
          sum(as.numeric(sf::st_length(s)))      # GEOS: punti 0, linee e collezioni = somma delle parti lineari
        }, 0)
        if (max(len) <= 0) next
        j <- hit[which(len == max(len))][1]
        u <- sf::st_union(sf::st_sfc(geoms[[j]]), sf::st_sfc(og))[[1]]
        if (inherits(u, "MULTIPOLYGON") && length(unclass(u)) > 1) {
          # confini non coincidenti all'ultima cifra (deldir): aggancio dell'orfano entro 1e-7 µm
          u <- sf::st_union(sf::st_sfc(geoms[[j]]), sf::st_snap(sf::st_sfc(og), sf::st_sfc(geoms[[j]]), 1e-7))[[1]]
          n_snap <- n_snap + 1L
        }
        geoms[[j]] <- u
        gained[j] <- gained[j] + 1L; n_re <- n_re + 1L
        orphans[[o]]$done <- TRUE; progress <- TRUE
      }
      pending <- Filter(function(o) is.null(orphans[[o]]$done), pending)
      if (!length(pending) || !progress) break
    }
    # orfani senza territorio confinante (solo regioni MULTIPOLYGON): restano al generatore
    for (o in pending) {
      i <- orphans[[o]]$from
      geoms[[i]] <- sf::st_union(sf::st_sfc(geoms[[i]]), sf::st_sfc(orphans[[o]]$g))[[1]]
      lost[i] <- lost[i] - 1L
    }
  }
  list(geoms = geoms, lost = lost, gained = gained, n_fragments = n_frag, n_reassigned = n_re, n_snapped = n_snap)
}

#' Chaikin (taglio 1/4–3/4) su un anello chiuso.
#' @keywords internal
.tv_chaikin_ring <- function(m, n_iter) {
  for (it in seq_len(n_iter)) {
    p <- m[-nrow(m), , drop = FALSE]; q <- rbind(p[-1, , drop = FALSE], p[1, ])
    a <- 0.75 * p + 0.25 * q; b <- 0.25 * p + 0.75 * q
    m <- matrix(t(cbind(a, b)), ncol = 2, byrow = TRUE)
    m <- rbind(m, m[1, ])
  }
  m
}

#' Smussatura (D-S1.3.3): Chaikin su ogni anello ∩ territorio non smussato; se
#' il risultato è in più pezzi resta quello con il generatore (gli altri → lacuna).
#' @keywords internal
.tv_smooth <- function(geoms, x, y, n_iter) {
  n_drop <- 0L
  out <- lapply(seq_along(geoms), function(i) {
    g <- geoms[[i]]
    sm <- if (inherits(g, "POLYGON")) sf::st_polygon(lapply(unclass(g), .tv_chaikin_ring, n_iter = n_iter)) else
      sf::st_multipolygon(lapply(unclass(g), function(p) lapply(p, .tv_chaikin_ring, n_iter = n_iter)))
    s <- sf::st_sfc(sm)
    if (!isTRUE(sf::st_is_valid(s))) s <- sf::st_make_valid(s)
    r <- .tv_polys_only(suppressWarnings(sf::st_intersection(s, sf::st_sfc(g))))
    if (inherits(r, "MULTIPOLYGON") && length(unclass(r)) > 1) {
      pcs <- .tv_pieces(r)
      own <- which(lengths(sf::st_intersects(sf::st_sfc(pcs), sf::st_sfc(sf::st_point(c(x[i], y[i]))))) > 0)
      if (!length(own)) own <- which.max(vapply(pcs, .tv_area, 0))
      n_drop <<- n_drop + length(pcs) - 1L
      r <- pcs[[own[1]]]
    }
    r
  })
  list(geoms = out, n_dropped = n_drop)
}

#' Geometria finale: POLYGON dove possibile, validità controllata (riparazione
#' solo come ultima risorsa, contata).
#' @keywords internal
.tv_finalise <- function(geoms) {
  geoms <- lapply(geoms, function(g) if (inherits(g, "MULTIPOLYGON") && length(unclass(g)) == 1) sf::st_polygon(unclass(g)[[1]]) else g)
  s <- sf::st_sfc(geoms)
  ok <- sf::st_is_valid(s)
  n_rep <- sum(!ok, na.rm = TRUE) + sum(is.na(ok))
  if (n_rep) s[!ok | is.na(ok)] <- sf::st_make_valid(s[!ok | is.na(ok)])
  n_multi <- sum(as.character(sf::st_geometry_type(s)) != "POLYGON")
  list(sfc = s, n_multipart = n_multi, n_repaired = n_rep)
}
