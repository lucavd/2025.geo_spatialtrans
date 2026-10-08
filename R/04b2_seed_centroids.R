#' Posiziona i centroidi cellulari nelle regioni (Step 1, S1.2)
#'
#' Per ogni regione di `extract_regions()`: densità target della miscela del
#' cluster (media **armonica** 1 / Σ f_i/ρ_i, decisione D3 del 2026-10-08),
#' numero di cellule con arrotondamento stocastico non distorto, ripartizione
#' per tipo con campionamento sistematico di Madow sui resti, piazzamento con
#' **RSA a n fissato** (adsorbimento sequenziale casuale: proposte uniformi nel
#' poligono, accettate se nessun centroide già posto è più vicino di
#' d_ij = (d_i + d_j)/2). La distanza minima è **globale** (vale anche fra
#' regioni confinanti). I tipi sono posti in ordine di densità decrescente.
#' Design: `docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md` §5.2,
#' con le deviazioni approvate in `results/S1.2/S1.2_preregistration.md`.
#'
#' @param regions output di `extract_regions()`: lista con `region_df`
#'   (`region_id`, `cluster_id`, `area_um2`) e `region_polygons` (sfc, µm).
#' @param cell_types data.frame con `cell_type`, `density` (cellule/mm²) e,
#'   opzionale, `min_dist_um` (NA → `dmin_factor × eq_radius_target`,
#'   eq_radius_target = sqrt(1e6 / (π density))).
#' @param region_composition data.frame con `cluster_id`, `cell_type`,
#'   `fraction` (frazioni per numero di cellule; normalizzate con warning se
#'   non sommano a 1 per cluster).
#' @param random_seed seme; lo stato globale del generatore è ripristinato
#'   all'uscita.
#' @param dmin_factor fattore della distanza minima di default (design: 2/3).
#' @param max_attempts_factor proposte **dentro il poligono** massime per
#'   (regione, tipo) = `max_attempts_factor × n + 1000`; il compito si ferma
#'   prima anche dopo `max_consecutive_rejections` proposte interne rifiutate
#'   di fila (saturazione, RA-11) o dopo 100 blocchi senza proposte interne.
#'   Le cellule mancanti sono registrate in `n_failed` con un warning.
#'   `info$n_attempts` conta le proposte interne valutate, non quelle generate
#'   nel riquadro.
#' @param max_consecutive_rejections vedi sopra (default 10 000: con
#'   probabilita' di accettazione >= 1e-3 l'arresto spurio ha probabilita' < e^-10).
#' @param verbose stampa un riepilogo.
#' @return lista con `centroids` (data.frame `cell_id, region_id, x, y,
#'   cell_type, cluster_id`), `region_df` (input + `target_density_weighted`,
#'   `lambda`, `n_target`, `n_cells`, `n_failed`, `achieved_density`),
#'   `type_df` (per tipo: densità, `min_dist_um`, ordine), `info`.
#' @export
seed_centroids <- function(regions,
                           cell_types,
                           region_composition,
                           random_seed = 42,
                           dmin_factor = 2 / 3,
                           max_attempts_factor = 1000,
                           max_consecutive_rejections = 10000L,
                           verbose = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  a <- .sc_check_args(regions, cell_types, region_composition, random_seed,
                      dmin_factor, max_attempts_factor)
  type_df <- .sc_type_table(a$cell_types, a$composition, dmin_factor)
  rdf <- regions$region_df
  res <- .sc_with_seed(random_seed, {
    plan <- .sc_plan(rdf, type_df, a$composition)
    pl   <- .sc_place_all(plan$tasks, regions$region_polygons, rdf, type_df,
                          max_attempts_factor, max_consecutive_rejections)
    list(plan = plan, pl = pl)
  })
  plan <- res$plan; pl <- res$pl

  n_cells  <- tabulate(match(pl$region_id, rdf$region_id), nbins = nrow(rdf))
  n_failed <- tabulate(match(rep(plan$tasks$region_id, pl$task_failed), rdf$region_id),
                       nbins = nrow(rdf))
  out_rdf <- rdf
  out_rdf$target_density_weighted <- plan$rho_mix
  out_rdf$lambda           <- plan$lambda
  out_rdf$n_target         <- plan$n_region
  out_rdf$n_cells          <- n_cells
  out_rdf$n_failed         <- n_failed
  out_rdf$achieved_density <- if (nrow(rdf)) n_cells / rdf$area_um2 * 1e6 else numeric(0)

  centroids <- data.frame(
    cell_id    = seq_along(pl$x),
    region_id  = pl$region_id,
    x          = pl$x,
    y          = pl$y,
    cell_type  = factor(type_df$cell_type[pl$type_idx], levels = type_df$cell_type),
    cluster_id = rdf$cluster_id[match(pl$region_id, rdf$region_id)]
  )
  if (sum(n_failed) > 0) {
    warning(sprintf("seed_centroids(): %d cellule non piazzate in %d regioni (RSA oltre il limite di proposte: regione satura o troppo stretta per la distanza minima).",
                    sum(n_failed), sum(n_failed > 0)), call. = FALSE)
  }
  info <- list(
    random_seed = random_seed, dmin_factor = dmin_factor,
    max_attempts_factor = max_attempts_factor,
    n_cells = nrow(centroids), n_target = sum(plan$n_region),
    n_failed = sum(n_failed), n_attempts = sum(pl$task_attempts),
    grid_cell_um = pl$grid_cell_um, hardcore = pl$hardcore,
    elapsed_s = proc.time()[["elapsed"]] - t0
  )
  if (verbose) {
    cat(sprintf("[seed_centroids] %d regioni, %d cellule (target %d, non piazzate %d), %d proposte, %.1f s\n",
                nrow(rdf), info$n_cells, info$n_target, info$n_failed,
                info$n_attempts, info$elapsed_s))
  }
  list(centroids = centroids, region_df = out_rdf, type_df = type_df, info = info)
}

# ---- helper interni ---------------------------------------------------------

#' Validazione degli argomenti; restituisce catalogo e composizione puliti.
#' @keywords internal
.sc_check_args <- function(regions, cell_types, region_composition, random_seed,
                           dmin_factor, max_attempts_factor) {
  if (!is.list(regions) || is.null(regions$region_df) || is.null(regions$region_polygons)) {
    stop("seed_centroids(): `regions` deve essere l'output di extract_regions() (region_df + region_polygons).", call. = FALSE)
  }
  rdf <- regions$region_df
  need_r <- c("region_id", "cluster_id", "area_um2")
  if (!all(need_r %in% names(rdf))) {
    stop("seed_centroids(): region_df senza colonne ", paste(setdiff(need_r, names(rdf)), collapse = ", "), call. = FALSE)
  }
  if (length(regions$region_polygons) != nrow(rdf)) {
    stop("seed_centroids(): region_polygons e region_df hanno lunghezze diverse.", call. = FALSE)
  }
  if (is.factor(rdf$region_id)) stop("seed_centroids(): region_id e' un factor: passare interi o caratteri (RA-07, S1.2).", call. = FALSE)
  if (anyDuplicated(rdf$region_id)) stop("seed_centroids(): region_id duplicati.", call. = FALSE)
  if (!is.numeric(rdf$area_um2)) stop("seed_centroids(): area_um2 deve essere numerica.", call. = FALSE)
  if (any(!is.finite(rdf$area_um2)) || any(rdf$area_um2 < 0)) {
    stop("seed_centroids(): area_um2 non valida.", call. = FALSE)
  }
  if (nrow(rdf)) {
    if (any(sf::st_is_empty(regions$region_polygons))) stop("seed_centroids(): region_polygons contiene geometrie vuote.", call. = FALSE)
    a_geom <- as.numeric(sf::st_area(regions$region_polygons))
    bad_a <- abs(a_geom - rdf$area_um2) > 0.01 * pmax(a_geom, 1e-9) + 1e-6
    if (any(bad_a)) warning(sprintf("seed_centroids(): area_um2 differisce > 1%% dall'area del poligono in %d regioni (RA-16): lambda usa area_um2.", sum(bad_a)), call. = FALSE)
  }
  if (!is.data.frame(cell_types) || !all(c("cell_type", "density") %in% names(cell_types))) {
    stop("seed_centroids(): `cell_types` deve avere colonne cell_type e density.", call. = FALSE)
  }
  for (nm in intersect(c("density", "min_dist_um"), names(cell_types))) {
    if (!is.numeric(cell_types[[nm]]) && !all(is.na(cell_types[[nm]]))) stop(sprintf("seed_centroids(): cell_types$%s deve essere numerica (factor/carattere non ammessi, RA-08).", nm), call. = FALSE)
  }
  ct <- data.frame(cell_type = as.character(cell_types$cell_type),
                   density = as.numeric(cell_types$density),
                   min_dist_um = if ("min_dist_um" %in% names(cell_types)) as.numeric(cell_types$min_dist_um) else NA_real_,
                   stringsAsFactors = FALSE)
  if (anyDuplicated(ct$cell_type)) stop("seed_centroids(): cell_type duplicati nel catalogo.", call. = FALSE)
  if (any(!is.finite(ct$density)) || any(ct$density < 0)) {
    stop("seed_centroids(): density non valida (deve essere >= 0 e finita).", call. = FALSE)
  }
  if (any(!is.na(ct$min_dist_um) & (!is.finite(ct$min_dist_um) | ct$min_dist_um < 0))) {
    stop("seed_centroids(): min_dist_um non valida.", call. = FALSE)
  }
  if (!is.data.frame(region_composition) || !all(c("cluster_id", "cell_type", "fraction") %in% names(region_composition))) {
    stop("seed_centroids(): `region_composition` deve avere colonne cluster_id, cell_type, fraction.", call. = FALSE)
  }
  if (!is.numeric(region_composition$fraction)) stop("seed_centroids(): fraction deve essere numerica (RA-08).", call. = FALSE)
  comp <- data.frame(cluster_id = as.character(region_composition$cluster_id),
                     cell_type = as.character(region_composition$cell_type),
                     fraction = as.numeric(region_composition$fraction),
                     stringsAsFactors = FALSE)
  if (any(!is.finite(comp$fraction)) || any(comp$fraction < 0)) {
    stop("seed_centroids(): fraction non valida.", call. = FALSE)
  }
  unk <- setdiff(comp$cell_type, ct$cell_type)
  if (length(unk)) stop("seed_centroids(): tipi non nel catalogo: ", paste(unk, collapse = ", "), call. = FALSE)
  comp <- comp[comp$fraction > 0, , drop = FALSE]
  if (anyDuplicated(paste(comp$cluster_id, comp$cell_type, sep = "\r"))) {
    stop("seed_centroids(): coppie (cluster_id, cell_type) duplicate nella composizione.", call. = FALSE)
  }
  zero_d <- comp$cell_type %in% ct$cell_type[ct$density <= 0]
  if (any(zero_d)) {
    stop("seed_centroids(): tipi a densità 0 con fraction > 0 nella composizione: ",
         paste(unique(comp$cell_type[zero_d]), collapse = ", "), call. = FALSE)
  }
  miss <- setdiff(unique(as.character(rdf$cluster_id)), comp$cluster_id)
  if (length(miss)) stop("seed_centroids(): cluster senza composizione: ", paste(miss, collapse = ", "), call. = FALSE)
  s <- tapply(comp$fraction, comp$cluster_id, sum)
  bad <- names(s)[abs(s - 1) > 1e-9]
  if (length(bad)) {
    warning(sprintf("seed_centroids(): frazioni normalizzate a 1 nei cluster %s.", paste(bad, collapse = ", ")), call. = FALSE)
    comp$fraction <- comp$fraction / s[comp$cluster_id]
  }
  if (length(random_seed) != 1 || !is.finite(random_seed)) stop("seed_centroids(): random_seed non valido.", call. = FALSE)
  if (length(dmin_factor) != 1 || !is.finite(dmin_factor) || dmin_factor < 0) stop("seed_centroids(): dmin_factor non valido.", call. = FALSE)
  if (length(max_attempts_factor) != 1 || !is.finite(max_attempts_factor) || max_attempts_factor < 1) {
    stop("seed_centroids(): max_attempts_factor non valido.", call. = FALSE)
  }
  list(cell_types = ct, composition = comp)
}

#' Tabella dei tipi usati: distanza minima e ordine di piazzamento (densità decrescente).
#' @keywords internal
.sc_type_table <- function(ct, comp, dmin_factor) {
  used <- ct[ct$cell_type %in% comp$cell_type, , drop = FALSE]
  used <- used[order(-used$density, used$cell_type, method = "radix"), , drop = FALSE]   # radix = ordine C, indipendente dalla locale (RA-09)
  eq_r <- sqrt(1e6 / (pi * used$density))
  used$eq_radius_target_um <- eq_r
  used$min_dist_source <- ifelse(is.na(used$min_dist_um), "factor", "user")
  used$min_dist_um <- ifelse(is.na(used$min_dist_um), dmin_factor * eq_r, used$min_dist_um)
  used$order <- seq_len(nrow(used))
  rownames(used) <- NULL
  used
}

#' Densità della miscela: media armonica (D3). Frazioni per numero di cellule.
#' @keywords internal
.sc_mix_density <- function(fraction, density) 1 / sum(fraction / density)

#' Numero di cellule per regione: parte intera + Bernoulli sulla parte frazionaria.
#' @keywords internal
.sc_round_regions <- function(lambda) {
  fl <- floor(lambda)
  as.integer(fl + (stats::runif(length(lambda)) < (lambda - fl)))
}

#' Ripartizione per tipo: floor(n f) + resti con campionamento sistematico di
#' Madow (probabilità di inclusione = parti frazionarie, somma esatta n).
#' @keywords internal
.sc_allocate_types <- function(n, fraction) {
  e  <- n * fraction
  fl <- floor(e + 1e-9)
  fr <- pmax(e - fl, 0)
  r  <- n - sum(fl)
  if (r <= 0) return(as.integer(fl))
  cs <- cumsum(fr); cs <- cs * (r / cs[length(cs)])
  u  <- stats::runif(1) + 0:(r - 1)
  hit <- tabulate(findInterval(u, c(0, cs), left.open = TRUE), nbins = length(fr))
  as.integer(fl + hit)
}

#' Piano: densità, lambda, n per regione e compiti (regione, tipo, n).
#' @keywords internal
.sc_plan <- function(rdf, type_df, comp) {
  nr <- nrow(rdf)
  rho_mix <- numeric(nr)
  comp_by <- split(comp, comp$cluster_id)
  cl <- as.character(rdf$cluster_id)
  dens <- type_df$density[match(comp$cell_type, type_df$cell_type)]
  rho_by <- vapply(comp_by, function(d) .sc_mix_density(d$fraction, type_df$density[match(d$cell_type, type_df$cell_type)]), numeric(1))
  rho_mix <- unname(rho_by[cl])
  lambda <- rho_mix * rdf$area_um2 / 1e6
  n_region <- .sc_round_regions(lambda)
  tasks <- vector("list", nr)
  for (i in seq_len(nr)) {
    if (n_region[i] == 0L) next
    d <- comp_by[[cl[i]]]
    ni <- .sc_allocate_types(n_region[i], d$fraction)
    k <- ni > 0
    if (any(k)) tasks[[i]] <- data.frame(region_idx = i, region_id = rdf$region_id[i],
                                          type_idx = match(d$cell_type[k], type_df$cell_type), n = ni[k])
  }
  tasks <- do.call(rbind, tasks)
  if (is.null(tasks)) tasks <- data.frame(region_idx = integer(0), region_id = integer(0), type_idx = integer(0), n = integer(0))
  tasks <- tasks[order(tasks$type_idx, tasks$region_idx), , drop = FALSE]
  rownames(tasks) <- NULL
  list(rho_mix = rho_mix, lambda = lambda, n_region = n_region, tasks = tasks)
}

#' Esegue expr con un seme fisso e ripristina lo stato globale del generatore.
#' @keywords internal
.sc_with_seed <- function(seed, expr) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  expr
}

#' Punti (x, y) dentro il poligono (GEOS, poligono preparato; bordo incluso).
#' @keywords internal
.sc_in_polygon <- function(poly, x, y) {
  if (!length(x)) return(logical(0))
  pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"))
  hit <- sf::st_intersects(poly, sf::st_geometry(pts))[[1]]
  out <- logical(length(x)); out[hit] <- TRUE
  out
}

#' Distanza minima fra due tipi: media delle distanze minime (dischi additivi).
#' @keywords internal
.sc_pair_dist <- function(di, dj) (di + dj) * 0.5

#' Scostamenti delle 9 celle del vicinato su una griglia con bordo di 1 cella.
#' @keywords internal
.sc_neighbour_offsets <- function(nxp) {
  as.integer(c(-nxp - 1L, -nxp, -nxp + 1L, -1L, 0L, 1L, nxp - 1L, nxp, nxp + 1L))
}

#' RSA globale su tutti i compiti (tipo per densità decrescente, poi regione).
#' @keywords internal
.sc_place_all <- function(tasks, polys, rdf, type_df, max_attempts_factor,
                          max_consecutive_rejections = 10000L) {
  ntot <- sum(tasks$n)
  X <- numeric(ntot); Y <- numeric(ntot); D <- numeric(ntot)
  TI <- integer(ntot); RI <- integer(ntot)
  task_failed <- integer(nrow(tasks)); task_attempts <- integer(nrow(tasks))
  dtype <- type_df$min_dist_um
  h <- if (nrow(tasks)) max(dtype[unique(tasks$type_idx)]) else 0
  hardcore <- h > 0
  k <- 0L
  grid_cell <- NA_real_
  if (nrow(tasks)) {
    bb_all <- sf::st_bbox(polys[unique(tasks$region_idx)])
    gx0 <- bb_all[["xmin"]]; gy0 <- bb_all[["ymin"]]
    ext <- max(bb_all[["xmax"]] - gx0, bb_all[["ymax"]] - gy0)
    if (hardcore) {
      grid_cell <- max(h, ext / 4000)                 # <= ~16 M celle; cella >= h basta per il vicinato 3x3
      if (grid_cell > 4 * h) warning(sprintf("seed_centroids(): estensione %.0f µm >> distanza minima %.2f µm: griglia grossolana (cella %.1f µm), memoria e tempo crescono (RA-14).", ext, h, grid_cell), call. = FALSE)
      nx <- as.integer(floor((bb_all[["xmax"]] - gx0) / grid_cell)) + 1L
      ny <- as.integer(floor((bb_all[["ymax"]] - gy0) / grid_cell)) + 1L
      nxp <- nx + 2L; nyp <- ny + 2L
      dmin_used <- min(dtype[unique(tasks$type_idx)])
      K <- if (dmin_used > 0) min(8L, as.integer(ceiling((grid_cell / dmin_used + 1)^2))) else 8L
      K <- max(K, 1L)
      cnt <- integer(nxp * nyp)
      cellpts <- matrix(0L, nrow = nxp * nyp, ncol = K)
      offs <- .sc_neighbour_offsets(nxp)
    }
  }
  for (t in seq_len(nrow(tasks))) {
    ri <- tasks$region_idx[t]; ti <- tasks$type_idx[t]; n <- tasks$n[t]
    dt <- dtype[ti]
    poly <- polys[ri]
    bb <- sf::st_bbox(poly)
    bw <- bb[["xmax"]] - bb[["xmin"]]; bh <- bb[["ymax"]] - bb[["ymin"]]
    p_acc <- min(1, max(rdf$area_um2[ri] / max(bw * bh, 1e-12), 1e-3))
    max_att <- as.integer(min(max_attempts_factor * n + 1000, .Machine$integer.max - 1))
    placed <- 0L; att <- 0L; empty_blocks <- 0L; rej_run <- 0L
    while (placed < n && att < max_att && empty_blocks < 100L && rej_run < max_consecutive_rejections) {
      need <- n - placed
      block <- as.integer(min(max(ceiling(2 * need / p_acc), 64), 1e5))
      cx <- stats::runif(block, bb[["xmin"]], bb[["xmax"]])
      cy <- stats::runif(block, bb[["ymin"]], bb[["ymax"]])
      ins <- .sc_in_polygon(poly, cx, cy)
      cx <- cx[ins]; cy <- cy[ins]
      empty_blocks <- if (length(cx)) 0L else empty_blocks + 1L   # poligono di area ~0: evita il ciclo infinito
      if (!hardcore) {
        m <- min(length(cx), need)
        if (m > 0) {
          idx <- (k + 1L):(k + m)
          X[idx] <- cx[1:m]; Y[idx] <- cy[1:m]; D[idx] <- dt; TI[idx] <- ti; RI[idx] <- tasks$region_id[t]
          k <- k + m; placed <- placed + m
        }
        att <- att + m
        next
      }
      for (j in seq_along(cx)) {
        att <- att + 1L
        x <- cx[j]; y <- cy[j]
        cid <- (as.integer((x - gx0) %/% grid_cell) + 2L) + (as.integer((y - gy0) %/% grid_cell) + 1L) * nxp
        nb <- cid + offs
        ok <- TRUE
        if (any(cnt[nb] > 0L)) {
          ids <- cellpts[nb, , drop = FALSE]
          ids <- ids[ids > 0L]
          dx <- X[ids] - x; dy <- Y[ids] - y
          thr <- .sc_pair_dist(dt, D[ids])
          if (any(dx * dx + dy * dy < thr * thr)) ok <- FALSE
        }
        if (!ok) rej_run <- rej_run + 1L
        if (ok) {
          rej_run <- 0L
          k <- k + 1L
          X[k] <- x; Y[k] <- y; D[k] <- dt; TI[k] <- ti; RI[k] <- tasks$region_id[t]
          c1 <- cnt[cid] + 1L
          if (c1 > ncol(cellpts)) cellpts <- cbind(cellpts, matrix(0L, nrow(cellpts), ncol(cellpts)))
          cellpts[cid, c1] <- k
          cnt[cid] <- c1
          placed <- placed + 1L
          if (placed == n) break
        }
        if (att >= max_att || rej_run >= max_consecutive_rejections) break
      }
    }
    task_failed[t] <- n - placed
    task_attempts[t] <- att
  }
  keep <- seq_len(k)
  list(x = X[keep], y = Y[keep], type_idx = TI[keep], region_id = RI[keep],
       task_failed = task_failed, task_attempts = task_attempts,
       grid_cell_um = grid_cell, hardcore = hardcore)
}
