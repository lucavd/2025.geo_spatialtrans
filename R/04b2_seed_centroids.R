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
#' @param type_order `"sequential"` (predefinito: un compito RSA per tipo, in ordine di densita' decrescente) o `"blind"` (un compito
#'   per regione, tipi assegnati dopo a caso con i conteggi ripartiti; richiede la stessa distanza minima per tutti i tipi; S1.5b).
#' @param spatial `"uniform"` (predefinito: proposte uniformi nella regione, comportamento di S1.2-S1.5b invariato) o
#'   `"compartments"` (S1.5c: campi latenti per tipo, tipo locale e densita' locale dal campo; definizione in docs/steps/S1.5c.md).
#' @param compartments con `spatial = "compartments"`: data.frame `cell_type`, `ell_um` (lunghezza di correlazione del campo, µm),
#'   `s` (ampiezza sul logaritmo delle probabilita'); i tipi assenti o con s = 0 non hanno campo.
#' @param labels con `"compartments"`: `"field"` (tipo estratto da p(x) al momento della proposta) o `"blind"` (stesse posizioni,
#'   etichette permutate dentro la regione; richiede la stessa distanza minima per tutti i tipi; controprova di S1.5c).
#' @param n_features numero di caratteristiche di Fourier casuali per campo (256).
#' @param heterogeneity con `"compartments"` (S1.5d): NULL (predefinito, comportamento di S1.5c bit per bit) o data.frame di una riga
#'   `void_frac` (frazione attesa dell'area senza cellule), `void_ell_um`, `dens_sigma` (ampiezza sul logaritmo della densita'),
#'   `dens_ell_um`: campo dei vuoti V e campo della densita' W, indipendenti dal tipo; definizione in docs/steps/S1.5d.md.
#'   Colonne opzionali (S1.5e, docs/steps/S1.5e.md; assenti = S1.5d bit per bit): `void_mode` (`"move"`: il numero di cellule della
#'   regione non cambia e le cellule si addensano fuori dai vuoti; `"remove"`: numero di cellule moltiplicato per 1 - v, v = frazione della
#'   griglia della regione nei vuoti, cioe' la densita' vale nella parte non vuota), `void_aspect` (kappa >= 1: campo dei vuoti anisotropo,
#'   lunghezze ell sqrt(kappa) e ell / sqrt(kappa), orientamento uniforme per regione), `dens_cap` (TRUE: fattore della densita' riscalato a
#'   media 1 fuori dai vuoti e troncato a lambda_cap / lambda medio, lambda_cap = 0.8 * 0.547 / (pi dbar^2 / 4), dbar = media pesata delle
#'   distanze minime dei tipi della regione; il tetto vale sul fattore della densita', non sulla densita' locale, che con i compartimenti dei
#'   tipi puo' superare lambda_cap di qualche per cento, revisione RC-4 di S1.5e).
#' @param quotas con `"compartments"`: FALSE (predefinito da S1.6: nessun rifiuto per quota, regola Q di S1.5e, DEC-025) o TRUE (quote per tipo
#'   di S1.5c: la proposta il cui tipo ha esaurito la quota e' rifiutata; comportamento di S1.5c-S1.5e). Predefinito di sessione: `options(geo.spatialtrans.quotas = )`.
#' @param thinning_bound con `"compartments"`: limite del diradamento delle proposte. `"exact"` (predefinito da S1.6, BL-161): per regione
#'   max_t rho_t dei tipi della regione (w = 1 / sum_t p_t / rho_t <= max_t rho_t), moltiplicato per il massimo del fattore della densita'
#'   (1 senza campo della densita'; il tetto H con `dens_cap`); con un campo della densita' senza tetto non c'e' un limite finito utile e si usa
#'   `"grid"` (registrato in `info$compartments$regions[[k]]$bound`). `"grid"`: 1.1 volte il massimo di w sulla griglia della regione
#'   (euristica di S1.5c-S1.5e, per riprodurre i banchi passati). Predefinito di sessione: `options(geo.spatialtrans.thinning_bound = )`.
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
                           type_order = c("sequential", "blind"),
                           spatial = c("uniform", "compartments"),
                           compartments = NULL,
                           labels = c("field", "blind"),
                           n_features = 256L,
                           heterogeneity = NULL,
                           quotas = getOption("geo.spatialtrans.quotas", FALSE),
                           thinning_bound = getOption("geo.spatialtrans.thinning_bound", "exact"),
                           verbose = TRUE) {
  type_order <- match.arg(type_order); spatial <- match.arg(spatial); labels <- match.arg(labels); thinning_bound <- match.arg(thinning_bound, c("exact", "grid"))
  t0 <- proc.time()[["elapsed"]]
  a <- .sc_check_args(regions, cell_types, region_composition, random_seed,
                      dmin_factor, max_attempts_factor)
  type_df <- .sc_type_table(a$cell_types, a$composition, dmin_factor)
  rdf <- regions$region_df
  res <- .sc_with_seed(random_seed, {
    plan <- .sc_plan(rdf, type_df, a$composition)
    if (spatial == "compartments") {
      pl <- .sc_place_compart(plan, regions$region_polygons, rdf, type_df, a$composition,
                              .sc_check_compart(compartments, type_df), labels, as.integer(n_features),
                              max_attempts_factor, max_consecutive_rejections,
                              het = .sc_check_het(heterogeneity), quotas = isTRUE(quotas), bound = thinning_bound)
    } else if (type_order == "blind" && nrow(plan$tasks)) {
      bl <- .sc_blind_tasks(plan$tasks, type_df)
      pl <- .sc_place_all(bl$tasks, regions$region_polygons, rdf, type_df,
                          max_attempts_factor, max_consecutive_rejections)
      pl$type_idx <- .sc_blind_labels(pl, plan$tasks)
      plan$tasks <- bl$tasks
    } else {
      pl <- .sc_place_all(plan$tasks, regions$region_polygons, rdf, type_df,
                          max_attempts_factor, max_consecutive_rejections)
    }
    list(plan = plan, pl = pl)
  })
  plan <- res$plan; pl <- res$pl

  n_cells  <- tabulate(match(pl$region_id, rdf$region_id), nbins = nrow(rdf))
  n_failed <- if (spatial == "compartments") pl$region_failed else
    tabulate(match(rep(plan$tasks$region_id, pl$task_failed), rdf$region_id), nbins = nrow(rdf))
  out_rdf <- rdf
  out_rdf$target_density_weighted <- plan$rho_mix
  out_rdf$lambda           <- plan$lambda
  out_rdf$n_target         <- plan$n_region
  out_rdf$n_cells          <- n_cells
  out_rdf$n_failed         <- n_failed
  out_rdf$achieved_density <- if (nrow(rdf)) n_cells / rdf$area_um2 * 1e6 else numeric(0)
  if (spatial == "compartments" && any(pl$n_removed > 0L)) {          # S1.5e, void_mode = "remove": bersaglio ridotto della parte nei vuoti
    out_rdf$n_target <- plan$n_region - pl$n_removed
    out_rdf$target_density_weighted <- plan$rho_mix * (1 - pl$void_v)
    out_rdf$void_removed_frac <- pl$void_v; out_rdf$n_removed_void <- pl$n_removed
  }

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
    random_seed = random_seed, dmin_factor = dmin_factor, type_order = type_order, spatial = spatial,
    labels = if (spatial == "compartments") labels else NA_character_, compartments = pl$compart_info,
    max_attempts_factor = max_attempts_factor,
    n_cells = nrow(centroids), n_target = sum(out_rdf$n_target),
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

#' Semina cieca al tipo (S1.5b, revisione RC-01): un solo compito RSA per regione con tutte le cellule della regione; i tipi si
#' assegnano dopo, permutando a caso le etichette ripartite (stessi conteggi per tipo). Richiede la stessa distanza minima per
#' tutti i tipi della regione, altrimenti l'ordine di semina conta (BL-055) e si usa `type_order = "sequential"`.
#' @keywords internal
.sc_blind_tasks <- function(tasks, type_df) {
  d <- type_df$min_dist_um[tasks$type_idx]
  rng <- tapply(d, tasks$region_idx, function(v) max(v) - min(v))
  if (any(rng > 0)) stop("seed_centroids(): type_order = \"blind\" richiede la stessa distanza minima per tutti i tipi di una regione.", call. = FALSE)
  ri <- sort(unique(tasks$region_idx))
  out <- data.frame(region_idx = ri, region_id = tasks$region_id[match(ri, tasks$region_idx)],
                    type_idx = tasks$type_idx[match(ri, tasks$region_idx)],
                    n = as.integer(tapply(tasks$n, tasks$region_idx, sum)[as.character(ri)]))
  list(tasks = out)
}

#' Etichette dei tipi per la semina cieca: per regione, permutazione casuale delle etichette ripartite (in ordine di regione).
#' Se la semina ha piazzato meno cellule del previsto, le etichette sono un campione senza reinserimento di quelle ripartite.
#' @keywords internal
.sc_blind_labels <- function(pl, tasks) {
  ti <- pl$type_idx
  for (rid in unique(tasks$region_id)) {
    k <- which(pl$region_id == rid); if (!length(k)) next
    tk <- tasks[tasks$region_id == rid, , drop = FALSE]
    labs <- rep(tk$type_idx, tk$n)
    ti[k] <- if (length(labs) == 1L) labs else sample(labs, length(k))
  }
  ti
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


# ---- semina a compartimenti (S1.5c) -----------------------------------------
# Definizione: docs/steps/S1.5c.md, «Cosa si implementa». Per regione: un campo gaussiano latente per tipo con s > 0
# (covarianza exp(-r^2 / (2 l^2)), caratteristiche di Fourier casuali), p_t(x) = softmax(log f_t + a_t + s_t Z_t(x)),
# densita' locale lambda(x) proporzionale a 1 / sum_t p_t(x) / rho_t; a_t sulla griglia della regione in modo che le frazioni
# attese per numero (pesate con lambda) siano f_t; numero di cellule della regione = quello del piano (densita' armonica).
# Proposte uniformi, diradate con lambda / lambda_max, tipo estratto da p(x), poi la distanza minima della coppia (RSA).
# Quote per tipo (deviazione D-1 dalla pre-registrazione, docs/steps/S1.5c.md): la proposta il cui tipo ha esaurito la quota e' rifiutata,
# cosi' i conteggi per tipo sono quelli ripartiti dal piano (senza quote la distanza minima sfavoriva i tipi densi: test T6).
# Da S1.6 le quote sono spente per predefinito (regola Q di S1.5e, DEC-025) e il diradamento usa il limite esatto max_t rho_t (BL-161).

#' Tabella dei compartimenti validata e allineata ai tipi usati (ell, s; s = 0 se assente).
#' @keywords internal
.sc_check_compart <- function(compartments, type_df) {
  out <- data.frame(cell_type = type_df$cell_type, ell_um = NA_real_, s = 0, stringsAsFactors = FALSE)
  if (is.null(compartments)) return(out)
  if (!is.data.frame(compartments) || !all(c("cell_type", "ell_um", "s") %in% names(compartments)))
    stop("seed_centroids(): `compartments` deve avere colonne cell_type, ell_um, s.", call. = FALSE)
  if (anyDuplicated(compartments$cell_type)) stop("seed_centroids(): tipi duplicati in `compartments`.", call. = FALSE)
  k <- match(out$cell_type, as.character(compartments$cell_type)); h <- !is.na(k)
  out$ell_um[h] <- as.numeric(compartments$ell_um[k[h]]); out$s[h] <- as.numeric(compartments$s[k[h]])
  bad <- out$s < 0 | !is.finite(out$s) | (out$s > 0 & !(is.finite(out$ell_um) & out$ell_um > 0))
  if (any(bad)) stop("seed_centroids(): `compartments` con s < 0 o ell_um non positiva per s > 0: ", paste(out$cell_type[bad], collapse = ", "), call. = FALSE)
  out
}

#' Parametri di un campo gaussiano a caratteristiche di Fourier casuali (omega ~ N(0, I / ell^2), b ~ U(0, 2 pi)).
#' @keywords internal
.sc_rff_draw <- function(ell, M) list(w = matrix(stats::rnorm(2L * M, 0, 1 / ell), ncol = 2L), b = stats::runif(M, 0, 2 * pi), M = M)

#' Valore del campo nei punti (x, y): sqrt(2 / M) * sum_k cos(w_k . x + b_k); media 0, varianza 1.
#' @keywords internal
.sc_rff_eval <- function(f, x, y) {
  if (!length(x)) return(numeric(0))
  ang <- outer(x, f$w[, 1]) + outer(y, f$w[, 2]) + rep(f$b, each = length(x))
  sqrt(2 / f$M) * rowSums(cos(ang))
}

#' Logit dei tipi nei punti: log f + s * Z (senza a); matrice n x K.
#' @keywords internal
.sc_eta <- function(fields, s, logf, x, y) {
  eta <- matrix(logf, nrow = length(x), ncol = length(logf), byrow = TRUE)
  for (t in seq_along(fields)) if (!is.null(fields[[t]])) eta[, t] <- eta[, t] + s[t] * .sc_rff_eval(fields[[t]], x, y)
  eta
}

#' Probabilita' locali (softmax) e densita' relativa w = 1 / sum p / rho dati i logit e lo spostamento a.
#' @keywords internal
.sc_prob <- function(eta, a, rho) {
  e <- sweep(eta, 2L, a, "+"); e <- e - apply(e, 1L, max); p <- exp(e); p <- p / rowSums(p)
  list(p = p, w = 1 / as.vector(p %*% (1 / rho)))
}

#' Spostamenti a_t: frazioni attese per numero (pesate con w) = f sui punti della griglia.
#' @keywords internal
.sc_solve_a <- function(eta, f, rho, iters = 500L, tol = 1e-12, h = NULL) {
  a <- numeric(length(f))
  for (it in seq_len(iters)) {
    pr <- .sc_prob(eta, a, rho)
    fr <- if (is.null(h)) colSums(pr$w * pr$p) / sum(pr$w) else colSums(h * pr$w * pr$p) / sum(h * pr$w)   # S1.5d: pesi dell'eterogeneita'
    step <- log(f) - log(fr); a <- a + step
    if (max(abs(step)) < tol) break
  }
  a - a[1]
}

#' Griglia di punti della regione per a_t e lambda_max: passo min(ell_min / 2, sqrt(area / 4000)), almeno 50 punti.
#' @keywords internal
.sc_grid_points <- function(poly, area, ell_min) {
  bb <- sf::st_bbox(poly)
  h <- min(if (is.finite(ell_min)) ell_min / 2 else Inf, sqrt(max(area, 1e-12) / 4000))
  for (k in 1:6) {
    gx <- seq(bb[["xmin"]] + h / 2, bb[["xmax"]], by = h); gy <- seq(bb[["ymin"]] + h / 2, bb[["ymax"]], by = h)
    g <- expand.grid(x = gx, y = gy); ins <- .sc_in_polygon(poly, g$x, g$y)
    if (sum(ins) >= 50L) return(g[ins, , drop = FALSE])
    h <- h / 2
  }
  if (any(ins)) return(g[ins, , drop = FALSE])
  xy <- sf::st_coordinates(sf::st_point_on_surface(poly))   # regione troppo sottile per la griglia (revisione RC-09): un punto interno
  data.frame(x = xy[1, 1], y = xy[1, 2])
}

#' Indice del tipo estratto da ogni riga di P (probabilita') con le uniformi u.
#' @keywords internal
.sc_draw_type <- function(P, u) {
  K <- ncol(P); if (K == 1L) return(rep(1L, nrow(P)))
  cp <- P[, 1]; out <- rep(1L, nrow(P))
  for (t in 2:K) { out <- out + (u >= cp); cp <- cp + P[, t] }
  pmin(out, K)
}

#' Tabella dell'eterogeneita' (S1.5d) validata; NULL se assente o con vuoti e densita' spenti (comportamento di S1.5c).
#' @keywords internal
.sc_check_het <- function(het) {
  if (is.null(het)) return(NULL)
  need <- c("void_frac", "void_ell_um", "dens_sigma", "dens_ell_um")
  if (!is.data.frame(het) || nrow(het) != 1L || !all(need %in% names(het)))
    stop("seed_centroids(): `heterogeneity` deve essere un data.frame di una riga con ", paste(need, collapse = ", "), ".", call. = FALSE)
  h <- lapply(het[need], as.numeric)
  if (!is.finite(h$void_frac) || h$void_frac < 0 || h$void_frac >= 1) stop("seed_centroids(): void_frac in [0, 1).", call. = FALSE)
  if (!is.finite(h$dens_sigma) || h$dens_sigma < 0) stop("seed_centroids(): dens_sigma >= 0.", call. = FALSE)
  if (h$void_frac > 0 && !(is.finite(h$void_ell_um) && h$void_ell_um > 0)) stop("seed_centroids(): void_ell_um > 0 con void_frac > 0.", call. = FALSE)
  if (h$dens_sigma > 0 && !(is.finite(h$dens_ell_um) && h$dens_ell_um > 0)) stop("seed_centroids(): dens_ell_um > 0 con dens_sigma > 0.", call. = FALSE)
  if (h$void_frac == 0 && h$dens_sigma == 0) return(NULL)
  h$tau <- if (h$void_frac > 0) stats::qnorm(1 - h$void_frac) else Inf
  # S1.5e: colonne opzionali; assenti o NA = comportamento di S1.5d
  h$void_mode <- if ("void_mode" %in% names(het) && !is.na(het$void_mode[1])) as.character(het$void_mode[1]) else "move"
  if (!h$void_mode %in% c("move", "remove")) stop("seed_centroids(): void_mode \"move\" o \"remove\".", call. = FALSE)
  h$void_aspect <- if ("void_aspect" %in% names(het) && !is.na(het$void_aspect[1])) as.numeric(het$void_aspect[1]) else 1
  if (!is.finite(h$void_aspect) || h$void_aspect < 1) stop("seed_centroids(): void_aspect >= 1.", call. = FALSE)
  h$dens_cap <- if ("dens_cap" %in% names(het) && !is.na(het$dens_cap[1])) isTRUE(as.logical(het$dens_cap[1])) else FALSE
  h
}

#' Campo anisotropo (S1.5e): omega = R(theta) diag(1 / (ell sqrt(k)), sqrt(k) / ell) z, z ~ N(0, I); stesse estrazioni di .sc_rff_draw
#' (z, poi b), poi theta ~ U(0, pi). Con k = 1 non si usa (si chiama .sc_rff_draw: nessuna estrazione in piu').
#' @keywords internal
.sc_rff_draw_aniso <- function(ell, k, M) {
  z <- matrix(stats::rnorm(2L * M), ncol = 2L); b <- stats::runif(M, 0, 2 * pi); th <- stats::runif(1L, 0, pi)
  u <- cbind(z[, 1] / (ell * sqrt(k)), z[, 2] * sqrt(k) / ell)
  R <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2L)                                 # colonne: asse lungo, asse corto
  list(w = u %*% t(R), b = b, M = M, theta = th, kappa = k)
}

#' Scala del fattore della densita' con il tetto (S1.5e): s tale che mean(min(s * h0, H)) = 1 sui punti non vuoti (bisezione).
#' H <= 1: tetto sotto la densita' media, fattore uniforme (s = 0, restituisce NA e il chiamante usa h = 1).
#' @keywords internal
.sc_cap_scale <- function(h0, H) {
  if (!length(h0) || !is.finite(H)) return(1)
  if (H <= 1) return(NA_real_)
  lo <- 0; hi <- 1 / max(min(h0), 1e-300)
  for (it in 1:200) { md <- (lo + hi) / 2; if (mean(pmin(md * h0, H)) < 1) lo <- md else hi <- md }
  (lo + hi) / 2
}

#' Peso dell'eterogeneita' nei punti: exp(dens_sigma * W) * 1[V <= tau] (S1.5d); con il tetto (S1.5e) min(s * exp(dens_sigma * W), H).
#' @keywords internal
.sc_het_weight <- function(het, fv, fw, x, y, cap = NULL) {
  h <- rep(1, length(x))
  if (!is.null(fw)) h <- exp(het$dens_sigma * .sc_rff_eval(fw, x, y))
  if (!is.null(fw) && !is.null(cap)) h <- if (is.na(cap$s)) rep(1, length(x)) else pmin(cap$s * h, cap$H)
  if (!is.null(fv)) h[.sc_rff_eval(fv, x, y) > het$tau] <- 0
  h
}

#' RSA a compartimenti su tutte le regioni (distanza minima globale come in .sc_place_all).
#' @keywords internal
.sc_place_compart <- function(plan, polys, rdf, type_df, comp, ctab, labels, M, max_attempts_factor, max_consecutive_rejections,
                              het = NULL, quotas = FALSE, bound = "exact") {
  nr <- nrow(rdf); ntot <- sum(plan$n_region)
  X <- numeric(ntot); Y <- numeric(ntot); D <- numeric(ntot); TI <- integer(ntot); RI <- integer(ntot); RX <- integer(ntot)
  region_failed <- integer(nr); attempts <- 0L; exceed <- 0L; info <- vector("list", nr); n_quota_rej <- 0L; n_thin_acc <- 0L
  n_void_disabled <- 0L; n_removed <- integer(nr); void_v <- numeric(nr)
  n_drawn <- integer(nrow(type_df))                     # tipi estratti da p(x) alle proposte accettate dal diradamento (prima di quote e RSA)
  dtype <- type_df$min_dist_um
  used <- sort(unique(match(comp$cell_type, type_df$cell_type)))
  h <- if (length(used)) max(dtype[used]) else 0; hardcore <- h > 0; k <- 0L; grid_cell <- NA_real_
  act <- which(plan$n_region > 0L)
  if (hardcore && length(act)) {
    bb_all <- sf::st_bbox(polys[act]); gx0 <- bb_all[["xmin"]]; gy0 <- bb_all[["ymin"]]
    ext <- max(bb_all[["xmax"]] - gx0, bb_all[["ymax"]] - gy0)
    grid_cell <- max(h, ext / 4000)
    nx <- as.integer(floor((bb_all[["xmax"]] - gx0) / grid_cell)) + 1L; ny <- as.integer(floor((bb_all[["ymax"]] - gy0) / grid_cell)) + 1L
    nxp <- nx + 2L; nyp <- ny + 2L
    cnt <- integer(nxp * nyp); cellpts <- matrix(0L, nrow = nxp * nyp, ncol = 8L); offs <- .sc_neighbour_offsets(nxp)
  }
  cl <- as.character(rdf$cluster_id); comp_by <- split(comp, comp$cluster_id)
  for (i in act) {
    n <- plan$n_region[i]; d <- comp_by[[cl[i]]]
    ti_r <- match(d$cell_type, type_df$cell_type); f <- d$fraction / sum(d$fraction); K <- length(ti_r)
    rho <- type_df$density[ti_r]; s <- ctab$s[match(type_df$cell_type[ti_r], ctab$cell_type)]
    ell <- ctab$ell_um[match(type_df$cell_type[ti_r], ctab$cell_type)]
    tk <- plan$tasks[plan$tasks$region_idx == i, , drop = FALSE]                 # quote per tipo (ripartizione del piano, come S1.5b)
    quota <- integer(K); quota[match(tk$type_idx, ti_r)] <- tk$n
    fields <- vector("list", K)
    for (t in seq_len(K)) if (s[t] > 0) fields[[t]] <- .sc_rff_draw(ell[t], M)
    fv <- fw <- NULL                                                              # S1.5d: campi dei vuoti e della densita', estratti solo se attivi
    if (!is.null(het) && het$void_frac > 0) fv <- if (het$void_aspect > 1) .sc_rff_draw_aniso(het$void_ell_um, het$void_aspect, M) else .sc_rff_draw(het$void_ell_um, M)
    if (!is.null(het) && het$dens_sigma > 0) fw <- .sc_rff_draw(het$dens_ell_um, M)
    poly <- polys[i]; bb <- sf::st_bbox(poly)
    ell_g <- if (any(s > 0)) min(ell[s > 0]) else Inf
    if (!is.null(fw)) ell_g <- min(ell_g, het$dens_ell_um)
    gp <- .sc_grid_points(poly, rdf$area_um2[i], ell_g)
    eta_g <- .sc_eta(fields, s, log(f), gp$x, gp$y)
    hmax <- 1                                                                      # massimo del fattore della densita' (S1.6, BL-161)
    if (is.null(het)) {
      hg <- NULL; a <- .sc_solve_a(eta_g, f, rho)
      wg <- .sc_prob(eta_g, a, rho)$w
    } else {
      hg <- .sc_het_weight(het, fv, fw, gp$x, gp$y)
      if (!any(hg > 0)) { fv <- NULL; n_void_disabled <- n_void_disabled + 1L; hg <- .sc_het_weight(het, NULL, fw, gp$x, gp$y) }   # regione tutta nei vuoti: senza vuoti
      k_in <- hg > 0
      if (het$void_mode == "remove" && !is.null(fv)) {                   # S1.5e: le cellule della parte nei vuoti si tolgono (per tipo, valore atteso esatto)
        v <- mean(!k_in); void_v[i] <- v
        qn <- quota * (1 - v); q2 <- floor(qn) + as.integer(stats::runif(K) < qn - floor(qn))
        n_removed[i] <- sum(quota) - sum(q2); n <- sum(q2)
        quota <- q2
      }
      cap <- NULL
      if (isTRUE(het$dens_cap) && !is.null(fw)) {                         # S1.5e: tetto della densita' locale (nota D-1 di docs/steps/S1.5e.md)
        dbar <- sum(f * dtype[ti_r]); lam_cap <- 0.8 * 0.547 / (pi * dbar^2 / 4)
        lam_bar <- n / (rdf$area_um2[i] * mean(k_in))
        h0 <- exp(het$dens_sigma * .sc_rff_eval(fw, gp$x[k_in], gp$y[k_in]))
        cap <- list(H = lam_cap / lam_bar); cap$s <- .sc_cap_scale(h0, cap$H)
        hg <- .sc_het_weight(het, fv, fw, gp$x, gp$y, cap = cap); k_in <- hg > 0
      }
      a <- .sc_solve_a(eta_g[k_in, , drop = FALSE], f, rho, h = hg[k_in])
      wg <- .sc_prob(eta_g, a, rho)$w * hg
      if (!is.null(fw)) hmax <- if (is.null(cap)) NA_real_ else if (is.na(cap$s)) 1 else cap$H   # senza tetto: nessun limite finito utile
    }
    b_used <- if (bound == "exact" && is.finite(hmax)) "exact" else "grid"
    wmax <- if (b_used == "exact") max(rho) * hmax * (1 + 1e-12) else 1.1 * max(wg)    # 1e-12: arrotondamento di 1 / sum p / rho
    info[[i]] <- list(region_id = rdf$region_id[i], a = a, wmax = wmax, w_mean = mean(wg), n_grid = nrow(gp), n_fields = sum(s > 0),
                      void_grid_frac = if (is.null(hg)) 0 else mean(hg == 0), field_void = fv, field_dens = fw,   # D-5: campi per i test
                      n_removed = n_removed[i], cap = if (is.null(het)) NULL else cap, n_plan = sum(tk$n),     # S1.5e
                      field_types = fields, ell_grid = ell_g, w_grid_max = max(wg), wmax_factor = wmax / max(wg),   # revisione: M3, M5, N6
                      bound = b_used)                                                                                # S1.6
    bw <- bb[["xmax"]] - bb[["xmin"]]; bh <- bb[["ymax"]] - bb[["ymin"]]
    p_in <- min(1, max(rdf$area_um2[i] / max(bw * bh, 1e-12), 1e-3)); p_thin <- max(mean(wg) / wmax, 1e-3)
    max_att <- as.integer(min(max_attempts_factor * n + 1000, .Machine$integer.max - 1))
    placed <- 0L; att <- 0L; empty_blocks <- 0L; rej_run <- 0L; dr <- dtype[ti_r]
    while (placed < n && att < max_att && empty_blocks < 100L && rej_run < max_consecutive_rejections) {
      need <- n - placed
      block <- as.integer(min(max(ceiling(2 * need / (p_in * p_thin)), 64), 1e5))
      cx <- stats::runif(block, bb[["xmin"]], bb[["xmax"]]); cy <- stats::runif(block, bb[["ymin"]], bb[["ymax"]])
      ins <- .sc_in_polygon(poly, cx, cy); cx <- cx[ins]; cy <- cy[ins]
      empty_blocks <- if (length(cx)) 0L else empty_blocks + 1L
      if (!length(cx)) next
      pr <- .sc_prob(.sc_eta(fields, s, log(f), cx, cy), a, rho)
      if (!is.null(het)) pr$w <- pr$w * .sc_het_weight(het, fv, fw, cx, cy, cap = cap)
      u_thin <- stats::runif(length(cx)); u_type <- stats::runif(length(cx))
      exceed <- exceed + sum(pr$w > wmax)
      keep <- u_thin < pr$w / wmax
      att <- att + length(cx)
      if (!any(keep)) next
      cx <- cx[keep]; cy <- cy[keep]
      tl <- .sc_draw_type(pr$p[keep, , drop = FALSE], u_type[keep])
      n_drawn <- n_drawn + tabulate(ti_r[tl], nbins = length(n_drawn))
      for (j in seq_along(cx)) {
        x <- cx[j]; y <- cy[j]; tj <- tl[j]; dt <- dr[tj]; ok <- TRUE
        n_thin_acc <- n_thin_acc + 1L
        if (quotas && quota[tj] == 0L) { n_quota_rej <- n_quota_rej + 1L; next }   # quota del tipo esaurita: proposta rifiutata (quotas = FALSE: S1.5d, braccio Q)
        if (hardcore) {
          cid <- (as.integer((x - gx0) %/% grid_cell) + 2L) + (as.integer((y - gy0) %/% grid_cell) + 1L) * nxp
          nb <- cid + offs
          if (any(cnt[nb] > 0L)) {
            ids <- cellpts[nb, , drop = FALSE]; ids <- ids[ids > 0L]
            dx <- X[ids] - x; dy <- Y[ids] - y; thr <- .sc_pair_dist(dt, D[ids])
            if (any(dx * dx + dy * dy < thr * thr)) ok <- FALSE
          }
        }
        if (!ok) { rej_run <- rej_run + 1L; if (rej_run >= max_consecutive_rejections) break; next }
        rej_run <- 0L; k <- k + 1L; quota[tj] <- quota[tj] - 1L
        X[k] <- x; Y[k] <- y; D[k] <- dt; TI[k] <- ti_r[tj]; RI[k] <- rdf$region_id[i]; RX[k] <- i
        if (hardcore) {
          c1 <- cnt[cid] + 1L
          if (c1 > ncol(cellpts)) cellpts <- cbind(cellpts, matrix(0L, nrow(cellpts), ncol(cellpts)))
          cellpts[cid, c1] <- k; cnt[cid] <- c1
        }
        placed <- placed + 1L
        if (placed == n) break
      }
    }
    region_failed[i] <- n - placed; attempts <- attempts + att
  }
  keep <- seq_len(k); TI <- TI[keep]; RX <- RX[keep]
  if (labels == "blind") {
    if (length(unique(dtype[used])) > 1L) stop("seed_centroids(): labels = \"blind\" richiede la stessa distanza minima per tutti i tipi.", call. = FALSE)
    for (i in unique(RX)) { w <- which(RX == i); if (length(w) > 1L) TI[w] <- TI[w][sample.int(length(w))] }
  }
  list(x = X[keep], y = Y[keep], type_idx = TI, region_id = RI[keep], region_failed = region_failed, n_removed = n_removed, void_v = void_v,
       task_attempts = attempts, grid_cell_um = grid_cell, hardcore = hardcore,
       compart_info = list(regions = Filter(Negate(is.null), info), n_exceed_wmax = exceed, n_quota_reject = n_quota_rej,
                           n_thin_accept = n_thin_acc, n_drawn = stats::setNames(n_drawn, type_df$cell_type), table = ctab,
                           heterogeneity = het, quotas = quotas, n_void_disabled = n_void_disabled, thinning_bound = bound))
}
