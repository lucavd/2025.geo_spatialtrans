# R/04b_cell_placement.R - Step 1, S1.5: orchestratore place_cells(), componenti piccole, figura del livello cellulare.
# Design §3, §4.1-4.7, §6.3 (v1.5); definizioni operative in docs/steps/S1.5.md.
# Catena: extract_regions -> seed_centroids -> tessellate_voronoi -> derive_cell_geometry -> validate_cell_layer.
# Dipendenze: sf, geos, ggplot2 (solo plot_cell_layer); file R/04b0...R/04b6.

.pc_mut <- function() Sys.getenv("S15_MUTANT", "none")

#' Semi indipendenti per gli stadi dell'orchestratore
#'
#' `set.seed(random_seed); sample.int(.Machine$integer.max, 3L)`: semina, taglia dei nuclei, nuclei pallidi. Lo stato
#' globale del generatore e' ripristinato all'uscita.
#' @export
stage_seeds <- function(random_seed) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(random_seed)
  s <- sample.int(.Machine$integer.max, 3L)
  if (.pc_mut() == "M4") s[2] <- s[1]
  c(seed = s[1], size = s[2], pale = s[3])
}

#' Componenti piu' piccole della soglia assegnate al cluster adiacente (BL-043)
#'
#' Ogni componente connessa (4-vicinato) di area < `min_region_area_um2` passa al cluster della componente grande con cui
#' condivide piu' lati di pixel (a pari merito il cluster con l'id minore); si ripete finche' nulla cambia. Le componenti
#' che toccano solo il fondo o solo componenti piccole restano come sono (le escludera' `extract_regions()`).
#' @param clust data.frame con `x`, `y` interi e la colonna dei cluster.
#' @param stride passo dei pixel (intero).
#' @return `clust` con la colonna dei cluster aggiornata e l'attributo `merge_info` (passate, pixel riassegnati).
#' @export
merge_small_components <- function(clust, pixel_size_um, min_region_area_um2, cluster_col = "intensity_cluster", stride = 1L) {
  x <- as.integer(round(clust$x)); y <- as.integer(round(clust$y)); s <- as.integer(stride)
  i <- (x - min(x)) %/% s + 1L; j <- (y - min(y)) %/% s + 1L; nx <- max(i); ny <- max(j)
  cl <- clust[[cluster_col]]; lev <- sort(unique(if (is.factor(cl)) as.character(cl) else cl))
  code <- match(if (is.factor(cl)) as.character(cl) else cl, lev)
  lab <- matrix(0L, nx, ny); lab[cbind(i, j)] <- code
  a_px <- (s * pixel_size_um)^2; n_moved <- 0L; passes <- 0L
  repeat {
    passes <- passes + 1L
    comp <- matrix(0L, nx, ny); off <- 0L; size <- integer(0); ccode <- integer(0)
    for (k in seq_along(lev)) {
      m <- lab == k; if (!any(m)) next
      lk <- .pc_label4(m); nk <- max(lk)
      comp[m] <- lk[m] + off; size <- c(size, tabulate(lk[m], nk)); ccode <- c(ccode, rep(k, nk)); off <- off + nk
    }
    small <- size * a_px < min_region_area_um2
    if (!any(small)) break
    # coppie di pixel adiacenti con componenti diverse (in orizzontale e in verticale)
    e1 <- c(as.vector(comp[-nx, ]), as.vector(comp[, -ny])); e2 <- c(as.vector(comp[-1, ]), as.vector(comp[, -1]))
    k <- e1 != e2 & e1 > 0 & e2 > 0
    a <- c(e1[k], e2[k]); b <- c(e2[k], e1[k])
    k2 <- small[a] & !small[b]
    if (!any(k2)) break
    tb <- stats::aggregate(list(n = rep(1L, sum(k2))), list(comp = a[k2], code = ccode[b[k2]]), length)
    tb <- tb[order(tb$comp, -tb$n, tb$code), ]; tb <- tb[!duplicated(tb$comp), ]
    mv <- comp %in% tb$comp
    newcode <- tb$code[match(comp[mv], tb$comp)]
    lab[mv] <- newcode; n_moved <- n_moved + sum(mv)
  }
  out <- clust
  newv <- lev[lab[cbind(i, j)]]
  out[[cluster_col]] <- if (is.factor(cl)) factor(newv, levels = levels(cl)) else newv
  attr(out, "merge_info") <- list(passes = passes, n_px_moved = n_moved)
  out
}

# componenti 4-connesse di una matrice logica (union-find sui run orizzontali, come extract_regions)
.pc_label4 <- function(m) {
  nx <- nrow(m); ny <- ncol(m); v <- as.vector(m)
  idx <- which(v); if (!length(idx)) return(matrix(0L, nx, ny))
  # run lungo x (prima dimensione) per ogni colonna y
  brk <- c(TRUE, diff(idx) != 1L | ((idx[-1] - 1L) %% nx) == 0L)
  rid <- cumsum(brk); n_r <- max(rid)
  run_of <- integer(length(v)); run_of[idx] <- rid
  # adiacenza verticale fra run: pixel (p) e (p + nx)
  up <- idx[idx + nx <= length(v)]; up <- up[v[up + nx]]
  pr <- unique(cbind(run_of[up], run_of[up + nx]))
  par <- seq_len(n_r)
  find <- function(a) { while (par[a] != a) a <- par[a]; a }
  if (nrow(pr)) for (q in seq_len(nrow(pr))) {
    ra <- find(pr[q, 1]); rb <- find(pr[q, 2]); if (ra != rb) par[max(ra, rb)] <- min(ra, rb)
  }
  root <- vapply(seq_len(n_r), find, 1L)
  lab <- integer(length(v)); lab[idx] <- match(root[rid], unique(root[rid]))
  matrix(lab, nx, ny)
}

#' Posiziona le cellule: il livello cellulare dello Step 1
#'
#' @param clust output di `cluster_image()`: pixel (x, y interi), colonna dei cluster e, per la stima della scala, `value`.
#' @param tissue_preset uno di `tissue_presets()`; ignorato se `cell_types` e' dato.
#' @param cell_types,region_composition catalogo e composizione dell'utente (design §4.4, §4.5); NULL = dal preset.
#' @param override_cell_types,override_composition righe che sostituiscono (per `cell_type` / `cluster_id` + `cell_type`)
#'   o si aggiungono a catalogo e composizione risolti.
#' @param cluster_to_dominant data.frame `cluster_id`, `cell_type`: tipo dominante per cluster quando la composizione e'
#'   costruita dal catalogo (un tipo per cluster); NULL con un catalogo di un tipo = quel tipo per tutti i cluster.
#' @param pixel_size_um obbligatorio: numero (µm/px), `"estimate"` o `list(scalefactors_json = <file>)`.
#' @param size_conversion conversione del modello di taglia dei preset (`tissue_preset_spec()`).
#' @param types,size_arm,composition preset per tipo (`"multi"`, predefinito da S1.5b) o a un tipo (`"single"`, S1.5), braccio del
#'   modello di taglia (NULL = `"recentred"`)
#'   e composizione (`tissue_preset_spec()`); con `"multi"` ogni cluster riceve la composizione del preset, la semina e' cieca al
#'   tipo (`seed_centroids(type_order = "blind")`) e `cell_df` riceve `size_group` ed `epithelium_label`.
#' @param spatial,min_dist_mode semina (`"uniform"` o `"compartments"`, NULL = predefinito del preset) e distanza minima
#'   (`"tissue"` o `"nucleus"`) dei preset per tipo (S1.5c, `tissue_preset_spec()`).
#' @param labels con i compartimenti: `"field"` o `"blind"` (etichette permutate, controprova di S1.5c).
#' @param override_compartments tabella `cell_type`, `ell_um`, `s` che sostituisce quella del preset (banco di S1.5c: pieghe).
#' @param heterogeneity,override_heterogeneity,quotas S1.5d, solo con i compartimenti: NULL = tabella del preset; `override_heterogeneity`
#'   (una riga `void_frac`, `void_ell_um`, `dens_sigma`, `dens_ell_um`) la sostituisce (banco: pieghe); `quotas` come `seed_centroids()`.
#'   Attenzione: con `spatial = "compartments"` e senza `override_heterogeneity` si usano i vuoti calibrati del preset (S1.5d, righe
#'   B-259...B-264; per esempio intestino 0.25 con 16 um), quindi la semina NON e' quella di S1.5c. Per riprodurre S1.5c passare
#'   `override_heterogeneity = data.frame(void_frac = 0, void_ell_um = 0, dens_sigma = 0, dens_ell_um = 0)` (test T7a, T19).
#' @param size_resid,resid_mode S1.5d: residui del modello di taglia (`"gaussian"` o `"empirical"`, NULL = predefinito del preset) e
#'   modo dei residui empirici (`"decile"`; `"pooled"` solo per la controprova CP-3).
#' @param small_components `"merge"` (BL-043) o `"exclude"`.
#' @param stride passo dei pixel; NULL = attributo `stride` di `clust`, altrimenti 1 se i dati sono contigui, errore se
#'   il passo rilevato e' > 1 (BL-042).
#' @param allow_estimate se TRUE accetta `pixel_size_um = "estimate"` (sconsigliato: controprova CP-1 di S1.5).
#' @param keep_tiles se TRUE restituisce anche i tile di Voronoi non ritagliati (`tiles`, diagnostica: regola d'interna di R3).
#' @return `cell_layer` (design §4.1): `cell_df`, `region_df`, `region_polygons`, `cell_territories`, `cell_nuclei`,
#'   `territory_df`, `region_check`, `region_composition`, `cell_types_resolved`, `metadata`.
#' @export
place_cells <- function(clust, tissue_preset = "epithelial", cell_types = NULL, region_composition = NULL,
                        override_cell_types = NULL, override_composition = NULL, cluster_to_dominant = NULL,
                        pixel_size_um, corner_smoothing = 0, min_region_area_um2 = 100, random_seed = 42,
                        cluster_col = "intensity_cluster", stride = NULL, simplify_tol_um = 0.5,
                        small_components = c("merge", "exclude"), nucleus_model = "clip",
                        size_conversion = c("preserve_area", "none"), types = c("multi", "single"), size_arm = NULL,
                        composition = "tissue", spatial = NULL, min_dist_mode = c("tissue", "nucleus"), labels = c("field", "blind"),
                        override_compartments = NULL, override_heterogeneity = NULL, quotas = TRUE, size_resid = NULL,
                        resid_mode = c("decile", "pooled"), validate = TRUE, strict = TRUE, keep_tiles = FALSE,
                        allow_estimate = FALSE, verbose = TRUE) {
  min_dist_mode <- match.arg(min_dist_mode); labels <- match.arg(labels); resid_mode <- match.arg(resid_mode)
  t0 <- proc.time()[["elapsed"]]
  small_components <- match.arg(small_components); size_conversion <- match.arg(size_conversion); types <- match.arg(types)
  if (missing(pixel_size_um))
    stop("place_cells(): `pixel_size_um` e' obbligatorio (decisione 10 / C-A): un numero, \"estimate\" o list(scalefactors_json = ...).", call. = FALSE)
  if (!is.data.frame(clust) || !all(c("x", "y", cluster_col) %in% names(clust)))
    stop("place_cells(): `clust` deve avere le colonne x, y e ", cluster_col, call. = FALSE)
  # ---- catalogo e composizione ----
  spec <- NULL
  if (is.null(cell_types)) {
    spec <- tissue_preset_spec(tissue_preset, size_conversion, types = types, size_arm = size_arm, composition = composition,
                               spatial = spatial, min_dist_mode = min_dist_mode, size_resid = size_resid)
    ct <- spec$cell_types; c8 <- spec$c8_interval
  } else { ct <- cell_types; c8 <- NULL }
  if (!is.null(override_cell_types)) {
    for (k in seq_len(nrow(override_cell_types))) {
      r <- override_cell_types[k, , drop = FALSE]; w <- match(r$cell_type, ct$cell_type)
      for (nm in setdiff(names(r), names(ct))) ct[[nm]] <- NA
      if (is.na(w)) { ct <- rbind(ct, r[, names(ct), drop = FALSE]) } else for (nm in names(r)) ct[w, nm] <- r[[nm]]
    }
  }
  cids <- sort(unique(if (is.factor(clust[[cluster_col]])) as.character(clust[[cluster_col]]) else clust[[cluster_col]]))
  if (is.null(region_composition)) {
    if (!is.null(cluster_to_dominant)) {
      comp <- data.frame(cluster_id = cluster_to_dominant$cluster_id, cell_type = cluster_to_dominant$cell_type, fraction = 1)
    } else if (nrow(ct) == 1L) {
      comp <- data.frame(cluster_id = cids, cell_type = ct$cell_type, fraction = 1)
    } else if (!is.null(spec) && !is.null(spec$composition)) {   # preset per tipo (S1.5b): stessa composizione in ogni cluster
      comp <- data.frame(cluster_id = rep(cids, each = nrow(spec$composition)), cell_type = rep(spec$composition$cell_type, length(cids)),
                         fraction = rep(spec$composition$fraction, length(cids)))
    } else stop("place_cells(): con piu' tipi nel catalogo servono `region_composition` o `cluster_to_dominant`.", call. = FALSE)
  } else comp <- region_composition
  if (!is.null(override_composition)) {
    key <- paste(comp$cluster_id, comp$cell_type); ko <- paste(override_composition$cluster_id, override_composition$cell_type)
    comp <- rbind(comp[!key %in% ko, c("cluster_id", "cell_type", "fraction")], override_composition[, c("cluster_id", "cell_type", "fraction")])
  }
  comp <- comp[comp$fraction > 0, , drop = FALSE]
  miss <- setdiff(comp$cell_type, ct$cell_type)
  if (length(miss)) stop("place_cells(): tipi nella composizione ma non nel catalogo: ", paste(miss, collapse = ", "), call. = FALSE)
  # ---- scala ----
  src <- "user"; est <- NULL
  if (identical(pixel_size_um, "estimate")) {
    if (!allow_estimate)
      stop("place_cells(): la stima della scala dall'immagine non e' raccomandata (S1.5, controprova CP-1: sbaglia in 5 tessuti su 6 e ",
           "puo' dare numeri sbagliati senza rifiuto). Passare pixel_size_um o list(scalefactors_json = ...); per usarla comunque allow_estimate = TRUE.", call. = FALSE)
    est <- estimate_pixel_size(clust, ct, cluster_to_dominant = cluster_to_dominant, method = "nuclear_spacing", cluster_col = cluster_col)
    if (isTRUE(est$refused)) stop("place_cells(): stima della scala rifiutata: ", est$reason, ". Passare pixel_size_um.", call. = FALSE)
    px <- est$pixel_size_um; src <- "estimated"
    warning(sprintf("place_cells(): scala stimata dall'immagine, %.4g um/px (IC 95%% %.4g-%.4g).", px, est$ci[1], est$ci[2]), call. = FALSE)
  } else if (is.list(pixel_size_um) && !is.null(pixel_size_um$scalefactors_json)) {
    est <- estimate_pixel_size(clust, ct, method = "metadata", scalefactors_json = pixel_size_um$scalefactors_json)
    px <- est$pixel_size_um; src <- "metadata"
  } else if (is.numeric(pixel_size_um) && length(pixel_size_um) == 1L && is.finite(pixel_size_um) && pixel_size_um > 0) {
    px <- pixel_size_um
  } else stop("place_cells(): `pixel_size_um` non valido.", call. = FALSE)
  px_used <- if (.pc_mut() == "M7") 1 else px
  # ---- passo dei pixel (BL-042) ----
  if (is.null(stride)) {
    if (!is.null(attr(clust, "stride"))) stride <- as.integer(attr(clust, "stride")) else {
      sx <- .er_axis_stride(as.integer(round(clust$x))); sy <- .er_axis_stride(as.integer(round(clust$y)))
      sdet <- max(c(sx, sy, 1L), na.rm = TRUE)
      if (sdet > 1L) stop(sprintf("place_cells(): passo dei pixel %d rilevato dai dati e non dichiarato: passare `stride` (BL-042).", sdet), call. = FALSE)
      stride <- 1L
    }
  }
  # ---- catena ----
  ss <- stage_seeds(random_seed)
  cl_in <- if (small_components == "merge" && min_region_area_um2 > 0)
    merge_small_components(clust, px_used, min_region_area_um2, cluster_col, stride) else clust
  reg <- extract_regions(cl_in, px_used, min_region_area_um2 = min_region_area_um2, cluster_col = cluster_col,
                         simplify_tol_um = simplify_tol_um, stride = stride, verbose = FALSE)
  comp_used <- comp[as.character(comp$cluster_id) %in% as.character(reg$region_df$cluster_id), , drop = FALSE]
  if (!nrow(reg$region_df)) stop("place_cells(): nessuna regione sopra la soglia di area.", call. = FALSE)
  comp_used$cluster_id <- reg$region_df$cluster_id[match(as.character(comp_used$cluster_id), as.character(reg$region_df$cluster_id))]
  # preset per tipo (S1.5b): semina cieca al tipo, cosi' il territorio non dipende dal tipo (revisione RC-01)
  t_order <- if (!is.null(spec) && identical(spec$types, "multi") && Sys.getenv("S15B_MUTANT", "none") != "M8") "blind" else "sequential"   # M8: mutante
  sp_used <- if (!is.null(spec)) spec$spatial else if (is.null(spatial)) "uniform" else spatial
  ctab <- if (!is.null(override_compartments)) override_compartments else if (!is.null(spec)) spec$compartments else NULL
  sc <- if (sp_used == "compartments") {
    if (is.null(ctab)) stop("place_cells(): spatial = \"compartments\" richiede una tabella dei compartimenti.", call. = FALSE)
    htab <- if (!is.null(override_heterogeneity)) override_heterogeneity else if (!is.null(spec)) spec$heterogeneity else NULL
    seed_centroids(reg, ct, comp_used, random_seed = ss[["seed"]], spatial = "compartments", compartments = ctab,
                   labels = labels, heterogeneity = htab, quotas = quotas, verbose = FALSE)
  } else seed_centroids(reg, ct, comp_used, random_seed = ss[["seed"]], type_order = t_order, verbose = FALSE)
  tv <- tessellate_voronoi(sc$centroids, reg, corner_smoothing = corner_smoothing, keep_tiles = keep_tiles, verbose = FALSE)
  rtab <- if ("nucleus_area_resid" %in% names(ct) && any(!is.na(ct$nucleus_area_resid))) .tp_resid_table() else NULL   # S1.5d
  dg <- derive_cell_geometry(sc$centroids, tv, ct, random_seed = ss[["size"]], nucleus_model = nucleus_model, verbose = FALSE,
                             resid_table = rtab, resid_mode = resid_mode)
  cd <- dg$cell_df
  # tipo dominante per regione e nuclei pallidi
  dom <- comp_used[order(comp_used$cluster_id, -comp_used$fraction), ]; dom <- dom[!duplicated(dom$cluster_id), ]
  rdf <- sc$region_df
  rdf$dominant_cell_type <- dom$cell_type[match(rdf$cluster_id, dom$cluster_id)]
  cd$is_minority_type <- as.character(cd$cell_type) != rdf$dominant_cell_type[match(cd$region_id, rdf$region_id)]
  for (nm in intersect(c("size_group", "epithelium_label"), names(ct)))      # preset per tipo (S1.5b, revisione RC-14)
    cd[[nm]] <- ct[[nm]][match(as.character(cd$cell_type), ct$cell_type)]
  pf <- if ("pale_fraction" %in% names(ct)) ct$pale_fraction[match(as.character(cd$cell_type), ct$cell_type)] else rep(0, nrow(cd))
  pf[is.na(pf)] <- 0
  u <- .pc_with_seed(ss[["pale"]], stats::runif(nrow(cd)))
  cd$nucleus_pale <- u < pf
  nuc <- dg$cell_nuclei
  if (.pc_mut() == "M5" && length(nuc) > 1) nuc <- nuc[c(2:length(nuc), 1)]
  cl_out <- list(cell_df = cd, region_df = rdf, region_polygons = reg$region_polygons, cell_territories = tv$cell_territories,
                 cell_nuclei = nuc, territory_df = tv$territory_df, region_check = tv$region_check,
                 region_composition = comp_used, cell_types_resolved = ct)
  if (keep_tiles) cl_out$tiles <- tv$tiles
  val <- NULL
  if (validate) {
    val <- validate_cell_layer(cl_out, c8_interval = c8, pixel_size_source = src, min_region_area_um2 = min_region_area_um2)
    if (strict && !val$ok) {
      bad <- val$table[val$table$status == "ERROR", ]
      stop("place_cells(): controlli falliti: ", paste(bad$check, bad$value, collapse = " | "), call. = FALSE)
    }
  }
  cl_out$metadata <- list(pixel_size_um = px, pixel_size_source = src, pixel_size_estimate = est,
                          n_cells = nrow(cd), n_regions = nrow(rdf), tissue_preset_used = if (is.null(spec)) NA_character_ else tissue_preset,
                          preset_sources = if (is.null(spec)) NULL else spec$sources, size_conversion = if (is.null(spec)) NA_character_ else size_conversion,
                          preset_types = if (is.null(spec)) NA_character_ else spec$types, size_arm = if (is.null(spec)) NA_character_ else spec$size_arm,
                          composition_name = if (is.null(spec)) NA_character_ else spec$composition_name,
                          type_sources = if (is.null(spec)) NULL else spec$type_sources,
                          spatial = sp_used, min_dist_mode = if (is.null(spec)) NA_character_ else spec$min_dist_mode,
                          labels = if (sp_used == "compartments") labels else NA_character_,
                          compartments = if (sp_used == "compartments") ctab else NULL,
                          heterogeneity = if (sp_used == "compartments") sc$info$compartments$heterogeneity else NULL, quotas = quotas,
                          size_resid = if (is.null(spec)) (if (is.null(rtab)) "gaussian" else "empirical") else spec$size_resid, resid_mode = resid_mode,
                          corner_smoothing = corner_smoothing, min_region_area_um2 = min_region_area_um2, small_components = small_components,
                          stride = stride, cell_types_resolved = ct, region_composition_resolved = comp_used,
                          random_seed = random_seed, stage_seeds = ss, nucleus_model = nucleus_model,
                          merge_info = attr(cl_in, "merge_info"), excluded_df = reg$excluded_df,
                          seed_info = sc$info, tessellation_info = tv$info, geometry_info = dg$info,
                          package_version = "geo.spatialtrans step1-cell-layer (S1.5d)", timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                          validation = val, elapsed_s = proc.time()[["elapsed"]] - t0)
  if (verbose) cat(sprintf("[place_cells] %s: %d regioni, %d cellule, scala %.4g um/px (%s), controlli %s, %.1f s\n",
                           if (is.null(spec)) "catalogo dell'utente" else tissue_preset, nrow(rdf), nrow(cd), px, src,
                           if (is.null(val)) "non eseguiti" else sprintf("%d ERROR / %d WARN", val$n_error, val$n_warn),
                           cl_out$metadata$elapsed_s))
  cl_out
}

.pc_with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed); force(expr)
}

#' Figura del livello cellulare (design §6.3)
#' @param cell_layer output di `place_cells()`.
#' @param mode `"voronoi"`, `"types"`, `"density"` o `"regions"`.
#' @return oggetto ggplot (sfondo bianco, y verso il basso come l'immagine).
#' @export
plot_cell_layer <- function(cell_layer, mode = c("voronoi", "types", "density", "regions")) {
  mode <- match.arg(mode)
  cl <- cell_layer; rd <- cl$region_df; cd <- cl$cell_df
  R <- sf::st_sf(cluster_id = factor(rd$cluster_id), ratio = 100 * rd$achieved_density / rd$target_density_weighted, geometry = cl$region_polygons)
  th <- ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                                                               panel.grid = ggplot2::element_blank())
  p <- ggplot2::ggplot() + ggplot2::scale_y_reverse() + th + ggplot2::labs(x = "x (um)", y = "y (um)")
  if (mode == "regions") {
    p <- p + ggplot2::geom_sf(data = R, ggplot2::aes(fill = cluster_id), colour = "grey30", linewidth = 0.2) +
      ggplot2::labs(title = "Regioni del tessuto", fill = "cluster")
  } else if (mode == "types") {
    P <- sf::st_as_sf(data.frame(cell_type = cd$cell_type, x = cd$x, y = cd$y), coords = c("x", "y"))
    p <- p + ggplot2::geom_sf(data = R, fill = NA, colour = "grey60", linewidth = 0.2) +
      ggplot2::geom_sf(data = P, ggplot2::aes(colour = cell_type), size = 0.3) + ggplot2::labs(title = "Centri e tipi cellulari", colour = "tipo")
  } else if (mode == "voronoi") {
    Tt <- sf::st_sf(cell_type = cd$cell_type, geometry = cl$cell_territories)
    N <- sf::st_sf(pale = cd$nucleus_pale, geometry = cl$cell_nuclei)
    p <- p + ggplot2::geom_sf(data = Tt, ggplot2::aes(fill = cell_type), colour = "black", linewidth = 0.1, alpha = 0.35) +
      ggplot2::geom_sf(data = N, ggplot2::aes(colour = pale), fill = "grey25", linewidth = 0.05) +
      ggplot2::scale_colour_manual(values = c(`FALSE` = "grey25", `TRUE` = "grey70")) +
      ggplot2::labs(title = "Territori e nuclei", fill = "tipo", colour = "nucleo pallido")
  } else {
    p <- p + ggplot2::geom_sf(data = R, ggplot2::aes(fill = ratio), colour = "grey30", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(midpoint = 100, low = "#2166ac", high = "#b2182b") +
      ggplot2::labs(title = "Densita' raggiunta / target (%)", fill = "%")
  }
  p
}
