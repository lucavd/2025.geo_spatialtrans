# R/05b5_simulate_visium_hd.R - Step 2, S2.5: simulatore completo (orchestratore dei passi 1-6).
# Documento di progetto dello Step 2 (docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md), sezioni 1, 7, 8, 11, 12, 13, 15;
# sessione S2.5 (docs/steps/S2.5.md). Mette in fila build_hd_grid() (S2.1), simulate_stripes() (S2.4), simulate_cell_profiles() (S2.2),
# place_and_move_transcripts() (S2.3, S2.4) e capture_and_bin() (S2.4); sceglie i blocchi di cellule dai trascritti (BL-227) e somma a 8 um
# anche la verita'. Dipendenze: Matrix, sf, geos; le funzioni di R/04b* e R/05b1...R/05b4.

#' Archetipo di un livello cellulare
#'
#' @param cell_layer uscita di `place_cells()`.
#' @return `"A1"` ... `"A6"`, dal preset usato (`metadata$tissue_preset_used`, tabella dei preset di R/04b0).
#' @export
visium_hd_archetype <- function(cell_layer) {
  p <- cell_layer$metadata$tissue_preset_used
  a <- if (is.null(p)) NA_character_ else .TP_TABLE$archetype[match(p, .TP_TABLE$preset)]
  if (length(a) != 1L || is.na(a)) stop("simulate_visium_hd(): archetipo non ricavabile dal preset del livello cellulare; indicare `archetype`.", call. = FALSE)
  a
}

#' Cellule per blocco dai trascritti (BL-227)
#'
#' Il picco di memoria di `place_and_move_transcripts()` cresce con i trascritti generati per blocco di cellule (S2.4, PH-1): il blocco
#' si sceglie perche' contenga in media `chunk_transcripts` trascritti generati.
#' @param n_transcripts trascritti generati per cellula (sovra-generazione compresa).
#' @param chunk_transcripts trascritti per blocco (predefinito 2e6: S2.5, regola R-1 con 4e6 e correzione dopo PH-1, area piu' carica del tumore
#'   8.9 GB con 4e6, 7.1 GB con 2e6).
#' @param max_cells tetto del blocco (cellule; 2000 = valore delle funzioni di S2.2 e S2.3).
#' @return intero >= 1.
#' @export
visium_hd_chunk_size <- function(n_transcripts, chunk_transcripts = 2e6, max_cells = 2000L) {
  if (!is.numeric(chunk_transcripts) || length(chunk_transcripts) != 1L || !is.finite(chunk_transcripts) || chunk_transcripts < 1)
    stop("visium_hd_chunk_size(): `chunk_transcripts` deve essere un numero >= 1.", call. = FALSE)
  m <- if (length(n_transcripts)) mean(as.numeric(n_transcripts)) else 0
  if (!is.finite(m) || m <= 0) return(as.integer(max_cells))
  as.integer(max(1, min(max_cells, floor(chunk_transcripts / m))))
}

#' Fattore di profondita' del preset dopo S2.5
#'
#' Fattore della regola R-2 di S2.5 (ritaratura sulla catena completa: f' = f (U_vero - b) / (U_sim - b), UMI per quadratino di 2 um in tessuto,
#' pannello delle 30 aree; `results/S2.5/S2.5_R2_depth.csv`). Sostituisce, per i profili dei preset, il fattore di S2.2.
#' @param archetype `"A1"` ... `"A6"`.
#' @return numero.
#' @export
visium_hd_depth_factor <- function(archetype) {
  tab <- c(A1 = 1.0127249977836763, A2 = 0.9540327133313311, A3 = 1.5614917408538553, A4 = 0.8191399388185701, A5 = 1.4278291184727523, A6 = 2.9386716762361265)
  if (!archetype %in% names(tab)) stop("visium_hd_depth_factor(): archetipo sconosciuto: ", archetype, call. = FALSE)
  tab[[archetype]]
}

# somma a 8 um di una tabella (quadratino, cellula) con chiavi intere (equivale a sum_to_8um(x, grid, cols, by = "cell_id"), piu' veloce)
.s25_to_8um <- function(df, grid, cols) {
  nc8 <- grid$n_col %/% grid$block
  b8 <- as.integer(((df$bin_id - 1L) %/% grid$n_col) %/% grid$block * nc8 + ((df$bin_id - 1L) %% grid$n_col) %/% grid$block + 1L)
  if (!nrow(df)) return(data.frame(bin8_id = integer(0), cell_id = df$cell_id[0], df[0, cols, drop = FALSE]))
  ucell <- sort(unique(df$cell_id)); ci <- match(df$cell_id, ucell)
  key <- (as.numeric(b8) - 1) * length(ucell) + (ci - 1)
  s <- rowsum(as.matrix(df[cols]), key, reorder = TRUE)
  k <- as.numeric(rownames(s))
  out <- data.frame(bin8_id = as.integer(k %/% length(ucell) + 1), cell_id = ucell[k %% length(ucell) + 1])
  for (cc in cols) out[[cc]] <- if (is.integer(df[[cc]])) as.integer(s[, cc]) else as.numeric(s[, cc])
  rownames(out) <- NULL
  out
}

# impronte e versione (S2.8, BL-247) ---------------------------------------------------------------------------------------------------
# impronta md5 del codice: testo senza commenti (deparse a 17 cifre) della funzione `root`, di tutte le funzioni dello stesso ambiente che
# chiama, seguite ricorsivamente (codetools::findGlobals), e degli oggetti dello stesso ambiente che leggono come variabili (per esempio
# la tabella .s23_ztab); una modifica del codice o di una costante cambia l'impronta, una modifica dei soli commenti no. Limiti (S2.8,
# revisione): non segue funzioni chiamate per nome con do.call("nome") o match.fun, metodi S3/S4 e chiusure create a runtime (oggi assenti).
.geo_code_parts <- function(root) {
  home <- environment(.geo_code_parts); seen <- character(0); todo_f <- root; todo_v <- character(0); parts <- list()
  ctl <- c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17")
  while (length(todo_f) || length(todo_v)) {
    is_f <- length(todo_f) > 0L
    nm <- if (is_f) todo_f[[1L]] else todo_v[[1L]]
    if (is_f) todo_f <- todo_f[-1L] else todo_v <- todo_v[-1L]
    key <- paste0(if (is_f) "f:" else "v:", nm)
    if (key %in% seen) next
    seen <- c(seen, key)
    if (is_f) {                                                    # chiamata: R cerca solo funzioni (un oggetto omonimo dell'utente non conta)
      fn <- get0(nm, envir = home, mode = "function", inherits = TRUE)
      if (is.null(fn) || is.primitive(fn) || !identical(environment(fn), home)) next
      parts[[paste0("## ", nm)]] <- deparse(fn, control = ctl)
      gl <- codetools::findGlobals(fn, merge = FALSE)
      todo_f <- c(todo_f, gl$functions); todo_v <- c(todo_v, gl$variables)
    } else if (exists(nm, envir = home, inherits = FALSE)) {      # variabile letta: conta l'oggetto dell'ambiente del codice
      v <- get(nm, envir = home, inherits = FALSE)
      if (is.function(v)) {
        if (!identical(environment(v), home)) next
        todo_f <- c(todo_f, nm)                                    # funzione passata come valore: si segue come una chiamata
      } else parts[[paste0("## <costante> ", nm)]] <- deparse(v, control = ctl)
    }
  }
  if (!length(parts)) stop(".geo_code_parts(): funzione sconosciuta: ", root, call. = FALSE)
  parts[order(names(parts), method = "radix")]                    # ordine indipendente dalla lingua del sistema
}
.geo_code_hash <- function(root) {
  parts <- .geo_code_parts(root)
  tf <- tempfile(); on.exit(unlink(tf))
  writeLines(unlist(Map(function(n, d) c(n, d), names(parts), parts), use.names = FALSE), tf, useBytes = TRUE)
  unname(tools::md5sum(tf))
}
# impronta md5 di un oggetto R (serializzazione versione 2: l'intestazione non contiene la codifica della lingua del sistema)
.geo_object_hash <- function(x) {
  tf <- tempfile(); on.exit(unlink(tf)); con <- file(tf, "wb"); serialize(x, con, version = 2L); close(con)
  unname(tools::md5sum(tf))
}
# versione del pacchetto: installato, altrimenti DESCRIPTION nella cartella `geo.spatialtrans.root` (predefinito: cartella di lavoro)
.geo_version <- function() {
  v <- tryCatch(as.character(utils::packageVersion("geo.spatialtrans")), error = function(e) NA_character_)
  if (is.na(v)) {
    f <- file.path(getOption("geo.spatialtrans.root", "."), "DESCRIPTION")
    if (file.exists(f)) { d <- read.dcf(f, fields = c("Package", "Version")); if (identical(unname(d[1, "Package"]), "geo.spatialtrans")) v <- unname(d[1, "Version"]) }
  }
  paste("geo.spatialtrans", if (is.na(v)) "(versione sconosciuta)" else v)
}

.s25_draw_seeds <- function(seed) {
  old_kind <- RNGkind(); had <- exists(".Random.seed", envir = globalenv())
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had) assign(".Random.seed", old, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }, add = TRUE)
  # semi dei quattro passi estratti con Mersenne-Twister dal seme dell'utente (niente semi consecutivi: S2.4 C-1, DEC-057, BL-224)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  s <- sample.int(1e8, 4L)
  c(stripes = s[1], profiles = s[2], place = s[3], capture = s[4])
}

#' Simulatore di Visium HD dal livello cellulare (Step 2 completo)
#'
#' Passi: griglia e verita' geometrica (`build_hd_grid()`), strisce (`simulate_stripes()`), trascritti per cellula con il fattore di
#' profondita' per la sovra-generazione M (`simulate_cell_profiles()`), posizione e spostamento con la seconda coda
#' (`place_and_move_transcripts()`), cattura, fondo e somma a 8 um (`capture_and_bin()`). Parametri dei preset per tessuto
#' (`transcript_params()`, `capture_params()`, profili di `load_cell_profiles()`). I blocchi di cellule si scelgono dai trascritti
#' (`visium_hd_chunk_size()`); i semi dei passi si estraggono dal seme dell'utente; lo stato del generatore dell'utente e' ripristinato.
#' @param cell_layer uscita di `place_cells()` (serve `metadata$pixel_size_um`).
#' @param seed intero.
#' @param archetype `"A1"` ... `"A6"`; predefinito: dal preset del livello.
#' @param layers strati accesi: `stripes`, `background`, `field`, `second_tail` (tutti `TRUE` per difetto; i metodi che correggono uno
#'   strato si provano con lo strato acceso).
#' @param frame_open `TRUE` per aree ritagliate da un tessuto che continua (il riquadro non e' bordo; seconda coda riflessa).
#' @param profiles uscita di `load_cell_profiles()`; predefinito: quella del tessuto.
#' @param depth_factor fattore di profondita' prima della sovra-generazione; predefinito: con i profili del preset `visium_hd_depth_factor()`
#'   (S2.5, R-2), con profili dati dall'utente il loro `depth_factor`.
#' @param chunk_size cellule per blocco della posizione e dello spostamento; predefinito (NULL): dai trascritti generati.
#' @param chunk_transcripts trascritti per blocco usati quando `chunk_size` e' NULL.
#' @param component_seeds semi dei quattro passi (`stripes`, `profiles`, `place`, `capture`); se dati, `seed` non si usa per estrarli (test).
#' @param hd griglia gia' costruita con `build_hd_grid()` sullo stesso livello (si controlla); predefinito: si costruisce.
#' @param workers processi per i passi che lavorano a blocchi (il risultato non dipende da `workers`).
#' @param verbose messaggi di avanzamento.
#' @return lista `counts_2um`, `counts_8um` (dgCMatrix geni x quadratini o blocchi), `bins`, `bins_8um`, `truth` (`overlap`, `overlap_8um`,
#'   `origin`, `origin_8um`, `background`, `background_8um`, `arrived_2um` (trascritti arrivati per quadratino prima della cattura), `cell_counts`, `cells`, `genes`, `layers`), `metadata` (design sezione 1).
#' @export
simulate_visium_hd <- function(cell_layer, seed, archetype = NULL,
                               layers = list(stripes = TRUE, background = TRUE, field = TRUE, second_tail = TRUE),
                               frame_open = FALSE, profiles = NULL, depth_factor = NULL, chunk_size = NULL, chunk_transcripts = 2e6,
                               component_seeds = NULL, hd = NULL, workers = 1L, verbose = FALSE) {
  t0 <- proc.time()[["elapsed"]]; tm <- c()
  # stato del generatore dell'utente registrato in testa e ripristinato all'uscita: alcune funzioni compilate chiamate dai passi
  # (imager, Rcpp) creano un .Random.seed se manca (revisione G-5 di S2.5)
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
            else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }, add = TRUE)
  px <- cell_layer$metadata$pixel_size_um
  if (is.null(px) || !is.numeric(px) || length(px) != 1L || !is.finite(px) || px <= 0)
    stop("simulate_visium_hd(): `cell_layer$metadata$pixel_size_um` manca o non e' valido (C-A): nessun default.", call. = FALSE)
  if (is.null(component_seeds) && (missing(seed) || length(seed) != 1L || !is.finite(seed))) stop("simulate_visium_hd(): serve un `seed`.", call. = FALSE)
  ly <- list(stripes = TRUE, background = TRUE, field = TRUE, second_tail = TRUE)
  if (!is.list(layers) || (length(layers) && (is.null(names(layers)) || any(!nzchar(names(layers))) || anyDuplicated(names(layers)))) ||
      length(setdiff(names(layers), names(ly))) || any(!vapply(layers, function(v) is.logical(v) && length(v) == 1L && !is.na(v), logical(1))))
    stop("simulate_visium_hd(): `layers` accetta solo stripes, background, field, second_tail (TRUE o FALSE).", call. = FALSE)
  ly[names(layers)] <- layers
  if (is.null(archetype)) archetype <- visium_hd_archetype(cell_layer)
  if (!archetype %in% paste0("A", 1:6)) stop("simulate_visium_hd(): archetipo sconosciuto: ", archetype, call. = FALSE)
  if (!is.null(chunk_size) && (!is.numeric(chunk_size) || length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1))
    stop("simulate_visium_hd(): `chunk_size` deve essere un intero >= 1.", call. = FALSE)
  sd_ <- if (!is.null(component_seeds)) {
    s <- component_seeds[c("stripes", "profiles", "place", "capture")]
    if (any(is.na(names(s))) || any(!is.finite(unlist(s)))) stop("simulate_visium_hd(): `component_seeds` deve avere stripes, profiles, place, capture.", call. = FALSE)
    unlist(s)
  } else .s25_draw_seeds(seed)
  # (1) griglia
  if (is.null(hd)) {
    hd <- build_hd_grid(cell_layer)
  } else {
    bb <- sf::st_bbox(c(sf::st_set_crs(cell_layer$cell_territories, NA), sf::st_set_crs(cell_layer$region_polygons, NA)))
    g0 <- hd_grid_geometry(bb, margin_um = hd$grid$margin_um, bin_um = hd$grid$bin_um)
    if (!isTRUE(all.equal(g0[c("x0", "y0", "n_row", "n_col")], hd$grid[c("x0", "y0", "n_row", "n_col")], tolerance = 0)) ||
        !identical(as.integer(hd$metadata$n_cells), nrow(cell_layer$cell_df)))
      stop("simulate_visium_hd(): `hd` non e' la griglia di questo livello cellulare.", call. = FALSE)
    # stesso livello anche dentro il riquadro: seme dello Step 1 e area del territorio di ogni cellula (revisione G-5 di S2.5)
    s1 <- cell_layer$metadata$random_seed
    aT <- geos::geos_area(geos::as_geos_geometry(sf::st_set_crs(cell_layer$cell_territories, NA)))
    sa <- rowsum(hd$truth$overlap$area_territory_um2, hd$truth$overlap$cell_id)
    aH <- sa[match(as.character(cell_layer$cell_df$cell_id), rownames(sa)), 1]
    if ((!is.null(s1) && !is.null(hd$metadata$step1_random_seed) && !identical(as.numeric(s1), as.numeric(hd$metadata$step1_random_seed))) ||
        anyNA(aH) || max(abs(aH - aT)) > 1e-3)
      stop("simulate_visium_hd(): `hd` ha territori diversi da quelli del livello cellulare (seme dello Step 1 o aree per cellula).", call. = FALSE)
  }
  tm["grid"] <- proc.time()[["elapsed"]] - t0
  cp <- capture_params(archetype); p <- transcript_params(archetype)
  p$w2 <- if (ly$second_tail) cp$w2 else 0; p$ell2 <- cp$ell2
  # (2) strisce
  t1 <- proc.time()[["elapsed"]]
  St <- if (ly$stripes) simulate_stripes(hd, cp, seed = sd_[["stripes"]]) else NULL
  M <- if (is.null(St)) 1 else St$M
  tm["stripes"] <- proc.time()[["elapsed"]] - t1
  # (3) profili, con il fattore di profondita' per M (sovra-generazione delle strisce)
  t1 <- proc.time()[["elapsed"]]
  prof_default <- is.null(profiles)
  if (prof_default) profiles <- load_cell_profiles(archetype)
  # profondita': argomento; altrimenti, con i profili del preset, il fattore ritarato in S2.5 (R-2); altrimenti quello dei profili dati
  df0 <- if (!is.null(depth_factor)) depth_factor else if (prof_default) visium_hd_depth_factor(archetype) else profiles$depth_factor
  prof <- simulate_cell_profiles(cell_layer$cell_df[, c("cell_id", "cell_type")], profiles, seed = sd_[["profiles"]], depth_factor = df0 * M, workers = workers)
  tm["profiles"] <- proc.time()[["elapsed"]] - t1
  ntx <- Matrix::colSums(prof$counts)
  cs <- if (is.null(chunk_size)) visium_hd_chunk_size(ntx, chunk_transcripts) else as.integer(chunk_size)
  if (verbose) message(sprintf("simulate_visium_hd(): %s, %d cellule, M %.3f, %.3g trascritti generati, blocchi di %d cellule", archetype, length(ntx), M, sum(ntx), cs))
  # (4) posizione e spostamento
  t1 <- proc.time()[["elapsed"]]
  tx <- place_and_move_transcripts(cell_layer, hd, prof, p, seed = sd_[["place"]], field = ly$field, frame_open = frame_open, chunk_size = cs, workers = workers)
  tm["place"] <- proc.time()[["elapsed"]] - t1
  # (5) cattura, fondo, somma a 8 um
  t1 <- proc.time()[["elapsed"]]
  cap <- capture_and_bin(tx, hd, St, list(b = cp$b), seed = sd_[["capture"]], overgeneration = M, background = ly$background, workers = workers)
  tm["capture"] <- proc.time()[["elapsed"]] - t1
  n_tx <- nrow(tx$transcripts)
  arr <- integer(hd$grid$n_row * hd$grid$n_col)                   # trascritti arrivati per quadratino, prima dell'assottigliamento
  if (n_tx) { a_ <- rowsum(tx$transcripts$n, tx$transcripts$bin_id); arr[as.integer(rownames(a_))] <- as.integer(a_[, 1]) }
  rm(tx); gc(FALSE)
  # (6) verita' a 8 um e uscite
  t1 <- proc.time()[["elapsed"]]
  g <- hd$grid
  ov <- hd$truth$overlap
  ov8 <- .s25_to_8um(ov[, c("bin_id", "cell_id", "area_territory_um2", "area_nucleus_um2")], g, c("area_territory_um2", "area_nucleus_um2"))
  ov8$frac_territory <- ov8$area_territory_um2 / (g$bin_um * g$block)^2; ov8$frac_nucleus <- ov8$area_nucleus_um2 / (g$bin_um * g$block)^2
  og8 <- .s25_to_8um(cap$truth$origin, g, c("n", "n_born_nucleus"))
  bg <- cap$truth$background
  bg8 <- if (nrow(bg)) { s <- sum_to_8um(bg, g, cols = "n"); s$n <- as.integer(s$n); s } else data.frame(bin8_id = integer(0), n = integer(0))
  cells <- merge(prof$cells[, c("cell_id", "cell_type", "subtype", "n_transcripts")], cap$truth$cells, by = "cell_id", sort = FALSE)
  cells <- cells[match(cell_layer$cell_df$cell_id, cells$cell_id), ]; rownames(cells) <- NULL
  names(cells)[names(cells) == "n_transcripts"] <- "n_transcripts_generated"
  cells$n_expected_captured <- cells$n_transcripts_generated / M
  tm["truth8"] <- proc.time()[["elapsed"]] - t1
  lay <- cap$truth$layers
  lay$diffusion <- p[c("r", "sigma_in", "sigma_e", "D", "w", "ell", "w2", "ell2")]
  lay$field <- list(on = ly$field, s2 = p$field_s2, L = p$field_L)
  lay$depth <- list(depth_factor = df0, overgeneration = M, profiles_winner = profiles$meta$winner, noise = profiles$noise)
  list(counts_2um = cap$counts_2um, counts_8um = cap$counts_8um, bins = cap$bins, bins_8um = hd$bins_8um,
       truth = list(overlap = ov, overlap_8um = ov8, origin = cap$truth$origin, origin_8um = og8, background = bg, background_8um = bg8, arrived_2um = arr,
                    cell_counts = prof$counts, cells = cells, genes = cap$truth$genes, layers = lay),
       metadata = list(step = "S2.5", archetype = archetype, seed = if (is.null(component_seeds)) seed else NA, component_seeds = sd_,
                       layers = ly, frame_open = frame_open, chunk_size = cs, chunk_transcripts = chunk_transcripts, overgeneration = M,
                       pixel_size_um = px, grid = g, n_cells = nrow(cell_layer$cell_df), n_transcript_triples = n_tx,
                       n_generated = sum(cells$n_generated), n_captured = sum(cells$n_captured), n_background = sum(bg$n), n_lost = sum(cells$n_lost),
                       step1 = cell_layer$metadata[c("tissue_preset_used", "random_seed", "package_version")],
                       package_version = .geo_version(), code_hash = .geo_code_hash("simulate_visium_hd"), elapsed_s = c(tm, total = proc.time()[["elapsed"]] - t0)))
}
