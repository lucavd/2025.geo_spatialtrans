# R/04b5_cell_layer_validation.R - Step 1, S1.5: controlli del livello cellulare (C1-C9) e helper unico per gli overlay.
# Design §6.1 (C1-C8) e §7.2 (C9); definizioni operative in docs/steps/S1.5.md. Overlay sempre su griglia fissa (F1 di
# S1.3d: geos_set_precision + geos_*_prec), mai sf::st_union / st_intersection sui territori (BL-069, BL-070, BL-081).
# Dipendenze: sf, geos.

.v5_mut <- function() Sys.getenv("S15_MUTANT", "none")

#' Overlay su griglia fissa (helper unico, F1 di S1.3d)
#'
#' Trasla x e y nella stessa origine locale (BL-094, BL-104), porta entrambi sulla griglia `grid_um` e applica
#' l'operazione elemento per elemento con `geos_*_prec`. Con `area = TRUE` restituisce le aree (numeric, µm²),
#' altrimenti le geometrie in coordinate originali come `sfc` via WKB (BL-116).
#' @param op `"intersection"`, `"difference"` o `"union"`.
#' @param x,y geometrie `sfc` (stessa lunghezza, oppure y di lunghezza 1).
#' @param grid_um passo della griglia (µm).
#' @param origin origine locale (x, y); NULL = angolo minimo del riquadro di x.
#' @param area se TRUE restituisce solo le aree.
#' @export
overlay_prec <- function(op = c("intersection", "difference", "union"), x, y, grid_um = 1e-6, origin = NULL, area = TRUE) {
  op <- match.arg(op)
  if (!inherits(x, "sfc") || !inherits(y, "sfc")) stop("overlay_prec(): x e y devono essere sfc", call. = FALSE)
  if (is.null(origin)) { bb <- sf::st_bbox(x); origin <- c(bb[["xmin"]], bb[["ymin"]]) }
  o <- as.numeric(origin)
  gx <- geos::geos_set_precision(geos::as_geos_geometry(sf::st_set_crs(x - o, NA)), grid_um)
  gy <- geos::geos_set_precision(geos::as_geos_geometry(sf::st_set_crs(y - o, NA)), grid_um)
  r <- switch(op, intersection = geos::geos_intersection_prec(gx, gy, grid_um),
              difference = geos::geos_difference_prec(gx, gy, grid_um), union = geos::geos_union_prec(gx, gy, grid_um))
  if (area) return(geos::geos_area(r))
  sf::st_as_sfc(geos::geos_write_wkb(r)) + o
}

# conteggio dei fallimenti: NA/NaN conta come fallimento (BL-095)
.nbad <- function(ok) sum(is.na(ok) | !ok)

.v5_row <- function(check, status, value, threshold, n_bad = NA_integer_, note = "")
  data.frame(check = check, status = status, value = as.character(value), threshold = as.character(threshold),
             n_bad = as.integer(n_bad), note = note, stringsAsFactors = FALSE)

# coppie (i < j) di geometrie che si toccano o si intersecano (predicato GEOS, albero STR di sf)
.v5_pairs <- function(g) {
  s <- sf::st_intersects(g)
  i <- rep(seq_along(s), lengths(s)); j <- unlist(s, use.names = FALSE)
  k <- i < j
  cbind(i[k], j[k])
}

#' Controlli del livello cellulare (C1-C9)
#'
#' @param cell_layer lista con `cell_df`, `region_df`, `region_polygons`, `cell_territories`, `cell_nuclei` e
#'   (facoltativi) `region_check` (lacune della smussatura), `territory_df` (`clipped`, `n_pieces_gained`: territori di
#'   bordo per il contenimento nella regione), `region_composition`, `cell_types_resolved`.
#' @param c8_interval intervallo plausibile della mediana di eq_r (µm); NULL = provvisorio [5, 25] del design.
#' @param pixel_size_source `"user"`, `"metadata"` o `"estimated"` (C8 e' ERROR solo con `"estimated"`).
#' @param grid_um griglia degli overlay (µm).
#' @param sum_rel_tol tolleranza relativa di C2 sulla somma delle aree per regione.
#' @param min_region_area_um2 soglia di C6.
#' @return lista con `table` (check, status PASS/WARN/ERROR, value, threshold, n_bad, note), `n_error`, `n_warn`, `ok`.
#' @export
validate_cell_layer <- function(cell_layer, c8_interval = NULL, pixel_size_source = "user", grid_um = 1e-6,
                                sum_rel_tol = 1e-9, min_region_area_um2 = 0) {
  cl <- cell_layer; cd <- cl$cell_df; rd <- cl$region_df
  Tt <- cl$cell_territories; Nn <- cl$cell_nuclei; R <- cl$region_polygons
  g <- grid_um; out <- list()
  # ---- C9 validita', ordine e legame righe-geometrie (FND-03, FND-04 della revisione di S1.5) ----
  n <- nrow(cd)
  ord_ok <- length(Tt) == n && length(Nn) == n && identical(as.numeric(cd$cell_id), as.numeric(seq_len(n)))
  rid_bad <- .nbad(cd$region_id %in% rd$region_id)
  cnt_bad <- if ("n_cells" %in% names(rd) && rid_bad == 0)   # region_df$n_cells = conteggio dalle cellule
    sum(rd$n_cells != as.numeric(table(factor(cd$region_id, levels = rd$region_id))), na.rm = TRUE) + sum(is.na(rd$n_cells)) else 0L
  gT <- geos::as_geos_geometry(sf::st_set_crs(Tt, NA)); gN <- geos::as_geos_geometry(sf::st_set_crs(Nn, NA))
  bad_valid <- .nbad(geos::geos_is_valid(gT)) + .nbad(geos::geos_is_valid(gN))
  bad_empty <- sum(geos::geos_is_empty(gT)) + sum(geos::geos_is_empty(gN))
  bad_link <- 0L
  if (ord_ok && bad_valid == 0 && bad_empty == 0) {
    gen <- geos::geos_make_point(cd$x, cd$y)
    a_t <- geos::geos_area(gT); a_n <- geos::geos_area(gN)
    bad_link <- .nbad(geos::geos_intersects(gen, gT)) +                                   # generatore nel proprio territorio
      .nbad(!is.finite(cd$territory_area) | abs(a_t - cd$territory_area) <= 1e-9 * pmax(a_t, 1)) +                    # area del territorio = cell_df
      .nbad(!is.finite(cd$nucleus_area) | abs(a_n - cd$nucleus_area) <= 1e-9 * pmax(a_n, 1) + grid_um * geos::geos_length(gN))   # area del nucleo = cell_df (NaN: C5)
  }
  nb9 <- (!ord_ok) + rid_bad + cnt_bad + bad_valid + bad_empty + bad_link
  out$c9 <- .v5_row("C9", if (nb9 == 0) "PASS" else "ERROR",
                    sprintf("ordine %s; region_id fuori o NA %d; n_cells di region_df diversi dalle cellule %d; non valide %d; vuote %d; righe non legate alle geometrie %d",
                            if (ord_ok) "ok" else "DIVERSO", rid_bad, cnt_bad, bad_valid, bad_empty, bad_link),
                    "ordine identico, region_id validi, n_cells coerenti, 0 non valide, 0 vuote, generatore nel territorio, aree = cell_df", nb9)
  if (nb9 > 0) return(.v5_finish(out))   # gli altri controlli non sono affidabili
  aT <- geos::geos_area(gT); pT <- geos::geos_length(gT); aN <- geos::geos_area(gN); pN <- geos::geos_length(gN)
  # una sola conversione: origine locale (BL-094/BL-104) e griglia g (F1), poi solo indici
  bb <- sf::st_bbox(Tt); o <- c(bb[["xmin"]], bb[["ymin"]])
  loc <- function(z) geos::geos_set_precision(wk::wk_transform(z, wk::wk_affine_translate(-o[1], -o[2])), g)
  pT_loc <- loc(gT); pN_loc <- loc(gN)
  gR <- geos::as_geos_geometry(sf::st_set_crs(R, NA))
  # ---- C1 densita' per regione ----
  if (all(c("achieved_density", "target_density_weighted") %in% names(rd))) {
    rel <- abs(rd$achieved_density - rd$target_density_weighted) / rd$target_density_weighted
    nb <- sum(rel > 0.2 | !is.finite(rel), na.rm = TRUE)
    out$c1 <- .v5_row("C1", if (nb == 0) "PASS" else "WARN", sprintf("max %.3g", suppressWarnings(max(rel, na.rm = TRUE))),
                      "|raggiunta - target| / target <= 0.20", nb)
  }
  # ---- C2 partizione ----
  # lacuna ammessa = sola perdita della smussatura dichiarata per territorio (FND-01: non la differenza calcolata dalla tassellazione)
  gap <- rep(0, nrow(rd)); td0 <- cl$territory_df
  if (!is.null(td0) && "smooth_loss" %in% names(td0) && nrow(td0) == n) {
    sl <- tapply(td0$smooth_loss, factor(cd$region_id, levels = rd$region_id), sum); sl[is.na(sl)] <- 0; gap <- as.numeric(sl)
  }
  if (!is.null(cl$region_check) && "gap_area" %in% names(cl$region_check)) {   # una lacuna NaN dichiarata resta un errore (BL-095)
    mi <- match(rd$region_id, cl$region_check$region_id); gq <- cl$region_check$gap_area[mi]
    gap[!is.na(mi) & is.nan(gq)] <- NaN
  }
  aR <- geos::geos_area(gR)
  sT <- tapply(aT, factor(cd$region_id, levels = rd$region_id), sum); sT[is.na(sT)] <- 0
  has_cells <- rd$region_id %in% cd$region_id            # le regioni senza cellule sono il caso di C6, non di C2
  rel_sum <- (abs(as.numeric(sT) + gap - aR) / aR)[has_cells]
  if (.v5_mut() == "M3") rel_sum[is.na(rel_sum)] <- 0
  nb_sum <- .nbad(rel_sum <= sum_rel_tol)                                    # NA/NaN -> errore (BL-095)
  pr <- .v5_pairs(Tt)
  if (nrow(pr)) {
    tau <- g * (pT[pr[, 1]] + pT[pr[, 2]]) + g^2
    a12 <- geos::geos_area(geos::geos_intersection_prec(pT_loc[pr[, 1]], pT_loc[pr[, 2]], g)); n_eval <- length(a12)
    a21 <- if (.v5_mut() == "M2") a12 else { z <- geos::geos_area(geos::geos_intersection_prec(pT_loc[pr[, 2]], pT_loc[pr[, 1]], g)); n_eval <- n_eval + length(z); z }
    if (.v5_mut() == "M3") { a12[is.na(a12)] <- 0; a21[is.na(a21)] <- 0 }
    nb_pair <- .nbad(a12 <= tau & a21 <= tau)
    max_pair <- suppressWarnings(max(c(a12, a21)))
  } else { nb_pair <- 0L; max_pair <- 0; n_eval <- 0L }
  # contenimento nella regione: territori di bordo (ritagliati o che hanno ricevuto frammenti)
  td <- cl$territory_df
  flag <- if (!is.null(td) && all(c("clipped", "n_pieces_gained") %in% names(td)) && nrow(td) == n) td$clipped | td$n_pieces_gained > 0 else rep(FALSE, n)
  # FND-02: anche dalla geometria (territorio non dentro la propria regione, predicato su regione preparata)
  geo_out <- rep(FALSE, n)
  for (k in seq_len(nrow(rd))) { idx <- which(cd$region_id == rd$region_id[k])
    if (length(idx)) geo_out[idx] <- !geos::geos_prepared_covers(gR[k], gT[idx]) }
  bnd <- which(flag | geo_out | is.na(geo_out))
  nb_out <- 0L; max_out <- 0
  if (length(bnd)) {
    ext <- geos::geos_extent(gT[bnd])
    rect <- list(ext$xmin - 1, ext$ymin - 1, ext$xmax + 1, ext$ymax + 1)
    Rloc <- loc(geos::geos_clip_by_rect(gR[match(cd$region_id[bnd], rd$region_id)], rect))
    ao <- geos::geos_area(geos::geos_difference_prec(pT_loc[bnd], Rloc, g))
    nb_out <- .nbad(ao <= g * pT[bnd] + g^2); max_out <- suppressWarnings(max(ao))
  }
  if (.v5_mut() == "M1") {   # mutante: partizione misurata con st_union (BL-070)
    U <- vapply(rd$region_id, function(r) as.numeric(sf::st_area(sf::st_union(Tt[cd$region_id == r]))), 0)
    nb_sum <- sum(abs(U - aR) / aR > 0.01); nb_pair <- 0L; nb_out <- 0L
  }
  nb2 <- nb_sum + nb_pair + nb_out
  out$c2 <- .v5_row("C2", if (nb2 == 0) "PASS" else "ERROR",
                    sprintf("somma: max rel %.3g; coppie %d (max %.3g um2); fuori regione max %.3g um2",
                            suppressWarnings(max(rel_sum)), nrow(pr), max_pair, max_out),
                    sprintf("somma <= %g rel; coppie e uscite <= g*P + g^2, g = %g", sum_rel_tol, g), nb2,
                    sprintf("regioni %d, coppie %d, uscite %d; intersezioni valutate %d", nb_sum, nb_pair, nb_out, n_eval))
  # ---- C3 contenimento del nucleo e nuclei disgiunti ----
  dn <- geos::geos_area(geos::geos_difference_prec(pN_loc, pT_loc, g))
  wi <- geos::geos_within(geos::geos_read_wkb(sf::st_as_binary(Nn)), geos::geos_read_wkb(sf::st_as_binary(Tt)))
  nb_c3 <- .nbad(dn <= g * pN + g^2) + .nbad(wi)
  pn <- .v5_pairs(Nn)
  nb_ov <- 0L; max_ov <- 0
  if (nrow(pn)) {
    an <- geos::geos_area(geos::geos_intersection_prec(pN_loc[pn[, 1]], pN_loc[pn[, 2]], g))
    nb_ov <- .nbad(an <= g * (pN[pn[, 1]] + pN[pn[, 2]]) + g^2); max_ov <- suppressWarnings(max(an))
  }
  out$c3 <- .v5_row("C3", if (nb_c3 + nb_ov == 0) "PASS" else "ERROR",
                    sprintf("fuori %d (max %.3g um2), non within %d; coppie di nuclei %d, sovrapposte %d (max %.3g um2)",
                            .nbad(dn <= g * pN + g^2), suppressWarnings(max(dn)), .nbad(wi), nrow(pn), nb_ov, max_ov),
                    "area(N \\ T) <= g*P_N; within; intersezioni fra nuclei <= g*(P_i+P_j) + g^2", nb_c3 + nb_ov)
  # ---- C4 frazioni per tipo ----
  if (!is.null(cl$region_composition) && nrow(cd)) {
    comp <- cl$region_composition; nb4 <- 0L; mx <- 0
    ncr <- table(factor(cd$region_id, levels = rd$region_id))
    for (r in rd$region_id[as.numeric(ncr) >= 10]) {
      cid <- rd$cluster_id[rd$region_id == r]; cc <- comp[comp$cluster_id == cid, ]
      obs <- table(factor(cd$cell_type[cd$region_id == r], levels = cc$cell_type)) / sum(cd$region_id == r)
      d <- max(abs(as.numeric(obs) - cc$fraction / sum(cc$fraction))); mx <- max(mx, d, na.rm = TRUE); nb4 <- nb4 + (is.na(d) || d > 0.10)
    }
    out$c4 <- .v5_row("C4", if (nb4 == 0) "PASS" else "WARN", sprintf("max %.3g", mx), "|osservata - target| <= 0.10 (regioni con >= 10 cellule)", nb4)
  }
  # ---- C5 area del nucleo ----
  nb5 <- .nbad(cd$nucleus_area > 0) + .nbad(aN > 0)
  out$c5 <- .v5_row("C5", if (nb5 == 0) "PASS" else "ERROR", sprintf("min %.3g um2", suppressWarnings(min(aN))), "> 0", nb5)
  # ---- C6 regioni senza cellule ----
  if (TRUE) {
    ncr6 <- as.numeric(table(factor(cd$region_id, levels = rd$region_id)))   # FND-04: dalle cellule, non da region_df$n_cells
    nb6 <- sum(rd$area_um2 >= min_region_area_um2 & ncr6 == 0)
    out$c6 <- .v5_row("C6", if (nb6 == 0) "PASS" else "WARN", nb6, "0 regioni senza cellule", nb6)
  }
  # ---- C7 eq_r per tipo ----
  if (!is.null(cl$cell_types_resolved)) {
    ct <- cl$cell_types_resolved; nb7 <- 0L; vals <- character(0)
    for (k in seq_len(nrow(ct))) {
      e <- cd$eq_radius[cd$cell_type == ct$cell_type[k]]; if (!length(e)) next
      r7 <- unname(stats::quantile(e, 0.9)) / sqrt(1e6 / (pi * ct$density[k]))
      vals <- c(vals, sprintf("%s %.3g", ct$cell_type[k], r7)); nb7 <- nb7 + !(r7 >= 0.5 && r7 <= 2)
    }
    out$c7 <- .v5_row("C7", if (nb7 == 0) "PASS" else "WARN", paste(vals, collapse = "; "), "p90 / atteso in [0.5, 2]", nb7)
  }
  # ---- C8 scala ----
  iv <- if (is.null(c8_interval)) c(5, 25) else c8_interval
  med <- stats::median(cd$eq_radius)
  ok8 <- is.finite(med) && med >= iv[1] && med <= iv[2]
  out$c8 <- .v5_row("C8", if (ok8) "PASS" else if (identical(pixel_size_source, "estimated")) "ERROR" else "WARN",
                    sprintf("mediana eq_r %.3g um", med), sprintf("[%g, %g] um%s", iv[1], iv[2], if (is.null(c8_interval)) " (provvisorio)" else ""),
                    as.integer(!ok8))
  .v5_finish(out)
}

.v5_finish <- function(out) {
  tab <- do.call(rbind, out); rownames(tab) <- NULL
  list(table = tab, n_error = sum(tab$status == "ERROR"), n_warn = sum(tab$status == "WARN"), ok = !any(tab$status == "ERROR"))
}
