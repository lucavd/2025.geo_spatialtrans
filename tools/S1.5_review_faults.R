# tools/S1.5_review_faults.R - S1.5, controprova dopo la revisione: i guasti indipendenti del revisore del codice
# (copia di /mnt/micron/geo_spatialtrans/S1.5/review_code/reviewer_faults.R, stessa costruzione) rieseguiti sul validatore corretto.
# Aggiunte: detected_any (qualsiasi ERROR; i layer di base hanno 0 ERROR) e by (controlli in ERROR). Uso: Rscript --vanilla tools/S1.5_review_faults.R
.libPaths("/home/user/2025.geo_spatialtrans/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

suppressPackageStartupMessages({ library(sf); library(geos); library(parallel) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
OUT <- "/mnt/micron/geo_spatialtrans/S1.5/postreview"; dir.create(OUT, showWarnings = FALSE); LAY <- "/mnt/micron/geo_spatialtrans/S1.5/layers"
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
sq <- function(p, delta) { h <- sqrt(delta) / 2; xy <- sf::st_coordinates(p)[1, 1:2]
  sf::st_polygon(list(rbind(xy + c(-h, -h), xy + c(h, -h), xy + c(h, h), xy + c(-h, h), xy + c(-h, -h)))) }
setg <- function(sfc, i, g) { sfc[[i]] <- g; sf::st_sfc(sfc) }
U <- function(a, b) sf::st_union(sf::st_sfc(a), sf::st_sfc(b))[[1]]
D <- function(a, b) sf::st_difference(sf::st_sfc(a), sf::st_sfc(b))[[1]]
drop_cells <- function(cl, keep) {
  cl$cell_df <- cl$cell_df[keep, ]; cl$cell_df$cell_id <- seq_len(nrow(cl$cell_df)); rownames(cl$cell_df) <- NULL
  cl$cell_territories <- cl$cell_territories[keep]; cl$cell_nuclei <- cl$cell_nuclei[keep]
  cl$territory_df <- cl$territory_df[keep, ]; cl$territory_df$cell_id <- seq_len(nrow(cl$territory_df)); cl }
inner_square <- function(cl, i, delta) sq(sf::st_point_on_surface(sf::st_sfc(D(cl$cell_territories[[i]], sf::st_buffer(cl$cell_nuclei[i], 1)[[1]]))), delta)
out_square <- function(cl, i, delta) { R <- cl$region_polygons[match(cl$cell_df$region_id[i], cl$region_df$region_id)]
  out <- sf::st_difference(sf::st_buffer(cl$cell_territories[i], 3), sf::st_buffer(R, 0.5))
  # lontano da ogni territorio (fondo)
  far <- sf::st_difference(out, sf::st_buffer(sf::st_union(cl$cell_territories[sf::st_intersects(sf::st_buffer(cl$cell_territories[i], 4), cl$cell_territories)[[1]]]), 0.3))
  if (all(sf::st_is_empty(far))) stop("nessun fondo vicino")
  sq(sf::st_point_on_surface(far), delta) }
FAULTS <- list(
  F1_exit_gap_unflagged = function(cl, d) { td <- cl$territory_df; i <- which(td$clipped & td$territory_area > 20 * d + 5)[1]
    S <- out_square(cl, i, d); G <- inner_square(cl, i, d)
    cl$cell_territories <- setg(cl$cell_territories, i, D(U(cl$cell_territories[[i]], S), G))
    cl$territory_df$clipped[i] <- FALSE; cl$territory_df$n_pieces_gained[i] <- 0L; list(cl = cl, expect = "C2") },
  F1b_exit_gap_flagged = function(cl, d) { td <- cl$territory_df; i <- which(td$clipped & td$territory_area > 20 * d + 5)[1]
    S <- out_square(cl, i, d); G <- inner_square(cl, i, d)
    cl$cell_territories <- setg(cl$cell_territories, i, D(U(cl$cell_territories[[i]], S), G)); list(cl = cl, expect = "C2") },
  F2_region_cells_missing = function(cl, d) { rd <- cl$region_df; cand <- rd$region_id[rd$n_cells > 0]
    if (length(cand) < 2) stop("una sola regione con cellule")
    r <- cand[which.min(rd$n_cells[match(cand, rd$region_id)])]
    list(cl = drop_cells(cl, cl$cell_df$region_id != r), expect = "C2|C6") },
  F3_swap_two_cells = function(cl, d) { td <- cl$territory_df; ok <- which(!td$clipped & td$n_pieces_gained == 0)
    nb <- sf::st_touches(cl$cell_territories[ok], cl$cell_territories); k <- which(lengths(nb) > 0)[1]; i <- ok[k]; j <- intersect(nb[[k]], ok)[1]
    Tt <- cl$cell_territories; Nn <- cl$cell_nuclei; ti <- Tt[[i]]; Tt[[i]] <- Tt[[j]]; Tt[[j]] <- ti; ni <- Nn[[i]]; Nn[[i]] <- Nn[[j]]; Nn[[j]] <- ni
    cl$cell_territories <- sf::st_sfc(Tt); cl$cell_nuclei <- sf::st_sfc(Nn); list(cl = cl, expect = "C9|C3|C2") },
  F4_duplicate_cell = function(cl, d) { n <- nrow(cl$cell_df); cl$cell_df <- rbind(cl$cell_df, cl$cell_df[n, ]); cl$cell_df$cell_id <- seq_len(n + 1); rownames(cl$cell_df) <- NULL
    cl$cell_territories <- c(cl$cell_territories, cl$cell_territories[n]); cl$cell_nuclei <- c(cl$cell_nuclei, cl$cell_nuclei[n])
    cl$territory_df <- rbind(cl$territory_df, cl$territory_df[n, ]); list(cl = cl, expect = "C2") },
  F5_region_id_NA = function(cl, d) { cl$cell_df$region_id[10] <- NA; list(cl = cl, expect = "C2") },
  F6_nucleus_in_background = function(cl, d) { td <- cl$territory_df; i <- which(td$clipped)[1]
    a <- as.numeric(sf::st_area(cl$cell_nuclei[i])); S <- out_square(cl, i, a)
    cl$cell_nuclei <- setg(cl$cell_nuclei, i, S); list(cl = cl, expect = "C3") },
  F7_empty_nucleus = function(cl, d) { cl$cell_nuclei <- setg(cl$cell_nuclei, 5, sf::st_polygon()); list(cl = cl, expect = "C9") },
  F9_small_gap = function(cl, d) { td <- cl$territory_df; i <- which(!td$clipped & td$territory_area > 10)[1]
    cl$cell_territories <- setg(cl$cell_territories, i, D(cl$cell_territories[[i]], inner_square(cl, i, d))); list(cl = cl, expect = "C2") },
  F10_nucleus_halved = function(cl, d) { i <- 7; g <- cl$cell_nuclei[i]; ctr <- sf::st_coordinates(sf::st_centroid(g))[1, 1:2]
    cl$cell_nuclei <- setg(cl$cell_nuclei, i, ((g - ctr) * sqrt(0.5) + ctr)[[1]]); list(cl = cl, expect = "C5|C3") }
)
DELTAS <- c(F1_exit_gap_unflagged = 1, F1b_exit_gap_flagged = 1, F2_region_cells_missing = NA, F3_swap_two_cells = NA, F4_duplicate_cell = NA,
            F5_region_id_NA = NA, F6_nucleus_in_background = NA, F7_empty_nucleus = NA, F9_small_gap = 0.001, F10_nucleus_halved = NA)
fl <- list.files(LAY, pattern = "_preserve_area\\.rds$", full.names = TRUE)
res <- mclapply(fl, function(f) {
  L <- readRDS(f); A <- sub("_.*", "", basename(f)); roi <- strsplit(basename(f), "_")[[1]][2]
  c8 <- tissue_preset_spec(PRESET[[A]])$c8_interval
  v0 <- validate_cell_layer(L, c8_interval = c8)
  rows <- list(data.frame(archetype = A, roi = roi, fault = "none", delta = NA, expect = "-", detected = NA, detected_any = NA, by = "", n_error = v0$n_error,
                          statuses = paste(v0$table$check, v0$table$status, sep = ":", collapse = " "), note = sprintf("aR max %.4g um2; regions %d", max(as.numeric(sf::st_area(L$region_polygons))), nrow(L$region_df))))
  for (nm in names(FAULTS)) {
    d <- DELTAS[[nm]]
    r <- tryCatch({ fz <- FAULTS[[nm]](L, d)
      v <- tryCatch(validate_cell_layer(fz$cl, c8_interval = c8), error = function(e) e)
      if (inherits(v, "error")) data.frame(archetype = A, roi = roi, fault = nm, delta = d, expect = fz$expect, detected = FALSE, detected_any = FALSE, by = "CRASH", n_error = NA, statuses = "CRASH", note = conditionMessage(v)) else {
        ex <- strsplit(fz$expect, "\\|")[[1]]; st <- v$table$status[match(ex, v$table$check)]
        det <- any(st == "ERROR", na.rm = TRUE) || ("C6" %in% ex && isTRUE(st[ex == "C6"] == "WARN"))
        det_any <- v$n_error > 0 || ("C6" %in% ex && isTRUE(v$table$status[v$table$check == "C6"] == "WARN"))
        data.frame(archetype = A, roi = roi, fault = nm, delta = d, expect = fz$expect, detected = det, detected_any = det_any, by = paste(v$table$check[v$table$status == "ERROR"], collapse = "+"), n_error = v$n_error,
                   statuses = paste(v$table$check, v$table$status, sep = ":", collapse = " "), note = paste(v$table$value[v$table$check %in% c("C2", "C3")], collapse = " || ")) }
    }, error = function(e) data.frame(archetype = A, roi = roi, fault = nm, delta = d, expect = "?", detected = NA, detected_any = NA, by = "", n_error = NA, statuses = "NOT BUILT", note = conditionMessage(e)))
    rows[[length(rows) + 1]] <- r
  }
  do.call(rbind, rows)
}, mc.cores = 6, mc.preschedule = FALSE)
bad <- !vapply(res, is.data.frame, TRUE); if (any(bad)) print(res[bad])
tab <- do.call(rbind, res[!bad]); write.csv(tab, "results/S1.5/S1.5_review_faults_after.csv", row.names = FALSE)
print(aggregate(cbind(built = !is.na(detected), detected = detected %in% TRUE, detected_any = detected_any %in% TRUE, crash = statuses == "CRASH") ~ fault, data = tab, FUN = sum))
cat("baseline layers with ERROR:", sum(tab$fault == "none" & tab$n_error > 0), "of", sum(tab$fault == "none"), "\n")
