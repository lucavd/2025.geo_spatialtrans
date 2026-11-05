# tools/S2.1_panel.R - S2.1: build_hd_grid() sulle 30 aree simulate dello Step 1 (livelli U6 di S1.6), invarianti G-2, descrittivi simulati.
# Definizioni: docs/steps/S2.1.md. Uso (root del repo): OMP_NUM_THREADS=1 Rscript --vanilla tools/S2.1_panel.R run | check | desc
# Stadi rilanciabili (un file per area). Processi dal file <OUT>/WORKERS (al piu' 6, server condiviso). Uscite: <M>/<S21_RUN>/{grid,check,desc}.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(Matrix); library(parallel) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
STAGE <- commandArgs(TRUE)[1]
M <- "/mnt/micron/geo_spatialtrans"; RUN <- Sys.getenv("S21_RUN", "S2.1/run1"); OUT <- file.path(M, RUN); RES <- Sys.getenv("S21_RES", "results/S2.1")
for (d in c("grid", "check", "desc")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
dir.create(RES, recursive = TRUE, showWarnings = FALSE)
wk <- file.path(OUT, "WORKERS"); if (!file.exists(wk)) writeLines("6", wk)
NW <- function() min(6L, as.integer(readLines(wk)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
tags <- paste(rois$archetype, rois$roi_id, sep = "_")
layer_of <- function(tag) readRDS(file.path(M, "S1.6/layers", paste0(tag, "_U6.rds")))
g <- 1e-6
# pezzi (area = FALSE) delle coppie, con la stessa origine di tessera di hd_overlap(): perimetri per le tolleranze di DEC-038
piece_perimeters <- function(grid, ii, jj, geoms, tile_um = 500) {
  nt <- as.integer(round(tile_um / grid$bin_um)); tile <- (ii %/% nt) * 1e6 + jj %/% nt; P <- numeric(length(ii))
  for (t in unique(tile)) {
    k <- which(tile == t); org <- c(grid$x0 + (jj[k[1]] %/% nt) * nt * grid$bin_um, grid$y0 + (ii[k[1]] %/% nt) * nt * grid$bin_um)
    pc <- overlay_prec("intersection", .hd_sfc(.hd_squares(grid, ii[k], jj[k])), geoms[k], grid_um = g, origin = org, area = FALSE)
    P[k] <- geos::geos_length(geos::as_geos_geometry(sf::st_set_crs(pc, NA)))
  }
  P
}
if (STAGE == "run") {
  todo <- tags[!file.exists(file.path(OUT, "grid", paste0(tags, ".rds")))]
  invisible(mclapply(todo, function(tag) {
    L <- layer_of(tag); invisible(gc(reset = TRUE)); t0 <- proc.time()[["elapsed"]]
    H <- build_hd_grid(L); el <- proc.time()[["elapsed"]] - t0; mx <- sum(gc()[, 6])
    saveRDS(H[c("grid", "bins", "bins_8um", "truth")], file.path(OUT, "grid", paste0(tag, ".rds")))
    write.csv(data.frame(tag = tag, n_cells = nrow(L$cell_df), n_bins = nrow(H$bins), n_pairs = nrow(H$truth$overlap), elapsed_s = el,
                         r_max_used_mb = mx, load1 = as.numeric(strsplit(readLines("/proc/loadavg"), " ")[[1]][1])),
              file.path(OUT, "grid", paste0(tag, "_time.csv")), row.names = FALSE)
    cat(sprintf("%s: %d cellule, %d quadratini, %d coppie, %.1f s\n", tag, nrow(L$cell_df), nrow(H$bins), nrow(H$truth$overlap), el))
  }, mc.cores = NW(), mc.preschedule = FALSE))
}
if (STAGE == "check") {
  rows <- mclapply(tags, function(tag) {
    f <- file.path(OUT, "check", paste0(tag, ".csv")); if (file.exists(f)) return(read.csv(f))
    L <- layer_of(tag); H <- readRDS(file.path(OUT, "grid", paste0(tag, ".rds"))); G <- H$grid; ov <- H$truth$overlap; b <- H$bins
    cid <- match(ov$cell_id, L$cell_df$cell_id); ii <- (ov$bin_id - 1L) %/% G$n_col; jj <- (ov$bin_id - 1L) %% G$n_col
    PT <- piece_perimeters(G, ii, jj, L$cell_territories[cid])
    kn <- ov$area_nucleus_um2 > 0; PN <- numeric(nrow(ov)); PN[kn] <- piece_perimeters(G, ii[kn], jj[kn], L$cell_nuclei[cid[kn]])
    n <- nrow(L$cell_df); fc <- factor(cid, levels = seq_len(n))
    nk <- tabulate(cid, n); sT <- as.numeric(tapply(ov$area_territory_um2, fc, sum)); sN <- as.numeric(tapply(ov$area_nucleus_um2, fc, sum))
    sT[is.na(sT)] <- 0; sN[is.na(sN)] <- 0
    tauT <- g * as.numeric(tapply(PT, fc, sum)) + nk * g^2; tauN <- g * as.numeric(tapply(PN, fc, sum)) + nk * g^2; tauT[is.na(tauT)] <- 0; tauN[is.na(tauN)] <- 0
    a_err <- max(abs(sT - L$cell_df$territory_area) / tauT, abs(sN - L$cell_df$nucleus_area) / pmax(tauN, g^2))
    # (b) ogni cellula ha quadratini; i territori stanno a 280 um dal bordo della griglia
    bbT <- sf::st_bbox(L$cell_territories)
    b_bad <- sum(nk == 0) + (bbT[["xmin"]] - G$x0 < 280 - 1e-9) + (bbT[["ymin"]] - G$y0 < 280 - 1e-9) +
      (G$x0 + 2 * G$n_col - bbT[["xmax"]] < 280 - 1e-9) + (G$y0 + 2 * G$n_row - bbT[["ymax"]] < 280 - 1e-9)
    # (c) per quadratino: somma dei territori = area di tessuto; <= 4; nucleo <= territorio per coppia
    bd <- which(b$tissue_area_um2 > 0 & b$tissue_area_um2 < 4); Ptis <- ifelse(b$tissue_area_um2 > 0, 8, 0)
    if (length(bd)) { pr <- numeric(length(bd))
      for (r in seq_along(L$region_polygons)) pr <- pr + piece_perimeters(G, b$i[bd], b$j[bd], rep(L$region_polygons[r], length(bd)))
      Ptis[bd] <- pr }
    fb <- factor(as.integer(ov$bin_id), levels = as.integer(b$bin_id)); sb <- as.numeric(tapply(ov$area_territory_um2, fb, sum)); sb[is.na(sb)] <- 0
    nb <- tabulate(ov$bin_id, nrow(b)); Pb <- as.numeric(tapply(PT, fb, sum)); Pb[is.na(Pb)] <- 0
    tau_b <- g * (Pb + Ptis) + (nb + 1) * g^2
    c_err <- max(abs(sb - b$tissue_area_um2) / tau_b, (sb - 4) / tau_b, 0)
    c_pair <- max((ov$area_nucleus_um2 - ov$area_territory_um2) / (g * (PT + PN) + 2 * g^2), 0)
    # (d) tessuto totale = somma delle aree delle regioni
    d_err <- abs(sum(b$tissue_area_um2) - sum(as.numeric(sf::st_area(L$region_polygons)))) / (g * sum(Ptis) + sum(Ptis > 0) * g^2)
    # (e) somma a 8 um = verita' sui quadrati di 8 um
    o8 <- sum_to_8um(data.frame(bin_id = ov$bin_id, cell_id = ov$cell_id, aT = ov$area_territory_um2, aN = ov$area_nucleus_um2, PT = PT, PN = PN), G, by = "cell_id")
    G8 <- modifyList(G, list(bin_um = 8, n_row = G$n_row %/% 4L, n_col = G$n_col %/% 4L, block = 1L))
    d8 <- hd_overlap(G8, L$cell_territories, L$cell_nuclei); d8$cell_id <- L$cell_df$cell_id[d8$cell_id]
    m8 <- merge(o8, d8, by.x = c("bin8_id", "cell_id"), by.y = c("bin_id", "cell_id"), all = TRUE)
    e_miss <- sum(is.na(m8$aT) & m8$area_territory_um2 > 0) + sum(is.na(m8$area_territory_um2) & m8$aT > 0)
    m8[is.na(m8)] <- 0
    e_err <- max(abs(m8$aT - m8$area_territory_um2) / (2 * g * m8$PT + 17 * g^2), abs(m8$aN - m8$area_nucleus_um2) / pmax(2 * g * m8$PN + 17 * g^2, 17 * g^2))
    out <- data.frame(tag = tag, archetype = sub("_.*", "", tag), n_cells = n, n_pairs = nrow(ov), n_bins_tissue = sum(b$in_tissue),
                      a_ratio = a_err, b_bad = b_bad, c_ratio = c_err, c_pair_ratio = c_pair, d_ratio = d_err, e_ratio = e_err, e_missing = e_miss,
                      n_valid_bad = sum(!geos::geos_is_valid(geos::as_geos_geometry(sf::st_set_crs(L$cell_territories, NA)))))
    out$pass <- out$a_ratio <= 1 && out$b_bad == 0 && out$c_ratio <= 1 && out$c_pair_ratio <= 1 && out$d_ratio <= 1 && out$e_ratio <= 1 && out$e_missing == 0 && out$n_valid_bad == 0
    write.csv(out, f, row.names = FALSE); out
  }, mc.cores = NW(), mc.preschedule = FALSE)
  tab <- do.call(rbind, rows); write.csv(tab, file.path(RES, "S2.1_G2.csv"), row.names = FALSE)
  cat(sprintf("%s G-2 invarianti sulle 30 aree: %d/%d aree PASS; rapporti errore/tolleranza max a %.3f c %.3f d %.3f e %.3f\n",
              if (all(tab$pass)) "PASS" else "FAIL", sum(tab$pass), nrow(tab), max(tab$a_ratio), max(tab$c_ratio), max(tab$d_ratio), max(tab$e_ratio)))
}
if (STAGE == "desc") {
  rows <- mclapply(tags, function(tag) {
    L <- layer_of(tag); H <- readRDS(file.path(OUT, "grid", paste0(tag, ".rds"))); G <- H$grid; ov <- H$truth$overlap; b <- H$bins
    # D-1: centro del quadratino in un nucleo (coppie con area di nucleo > 0)
    kn <- which(ov$area_nucleus_um2 > 0); cid <- match(ov$cell_id[kn], L$cell_df$cell_id)
    pts <- geos::geos_make_point(b$x[ov$bin_id[kn]], b$y[ov$bin_id[kn]])
    inn <- geos::geos_intersects(pts, geos::as_geos_geometry(sf::st_set_crs(L$cell_nuclei[cid], NA)))
    cen_nuc <- logical(nrow(b)); cen_nuc[ov$bin_id[kn][inn]] <- TRUE
    d1 <- mean(cen_nuc[b$in_tissue])
    # D-3b: quadratini puri (un territorio >= 99 % dell'area di tessuto)
    mx <- tapply(ov$area_territory_um2, factor(as.integer(ov$bin_id), levels = as.integer(b$bin_id)), max); mx[is.na(mx)] <- 0
    pure <- mean((mx >= 0.99 * b$tissue_area_um2)[b$in_tissue])
    # D-2 e D-3a sui blocchi di 8 um in tessuto (almeno un figlio in tessuto, come Space Ranger)
    b8 <- H$bins_8um; inb8 <- b8$bin8_id[b8$in_tissue]; nc8 <- G$n_col %/% 4L
    b8_of <- function(id) { i <- (id - 1L) %/% G$n_col; j <- (id - 1L) %% G$n_col; (i %/% 4L) * nc8 + j %/% 4L + 1L }
    q <- data.frame(b8 = b8_of(ov$bin_id), cell = ov$cell_id, n = ov$area_nucleus_um2 > 0)
    nT <- tabulate(unique(q[, c("b8", "cell")])$b8, nrow(b8))[inb8]
    qn <- unique(q[q$n, c("b8", "cell")]); nN <- tabulate(qn$b8, nrow(b8))[inb8]
    cls <- function(v) c(mean = mean(v), p0 = mean(v == 0), p1 = mean(v == 1), p2 = mean(v == 2), p3 = mean(v == 3), p4 = mean(v >= 4))
    data.frame(tag = tag, archetype = sub("_.*", "", tag), source = "simulated", n_bins_tissue = sum(b$in_tissue), n_8um_tissue = length(inb8),
               d1_centre_in_nucleus = d1, d3_pure_bins = pure, t(setNames(cls(nN), paste0("d2_nuclei_", names(cls(nN))))),
               t(setNames(cls(nT), paste0("d3_cells_", names(cls(nT))))))
  }, mc.cores = NW(), mc.preschedule = FALSE)
  tab <- do.call(rbind, rows); write.csv(tab, file.path(RES, "S2.1_desc_sim_roi.csv"), row.names = FALSE); print(head(tab, 3))
}

if (STAGE == "ph1") {   # post hoc PH-1 (dopo il calcolo): G-2 (c) con il tessuto delle sole regioni con cellule
  rows <- mclapply(tags, function(tag) {
    L <- layer_of(tag); H <- readRDS(file.path(OUT, "grid", paste0(tag, ".rds"))); G <- H$grid; ov <- H$truth$overlap; b <- H$bins
    nc <- as.numeric(table(factor(L$cell_df$region_id, levels = L$region_df$region_id))); keepR <- nc > 0
    tis <- hd_tissue_area(G, L$region_polygons[keepR])
    cid <- match(ov$cell_id, L$cell_df$cell_id); ii <- (ov$bin_id - 1L) %/% G$n_col; jj <- (ov$bin_id - 1L) %% G$n_col
    PT <- piece_perimeters(G, ii, jj, L$cell_territories[cid])
    bd <- which(tis > 0 & tis < 4); Ptis <- ifelse(tis > 0, 8, 0)
    if (length(bd)) { pr <- numeric(length(bd))
      for (r in which(keepR)) pr <- pr + piece_perimeters(G, b$i[bd], b$j[bd], rep(L$region_polygons[r], length(bd)))
      Ptis[bd] <- pr }
    fb <- factor(as.integer(ov$bin_id), levels = as.integer(b$bin_id)); sb <- as.numeric(tapply(ov$area_territory_um2, fb, sum)); sb[is.na(sb)] <- 0
    nb <- tabulate(ov$bin_id, nrow(b)); Pb <- as.numeric(tapply(PT, fb, sum)); Pb[is.na(Pb)] <- 0
    tau_b <- g * (Pb + Ptis) + (nb + 1) * g^2
    data.frame(tag = tag, archetype = sub("_.*", "", tag), n_regions = length(keepR), n_regions_no_cells = sum(!keepR),
               area_no_cells_um2 = sum(as.numeric(sf::st_area(L$region_polygons[!keepR]))), area_tissue_um2 = sum(as.numeric(sf::st_area(L$region_polygons))),
               n_bins_touching_no_cell_regions = sum(abs(b$tissue_area_um2 - tis) > 0),
               c_ratio_cells_only = max(abs(sb - tis) / tau_b, (sb - 4) / tau_b, 0))
  }, mc.cores = NW(), mc.preschedule = FALSE)
  tab <- do.call(rbind, rows); tab$pass <- tab$c_ratio_cells_only <= 1
  write.csv(tab, file.path(RES, "S2.1_PH1_G2c_cells_only.csv"), row.names = FALSE)
  cat(sprintf("PH-1 G-2 (c) con le sole regioni con cellule: %d/%d aree entro la tolleranza (rapporto max %.3f); regioni senza cellule in %d aree, area totale %.1f um2\n",
              sum(tab$pass), nrow(tab), max(tab$c_ratio_cells_only), sum(tab$n_regions_no_cells > 0), sum(tab$area_no_cells_um2)))
}
