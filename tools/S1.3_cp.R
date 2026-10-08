# tools/S1.3_cp.R — S1.3: CP-1 (testo del design vs D-S1.3.1), CP-2 (potenza del check di forma), D-1 (BL-056), D-3.
# Pre-registrazione: results/S1.3/S1.3_preregistration.md (commit e2f32fc). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.3_cp.R <cp1|cp2|d1|c10|all>
# Definizioni operative: motore geos; CP-1 e D-3 su I1–I5 (extract_regions di default) con catalogo/composizione di test
# di R/testing/test_S1.2.R, seed 1..5; CP-2 interne = territori non ritagliati (regione = quadrato = finestra);
# eccentricita'/orientazione dai momenti del territorio (r3_poly_moments/r3_shape); D-1 seed 20261008 + 60000 + replica.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
source("tools/R3_voronoi_metrics.R"); source("tools/S1.3_mutants.R")
STAGE <- commandArgs(TRUE)[1]; RES <- "results/S1.3"; dir.create(RES, showWarnings = FALSE)
IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"; NCORES <- 25L; BASE_SEED <- 20261008
ct_test <- data.frame(cell_type = c("T1", "T2", "T3", "T4"), density = c(8000, 3000, 20000, 1000), min_dist_um = c(NA, NA, NA, 15))
comp_for <- function(cluster_ids) {
  base <- list(`1` = data.frame(cell_type = c("T1", "T2", "T3"), fraction = c(.92, .05, .03)),
               `2` = data.frame(cell_type = c("T2", "T3"), fraction = c(.95, .05)),
               `3` = data.frame(cell_type = c("T4", "T2"), fraction = c(.85, .15)),
               `4` = data.frame(cell_type = "T3", fraction = 1))
  ids <- sort(unique(as.character(cluster_ids)))
  do.call(rbind, lapply(seq_along(ids), function(k) {
    b <- if (k <= 4) base[[k]] else data.frame(cell_type = "T2", fraction = 1)
    data.frame(cluster_id = ids[k], b)
  }))
}
rect_region <- function(x0, y0, x1, y1, id = 1L) list(
  region_df = data.frame(region_id = id, cluster_id = 1L, area_um2 = (x1 - x0) * (y1 - y0)),
  region_polygons = st_sfc(st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))))
shape_of <- function(terr) {
  mom <- t(vapply(terr, function(p) { m <- unclass(p)[[1]]; m <- m[-nrow(m), , drop = FALSE]; r3_poly_moments(m[, 1], m[, 2]) }, numeric(6)))
  r3_shape(mom[, "mxx"], mom[, "myy"], mom[, "mxy"])
}

if (STAGE %in% c("cp1", "all")) {
  E3 <- tv_mutate(tv_env(), "M3")
  syn <- c("I1_syn600_c1", "I2_syn600_c2", "I3_syn600_c3", "I4_syn600_c4", "I5_syn600_c2_labels")
  jobs <- expand.grid(id = syn, seed = 1:5, stringsAsFactors = FALSE)
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    inp <- readRDS(file.path(IN, paste0(j$id, ".rds")))
    reg <- extract_regions(inp$clust, inp$pixel_size_um, cluster_col = inp$cluster_col, verbose = FALSE)
    cen <- suppressWarnings(seed_centroids(reg, ct_test, comp_for(reg$region_df$cluster_id), random_seed = j$seed, verbose = FALSE))$centroids
    a <- tessellate_voronoi(cen, reg, verbose = FALSE); b <- E3$tessellate_voronoi(cen, reg, verbose = FALSE)
    clip_frac <- tapply(a$territory_df$clipped, factor(a$territory_df$region_id, levels = reg$region_df$region_id), mean)
    data.frame(input = j$id, seed = j$seed, region_id = a$region_check$region_id, area_um2 = a$region_check$area_um2,
               n_cells = a$region_check$n_cells,
               uncovered_pkg = (a$region_check$area_um2 - a$region_check$sum_territory) / a$region_check$area_um2,
               uncovered_design = (b$region_check$area_um2 - b$region_check$sum_territory) / b$region_check$area_um2,
               clipped_frac = as.numeric(clip_frac))
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) !is.data.frame(r), TRUE); if (any(bad)) stop("cp1: errori ", sum(bad), ": ", as.character(res[[which(bad)[1]]]))
  out <- do.call(rbind, res); write.csv(out, file.path(RES, "S1.3_cp1_regions.csv"), row.names = FALSE)
  k <- out$n_cells >= 3 & out$n_cells <= 1000
  cat(sprintf("CP-1: regioni 3<=n<=1000: %d; design > 1%%: %.3f; pacchetto > 1e-9: %d\n", sum(k),
              mean(out$uncovered_design[k] > 0.01), sum(abs(out$uncovered_pkg[out$n_cells > 0]) > 1e-9)))
}

if (STAGE %in% c("cp2", "all")) {
  L <- 1000; rho <- 1e6 / 1012.1; ks <- c(1, 1.25, 1.5, 2, 3)
  jobs <- expand.grid(ki = seq_along(ks), rp = 1:20)
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    k <- ks[j$ki]; seed <- BASE_SEED + 50000 + 100 * j$ki + j$rp
    o <- suppressWarnings(seed_centroids(rect_region(0, 0, k * L, L), data.frame(cell_type = "c", density = rho / k, min_dist_um = 2.5),
                                         data.frame(cluster_id = 1L, cell_type = "c", fraction = 1), random_seed = seed, verbose = FALSE))
    cen <- o$centroids; cen$x <- cen$x / k
    tv <- tessellate_voronoi(cen, rect_region(0, 0, L, L), verbose = FALSE)
    i <- !tv$territory_df$clipped; sh <- shape_of(tv$cell_territories[i])
    data.frame(k = k, rep = j$rp, seed = seed, n = nrow(cen), n_interior = sum(i), median_ecc_T = median(sh$ecc),
               frac_near_y = mean(abs(sh$theta) > 3 * pi / 8))
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) !is.data.frame(r), TRUE); if (any(bad)) stop("cp2: errori ", sum(bad), ": ", as.character(res[[which(bad)[1]]]))
  out <- do.call(rbind, res); write.csv(out, file.path(RES, "S1.3_cp2_anisotropy.csv"), row.names = FALSE)
  print(aggregate(cbind(median_ecc_T, frac_near_y) ~ k, out, median))
}

if (STAGE %in% c("d1", "all")) {
  rho <- 12419; s <- 1e3 / sqrt(rho)
  reg <- list(region_df = data.frame(region_id = 1:2, cluster_id = 1L, area_um2 = c(250000, 250000)),
              region_polygons = st_sfc(st_polygon(list(rbind(c(0, 0), c(500, 0), c(500, 500), c(0, 500), c(0, 0)))),
                                       st_polygon(list(rbind(c(500, 0), c(1000, 0), c(1000, 500), c(500, 500), c(500, 0))))))
  res <- mclapply(1:20, function(rp) {
    seed <- BASE_SEED + 60000 + rp
    cen <- suppressWarnings(seed_centroids(reg, data.frame(cell_type = "c", density = rho), data.frame(cluster_id = 1L, cell_type = "c", fraction = 1),
                                           random_seed = seed, verbose = FALSE))$centroids
    tv <- tessellate_voronoi(cen, reg, verbose = FALSE); td <- tv$territory_df; td$x <- cen$x[order(cen$cell_id)]
    band <- abs(td$x - 500) <= s; core <- abs(td$x - 250) <= 100 | abs(td$x - 750) <= 100
    m <- function(r, sel) mean(td$territory_area[td$region_id == r & sel])
    data.frame(rep = rp, seed = seed, band_r1 = m(1, band), band_r2 = m(2, band), core_r1 = m(1, core), core_r2 = m(2, core),
               ratio_band = m(2, band) / m(1, band))
  }, mc.cores = 20L)
  out <- do.call(rbind, res); write.csv(out, file.path(RES, "S1.3_d1_boundary.csv"), row.names = FALSE)
  cat(sprintf("D-1: rapporto mediano (regione 2 / regione 1, fascia di %.2f um) = %.3f\n", s, median(out$ratio_band)))
}

if (STAGE %in% c("c10", "all")) {          # C-10a sugli ingressi del test (sintetici seed 1-3, avversari seed 1), cella per cella
  source("tools/S1.3_variants.R"); source("tools/S1.3_arbiter.R")
  EG <- tv_variant(tv_env(), "geos"); ED <- tv_variant(tv_env(), "deldir")
  syn <- c("I1_syn600_c1", "I2_syn600_c2", "I3_syn600_c3", "I4_syn600_c4", "I5_syn600_c2_labels")
  adv <- sub("\\\\.rds$", "", list.files(IN, pattern = "^adv_.*\\\\.rds$"))
  jobs <- rbind(expand.grid(id = syn, seed = 1:3, stringsAsFactors = FALSE), data.frame(id = adv, seed = 1))
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    inp <- readRDS(file.path(IN, paste0(j$id, ".rds"))); isadv <- grepl("^adv_", j$id)
    reg <- if (isadv) extract_regions(inp$clust, inp$pixel_size_um, min_region_area_um2 = 0, cluster_col = inp$cluster_col, simplify_tol_um = 0, stride = 1L, verbose = FALSE) else
      extract_regions(inp$clust, inp$pixel_size_um, cluster_col = inp$cluster_col, verbose = FALSE)
    cen <- suppressWarnings(seed_centroids(reg, ct_test, comp_for(reg$region_df$cluster_id), random_seed = j$seed, verbose = FALSE))$centroids[, c("cell_id", "region_id", "x", "y")]
    if (isadv) { empty <- which(!reg$region_df$region_id %in% cen$region_id)       # come R/testing/test_S1.3.R
      if (length(empty)) { pos <- st_coordinates(st_point_on_surface(reg$region_polygons[empty]))
        cen <- rbind(cen, data.frame(cell_id = nrow(cen) + seq_along(empty), region_id = reg$region_df$region_id[empty], x = pos[, 1], y = pos[, 2])) } }
    g <- EG$tessellate_voronoi(cen, reg, keep_tiles = TRUE, verbose = FALSE); d <- ED$tessellate_voronoi(cen, reg, keep_tiles = TRUE, verbose = FALSE)
    ga <- g$territory_df; da <- d$territory_df
    nsd <- function(t) vapply(t, function(p) nrow(unclass(p)[[1]]) - 1L, 1L)
    k <- !ga$clipped & !da$clipped
    ab <- arbitrate(g$tiles, d$tiles, cen$x[order(cen$cell_id)], cen$y[order(cen$cell_id)], k)
    data.frame(set = "syn", archetype = "-", roi_id = sprintf("%s_s%d", j$id, j$seed), model = "test",
               rel_area = max(abs(ga$territory_area - da$territory_area) / da$territory_area),
               interior_identical = identical(ga$clipped, da$clipped),
               frag_identical = identical(ga$n_pieces_lost, da$n_pieces_lost) && identical(ga$n_pieces_gained, da$n_pieces_gained),
               nsides_interior_identical = identical(nsd(g$tiles)[k], nsd(d$tiles)[k]),
               t_geos = g$info$elapsed_s, t_deldir = d$info$elapsed_s, snapped_geos = g$info$n_snapped, snapped_deldir = d$info$n_snapped,
               n_discord = nrow(ab), n_geos_ok = sum(ab$geos_ok), n_deldir_ok = sum(ab$deldir_ok))
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) !is.data.frame(r), TRUE); if (any(bad)) stop("c10: errori ", sum(bad), ": ", as.character(res[[which(bad)[1]]]))
  out <- do.call(rbind, res); write.csv(out, file.path(RES, "S1.3_c10a_syn.csv"), row.names = FALSE)
  cat(sprintf("C-10a sintetici: %d insiemi, max rel area %.2e\n", nrow(out), max(out$rel_area)))
}
