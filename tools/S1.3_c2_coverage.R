# tools/S1.3_c2_coverage.R — POST HOC dichiarato (2026-10-08): verifica della partizione senza operazioni d'insieme GEOS
# sui casi in cui C-2 (misurato con st_union) e' fallito. Diagnosi: l'unione GEOS perde area su maschere con buchi e su
# geometrie non annodate (A5 r5 RSA rep 15: 0 coppie sovrapposte, Σ = regione a 1.5e-9 µm², unione −1 359 µm² a cascata,
# −1 443 µm² incrementale). Misure: 1e5 punti uniformi nella regione → territori che coprono ogni punto (atteso 1);
# coppie con interni sovrapposti (st_relate "2********") e loro area totale.
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.3_c2_coverage.R  → results/S1.3/S1.3_c2_coverage.csv
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
source("tools/S1.3_roi.R"); source("R/04b2_seed_centroids.R")
measure <- function(terr, region, seed = 1) {
  set.seed(seed); p <- st_sample(region, 1e5); h <- lengths(st_intersects(p, terr))
  pr <- st_relate(terr, terr, pattern = "2********")
  pairs <- do.call(rbind, lapply(seq_along(pr), function(q) { j <- pr[[q]]; j <- j[j > q]; if (length(j)) cbind(q, j) }))
  ar <- if (is.null(pairs)) 0 else vapply(seq_len(nrow(pairs)), function(r) as.numeric(st_area(st_sfc(st_intersection(terr[[pairs[r, 1]]], terr[[pairs[r, 2]]])))), 0)
  c(n_pts = length(p), cov0 = sum(h == 0), cov1 = sum(h == 1), cov2 = sum(h >= 2), n_pairs = if (is.null(pairs)) 0 else nrow(pairs), ov_area = sum(ar))
}
rows <- list()
ns <- read.csv("results/S1.3/S1.3_null_summary.csv"); bn <- ns[ns$c2_overlap > TOL | ns$c2_symdiff > TOL, ]
nf <- list.files(file.path(OUT, "null"), "_01\\.rds$", full.names = TRUE)
dn <- do.call(rbind, lapply(nf, function(f) { o <- readRDS(f); if (!is.null(o$c123_deldir)) data.frame(archetype = o$archetype, roi_id = o$roi_id, model = o$model, rep = 1L, t(o$c123_deldir)) }))
dn <- dn[dn$c2_overlap > TOL | dn$c2_symdiff > TOL, ]
rr <- read.csv("results/S1.3/S1.3_roi_real.csv"); br <- rr[rr$c2_overlap > TOL | rr$c2_symdiff > TOL, ]
jobs <- c(lapply(seq_len(nrow(bn)), function(k) list(kind = "null", engine = "geos", r = bn[k, ])),
          lapply(seq_len(nrow(dn)), function(k) list(kind = "null", engine = "deldir", r = dn[k, ])),
          lapply(seq_len(nrow(br)), function(k) list(kind = "real", engine = br$engine[k], r = br[k, ])))
res <- mclapply(jobs, function(j) {
  A <- j$r$archetype; roi <- j$r$roi_id; i <- which(rois$archetype == A & rois$roi_id == roi); a <- ARCH[ARCH$archetype == A, ]
  w <- r3_roi_window(A, roi, rois)
  if (j$kind == "null") {
    n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
    p <- gen_points(i, j$r$model, j$r$rep, w, n_obs, a); x <- p$x; y <- p$y; lab <- sprintf("%s rep %d", j$r$model, j$r$rep)
  } else { g <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", A, roi, j$r$method)))); x <- g$x; y <- g$y; lab <- j$r$method }
  tv <- ENG[[j$engine]]$tessellate_voronoi(data.frame(cell_id = seq_along(x), region_id = 1L, x = x, y = y), one_region(w), verbose = FALSE)
  data.frame(set = j$kind, engine = j$engine, archetype = A, roi_id = roi, case = lab, c2_union = j$r$c2_overlap, t(measure(tv$cell_territories, st_sfc(w$poly[[1]]))))
}, mc.cores = 30L, mc.preschedule = FALSE)
bad <- vapply(res, function(r) !is.data.frame(r), TRUE); if (any(bad)) stop("errori ", sum(bad), ": ", as.character(res[[which(bad)[1]]]))
out <- do.call(rbind, res)
# casi sintetici del test deldir (I4 seed 2, 3)
src <- readLines("R/testing/test_S1.3.R"); source("tools/S1.3_variants.R")
ct_test <- data.frame(cell_type = c("T1", "T2", "T3", "T4"), density = c(8000, 3000, 20000, 1000), min_dist_um = c(NA, NA, NA, 15))
eval(parse(text = src[grep("^comp_for <- function", src):(grep("^load_case <- function", src) - 1)]))
inp <- readRDS("/mnt/micron/geo_spatialtrans/S1.1/inputs/I4_syn600_c4.rds"); reg <- extract_regions(inp$clust, inp$pixel_size_um, cluster_col = inp$cluster_col, verbose = FALSE)
for (e in c("geos", "deldir")) for (sd in 2:3) {
  cen <- suppressWarnings(seed_centroids(reg, ct_test, comp_for(reg$region_df$cluster_id), random_seed = sd, verbose = FALSE))$centroids
  tt <- ENG[[e]]$tessellate_voronoi(cen, reg, verbose = FALSE)
  for (r in seq_len(nrow(reg$region_df))) { kk <- which(tt$territory_df$region_id == reg$region_df$region_id[r]); if (length(kk) < 2) next
    out <- rbind(out, data.frame(set = "syn", engine = e, archetype = "-", roi_id = sprintf("I4_s%d_reg%d", sd, reg$region_df$region_id[r]), case = "test", c2_union = NA,
                                 t(measure(tt$cell_territories[kk], reg$region_polygons[r])))) }
}
write.csv(out, "results/S1.3/S1.3_c2_coverage.csv", row.names = FALSE)
print(aggregate(cbind(n = 1, cov0, cov1, cov2, n_pairs, ov_area) ~ set + engine, out, sum))
