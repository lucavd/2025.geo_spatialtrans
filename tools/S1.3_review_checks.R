# tools/S1.3_review_checks.R — POST HOC dichiarati, richiesti dalla revisione avversariale (2026-10-09)
#  (1) RA-claims-11: attribuzione di C-6 — repliche nulle ritassellate con la regola dei frammenti di R3 (variante r3rule):
#      se la causa e' D-S1.3.2, tornano identiche a R3 (resta solo l'effetto dei lati cortissimi di deldir in R3)
#  (2) RA-claims-12: A1 r5 — perche' anche la variante deldir differisce da R3 sulle celle senza frammenti (riquadro rw?)
#  (3) RA-claims-18: la coppia deldir di I4 (36.6 µm² di "sovrapposizione") verificata con i predicati
# Uso: R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.3_review_checks.R [outdir_prereview]
args <- c("none", commandArgs(TRUE)[1]); if (is.na(args[2])) args[2] <- "/mnt/micron/geo_spatialtrans/S1.3_prereview"
OUT0 <- args[2]; commandArgs <- function(trailingOnly = TRUE) args
source("tools/S1.3_roi.R"); source("R/04b2_seed_centroids.R")
ENG$r3rule <- tv_variant(tv_env(), "r3rule")
RES <- "results/S1.3"
# (1)
c6 <- read.csv(file.path(RES, "S1.3_c6_prereview.csv")); nsm <- read.csv(file.path(RES, "S1.3_null_summary_prereview.csv"))
m <- merge(c6, nsm[, c("archetype", "roi_id", "model", "rep", "n_fragments")])
set.seed(20261009)
pick <- rbind(m[m$n_fragments > 0, ][sample(sum(m$n_fragments > 0), 40), ], m[m$n_fragments == 0 & m$rel > TOL, ])
r3n <- read.csv("results/R3/R3_null_summary.csv")
cols <- c("n", "n_interior", "frac_interior", "mean_area", "median_area", "cv", "median_eq_r", "median_ecc_T", "mean_nsides", "var_nsides", "var_norm_area", "cv_loc")
one <- function(k) {
  b <- pick[k, ]; A <- b$archetype; roi <- b$roi_id; i <- which(rois$archetype == A & rois$roi_id == roi); a <- ARCH[ARCH$archetype == A, ]
  w <- r3_roi_window(A, roi, rois); n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  p <- gen_points(i, b$model, b$rep, w, n_obs, a); li <- r3_local_intensity(p$x, p$y, w)
  res <- lapply(c("r3rule", "deldir_r3"), function(v) {
    if (v == "r3rule") { s <- s13_tessellate(p$x, p$y, w, backend = "r3rule"); d <- s$df } else d <- r3_tessellate(p$x, p$y, w$poly, w$frame)
    sm <- r3_summary(d, li$lambda_loc); ref <- unlist(r3n[r3n$archetype == A & r3n$roi_id == roi & r3n$model == b$model & r3n$rep == b$rep, cols])
    rel <- abs(sm[cols] - ref) / pmax(abs(ref), 1e-300)
    data.frame(variant = v, max_rel = max(rel), worst = cols[which.max(rel)])
  })
  data.frame(archetype = A, roi_id = roi, model = b$model, rep = b$rep, n_fragments = b$n_fragments, rel_s13 = b$rel,
             r3rule_max_rel = res[[1]]$max_rel, r3rule_worst = res[[1]]$worst, r3code_max_rel = res[[2]]$max_rel)
}
ENG$r3rule$.tv_engine <- function() "r3rule"
o1 <- do.call(rbind, mclapply(seq_len(nrow(pick)), function(k) tryCatch(one(k), error = function(e) data.frame(archetype = conditionMessage(e))), mc.cores = 40L, mc.preschedule = FALSE))
write.csv(o1, file.path(RES, "S1.3_review_c6_attribution.csv"), row.names = FALSE)
cat("(1) repliche:", nrow(o1), "; identiche a R3 con regola R3 (geos):", sum(o1$r3rule_max_rel <= TOL, na.rm = TRUE), "; con codice R3 (deldir):", sum(o1$r3code_max_rel <= TOL, na.rm = TRUE), "\n")
print(table(o1$r3rule_worst[o1$r3rule_max_rel > TOL]))
# (2)
w <- r3_roi_window("A1", "r5", rois); g <- as.data.frame(read_parquet(file.path(R3, "real", "A1_r5_cellpose_rgb_gen.parquet")))
ref <- as.data.frame(read_parquet(file.path(R3, "real", "A1_r5_cellpose_rgb_cells.parquet")))
sd <- s13_tessellate(g$x, g$y, w, backend = "deldir")$df; sg <- s13_tessellate(g$x, g$y, w, backend = "geos")$df
mono <- sd$n_lost == 0 & sd$n_gained == 0
rel_d <- abs(sd$area - ref$area) / ref$area; k <- which(mono & rel_d > TOL)
rwS <- ENG$geos$.tv_window(w$poly); bb <- st_bbox(w$frame); rwR <- c(bb[["xmin"]] - 50, bb[["xmax"]] + 50, bb[["ymin"]] - 50, bb[["ymax"]] + 50)
t_S <- deldir::tile.list(deldir::deldir(g$x, g$y, rw = rwS, round = FALSE)); t_R <- deldir::tile.list(deldir::deldir(g$x, g$y, rw = rwR, round = FALSE))
aS <- vapply(t_S, function(t) t$area, 0)[order(vapply(t_S, function(t) t$ptNum, 1L))]; aR <- vapply(t_R, function(t) t$area, 0)[order(vapply(t_R, function(t) t$ptNum, 1L))]
o2 <- data.frame(cell = k, area_R3 = ref$area[k], area_deldir_S13 = sd$area[k], area_geos = sg$area[k], tile_deldir_rwS13 = aS[k], tile_deldir_rwR3 = aR[k], interior = ref$interior[k])
write.csv(o2, file.path(RES, "S1.3_review_A1r5_deldir_rw.csv"), row.names = FALSE)
cat("(2) A1 r5: celle senza frammenti con deldir S1.3 != R3:", length(k), "; tile deldir diversi fra i due riquadri (rel > 1e-9):", sum(abs(aS - aR) / aR > TOL), "\n"); print(o2, digits = 10)
# (3)
cv <- read.csv(file.path(RES, "S1.3_c2_coverage_prereview.csv")); cs <- cv[cv$set == "syn" & cv$engine == "deldir" & cv$ov_area > 1e-6, ]
ct_test <- data.frame(cell_type = c("T1", "T2", "T3", "T4"), density = c(8000, 3000, 20000, 1000), min_dist_um = c(NA, NA, NA, 15))
src <- readLines("R/testing/test_S1.3.R"); eval(parse(text = src[grep("^comp_for <- function", src):(grep("^load_case <- function", src) - 1)]))
inp <- readRDS("/mnt/micron/geo_spatialtrans/S1.1/inputs/I4_syn600_c4.rds"); reg <- extract_regions(inp$clust, inp$pixel_size_um, cluster_col = inp$cluster_col, verbose = FALSE)
o3 <- do.call(rbind, lapply(seq_len(nrow(cs)), function(q) {
  sdd <- as.integer(sub("I4_s(\\d+)_reg.*", "\\1", cs$roi_id[q])); rg <- as.integer(sub(".*_reg(\\d+)", "\\1", cs$roi_id[q]))
  cen <- suppressWarnings(seed_centroids(reg, ct_test, comp_for(reg$region_df$cluster_id), random_seed = sdd, verbose = FALSE))$centroids
  tt <- ENG$deldir$tessellate_voronoi(cen, reg, verbose = FALSE); kk <- which(tt$territory_df$region_id == rg); t <- tt$cell_territories[kk]
  pr <- st_relate(t, t, pattern = "2********"); pairs <- do.call(rbind, lapply(seq_along(pr), function(z) { j <- pr[[z]]; j <- j[j > z]; if (length(j)) cbind(z, j) }))
  ar <- vapply(seq_len(nrow(pairs)), function(r) sum(as.numeric(st_area(st_intersection(t[pairs[r, 1]], t[pairs[r, 2]])))), 0)
  big <- which(ar > 1e-6)
  do.call(rbind, lapply(big, function(r) { gi <- st_collection_extract(st_intersection(t[pairs[r, 1]], t[pairs[r, 2]]), "POLYGON"); set.seed(1); s <- st_sample(gi, 2000)
    data.frame(case = cs$roi_id[q], i = pairs[r, 1], j = pairs[r, 2], area = ar[r], n_pts = length(s),
               in_both = sum(lengths(st_intersects(s, t[pairs[r, 1]])) > 0 & lengths(st_intersects(s, t[pairs[r, 2]])) > 0)) }))
}))
write.csv(o3, file.path(RES, "S1.3_review_deldir_pair_check.csv"), row.names = FALSE); cat("(3)\n"); print(o3)
