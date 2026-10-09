# tools/S1.3_c2_pairs_check.R — POST HOC dichiarato: le "sovrapposizioni" fra coppie misurate con st_intersection nei 15 nulli GEOS
# con C-2 fallito sono reali? 2000 punti campionati nell'intersezione dichiarata, appartenenza verificata con il predicato st_intersects.

args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
source("tools/S1.3_roi.R")
ns <- read.csv("results/S1.3/S1.3_null_summary.csv"); bn <- ns[ns$c2_overlap > TOL, ]
out <- do.call(rbind, mclapply(seq_len(nrow(bn)), function(k) {
  b <- bn[k, ]; A <- b$archetype; roi <- b$roi_id; i <- which(rois$archetype == A & rois$roi_id == roi); a <- ARCH[ARCH$archetype == A, ]
  w <- r3_roi_window(A, roi, rois); n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  p <- gen_points(i, b$model, b$rep, w, n_obs, a)
  tv <- tessellate_voronoi(data.frame(cell_id = seq_along(p$x), region_id = 1L, x = p$x, y = p$y), one_region(w), verbose = FALSE); t <- tv$cell_territories
  pr <- st_relate(t, t, pattern = "2********")
  pairs <- do.call(rbind, lapply(seq_along(pr), function(q) { j <- pr[[q]]; j <- j[j > q]; if (length(j)) cbind(q, j) }))
  if (is.null(pairs)) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(pairs)), function(r) {
    g <- st_intersection(t[pairs[r, 1]], t[pairs[r, 2]]); ar <- sum(as.numeric(st_area(g)))
    if (ar < 1e-6) return(data.frame(case = paste(A, roi, b$model, b$rep), i = pairs[r, 1], j = pairs[r, 2], area = ar, n_pts = NA, in_both = NA, in_i = NA, in_j = NA,
                                     multi_i = as.character(st_geometry_type(t[pairs[r, 1]])), multi_j = as.character(st_geometry_type(t[pairs[r, 2]]))))
    set.seed(1); s <- st_sample(st_collection_extract(g, "POLYGON"), 2000)
    hi <- lengths(st_intersects(s, t[pairs[r, 1]])) > 0; hj <- lengths(st_intersects(s, t[pairs[r, 2]])) > 0
    data.frame(case = paste(A, roi, b$model, b$rep), i = pairs[r, 1], j = pairs[r, 2], area = ar, n_pts = length(s), in_both = sum(hi & hj), in_i = sum(hi), in_j = sum(hj),
               multi_i = as.character(st_geometry_type(t[pairs[r, 1]])), multi_j = as.character(st_geometry_type(t[pairs[r, 2]])))
  }))
}, mc.cores = 15L, mc.preschedule = FALSE))
print(out)
write.csv(out, "results/S1.3/S1.3_c2_overlap_pairs_check.csv", row.names = FALSE)
