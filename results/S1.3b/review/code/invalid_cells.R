
RV <- "/mnt/micron/geo_spatialtrans/S1.3b/review_code"
args <- c("none", file.path(RV, "s13out")); commandArgs <- function(trailingOnly = TRUE) args
suppressMessages(source("tools/S1.3_roi.R")); rm(commandArgs); source("tools/S1.3b_overlay.R")
src <- readLines("tools/S1.3b_run.R"); eval(parse(text = src[grep("^regen <- function", src):(grep("^reg_features <- function", src) - 1)]))
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)
for (cs in list(list("A5", "r1", "CSR", 12L, 6L), list("A1", "r4", "RSArule", 13L, 1655L))) {
  G <- regen(roi_i(cs[[1]], cs[[2]]), cs[[3]], cs[[4]]); t <- G$terr[cs[[5]]]
  v <- s3b_geoms_df(list(t[[1]])); cat(cs[[1]], cs[[2]], cs[[3]], cs[[4]], cs[[5]], "type", class(t[[1]])[2], "rings", nrow(unique(v[, c("p", "r")])), "nvert", nrow(v), "\n")
  best <- c(Inf, NA, NA)
  for (a in 1:(nrow(v) - 1)) for (b in (a + 1):nrow(v)) { if (v$p[a] == v$p[b] && v$r[a] == v$r[b]) { nr <- sum(v$p == v$p[a] & v$r == v$r[a]); if (abs(v$i[a] - v$i[b]) %in% c(1, nr - 1)) next }
    d <- sqrt((v$x[a] - v$x[b])^2 + (v$y[a] - v$y[b])^2); if (d < best[1]) best <- c(d, a, b) }
  a <- best[2]; b <- best[3]
  cat(sprintf("  min non-adjacent vertex distance %.3g: (%.17g, %.17g) ring %d/%d vs (%.17g, %.17g) ring %d/%d; dist to region boundary %.3g\n", best[1], v$x[a], v$y[a], v$p[a], v$r[a], v$x[b], v$y[b], v$p[b], v$r[b],
      as.numeric(st_distance(st_sfc(st_point(c(v$x[a], v$y[a]))), st_boundary(G$reg)))))
  V1 <- s3b_variant(G$terr, G$reg, "V1"); cat("  V1 reason:", st_is_valid(V1$terr[cs[[5]]], reason = TRUE), "\n")
}
