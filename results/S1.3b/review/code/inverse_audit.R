
# review/code/inverse_audit.R — (1) dallo stato annodato, 1 ulp su vertici condivisi: l'errore ricompare? (2) celle invalide V1/V3 (P2)
RV <- "/mnt/micron/geo_spatialtrans/S1.3b/review_code"
args <- c("none", file.path(RV, "s13out")); commandArgs <- function(trailingOnly = TRUE) args
suppressMessages(source("tools/S1.3_roi.R")); rm(commandArgs); source("tools/S1.3b_overlay.R")
src <- readLines("tools/S1.3b_run.R"); eval(parse(text = src[grep("^regen <- function", src):(grep("^reg_features <- function", src) - 1)]))
suppressPackageStartupMessages(library(parallel))
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)
ulp <- function(x) 2^(floor(log2(pmax(abs(x), 2^-1022))) - 52)
inter_area <- function(a, b) sum(as.numeric(st_area(st_intersection(st_sfc(a), st_sfc(b)))))
setv <- function(g, df) s3b_rebuild(df, class(g)[2])$geoms[[1]]
P1 <- read.csv("results/S1.3b/S1.3b_posthoc_P1_pairs_v0.csv")
one <- function(q) { cs <- P1[q, ]; G <- regen(roi_i(cs$archetype, cs$roi_id), cs$model, cs$rep); ti <- G$terr[[cs$i]]; tj <- G$terr[[cs$j]]
  vi <- s3b_geoms_df(list(ti)); vj <- s3b_geoms_df(list(tj)); D <- sqrt(outer(vi$x, vj$x, "-")^2 + outer(vi$y, vj$y, "-")^2)
  w <- which(D <= 1e-6, arr.ind = TRUE); dif <- w[D[w] > 0, , drop = FALSE]; idn <- w[D[w] == 0, , drop = FALSE]
  vf <- vj; vf$x[dif[, 2]] <- vi$x[dif[, 1]]; vf$y[dif[, 2]] <- vi$y[dif[, 1]]; tf <- setv(tj, vf)      # stato annodato
  rb <- st_boundary(G$reg); rows <- list()
  for (kind in c("boundary", "interior")) { b <- if (kind == "boundary") dif[1, 2] else idn[1, 2]
    for (ax in c("x", "y")) for (s in c(-1, 1)) for (nu in c(1, 4)) { vv <- vf; vv[[ax]][b] <- vv[[ax]][b] + s * nu * ulp(vv[[ax]][b])
      g2 <- setv(tj, vv); rows[[length(rows) + 1]] <- data.frame(pair = paste(cs$archetype, cs$roi_id, cs$model, cs$rep, cs$i, cs$j), kind = kind, axis = ax, sign = s, nulp = nu,
        dist_reg = as.numeric(st_distance(st_sfc(st_point(c(vf$x[b], vf$y[b]))), rb)), area = inter_area(ti, g2), area_i = as.numeric(st_area(st_sfc(ti)))) } }
  r <- do.call(rbind, rows); r$a_fixed <- inter_area(ti, tf); r }
R <- do.call(rbind, mclapply(seq_len(nrow(P1)), one, mc.cores = 10, mc.preschedule = FALSE))
write.csv(R, "results/S1.3b/review/code/inverse_audit.csv", row.names = FALSE)
R$err <- R$area > 1e-6
cat("a_fixed max:", max(R$a_fixed), "\n")
print(aggregate(err ~ kind + nulp, R, function(x) sprintf("%d/%d", sum(x), length(x))))
print(aggregate(err ~ pair + kind, R, sum))
cat("errors reproduce full cell i area:", sum(R$err & abs(R$area - R$area_i) < 1e-6), "of", sum(R$err), "\n")
# (2) celle invalide di P2
P2 <- read.csv("results/S1.3b/S1.3b_posthoc_P2_invalid.csv"); P2 <- unique(P2[, c("archetype", "roi_id", "model", "rep", "cell")])
for (q in seq_len(nrow(P2))) { cs <- P2[q, ]; G <- regen(roi_i(cs$archetype, cs$roi_id), cs$model, cs$rep); t <- G$terr[cs$cell]
  m <- st_coordinates(t)[, 1:2]; m <- m[-nrow(m), , drop = FALSE]; n <- nrow(m)
  Dm <- as.matrix(dist(m)); for (k in 1:n) { Dm[k, k] <- Inf; Dm[k, (k %% n) + 1] <- Inf; Dm[(k %% n) + 1, k] <- Inf }
  kk <- which(Dm == min(Dm), arr.ind = TRUE)[1, ]
  V1 <- s3b_variant(G$terr, G$reg, "V1"); V2 <- s3b_variant(G$terr, G$reg, "V2"); V3 <- s3b_variant(G$terr, G$reg, "V3")
  t2 <- st_as_sfc(as_wkb(V2$terr[cs$cell]))
  cat(sprintf("\n%s %s %s %d cell %d: n vert %d; min non-consecutive vertex distance %.3g (vertices %d,%d: %.17g %.17g / %.17g %.17g); dist to region boundary %.3g\n",
    cs$archetype, cs$roi_id, cs$model, cs$rep, cs$cell, n, min(Dm), kk[1], kk[2], m[kk[1], 1], m[kk[1], 2], m[kk[2], 1], m[kk[2], 2],
    as.numeric(st_distance(st_sfc(st_point(m[kk[1], ])), st_boundary(G$reg)))))
  cat(sprintf("  V0 valid %s; V1 valid %s (%s); V3 valid %s; V2 valid %s type %s n parts %d; area V0 %.6f V2 %.6f; V2 n vert %d\n",
    st_is_valid(t), st_is_valid(V1$terr[cs$cell]), st_is_valid(V1$terr[cs$cell], reason = TRUE), st_is_valid(V3$terr[cs$cell]), st_is_valid(t2),
    as.character(st_geometry_type(t2)), length(st_cast(t2, "POLYGON")), as.numeric(st_area(t)), as.numeric(st_area(t2)), nrow(st_coordinates(t2))))
}
