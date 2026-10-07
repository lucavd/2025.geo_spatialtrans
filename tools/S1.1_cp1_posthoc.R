# S1.1 — CP-1 ESPLORATIVO (NON pre-registrato, aggiunto dopo il FAIL di CP-1 su I1-I4)
# Domanda: la soglia del 2 % distingue il contorno di Moore sui centri quando le regioni sono
# piccole (mappe reali a 8 µm) o solo quando sono grandi? Uso:
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.1_cp1_posthoc.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(imager); library(parallel) })
source("tools/S1.1_moore_variant.R")
IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
idx <- read.csv("results/S1.1/S1.1_inputs.csv", stringsAsFactors = FALSE)
ids <- idx$id[idx$group == "real_roi" | grepl("^I[1-4]_", idx$id)]
one <- function(id) {
  o <- readRDS(file.path(IN, paste0(id, ".rds")))
  x <- round(o$clust$x); y <- round(o$clust$y); i <- x - min(x) + 1; j <- y - min(y) + 1
  lv <- sort(unique(as.character(o$clust[[o$cluster_col]])))
  G <- matrix(0L, max(j) + 2, max(i) + 2); G[cbind(j + 1, i + 1)] <- match(as.character(o$clust[[o$cluster_col]]), lv)
  L <- as.matrix(imager::label(as.cimg(G), high_connectivity = FALSE)); fg <- G != 0L
  pix <- split(which(fg), L[fg])
  do.call(rbind, lapply(pix, function(p) {
    w <- arrayInd(p, dim(G))
    B <- matrix(FALSE, diff(range(w[, 1])) + 1, diff(range(w[, 2])) + 1)
    B[cbind(w[, 1] - min(w[, 1]) + 1, w[, 2] - min(w[, 2]) + 1)] <- TRUE
    data.frame(id = id, px_um = o$pixel_size_um, n_px = nrow(w), area_moore_px = moore_centre_area(B))
  }))
}
res <- do.call(rbind, mclapply(ids, one, mc.cores = 12))
res$area_um2 <- res$n_px * res$px_um^2
summ <- do.call(rbind, lapply(split(res, res$id), function(r) {
  big <- r$area_um2 >= 100
  data.frame(id = r$id[1], n_comp = nrow(r), relerr_all = sum(r$area_moore_px - r$n_px) / sum(r$n_px),
             n_comp_ge100 = sum(big), relerr_ge100 = sum(r$area_moore_px[big] - r$n_px[big]) / sum(r$n_px[big]))
}))
write.csv(summ, "results/S1.1/S1.1_cp1_posthoc.csv", row.names = FALSE)
print(summ, digits = 3)
