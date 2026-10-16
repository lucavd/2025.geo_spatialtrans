# tools/S1.4c_bio.R — S1.4c: valori per docs/BIO_REFERENCES.md (B-121...), da grezzi di S1.4c. Uso: Rscript --vanilla tools/S1.4c_bio.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(arrow))
OUT <- "/mnt/micron/geo_spatialtrans/S1.4c"; RES <- "results/S1.4c"
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
cp0 <- read.csv(file.path(RES, "S1.4c_cp0.csv")); rf <- read.csv(file.path(RES, "S1.4c_real_refs.csv"))
prim <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
out <- do.call(rbind, lapply(paste0("A", 1:6), function(A) {
  rr <- rois$roi_id[rois$archetype == A]
  v <- sapply(c("cellpose_rgb", "stardist_he"), function(sg) { z <- lapply(rr, function(r) { d <- as.data.frame(read_parquet(file.path(OUT, "cp2", sprintf("%s_%s_%s_cp2.parquet", A, r, sg)))); d[d$ev, ] })
    c(s1 = mean(vapply(z, function(d) mean(d$n_c2 > 0), 0)), s1_min = min(vapply(z, function(d) mean(d$n_c2 > 0), 0)), s1_max = max(vapply(z, function(d) mean(d$n_c2 > 0), 0)),
      sol = median(unlist(lapply(z, function(d) d$sol)))) })
  nc <- unlist(lapply(rr, function(r) readRDS(file.path(OUT, "prep", sprintf("%s_%s_%s.rds", A, r, prim[[A]])))$nc_ev))
  data.frame(archetype = A, s1_cellpose = v["s1", 1], s1_cellpose_range = sprintf("%.3f–%.3f", v["s1_min", 1], v["s1_max", 1]),
             s1_stardist = v["s1", 2], s1_stardist_range = sprintf("%.3f–%.3f", v["s1_min", 2], v["s1_max", 2]),
             sol_cellpose = v["sol", 1], sol_stardist = v["sol", 2], s1_shape_ev = mean(rf$s1[rf$archetype == A]),
             auc_cp0 = cp0$auc[cp0$archetype == A], nc_q99 = unname(quantile(nc, 0.99, type = 7)), nc_n = length(nc)) }))
write.csv(out, file.path(RES, "S1.4c_bio_refs.csv"), row.names = FALSE); print(format(out, digits = 3), row.names = FALSE)
