# rv_export.R — revisione avversariale S1.2, track checkB.
# Solo ESPORTAZIONE di dati (nessun calcolo di metriche): rds del progetto -> CSV leggibili da Python.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_export.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
RV <- "results/S1.2/review/checkB/data"; dir.create(RV, recursive = TRUE, showWarnings = FALSE)
OUT <- "/mnt/micron/geo_spatialtrans/S1.2/checkB"
R_GRID <- seq(0, 30, by = 0.5)
real <- readRDS(file.path(OUT, "real_pcf.rds"))
rg <- do.call(rbind, lapply(real, function(z) do.call(rbind, lapply(c("primary", "secondary"), function(role)
  data.frame(archetype = z$archetype, roi_id = z$roi_id, role = role, method = z[[role]]$method, n = z[[role]]$n,
             n_lost = z[[role]]$n_lost, lambda = z[[role]]$lambda, bw = z$bw, area_mm2 = z$area_mm2, n_regions = z$n_regions,
             r = R_GRID, g = z[[role]]$g)))))
write.csv(rg, file.path(RV, "real_pcf_export.csv"), row.names = FALSE)
fs <- list.files(file.path(OUT, "sim"), full.names = TRUE)
sim <- do.call(rbind, lapply(fs, function(f) { s <- readRDS(f); s$file <- basename(f); s }))
write.csv(sim, file.path(RV, "sim_export.csv"), row.names = FALSE)
IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
rois <- read.csv("results/R2/R2_rois_checked.csv")
meta <- list()
for (i in seq_len(nrow(rois))) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]
  o <- readRDS(file.path(IN, sprintf("roi_%s_%s.rds", A, roi)))
  cl <- o$clust; cc <- o$cluster_col
  write.csv(data.frame(x = cl$x, y = cl$y, cluster = as.character(cl[[cc]])), file.path(RV, sprintf("clust_%s_%s.csv", A, roi)), row.names = FALSE)
  meta[[i]] <- data.frame(archetype = A, roi_id = roi, pixel_size_um = o$pixel_size_um, cluster_col = cc, n_px = nrow(cl),
                          names = paste(names(o), collapse = ";"))
}
write.csv(do.call(rbind, meta), file.path(RV, "clust_meta.csv"), row.names = FALSE)
cat("export ok:", nrow(rg), "righe reale,", nrow(sim), "righe sim,", length(fs), "file sim\n")
