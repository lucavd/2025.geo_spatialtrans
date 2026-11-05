# tools/S2.1_ph2.R - S2.1, post hoc PH-2 (dopo il calcolo): le differenze di D-1 fra simulato e segmentatore dello Step 1 si spiegano con
# densita' e area dei nuclei dello Step 1? Per area: nuclei per mm2 di tessuto e area media del nucleo, simulato (livello U6) e segmentatore
# (nuclei keep di R2 / consenso R3b, come R4). Uso (root del repo): Rscript --vanilla tools/S2.1_ph2.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(arrow) })
M <- "/mnt/micron/geo_spatialtrans"
rois <- read.csv("results/R2/R2_rois_checked.csv"); SRC <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "consensus", A5 = "cellpose_rgb", A6 = "consensus")
d1 <- read.csv("results/S2.1/S2.1_desc_sim_roi.csv"); d1r <- read.csv("results/S2.1/S2.1_desc_real_roi.csv")
out <- do.call(rbind, lapply(seq_len(nrow(rois)), function(r) {
  A <- rois$archetype[r]; tag <- paste(A, rois$roi_id[r], sep = "_")
  L <- readRDS(file.path(M, "S1.6/layers", paste0(tag, "_U6.rds"))); ta <- sum(as.numeric(st_area(L$region_polygons)))
  nu <- as.data.frame(read_parquet(if (SRC[[A]] == "consensus") file.path(M, "R3b/nuclei", paste0(tag, "_consensus_native.parquet")) else
    file.path("results/R2/nuclei", paste0(tag, "_cellpose_rgb_native.parquet")))); nu <- nu[nu$keep, ]
  data.frame(tag = tag, archetype = A, tissue_mm2 = ta / 1e6, sim_n_per_mm2 = nrow(L$cell_df) / (ta / 1e6), seg_n_per_mm2 = nrow(nu) / (ta / 1e6),
             sim_nuc_area = mean(L$cell_df$nucleus_area), seg_nuc_area = mean(nu$area_um2),
             sim_cover = sum(L$cell_df$nucleus_area) / ta, seg_cover = sum(nu$area_um2) / ta,
             sim_d1 = d1$d1_centre_in_nucleus[d1$tag == tag], seg_d1 = d1r$seg_d1_centre_in_nucleus[d1r$tag == tag])
}))
out$ratio_d1 <- out$sim_d1 / out$seg_d1; out$ratio_density <- out$sim_n_per_mm2 / out$seg_n_per_mm2; out$ratio_area <- out$sim_nuc_area / out$seg_nuc_area
out$ratio_cover <- out$sim_cover / out$seg_cover
write.csv(out, "results/S2.1/S2.1_PH2_density_area.csv", row.names = FALSE)
agg <- aggregate(out[, c("sim_n_per_mm2", "seg_n_per_mm2", "sim_nuc_area", "seg_nuc_area", "ratio_d1", "ratio_density", "ratio_area", "ratio_cover")], list(archetype = out$archetype), mean)
print(format(agg, digits = 3)); cat(sprintf("correlazione fra rapporto di D-1 e rapporto della copertura nucleare (30 aree): %.3f\n", cor(out$ratio_d1, out$ratio_cover)))
