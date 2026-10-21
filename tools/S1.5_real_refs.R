# tools/S1.5_real_refs.R - S1.5: riferimenti reali per ROI usati per gli intervalli delle previsioni (prima di ogni simulazione).
# Stesse definizioni di S1.3 (CV globale ed eccentricita' dei territori delle celle interne, r3_summary/r3_shape) e di
# S1.4d (insieme di valutazione interior & matched & a_n > 0; rapporto a_n / a_t; quasi senza citoplasma = a_n / a_t >= 0.95).
# Ingressi: territori del pacchetto sui nuclei reali (S1.4/real: cellpose_rgb in A1, A2, A3, A5) e maschera unica di R3b
# per A4, A6 (R3b/S1.4/real *_consensus.rds; nuclei di forma StarDist per taglia e N/C, DEC-006).
# Uscita: results/S1.5/S1.5_real_refs_roi.csv.  Uso: Rscript --vanilla tools/S1.5_real_refs.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(arrow) })
source("tools/R3_voronoi_metrics.R")
M <- "/mnt/micron/geo_spatialtrans"
rois <- read.csv("results/R2/R2_rois_checked.csv")
src <- data.frame(archetype = paste0("A", 1:6), method = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "consensus", "cellpose_rgb", "consensus"),
                  base = file.path(M, c("S1.4", "S1.4", "S1.4", "R3b/S1.4", "S1.4", "R3b/S1.4")), stringsAsFactors = FALSE)
ecc_of <- function(g) vapply(seq_along(g), function(i) { xy <- st_coordinates(g[[i]])[, 1:2]
  m <- r3_poly_moments(xy[, 1], xy[, 2]); r3_shape(m[["mxx"]], m[["myy"]], m[["mxy"]])$ecc }, 0)
out <- list()
for (k in seq_len(nrow(src))) {
  A <- src$archetype[k]; ids <- rois$roi_id[rois$archetype == A]
  for (r in ids) {
    f <- file.path(src$base[k], "real", sprintf("%s_%s_%s.rds", A, r, src$method[k])); z <- readRDS(f); cc <- z$cells
    if (src$method[k] == "consensus") {
      g <- as.data.frame(read_parquet(file.path(M, "R3b", "R3", "real", sprintf("%s_%s_consensus_gen.parquet", A, r))))
      mp <- as.data.frame(read_parquet(file.path(M, "R3b", "masks", sprintf("%s_%s_consensus_map.parquet", A, r))))
      cc$matched <- mp$src[match(g$label, mp$label)] == "sd"
    }
    ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
    it <- cc$interior; a <- cc$a_t[it]; nc <- cc$a_n[ev] / cc$a_t[ev]
    ecc <- ecc_of(z$terr[which(it)])
    out[[length(out) + 1]] <- data.frame(archetype = A, roi_id = r, method = src$method[k], n = z$n, n_interior = sum(it), n_eval = sum(ev),
      density_seg_mm2 = z$n / z$area_window * 1e6, median_area_nuc = median(cc$a_n[ev]), median_nc = median(nc),
      frac_nc095 = mean(nc >= 0.95), cv_terr = sd(a) / mean(a), median_ecc_T = median(ecc), median_eq_r = median(sqrt(a / pi)),
      md5 = unname(tools::md5sum(f)))
  }
}
res <- do.call(rbind, out)
dir.create("results/S1.5", showWarnings = FALSE)
write.csv(res, "results/S1.5/S1.5_real_refs_roi.csv", row.names = FALSE)
print(res[, setdiff(names(res), "md5")], digits = 4)
