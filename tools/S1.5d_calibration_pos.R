# tools/S1.5d_calibration_pos.R - S1.5d: semina per la calibrazione dei campi dei vuoti e della densita' (docs/steps/S1.5d.md).
# Per ROI i (30), combinazione della griglia e replica r (2): seed_centroids() a compartimenti, compartimenti ricalibrati della piega
# lascia-fuori della ROI, distanza per nucleo, modello del gruppo (DEC-020; la taglia non entra nella semina), eterogeneita' della combinazione;
# seme 20261025 + 1000 * i + r (diverso dal banco). Uscita: un parquet per unita' (x, y) in /mnt/micron/geo_spatialtrans/S1.5d/calpos/<combo>/,
# rilanciabile; n_failed in un csv per unita'. Le metriche: tools/S1.5d_calibration_pos.py.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5d_calibration_pos.R <file delle combinazioni>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5d_bench_lib.R")
options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")   # S1.6: predefiniti di S1.5c-S1.5e, per riprodurre i numeri di questo step
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); combos <- read.csv(args[1], stringsAsFactors = FALSE)
OUT <- "/mnt/micron/geo_spatialtrans/S1.5d/calpos"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
wk_file <- "/mnt/micron/geo_spatialtrans/S1.5d/WORKERS"; if (!file.exists(wk_file)) writeLines("20", wk_file)
NW <- function() min(28L, as.integer(readLines(wk_file)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
jobs <- expand.grid(i = rois$roi_idx, r = 1:2, k = seq_len(nrow(combos)), stringsAsFactors = FALSE)
if ("archetypes" %in% names(combos)) {                              # estensione della griglia solo per i tessuti indicati («A1|A3»)
  ok <- mapply(function(i, k) { a <- combos$archetypes[k]; is.na(a) || !nzchar(a) || rois$archetype[i] %in% strsplit(a, "|", fixed = TRUE)[[1]] }, jobs$i, jobs$k)
  jobs <- jobs[ok, ]
}
jobs <- jobs[order(-(rois$archetype[jobs$i] %in% c("A4", "A1")), jobs$k, jobs$i, jobs$r), ]
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(q) tryCatch({
  i <- jobs$i[q]; r <- jobs$r[q]; cb <- combos[jobs$k[q], ]; A <- rois$archetype[i]; roi <- rois$roi_id[i]
  d <- file.path(OUT, cb$combo); f <- file.path(d, sprintf("%s_%s_r%d.parquet", A, roi, r)); if (file.exists(f)) return("ok")
  dir.create(d, showWarnings = FALSE)
  w <- r3_roi_window(A, roi, rois); reg <- w$reg
  sp <- tissue_preset_spec(PRESET[[A]], size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus")
  comp <- data.frame(cluster_id = reg$region_df$cluster_id[1], cell_type = sp$composition$cell_type, fraction = sp$composition$fraction)
  h <- data.frame(void_frac = cb$void_frac, void_ell_um = cb$void_ell_um, dens_sigma = cb$dens_sigma, dens_ell_um = cb$dens_ell_um)
  sc <- suppressWarnings(seed_centroids(reg, sp$cell_types, comp, random_seed = 20261025 + 1000 * i + r, spatial = "compartments",
                                        compartments = compart_for_d(A, roi), heterogeneity = h, verbose = FALSE))
  write_parquet(data.frame(x = sc$centroids$x, y = sc$centroids$y, n_failed = sc$info$n_failed,
                           n_exceed = sc$info$compartments$n_exceed_wmax), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f); "ok"
}, error = function(e) paste("ERR", conditionMessage(e))), mc.cores = NW(), mc.preschedule = FALSE)
res <- unlist(res); cat(sprintf("[calpos] %d unita', %d errori, %.1f min\n", length(res), sum(grepl("^ERR", res)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (any(grepl("^ERR", res))) print(head(cbind(jobs[grepl("^ERR", res), ], err = res[grepl("^ERR", res)]), 10))
