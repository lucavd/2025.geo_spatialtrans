# tools/S1.5e_calibration_pos.R - S1.5e: semina per la calibrazione dei vuoti (X, XA) e della densita' con il tetto (XAD) (docs/steps/S1.5e.md).
# Come tools/S1.5d_calibration_pos.R: per ROI i (30), combinazione e replica r (2): seed_centroids() a compartimenti, compartimenti ricalibrati della
# piega (S1.5d), distanza per nucleo, eterogeneita' della combinazione con le colonne void_mode, void_aspect, dens_cap; seme 20261026 + 1000 * i + r
# (diverso dal banco e da S1.5d). Colonne facoltative del file delle combinazioni: archetypes («A1|A3»), rois («r1|r2»: solo queste ROI; XAD,
# vuoti della piega). Uscita: un parquet per unita' (x, y, n_failed, n_exceed) in /mnt/micron/geo_spatialtrans/S1.5e/calpos/<combo>/, rilanciabile.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5e_calibration_pos.R <file delle combinazioni>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5d_bench_lib.R")
options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")   # S1.6: predefiniti di S1.5c-S1.5e, per riprodurre i numeri di questo step
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); combos <- read.csv(args[1], stringsAsFactors = FALSE)
OUT <- "/mnt/micron/geo_spatialtrans/S1.5e/calpos"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
wk_file <- "/mnt/micron/geo_spatialtrans/S1.5e/WORKERS"; if (!file.exists(wk_file)) writeLines("22", wk_file)
NW <- function() min(28L, as.integer(readLines(wk_file)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
jobs <- expand.grid(i = rois$roi_idx, r = 1:2, k = seq_len(nrow(combos)), stringsAsFactors = FALSE)
flt <- function(col, val) { a <- if (col %in% names(combos)) combos[[col]] else rep(NA, nrow(combos))
  mapply(function(i, k) { z <- a[k]; is.na(z) || !nzchar(z) || val(i) %in% strsplit(z, "|", fixed = TRUE)[[1]] }, jobs$i, jobs$k) }
jobs <- jobs[flt("archetypes", function(i) rois$archetype[i]) & flt("rois", function(i) rois$roi_id[i]), ]
jobs <- jobs[order(-(rois$archetype[jobs$i] %in% c("A4", "A1")), jobs$k, jobs$i, jobs$r), ]
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(q) tryCatch({
  i <- jobs$i[q]; r <- jobs$r[q]; cb <- combos[jobs$k[q], ]; A <- rois$archetype[i]; roi <- rois$roi_id[i]
  d <- file.path(OUT, cb$combo); f <- file.path(d, sprintf("%s_%s_r%d.parquet", A, roi, r)); if (file.exists(f)) return("ok")
  dir.create(d, showWarnings = FALSE)
  w <- r3_roi_window(A, roi, rois); reg <- w$reg
  sp <- tissue_preset_spec(PRESET[[A]], size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus")
  comp <- data.frame(cluster_id = reg$region_df$cluster_id[1], cell_type = sp$composition$cell_type, fraction = sp$composition$fraction)
  h <- data.frame(void_frac = cb$void_frac, void_ell_um = cb$void_ell_um, dens_sigma = cb$dens_sigma, dens_ell_um = cb$dens_ell_um,
                  void_mode = cb$void_mode, void_aspect = cb$void_aspect, dens_cap = as.logical(cb$dens_cap))
  sc <- suppressWarnings(seed_centroids(reg, sp$cell_types, comp, random_seed = 20261026 + 1000 * i + r, spatial = "compartments",
                                        compartments = compart_for_d(A, roi), heterogeneity = h, verbose = FALSE))
  write_parquet(data.frame(x = sc$centroids$x, y = sc$centroids$y, n_failed = sc$info$n_failed,
                           n_exceed = sc$info$compartments$n_exceed_wmax), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f); "ok"
}, error = function(e) paste("ERR", conditionMessage(e))), mc.cores = NW(), mc.preschedule = FALSE)
res <- unlist(res); cat(sprintf("[calpos] %d unita', %d errori, %.1f min\n", length(res), sum(grepl("^ERR", res)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (any(grepl("^ERR", res))) print(head(cbind(jobs[grepl("^ERR", res), ], err = res[grepl("^ERR", res)]), 10))
