
# review/code/export_targets.R — revisore "code": per 30 fold x 10 semi esporta le aree estratte (stesso flusso del banco)
# e i dati reali necessari a ricalcolare B-1, B-2, B-3, B-4, B-7 di V-a/V-b con algoritmi di non sovrapposizione alternativi.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R"); source("tools/S1.4_mutants.R")
E <- dcg_env(); OD <- "/tmp/code_tgt"; dir.create(OD, showWarnings = FALSE)
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"; S14B <- "/mnt/micron/geo_spatialtrans/S1.4b"
ARCH <- data.frame(archetype = paste0("A", 1:6), primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
one <- function(A, k) { method <- ARCH$primary[ARCH$archetype == A]
  o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, k))); b <- readRDS(file.path(S14B, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  te <- readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, o$roi_test, method))); cc <- te$cells; n <- te$n
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  w <- r3_roi_window(A, o$roi_test, rois); reg <- geos_set_precision(as_geos_geometry(w$poly), 1e-6)
  f2 <- o$fits$M2
  ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]], nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]])
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c")
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t))
  df <- data.frame(x = cc$x, y = cc$y, d_reg = geos_distance(geos_make_point(cc$x, cc$y), geos_boundary(reg)), ev = ev, a_t = cc$a_t, a_n = cc$a_n, d_gen = cc$d_gen)
  for (s in 1:10) df[[paste0("tgt", s)]] <- E$derive_cell_geometry(cen, tess, ct, random_seed = o$seed + 100000 * s, verbose = FALSE)$cell_df$nucleus_area_target
  write_parquet(df, file.path(OD, sprintf("%s_fold%d.parquet", A, k)))
  data.frame(A = A, k = k, y_nc_med = median(o$yard$ks_nc), y_nc_max = max(o$yard$ks_nc), y_an_med = median(o$yard$ks_an), y_an_max = max(o$yard$ks_an),
             y_sp_med = median(o$yard$d_sp), y_sp_max = max(o$yard$d_sp), y_b7_med = median(b$yard_b7), y_b7_max = max(b$yard_b7), b7_real = b$b7_real, real_sp = o$real_sp) }
jobs <- expand.grid(A = ARCH$archetype, k = 1:5, stringsAsFactors = FALSE)
r <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(one(jobs$A[i], jobs$k[i]), error = function(e) conditionMessage(e)), mc.cores = 30, mc.preschedule = FALSE)
bad <- !vapply(r, is.data.frame, TRUE); print(sum(bad)); print(unlist(r[bad]))
write.csv(do.call(rbind, r[!bad]), file.path(OD, "yards.csv"), row.names = FALSE)
