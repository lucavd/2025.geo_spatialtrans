# tools/S2.6_explore.R - S2.6: esplorazioni dichiarate (prima della pre-registrazione). Definizioni: docs/steps/S2.6.md.
# Uso (root del repo), un processo per chiamata (picco di memoria del processo):
#   Rscript --vanilla tools/S2.6_explore.R pc  <A> <side_um> <px_um> <seed> <outdir>   livello cellulare su un quadrato pieno
#   Rscript --vanilla tools/S2.6_explore.R sim <layer.rds> <seed> <outdir> [workers]  simulate_visium_hd() sul livello, frame_open = TRUE
# Scrive <outdir>/<tag>.json (tempo, picco VmHWM, profilo Rprof per funzione) e, per pc, <outdir>/<tag>.rds (livello).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
peak_mb <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
a <- commandArgs(TRUE); MODE <- a[1]
prof_top <- function(pf, keys) {
  s <- tryCatch(summaryRprof(pf)$by.total, error = function(e) NULL)
  if (is.null(s)) return(list())
  rn <- gsub('"', "", rownames(s)); out <- list()
  for (k in keys) { w <- which(rn == k); out[[k]] <- if (length(w)) s$total.time[w[1]] else 0 }
  out
}
if (MODE == "pc") {
  A <- a[2]; side <- as.numeric(a[3]); px <- as.numeric(a[4]); seed <- as.integer(a[5]); OUT <- a[6]
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  tag <- sprintf("pc_%s_%gum_px%g_s%d", A, side, px, seed)
  n <- as.integer(round(side / px))
  clust <- data.frame(x = rep(seq_len(n), times = n), y = rep(seq_len(n), each = n), intensity_cluster = 1L)
  pf <- tempfile(fileext = ".rprof"); gc(); m0 <- peak_mb()
  Rprof(pf, interval = 0.2)
  t0 <- proc.time()[["elapsed"]]
  L <- suppressWarnings(place_cells(clust, PRESET[[A]], pixel_size_um = px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                    small_components = "exclude", random_seed = seed, types = "multi", size_arm = "recentred",
                    spatial = "uniform", size_resid = "gaussian", strict = FALSE, verbose = FALSE))
  te <- proc.time()[["elapsed"]] - t0
  Rprof(NULL)
  keys <- c("extract_regions", "seed_centroids", "tessellate_voronoi", "derive_cell_geometry", "validate_cell_layer", "merge_small_components")
  bb <- sf::st_bbox(L$region_polygons)
  chk <- list(tag = tag, archetype = A, side_um = side, px_um = px, n_px = n * n, seed = seed, n_cells = nrow(L$cell_df),
              area_mm2 = as.numeric(sum(sf::st_area(L$region_polygons))) / 1e6, bbox = as.list(bb),
              cells_per_mm2 = nrow(L$cell_df) / (as.numeric(sum(sf::st_area(L$region_polygons))) / 1e6),
              mean_territory_um2 = mean(L$cell_df$territory_area), mean_nucleus_um2 = mean(L$cell_df$nucleus_area),
              elapsed_s = te, peak_mb = peak_mb(), peak_mb_before = m0, rprof_total_s = prof_top(pf, keys),
              load = readLines("/proc/loadavg"))
  saveRDS(L, file.path(OUT, paste0(tag, ".rds")))
  writeLines(jsonlite::toJSON(chk, auto_unbox = TRUE, digits = NA), file.path(OUT, paste0(tag, ".json")))
  cat(sprintf("[%s] %d cellule, %.1f s, %.0f MB\n", tag, chk$n_cells, te, chk$peak_mb))
} else if (MODE == "sim") {
  f <- a[2]; seed <- as.integer(a[3]); OUT <- a[4]; W <- if (length(a) >= 5) as.integer(a[5]) else 1L
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  L <- readRDS(f); A <- visium_hd_archetype(L)
  tag <- sprintf("sim_%s_s%d_w%d", sub("\\.rds$", "", basename(f)), seed, W)
  gc(); m0 <- peak_mb(); t0 <- proc.time()[["elapsed"]]
  out <- simulate_visium_hd(L, seed = seed, frame_open = TRUE, workers = W)
  te <- proc.time()[["elapsed"]] - t0
  g <- out$metadata$grid
  chk <- list(tag = tag, layer = basename(f), archetype = A, seed = seed, workers = W, n_cells = nrow(L$cell_df),
              n_bins = g$n_row * g$n_col, n_generated = out$metadata$n_generated, n_captured = out$metadata$n_captured,
              chunk_size = out$metadata$chunk_size, M = out$metadata$overgeneration,
              elapsed = as.list(out$metadata$elapsed_s), elapsed_s = te, peak_mb = peak_mb(), peak_mb_before = m0,
              size_out_mb = as.numeric(object.size(out)) / 2^20, load = readLines("/proc/loadavg"))
  writeLines(jsonlite::toJSON(chk, auto_unbox = TRUE, digits = NA), file.path(OUT, paste0(tag, ".json")))
  cat(sprintf("[%s] %.3g generati, %.1f s, %.0f MB\n", tag, chk$n_generated, te, chk$peak_mb))
} else stop("modo sconosciuto")
