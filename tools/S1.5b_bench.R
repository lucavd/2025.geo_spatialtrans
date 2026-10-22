# tools/S1.5b_bench.R - S1.5b: banco del pannello con i preset per tipo (G-3, P-1...P-8, CP-1, descrittive). Definizioni: docs/steps/S1.5b.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5b_bench.R <stadio>
#   md5before | gen | check | md5after   (stadi rilanciabili: saltano le unita' fatte)
# Bracci: "group", "recentred" (preset per tipo) e "single" (preset di S1.5: riferimento nullo per CP-1 e compatibilita' per G-3).
# Uscite stabili in /mnt/micron/geo_spatialtrans/S1.5b (un file per unita'); tabelle in results/S1.5b.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R")
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); STAGE <- if (length(args)) args[1] else "none"
M <- "/mnt/micron/geo_spatialtrans"; OUT <- file.path(M, "S1.5b"); RES <- "results/S1.5b"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
wk_file <- file.path(OUT, "WORKERS"); if (!file.exists(wk_file)) writeLines("6", wk_file)
NW <- function() min(28L, as.integer(readLines(wk_file)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
ARMS <- c("group", "recentred", "single"); SEEDS <- 1:10; BASE <- 20261021
seed_of <- function(i, s) BASE + 1000 * i + s                      # come S1.5
THR <- read.csv(file.path(RES, "S1.5b_real_refs_roi.csv"))$small_threshold_um2; THR <- THR[is.finite(THR)][1]   # 59.1 um2
clust_of <- function(w) { idx <- which(w$m, arr.ind = TRUE); data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L) }
layer_md5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  on.exit(unlink(f)); unname(tools::md5sum(f)) }
run_pc <- function(i, s, arm, w) {
  A <- rois$archetype[i]
  if (arm == "single") return(suppressWarnings(place_cells(clust_of(w), PRESET[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0,
                               stride = 1L, small_components = "exclude", random_seed = seed_of(i, s), size_conversion = "preserve_area",
                               keep_tiles = TRUE, strict = FALSE, verbose = FALSE)))
  suppressWarnings(place_cells(clust_of(w), PRESET[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                               small_components = "exclude", random_seed = seed_of(i, s), types = "multi", size_arm = arm,
                               keep_tiles = TRUE, strict = FALSE, verbose = FALSE))
}
# etichette del nullo (CP-1): tipi estratti a caso secondo la composizione del preset, seme dell'unita' + 7
null_types <- function(A, n, seed) { sp <- tissue_preset_spec(PRESET[[A]], types = "multi", size_arm = "group")
  set.seed(seed + 7L); sample(sp$composition$cell_type, n, replace = TRUE, prob = sp$composition$fraction) }
metrics_of <- function(L, w, A, arm, seed) {
  cd <- L$cell_df; td <- L$territory_df
  in_frame <- geos::geos_within(geos::as_geos_geometry(L$tiles), geos::as_geos_geometry(w$frame))
  it <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
  a <- cd$territory_area[it]; nc <- cd$nucleus_area / cd$territory_area; vt <- L$metadata$validation$table
  ty <- if (arm == "single") null_types(A, nrow(cd), seed) else as.character(cd$cell_type)
  grp <- if (arm == "single") tissue_preset_spec(PRESET[[A]], types = "multi", size_arm = "group")$cell_types else L$cell_types_resolved
  gi <- grp$size_group[match(ty, grp$cell_type)]
  u <- data.frame(n = nrow(cd), n_interior = sum(it), median_area_nuc = median(cd$nucleus_area[it]), median_nc = median(nc[it]),
                  frac_nc095 = mean(nc[it] >= 0.95), frac_trunc = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target),
                  cv_terr = sd(a) / mean(a), n_failed = sum(L$metadata$seed_info$n_failed), n_error = L$metadata$validation$n_error,
                  n_warn = L$metadata$validation$n_warn, statuses = paste(vt$check, vt$status, sep = ":", collapse = " "),
                  elapsed_s = L$metadata$elapsed_s, md5 = layer_md5(L), frac_band = NA_real_, ratio_neuron_glia = NA_real_)
  if (A == "A5") {
    an <- cd$nucleus_area[it]; u$frac_band <- mean(an > THR / sqrt(2) & an < THR * sqrt(2))
    gN <- grp$size_group[grp$cell_type == "neuroni"]; gg <- gi[it]
    u$ratio_neuron_glia <- median(an[gg == gN]) / median(an[gg != gN])
  }
  tt <- do.call(rbind, lapply(sort(unique(ty)), function(t) { k <- it & ty == t
    data.frame(cell_type = t, n_type = sum(ty == t), n_interior = sum(k), median_area_nuc = if (any(k)) median(cd$nucleus_area[k]) else NA_real_,
               median_nc = if (any(k)) median(nc[k]) else NA_real_) }))
  list(unit = u, type = tt, it = it, ty = ty)
}
par_run <- function(jobs, fun, label) {
  t0 <- Sys.time(); r <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch({ fun(jobs[k, ]); "ok" }, error = function(e) paste("ERR", conditionMessage(e))),
                                  mc.cores = NW(), mc.preschedule = FALSE)
  r <- unlist(r); cat(sprintf("[%s] %d unita', %d errori, %.1f min\n", label, length(r), sum(grepl("^ERR", r)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (any(grepl("^ERR", r))) print(head(cbind(jobs[grepl("^ERR", r), ], err = r[grepl("^ERR", r)]), 10))
}
md5_inputs <- function() {
  f <- c(list.files("results/R2/tissue_masks", full.names = TRUE), list.files(file.path(M, "R2/rois"), full.names = TRUE),
         list.files(file.path(M, "R6/labels"), pattern = "exp_umi20_seed20261022_nuclei_ph1.parquet$", full.names = TRUE),
         list.files(file.path(M, "R6/terr"), full.names = TRUE), "results/R2b/R2b_points.csv", "results/R2b/R2b_windows.csv",
         "results/S1.5/S1.5_panel_units.csv", "results/S1.5/S1.5_real_refs_roi.csv", "docs/BIO_REFERENCES.md")
  data.frame(file = f, md5 = unname(tools::md5sum(f)))
}
if (STAGE == "md5before") { write.csv(md5_inputs(), file.path(OUT, "md5_before.csv"), row.names = FALSE); cat("md5 prima:", nrow(md5_inputs()), "file\n") }
if (STAGE == "md5after") {
  b <- read.csv(file.path(OUT, "md5_before.csv")); a <- md5_inputs(); write.csv(a, file.path(OUT, "md5_after.csv"), row.names = FALSE)
  same <- nrow(a) == nrow(b) && identical(a$file, b$file) && identical(a$md5, b$md5)
  cat(sprintf("%s md5 dei dati letti invariati: %d/%d\n", if (same) "PASS" else "FAIL", sum(a$md5 == b$md5[match(a$file, b$file)], na.rm = TRUE), nrow(b)))
}
if (STAGE == "gen") {
  for (arm in ARMS) dir.create(file.path(OUT, "gen", arm), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUT, "layers"), showWarnings = FALSE)
  jobs <- expand.grid(i = rois$roi_idx, s = SEEDS, arm = ARMS, stringsAsFactors = FALSE)
  jobs <- jobs[order(jobs$s, -(rois$archetype[jobs$i] == "A4")), ]
  par_run(jobs, function(j) {
    i <- j$i; s <- j$s; arm <- j$arm; A <- rois$archetype[i]; roi <- rois$roi_id[i]
    f <- file.path(OUT, "gen", arm, sprintf("%s_%s_s%02d.rds", A, roi, s)); if (file.exists(f)) return(invisible())
    w <- r3_roi_window(A, roi, rois); L <- run_pc(i, s, arm, w); m <- metrics_of(L, w, A, arm, seed_of(i, s))
    key <- data.frame(archetype = A, roi_id = roi, roi_idx = i, seed = s, arm = arm, random_seed = seed_of(i, s))
    if (s == 1 && arm != "single") saveRDS(L, file.path(OUT, "layers", sprintf("%s_%s_%s.rds", A, roi, arm)))
    saveRDS(list(unit = cbind(key, m$unit), type = cbind(key, m$type)), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  }, "gen")
  z <- lapply(list.files(file.path(OUT, "gen"), recursive = TRUE, pattern = "rds$", full.names = TRUE), readRDS)
  write.csv(do.call(rbind, lapply(z, `[[`, "unit")), file.path(RES, "S1.5b_panel_units.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(z, `[[`, "type")), file.path(RES, "S1.5b_panel_types.csv"), row.names = FALSE)
  cat("unita':", length(z), "\n")
}
if (STAGE == "check") {
  U <- read.csv(file.path(RES, "S1.5b_panel_units.csv"), stringsAsFactors = FALSE)
  multi <- U[U$arm != "single", ]
  e1 <- sum(grepl("C2:ERROR|C3:ERROR|C5:ERROR|C9:ERROR", U$statuses)); nf <- sum(U$n_failed)
  # compatibilita': il braccio "single" deve riprodurre le unita' di S1.5 (braccio preserve_area, stessi semi)
  S <- read.csv("results/S1.5/S1.5_panel_units.csv", stringsAsFactors = FALSE); S <- S[S$arm == "preserve_area", ]
  sg <- U[U$arm == "single", ]; k <- match(paste(sg$archetype, sg$roi_id, sg$seed), paste(S$archetype, S$roi_id, S$seed))
  compat <- sum(sg$md5 == S$md5[k], na.rm = TRUE)
  # riproducibilita': un'unita' per tessuto e braccio per tipo, rigenerata
  pick <- do.call(rbind, lapply(split(rois, rois$archetype), function(d) d[1, ]))
  rep_ok <- vapply(seq_len(nrow(pick)), function(r) { i <- pick$roi_idx[r]; arm <- if (r %% 2) "group" else "recentred"
    w <- r3_roi_window(pick$archetype[r], pick$roi_id[r], rois); L <- run_pc(i, 3, arm, w)
    identical(layer_md5(L), U$md5[U$roi_idx == i & U$seed == 3 & U$arm == arm]) }, TRUE)
  out <- data.frame(check = c("units", "no_error_C2_C3_C5_C9", "n_failed_total", "single_equals_S1.5_md5", "rerun_identical"),
                    value = c(nrow(U), nrow(U) - e1, nf, compat, sum(rep_ok)), expected = c(900, 900, 0, 300, 6))
  out$pass <- out$value == out$expected; write.csv(out, file.path(RES, "S1.5b_G3_check.csv"), row.names = FALSE); print(out)
  cat(sprintf("%s G-3 (md5 dei dati: stadio md5after)\n", if (all(out$pass)) "PASS" else "FAIL"))
}
