# tools/S1.5_bench.R - S1.5: banco del pannello (G-3, G-4, P-1...P-8, descrittive). Definizioni: docs/steps/S1.5.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.5_bench.R <stadio> [arg]
#   md5before | gen | faults | nulls | check | perf <preset> | cp1 | md5after   (stadi rilanciabili: saltano i file fatti)
# Uscite stabili in /mnt/micron/geo_spatialtrans/S1.5 (un file per unita'); tabelle riassuntive in results/S1.5.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5_faults.R")
data.table_ok <- requireNamespace("data.table", quietly = TRUE); if (data.table_ok) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); STAGE <- if (length(args)) args[1] else "none"
M <- "/mnt/micron/geo_spatialtrans"; OUT <- file.path(M, "S1.5"); RES <- "results/S1.5"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(RES, showWarnings = FALSE)
wk_file <- file.path(OUT, "WORKERS"); if (!file.exists(wk_file)) writeLines("20", wk_file)
NW <- function() min(28L, as.integer(readLines(wk_file)[1]))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
ARMS <- c("preserve_area", "none"); SEEDS <- 1:10; BASE <- 20261021
WIN <- read.csv("results/R2b/R2b_windows.csv"); WIN_ROI <- unique(paste(WIN$archetype, WIN$roi_id))
seed_of <- function(i, s) BASE + 1000 * i + s
clust_of <- function(w) { idx <- which(w$m, arr.ind = TRUE); data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = 1L) }
layer_md5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  on.exit(unlink(f)); unname(tools::md5sum(f)) }
run_pc <- function(i, s, arm, w = NULL) {
  A <- rois$archetype[i]; if (is.null(w)) w <- r3_roi_window(A, rois$roi_id[i], rois)
  suppressWarnings(place_cells(clust_of(w), PRESET[[A]], pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                               small_components = "exclude", random_seed = seed_of(i, s), size_conversion = arm, keep_tiles = TRUE,
                               strict = FALSE, verbose = FALSE))
}
ecc_of <- function(g) vapply(seq_along(g), function(k) { xy <- sf::st_coordinates(g[[k]])[, 1:2]
  m <- r3_poly_moments(xy[, 1], xy[, 2]); r3_shape(m[["mxx"]], m[["myy"]], m[["mxy"]])$ecc }, 0)
metrics_of <- function(L, w) {
  cd <- L$cell_df; td <- L$territory_df
  in_frame <- geos::geos_within(geos::as_geos_geometry(L$tiles), geos::as_geos_geometry(w$frame))
  it <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999)
  a <- cd$territory_area[it]; nc <- cd$nucleus_area / cd$territory_area; vt <- L$metadata$validation$table
  data.frame(n = nrow(cd), n_interior = sum(it), median_area_nuc = median(cd$nucleus_area[it]), median_nc = median(nc[it]),
             frac_nc095 = mean(nc[it] >= 0.95), frac_trunc = mean(cd$nucleus_area < (1 - 1e-6) * cd$nucleus_area_target),
             cv_terr = sd(a) / mean(a), median_ecc_T = median(ecc_of(L$cell_territories[which(it)])), median_eq_r = median(sqrt(a / pi)),
             frac_pale = mean(cd$nucleus_pale), n_failed = sum(L$metadata$seed_info$n_failed),
             desc_border_ratio = median(cd$territory_area[td$clipped]) / median(a),
             n_error = L$metadata$validation$n_error, n_warn = L$metadata$validation$n_warn,
             statuses = paste(vt$check, vt$status, sep = ":", collapse = " "), elapsed_s = L$metadata$elapsed_s, md5 = layer_md5(L))
}
par_run <- function(jobs, fun, label) {
  t0 <- Sys.time(); r <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch({ fun(jobs[k, ]); "ok" }, error = function(e) paste("ERR", conditionMessage(e))),
                                  mc.cores = NW(), mc.preschedule = FALSE)
  r <- unlist(r); cat(sprintf("[%s] %d unita', %d errori, %.1f min\n", label, length(r), sum(grepl("^ERR", r)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (any(grepl("^ERR", r))) print(head(cbind(jobs[grepl("^ERR", r), ], err = r[grepl("^ERR", r)]), 10))
}
md5_inputs <- function() {
  f <- c(list.files("results/R2/tissue_masks", full.names = TRUE), list.files(file.path(M, "R2/rois"), full.names = TRUE),
         list.files(file.path(M, "S1.4/real"), full.names = TRUE), list.files(file.path(M, "R3b/S1.4/real"), full.names = TRUE),
         "results/R2b/R2b_windows.csv", "results/R2c/R2c_dmax.csv", "results/R2b/annotations/R2b_annotations_Luca.json",
         "results/S1.4/S1.4_calibration.csv", "results/R3b/R3b_refs.csv", "docs/BIO_REFERENCES.md")
  data.frame(file = f, md5 = unname(tools::md5sum(f)))
}

if (STAGE == "md5before") { write.csv(md5_inputs(), file.path(OUT, "md5_before.csv"), row.names = FALSE); cat("md5 prima:", nrow(md5_inputs()), "file\n") }
if (STAGE == "md5after") {
  b <- read.csv(file.path(OUT, "md5_before.csv")); a <- md5_inputs(); write.csv(a, file.path(OUT, "md5_after.csv"), row.names = FALSE)
  same <- nrow(a) == nrow(b) && identical(a$file, b$file) && identical(a$md5, b$md5)
  cat(sprintf("%s md5 dei dati letti invariati: %d/%d\n", if (same) "PASS" else "FAIL", sum(a$md5 == b$md5[match(a$file, b$file)], na.rm = TRUE), nrow(b)))
}

# ---- gen: 30 ROI x 10 semi x 2 bracci ----
if (STAGE == "gen") {
  for (arm in ARMS) dir.create(file.path(OUT, "gen", arm), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUT, "layers"), showWarnings = FALSE); dir.create(file.path(OUT, "nuc"), showWarnings = FALSE)
  jobs <- expand.grid(i = rois$roi_idx, s = SEEDS, arm = ARMS, stringsAsFactors = FALSE)
  jobs <- jobs[order(jobs$s, -(rois$archetype[jobs$i] == "A4")), ]
  par_run(jobs, function(j) {
    i <- j$i; s <- j$s; arm <- j$arm; A <- rois$archetype[i]; roi <- rois$roi_id[i]
    f <- file.path(OUT, "gen", arm, sprintf("%s_%s_s%02d.rds", A, roi, s)); if (file.exists(f)) return(invisible())
    w <- r3_roi_window(A, roi, rois); L <- run_pc(i, s, arm, w)
    m <- cbind(data.frame(archetype = A, roi_id = roi, roi_idx = i, seed = s, arm = arm, random_seed = seed_of(i, s)), metrics_of(L, w))
    if (paste(A, roi) %in% WIN_ROI) {
      dir.create(file.path(OUT, "nuc", arm), showWarnings = FALSE)
      g <- geos::as_geos_geometry(L$cell_nuclei)
      write_parquet(data.frame(idx = seq_along(g), wkb = I(geos::geos_write_wkb(g))), file.path(OUT, "nuc", arm, sprintf("%s_%s_s%02d.parquet", A, roi, s)))
    }
    if (s == 1) saveRDS(L, file.path(OUT, "layers", sprintf("%s_%s_%s.rds", A, roi, arm)))
    saveRDS(m, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  }, "gen")
  tab <- do.call(rbind, lapply(list.files(file.path(OUT, "gen"), recursive = TRUE, pattern = "rds$", full.names = TRUE), readRDS))
  write.csv(tab, file.path(RES, "S1.5_panel_units.csv"), row.names = FALSE); cat("unita':", nrow(tab), "\n")
}

# ---- faults: guasti costruiti sulle uscite del seme 1 (braccio preserve_area) ----
if (STAGE == "faults") {
  FDIR <- Sys.getenv("S15_FAULTS_DIR", "faults"); dir.create(file.path(OUT, FDIR), showWarnings = FALSE)
  jobs <- expand.grid(i = rois$roi_idx, type = FAULT_TYPES, delta = c(0.01, 0.1, 1), stringsAsFactors = FALSE)
  par_run(jobs, function(j) {
    A <- rois$archetype[j$i]; roi <- rois$roi_id[j$i]
    f <- file.path(OUT, FDIR, sprintf("%s_%s_%s_%g.rds", A, roi, j$type, j$delta)); if (file.exists(f)) return(invisible())
    L <- readRDS(file.path(OUT, "layers", sprintf("%s_%s_preserve_area.rds", A, roi)))
    fz <- make_fault(L, j$type, j$delta); v <- validate_cell_layer(fz$cl, c8_interval = tissue_preset_spec(PRESET[[A]])$c8_interval)
    r <- data.frame(archetype = A, roi_id = roi, type = j$type, delta = j$delta, expect = fz$expect, cell = fz$cell,
                    status_expected = v$table$status[v$table$check == fz$expect], n_error = v$n_error,
                    statuses = paste(v$table$check, v$table$status, sep = ":", collapse = " "))
    saveRDS(r, f)
  }, "faults")
  tab <- do.call(rbind, lapply(list.files(file.path(OUT, FDIR), full.names = TRUE), readRDS))
  write.csv(tab, file.path(RES, if (FDIR == "faults") "S1.5_G4_faults.csv" else paste0("S1.5_G4_", FDIR, ".csv")), row.names = FALSE)
  cat(sprintf("%s G-4 guasti rilevati dal controllo atteso: %d/%d\n", if (all(tab$status_expected == "ERROR")) "PASS" else "FAIL", sum(tab$status_expected == "ERROR"), nrow(tab)))
}

# ---- nulls: i 1 800 nulli di S1.3 (codice e semi di tools/S1.3_roi.R, gen_points) con i nuclei del preset ----
if (STAGE == "nulls") {
  NDIR <- Sys.getenv("S15_NULLS_DIR", "nulls"); dir.create(file.path(OUT, NDIR), showWarnings = FALSE)
  DSTAR <- c(A1 = 2.5, A2 = 3.5, A3 = 3.5, A4 = 4, A5 = 7, A6 = 2.5)
  PRIM <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
  gen_points <- function(i, model, rp, w, n_obs, dstar) {      # copia di tools/S1.3_roi.R (BASE_SEED 20261008)
    if (model == "CSR") { seed <- 20261008 + 1000 * i + 100 + rp
      set.seed(seed); idx <- which(w$m, arr.ind = TRUE); s <- sample.int(nrow(idx), n_obs, replace = TRUE)
      return(list(x = (idx[s, 2] - 1 + runif(n_obs)) * w$px, y = (idx[s, 1] - 1 + runif(n_obs)) * w$px, seed = seed)) }
    comp <- data.frame(cluster_id = unique(as.character(w$reg$region_df$cluster_id)), cell_type = "c", fraction = 1)
    if (model == "RSA") { seed <- 20261008 + 1000 * i + 200 + rp; ct <- data.frame(cell_type = "c", density = n_obs / w$area_um2 * 1e6, min_dist_um = dstar)
    } else { seed <- 20261008 + 1000 * i + 300 + rp; ct <- data.frame(cell_type = "c", density = n_obs / w$area_um2 * 1e6) }
    o <- suppressWarnings(seed_centroids(w$reg, ct, comp, random_seed = seed, verbose = FALSE)); list(x = o$centroids$x, y = o$centroids$y, seed = seed)
  }
  jobs <- expand.grid(i = rois$roi_idx, model = c("CSR", "RSA", "RSArule"), rp = 1:20, stringsAsFactors = FALSE)
  jobs <- jobs[order(-(rois$archetype[jobs$i] == "A4")), ]
  par_run(jobs, function(j) {
    A <- rois$archetype[j$i]; roi <- rois$roi_id[j$i]
    f <- file.path(OUT, NDIR, sprintf("%s_%s_%s_%02d.rds", A, roi, j$model, j$rp)); if (file.exists(f)) return(invisible())
    w <- r3_roi_window(A, roi, rois)
    n_obs <- readRDS(file.path(M, "R3", "real", sprintf("%s_%s_%s_roi.rds", A, roi, PRIM[[A]])))$n
    p <- gen_points(j$i, j$model, j$rp, w, n_obs, DSTAR[[A]])
    reg1 <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2), region_polygons = sf::st_sfc(w$poly[[1]]))
    cen <- data.frame(cell_id = seq_along(p$x), region_id = 1L, x = p$x, y = p$y, cell_type = paste0(PRESET[[A]], "_all"), cluster_id = 1L)
    tv <- tessellate_voronoi(cen, reg1, verbose = FALSE)
    sp <- tissue_preset_spec(PRESET[[A]]); dg <- derive_cell_geometry(cen, tv, sp$cell_types, random_seed = stage_seeds(p$seed)[["size"]], strict = FALSE, verbose = FALSE)   # FND-14
    rd <- data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2, n_cells = nrow(cen))
    cl <- list(cell_df = dg$cell_df, region_df = rd, region_polygons = reg1$region_polygons, cell_territories = tv$cell_territories,
               cell_nuclei = dg$cell_nuclei, territory_df = tv$territory_df, region_check = tv$region_check)
    v <- validate_cell_layer(cl, c8_interval = sp$c8_interval)
    saveRDS(data.frame(archetype = A, roi_id = roi, model = j$model, rep = j$rp, seed = p$seed, n = nrow(cen), n_error = v$n_error,
                       statuses = paste(v$table$check, v$table$status, sep = ":", collapse = " "),
                       c2 = v$table$value[v$table$check == "C2"]), f)
  }, "nulls")
  tab <- do.call(rbind, lapply(list.files(file.path(OUT, NDIR), full.names = TRUE), readRDS))
  write.csv(tab, file.path(RES, if (NDIR == "nulls") "S1.5_G4_nulls.csv" else paste0("S1.5_G4_", NDIR, ".csv")), row.names = FALSE)
  cat(sprintf("%s G-4 nulli di S1.3 senza ERROR: %d/%d\n", if (all(tab$n_error == 0)) "PASS" else "FAIL", sum(tab$n_error == 0), nrow(tab)))
}

# ---- recheck: dopo la revisione, tutte le 600 unita' con il validatore corretto; md5 delle geometrie uguale al pannello? ----
if (STAGE == "recheck") {
  dir.create(file.path(OUT, "recheck"), showWarnings = FALSE)
  tab <- read.csv(file.path(RES, "S1.5_panel_units.csv"))
  jobs <- expand.grid(i = rois$roi_idx, s = SEEDS, arm = ARMS, stringsAsFactors = FALSE)
  par_run(jobs, function(j) {
    A <- rois$archetype[j$i]; roi <- rois$roi_id[j$i]; f <- file.path(OUT, "recheck", sprintf("%s_%s_s%02d_%s.rds", A, roi, j$s, j$arm)); if (file.exists(f)) return(invisible())
    L <- run_pc(j$i, j$s, j$arm); vt <- L$metadata$validation$table
    saveRDS(data.frame(archetype = A, roi_id = roi, seed = j$s, arm = j$arm, md5 = layer_md5(L), md5_panel = tab$md5[tab$roi_idx == j$i & tab$seed == j$s & tab$arm == j$arm],
                       n_error = L$metadata$validation$n_error, statuses = paste(vt$check, vt$status, sep = ":", collapse = " ")), f)
  }, "recheck")
  r <- do.call(rbind, lapply(list.files(file.path(OUT, "recheck"), full.names = TRUE), readRDS)); write.csv(r, file.path(RES, "S1.5_G3_recheck.csv"), row.names = FALSE)
  cat(sprintf("%s ricontrollo: %d/%d unita' senza ERROR, md5 uguale al pannello %d/%d\n", if (all(r$n_error == 0) && all(r$md5 == r$md5_panel)) "PASS" else "FAIL",
              sum(r$n_error == 0), nrow(r), sum(r$md5 == r$md5_panel), nrow(r)))
}

# ---- check: G-3 (ERROR, riproducibilita', equivalenza con la catena manuale) ----
if (STAGE == "check") {
  tab <- read.csv(file.path(RES, "S1.5_panel_units.csv"))
  cat(sprintf("%s G-3 unita' del pannello senza ERROR: %d/%d\n", if (all(tab$n_error == 0)) "PASS" else "FAIL", sum(tab$n_error == 0), nrow(tab)))
  pick <- sapply(names(PRESET), function(A) rois$roi_idx[rois$archetype == A][1])
  rr <- do.call(rbind, mclapply(pick, function(i) {
    A <- rois$archetype[i]; roi <- rois$roi_id[i]; w <- r3_roi_window(A, roi, rois)
    L <- run_pc(i, 1, "preserve_area", w)
    m0 <- tab$md5[tab$roi_idx == i & tab$seed == 1 & tab$arm == "preserve_area"]
    sp <- tissue_preset_spec(PRESET[[A]]); s <- stage_seeds(seed_of(i, 1)); cl <- clust_of(w)
    rg <- extract_regions(cl, w$px, 0, cluster_col = "intensity_cluster", simplify_tol_um = 0, stride = 1L, verbose = FALSE)
    comp <- unique(data.frame(cluster_id = rg$region_df$cluster_id, cell_type = sp$cell_types$cell_type, fraction = 1))
    sc <- suppressWarnings(seed_centroids(rg, sp$cell_types, comp, random_seed = s[["seed"]], verbose = FALSE))
    tv <- tessellate_voronoi(sc$centroids, rg, verbose = FALSE)
    dg <- derive_cell_geometry(sc$centroids, tv, sp$cell_types, random_seed = s[["size"]], verbose = FALSE)
    eq <- identical(dg$cell_df, L$cell_df[, names(dg$cell_df)]) && identical(sf::st_as_binary(dg$cell_nuclei), sf::st_as_binary(L$cell_nuclei)) &&
      identical(sf::st_as_binary(tv$cell_territories), sf::st_as_binary(L$cell_territories))
    data.frame(archetype = A, roi_id = roi, md5_rerun_identical = identical(layer_md5(L), m0), manual_chain_identical = eq)
  }, mc.cores = 6))
  write.csv(rr, file.path(RES, "S1.5_G3_repro.csv"), row.names = FALSE); print(rr)
  cat(sprintf("%s G-3 rerun identico %d/6, catena manuale identica %d/6\n", if (all(rr$md5_rerun_identical & rr$manual_chain_identical)) "PASS" else "FAIL",
              sum(rr$md5_rerun_identical), sum(rr$manual_chain_identical)))
}

# ---- perf: P-8, 4000 x 4000 px a 1 µm/px, un cluster, riscaldamento + 3 ripetizioni nello stesso processo ----
if (STAGE == "perf") {
  pr <- args[2]; dir.create(file.path(OUT, "perf"), showWarnings = FALSE)
  env <- function() c(load = readLines("/proc/loadavg"), tctl = paste(system("sensors 2>/dev/null | grep -i tctl", intern = TRUE), collapse = " "))
  big <- expand.grid(x = 1:4000, y = 1:4000); big$intensity_cluster <- 1L
  rows <- list()
  for (k in 0:3) {
    e0 <- env(); t <- system.time(L <- place_cells(big, pr, pixel_size_um = 1, random_seed = BASE + k, verbose = FALSE))[["elapsed"]]
    rows[[k + 1]] <- data.frame(preset = pr, rep = k, warmup = k == 0, elapsed_s = t, n_cells = nrow(L$cell_df), n_error = L$metadata$validation$n_error,
                                t_seed = L$metadata$seed_info$elapsed_s, t_tess = L$metadata$tessellation_info$elapsed_s,
                                t_geom = L$metadata$geometry_info$elapsed_s, load_before = e0[["load"]], tctl_before = e0[["tctl"]])
    rm(L); invisible(gc())
  }
  r <- do.call(rbind, rows); write.csv(r, file.path(OUT, "perf", paste0(pr, ".csv")), row.names = FALSE); print(r)
}

# ---- cp1: stima della scala sulle immagini H&E ricampionate (tools/S1.5_cp1_resample.py) ----
if (STAGE == "cp1") {
  TAG <- Sys.getenv("S15_CP1_TAG", "est"); SC <- Sys.getenv("S15_CP1_SIZECHECK", "FALSE") == "TRUE"
  dir.create(file.path(OUT, "cp1", TAG), recursive = TRUE, showWarnings = FALSE)
  fs <- list.files(file.path(OUT, "cp1", "img"), pattern = "\\.npz$", full.names = TRUE)
  jobs <- expand.grid(f = fs, variant = c("pipeline", "full_mask"), stringsAsFactors = FALSE)
  par_run(jobs, function(j) {
    b <- sub("\\.npz$", "", basename(j$f)); o <- file.path(OUT, "cp1", TAG, paste0(b, "_", j$variant, ".rds")); if (file.exists(o)) return(invisible())
    pz <- file.path(OUT, "cp1", "img", paste0(b, ".parquet")); d <- as.data.frame(read_parquet(pz))
    tk <- strsplit(b, "_")[[1]]; A <- tk[1]; roi <- tk[2]; s_true <- as.numeric(sub("um$", "", tk[3]))
    cl <- if (j$variant == "pipeline") d[d$value < 0.7, ] else d[d$valid, ]
    cl$intensity_cluster <- 1L
    e <- tryCatch(estimate_pixel_size(cl, tissue_preset_spec(PRESET[[A]])$cell_types, size_check = SC), error = function(e) list(refused = NA, reason = conditionMessage(e), pixel_size_um = NA, ci = c(NA, NA), n_features = NA, sigma_px = NA))
    saveRDS(data.frame(archetype = A, roi_id = roi, scale_true = s_true, variant = j$variant, n_px = nrow(cl), refused = e$refused, reason = e$reason,
                       scale_est = e$pixel_size_um, ratio = e$pixel_size_um / s_true, ci_lo = e$ci[1], ci_hi = e$ci[2], n_features = e$n_features, sigma_px = e$sigma_px), o)
  }, "cp1")
  tab <- do.call(rbind, lapply(list.files(file.path(OUT, "cp1", TAG), full.names = TRUE), readRDS))
  write.csv(tab, file.path(RES, if (TAG == "est") "S1.5_CP1_scale.csv" else paste0("S1.5_CP1_scale_", TAG, ".csv")), row.names = FALSE); print(table(tab$variant, tab$scale_true, tab$refused, useNA = "ifany"))
}
