# tools/R3_run.R — R3: Voronoi sui nuclei reali (stadi real, null, calib). Pre-registrazione:
# results/R3/R3_preregistration.md (commit 3fac4e7). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/R3_run.R <stage> [outdir]
#   stage = real | null | calib ; outdir default /mnt/micron/geo_spatialtrans/R3 (C-R3.6: rerun in un'altra cartella)
# Stadi rilanciabili: saltano i file gia' scritti. Dopo `real`: .venv/bin/python tools/R3_containment.py all
# Definizioni operative aggiunte qui (piu' specifiche del testo, non divergenti):
#  - intensita' locale con leaveoneout = FALSE (il generatore contribuisce alla propria stima);
#  - CSR: pixel validi estratti con reinserimento, posizione uniforme nel pixel (convenzione dei lati);
#  - RSA: seed_centroids() con densita' = n_obs / area valida; max_attempts_factor = 1000 (default);
#  - nulli: riassunto per replica + tabella per cellula della sola replica 1;
#  - D-2: nndist con meno-campionamento (punti con nndist < distanza dal bordo della finestra), clarkevans "cdf";
#  - calibrazione: regione = quadrato della finestra ∩ maschera valida − esclusioni di Luca (A6_w2–w4).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel); library(arrow); library(spatstat.geom)
                                 library(spatstat.explore) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("tools/R3_voronoi_metrics.R")
args <- commandArgs(TRUE); STAGE <- args[1]
OUT <- if (length(args) > 1) args[2] else "/mnt/micron/geo_spatialtrans/R3"
for (d in c("real", "null", "calib")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
NCORES <- 40L; NREP <- 20L; BASE_SEED <- 20261008
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   sens = c(NA, NA, NA, "stardist_he", NA, "stardist_he"),
                   dstar = c(2.5, 3.5, 3.5, 4, 7, 2.5), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))

load_nuclei <- function() {
  nuc <- read_parquet("results/R2/R2_nuclei_all.parquet",
                      col_select = c("archetype", "roi_id", "method", "scale", "keep", "label", "x_um", "y_um",
                                     "area_um2", "eccentricity", "orientation", "um_per_px"))
  as.data.frame(nuc[nuc$scale == 1 & nuc$keep, ])
}

# ---- stadio real ---------------------------------------------------------------------------------
real_one <- function(A, roi, method, nuc) {
  f_cells <- file.path(OUT, "real", sprintf("%s_%s_%s_cells.parquet", A, roi, method))
  if (file.exists(f_cells)) return(NULL)
  w <- r3_roi_window(A, roi, rois)
  d <- nuc[nuc$archetype == A & nuc$roi_id == roi & nuc$method == method, ]
  n_keep <- nrow(d)
  d <- d[r3_in_mask(d$x_um, d$y_um, w$m, w$px), ]
  gen <- data.frame(idx = seq_len(nrow(d)), label = d$label, x = d$x_um, y = d$y_um, area_nuc = d$area_um2,
                    ecc_parquet = d$eccentricity, orient_parquet = d$orientation, um_per_px = d$um_per_px)
  write_parquet(gen, file.path(OUT, "real", sprintf("%s_%s_%s_gen.parquet", A, roi, method)))
  tv <- r3_tessellate(gen$x, gen$y, w$poly, w$frame)
  li <- r3_local_intensity(gen$x, gen$y, w)
  tv$lambda_loc <- li$lambda_loc; tv$label <- gen$label; tv$area_nuc <- gen$area_nuc
  X <- li$X; nn <- nndist(X); bd <- bdist.points(X)
  ce <- tryCatch(clarkevans(X, correction = "cdf"), error = function(e) NA_real_)
  dip <- tryCatch(diptest::dip.test(tv$area[tv$interior])$p.value, error = function(e) NA_real_)
  write_parquet(tv, f_cells)
  roi_row <- data.frame(archetype = A, roi_id = roi, method = method, n_keep = n_keep, n = nrow(gen),
                        area_poly_um2 = w$area_um2, area_px_um2 = w$n_valid_px * w$px^2, sum_area = sum(tv$area),
                        gen_in_own = mean(tv$gen_in_own), sigma_loc_um = li$sigma, lambda_mm2 = li$lambda * 1e6,
                        nn_median_minus = median(nn[nn < bd]), nn_n_minus = sum(nn < bd), clark_evans_cdf = as.numeric(ce),
                        dip_p_area = dip, t(r3_summary(tv, tv$lambda_loc)))
  saveRDS(roi_row, file.path(OUT, "real", sprintf("%s_%s_%s_roi.rds", A, roi, method)))
  roi_row
}

if (STAGE == "real") {
  nuc <- load_nuclei()
  jobs <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) {
    a <- ARCH[ARCH$archetype == rois$archetype[i], ]
    m <- c(a$primary, if (!is.na(a$sens)) a$sens)
    data.frame(archetype = rois$archetype[i], roi_id = rois$roi_id[i], method = m, stringsAsFactors = FALSE)
  }))
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(real_one(j$archetype, j$roi_id, j$method, nuc), error = function(e) conditionMessage(e))
    cat(sprintf("[real] %s %s %s %.0fs %s\n", j$archetype, j$roi_id, j$method, proc.time()[["elapsed"]] - t0,
                if (is.character(r)) r else "ok")); r
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  errs <- vapply(res, is.character, TRUE); if (any(errs)) stop("errori nello stadio real: ", sum(errs))
}

# ---- stadio null ---------------------------------------------------------------------------------
null_one <- function(i, model, rp) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]
  f <- file.path(OUT, "null", sprintf("%s_%s_%s_%02d.rds", A, roi, model, rp))
  if (file.exists(f)) return(invisible(NULL))
  w <- r3_roi_window(A, roi, rois)
  n_obs <- readRDS(file.path(OUT, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  seed <- BASE_SEED + 1000 * i + 100 * (if (model == "CSR") 1 else 2) + rp
  if (model == "CSR") {
    set.seed(seed); idx <- which(w$m, arr.ind = TRUE); s <- sample.int(nrow(idx), n_obs, replace = TRUE)
    x <- (idx[s, 2] - 1 + runif(n_obs)) * w$px; y <- (idx[s, 1] - 1 + runif(n_obs)) * w$px
    n_failed <- 0L
  } else {
    comp <- data.frame(cluster_id = unique(as.character(w$reg$region_df$cluster_id)), cell_type = "c", fraction = 1)
    o <- seed_centroids(w$reg, data.frame(cell_type = "c", density = n_obs / w$area_um2 * 1e6, min_dist_um = a$dstar),
                        comp, random_seed = seed, verbose = FALSE)
    x <- o$centroids$x; y <- o$centroids$y; n_failed <- o$info$n_failed
  }
  tv <- r3_tessellate(x, y, w$poly, w$frame)
  li <- r3_local_intensity(x, y, w)
  s <- r3_summary(tv, li$lambda_loc)
  ns <- tabulate(pmin(tv$nsides[tv$interior], 15), 15)
  out <- list(archetype = A, roi_id = roi, model = model, rep = rp, seed = seed, n_failed = n_failed,
              summary = s, nsides = ns, sum_area = sum(tv$area), area_poly = w$area_um2,
              cells = if (rp == 1) tv[, c("x", "y", "area", "interior", "nsides", "ecc_T", "theta_T")] else NULL)
  saveRDS(out, f); invisible(NULL)
}
if (STAGE == "null") {
  jobs <- expand.grid(i = seq_len(nrow(rois)), model = c("CSR", "RSA"), rp = seq_len(NREP), stringsAsFactors = FALSE)
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
    r <- tryCatch(null_one(j$i, j$model, j$rp), error = function(e) conditionMessage(e))
    if (is.character(r)) cat(sprintf("[null] ERRORE %d %s %d: %s\n", j$i, j$model, j$rp, r)); r
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  errs <- vapply(res, is.character, TRUE); cat("null: errori", sum(errs), "su", length(res), "\n")
  if (any(errs)) quit(status = 1)
}

# ---- stadio calib (CP-4) --------------------------------------------------------------------------
if (STAGE == "calib") {
  nuc <- load_nuclei()
  win <- read.csv("results/R2b/R2b_windows.csv")
  pts <- read.csv("results/R2b/R2b_points.csv"); pts <- pts[pts$rater == "Luca", ]
  ann <- jsonlite::fromJSON("results/R2b/annotations/R2b_annotations_Luca.json", simplifyVector = FALSE)$windows
  names(ann) <- vapply(ann, function(z) z$win_id, "")
  res <- mclapply(seq_len(nrow(win)), function(k) {
    wv <- win[k, ]; A <- wv$archetype; roi <- wv$roi_id; a <- ARCH[ARCH$archetype == A, ]
    w <- r3_roi_window(A, roi, rois)
    sqw <- st_sfc(st_polygon(list(rbind(c(wv$x0_um, wv$y0_um), c(wv$x0_um + wv$side_um, wv$y0_um),
                                         c(wv$x0_um + wv$side_um, wv$y0_um + wv$side_um), c(wv$x0_um, wv$y0_um + wv$side_um),
                                         c(wv$x0_um, wv$y0_um)))))
    reg <- st_intersection(w$poly, sqw)
    aj <- ann[[wv$win_id]]; m0 <- aj$margin_px; upp <- aj$um_per_px
    if (length(aj$exclusions)) {
      ex <- st_union(st_sfc(lapply(aj$exclusions, function(e) {
        p <- do.call(rbind, lapply(e$pts, unlist)); p <- cbind(wv$x0_um + (p[, 1] - m0) * upp, wv$y0_um + (p[, 2] - m0) * upp)
        st_make_valid(st_polygon(list(rbind(p, p[1, ]))))
      })))
      reg <- st_difference(reg, ex)
    }
    reg <- st_union(reg)
    inreg <- function(x, y) lengths(st_within(st_as_sf(data.frame(x = x, y = y), coords = c("x", "y")), reg)) > 0
    pm <- pts[pts$win_id == wv$win_id, ]; pm <- pm[inreg(pm$x_um, pm$y_um), ]
    sg <- nuc[nuc$archetype == A & nuc$roi_id == roi & nuc$method == a$primary, ]
    sg <- sg[sg$x_um >= wv$x0_um - 1 & sg$x_um <= wv$x0_um + wv$side_um + 1 &
             sg$y_um >= wv$y0_um - 1 & sg$y_um <= wv$y0_um + wv$side_um + 1, ]
    sg <- sg[inreg(sg$x_um, sg$y_um), ]
    one <- function(x, y, src) {
      tv <- r3_tessellate(x, y, reg, sqw, pad = 30)
      write_parquet(cbind(win_id = wv$win_id, source = src, tv), file.path(OUT, "calib", sprintf("%s_%s.parquet", wv$win_id, src)))
      data.frame(win_id = wv$win_id, archetype = A, source = src, area_region = as.numeric(st_area(reg)), sum_area = sum(tv$area),
                 t(r3_summary(tv)))
    }
    rbind(one(pm$x_um, pm$y_um, "manual_Luca"), one(sg$x_um, sg$y_um, a$primary))
  }, mc.cores = 24L)
  errs <- vapply(res, function(z) inherits(z, "try-error"), TRUE); if (any(errs)) stop(res[errs][[1]])
  dir.create("results/R3", showWarnings = FALSE)
  write.csv(do.call(rbind, res), file.path(if (OUT == "/mnt/micron/geo_spatialtrans/R3") "results/R3" else OUT, "R3_calibration_windows.csv"),
            row.names = FALSE)
}
