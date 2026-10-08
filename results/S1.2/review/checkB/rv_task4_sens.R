# rv_task4_sens.R — revisione avversariale S1.2 / checkB, compito 4: sensibilita' alla banda.
# Rigenera le simulazioni con seed_centroids() del progetto (stessi seed, stessi modelli, stesse regole di arresto
# di tools/S1.2_checkB.R) e stima g(r) per ogni pattern con 4 bande: bw_ROI x1 (come il progetto), x0.5, x2, e
# "propria" (regola di Stoyan sul pattern stesso, default spatstat). Lo stimatore e' scritto qui: coppie da
# spatstat.geom::closepairs, pesi di traslazione da spatstat.explore::edge.Trans, kernel di Epanechnikov
# (sd = bw, semi-ampiezza sqrt(5) bw), divisore d, rinormalizzazione della massa del kernel su [0, inf) — verificato
# contro spatstat::pcf in rv_task1 (|diff| <= 0.002) e qui sul primo pattern di ogni unita' (colonna chk_maxabs).
# Le stesse 4 bande si applicano al riferimento reale (primario e secondario).
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task4_sens.R [test|run] [ncores]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel); library(spatstat.geom); library(spatstat.explore); library(arrow); library(png) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R")
args <- commandArgs(TRUE); MODE <- if (length(args)) args[1] else "run"; NC <- if (length(args) > 1) as.integer(args[2]) else 30L
SO <- "results/S1.2/review/checkB/sens"; dir.create(SO, recursive = TRUE, showWarnings = FALSE)
R_GRID <- seq(0, 30, by = 0.5); VAR <- c(x1 = 1, x0.5 = 0.5, x2 = 2, own = NA)
ARCH <- data.frame(archetype = paste0("A", 1:6), tot = c(12419, 8924, 3100, 28096, 1185, 988),
                   d_nuc = c(4.1, 6.3, 6.3, 4.2, 11.3, 4.8),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"), stringsAsFactors = FALSE)
ARCH$d_rule <- 2 / 3 * sqrt(1e6 / (pi * ARCH$tot))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px
real_ref <- readRDS("/mnt/micron/geo_spatialtrans/S1.2/checkB/real_pcf.rds")
bw_of <- setNames(vapply(real_ref, function(z) z$bw, 0), vapply(real_ref, function(z) paste(z$archetype, z$roi_id), ""))

mkwin <- function(i) {
  m <- readPNG(sprintf("results/R2/tissue_masks/%s_%s_valid_ds4.png", rois$archetype[i], rois$roi_id[i]))
  if (length(dim(m)) == 3) m <- m[, , 1]
  m <- m > 0.5; s <- rois$side_um[i]; px <- s / ncol(m)
  W <- owin(mask = m[nrow(m):1, ], xrange = c(0, s), yrange = c(0, s))
  idx <- which(m, arr.ind = TRUE)
  reg <- extract_regions(data.frame(x = idx[, 2], y = idx[, 1], cl = 1L), pixel_size_um = px, min_region_area_um2 = 0,
                         cluster_col = "cl", simplify_tol_um = 0, stride = 1L, verbose = FALSE)
  list(W = W, reg = reg, s = s, A = area(W))
}
pp <- function(x, y_down, w) { y <- w$s - y_down; ok <- inside.owin(x, y, w$W); ppp(x[ok], y[ok], window = w$W, check = FALSE) }
gmulti <- function(X, bw0) {
  n <- npoints(X); lam <- n / area(Window(X)); bws <- ifelse(is.na(VAR), 0.15 / sqrt(lam) / sqrt(5), VAR * bw0)
  cp <- closepairs(X, rmax = 30 + sqrt(5) * max(bws) + 0.5, what = "all", twice = TRUE)
  e <- edge.Trans(dx = cp$dx, dy = cp$dy, W = Window(X), paired = TRUE)
  l2a <- n * (n - 1) / area(Window(X)); wd <- e / (2 * pi * cp$d)
  out <- sapply(bws, function(bw) { h <- sqrt(5) * bw
    vapply(R_GRID, function(r) { u <- (r - cp$d) / h; s <- abs(u) < 1
      lo <- max(-1, -r / h); mass <- 0.75 * ((1 - lo) - (1 - lo^3) / 3)
      sum(0.75 / h * (1 - u[s]^2) * wd[s]) / l2a / mass }, 0) })
  colnames(out) <- names(VAR); out
}
models <- function(a) rbind(data.frame(model = c("CSR", "PD", "NUC"), d = c(0, a$d_rule, a$d_nuc), grid = FALSE),
                            data.frame(model = sprintf("G%04.1f", seq(0, 15, 0.5)), d = seq(0, 15, 0.5), grid = TRUE))
run_unit <- function(i, block, nrep = 20L, mods = NULL) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]; w <- mkwin(i); bw0 <- bw_of[[paste(A, roi)]]
  m <- models(a); m$model_idx <- seq_len(nrow(m)); m <- m[m$grid == (block == "grid"), ]
  if (!is.null(mods)) m <- m[m$model %in% mods, ]
  res <- list(); infeasible <- FALSE; chk <- NA_real_
  if (block == "main") {   # riferimento reale, stesse 4 bande
    nuc <- read_parquet("results/R2/R2_nuclei_all.parquet", col_select = c("archetype", "roi_id", "method", "scale", "keep", "x_um", "y_um"))
    nuc <- nuc[nuc$scale == 1 & nuc$keep & nuc$archetype == A & nuc$roi_id == roi, ]
    for (role in c("primary", "secondary")) { meth <- if (role == "primary") a$primary else "stardist_he"
      d <- nuc[nuc$method == meth, ]; X <- pp(d$x_um, d$y_um, w); G <- gmulti(X, bw0)
      res[[length(res) + 1]] <- data.frame(archetype = A, roi_id = roi, source = paste0("real_", role), d = NA, feasible = TRUE, reps = 1L, n_mean = npoints(X),
                                           variant = rep(colnames(G), each = length(R_GRID)), r = R_GRID, g = as.vector(G)) }
  }
  for (k in seq_len(nrow(m))) {
    if (m$grid[k] && infeasible) { res[[length(res) + 1]] <- data.frame(archetype = A, roi_id = roi, source = m$model[k], d = m$d[k], feasible = FALSE, reps = 0L,
                                                                       n_mean = NA, variant = rep(names(VAR), each = length(R_GRID)), r = R_GRID, g = NA); next }
    S <- array(NA_real_, c(nrep, length(R_GRID), length(VAR))); feas <- TRUE; reps <- 0L; nn <- numeric(0)
    for (rp in seq_len(nrep)) {
      seed <- 1e6 * i + 1e3 * m$model_idx[k] + rp
      o <- seed_centroids(w$reg, data.frame(cell_type = "c", density = a$tot, min_dist_um = m$d[k]), data.frame(cluster_id = 1L, cell_type = "c", fraction = 1),
                          random_seed = seed, max_attempts_factor = if (m$grid[k]) 50 else 1000, verbose = FALSE)
      if (o$info$n_failed > 0) { feas <- FALSE; if (m$grid[k]) { infeasible <- TRUE; break } }
      X <- pp(o$centroids$x, o$centroids$y, w); S[rp, , ] <- gmulti(X, bw0); reps <- reps + 1L; nn <- c(nn, npoints(X))
      if (is.na(chk)) chk <- max(abs(pcf(X, r = R_GRID, correction = "translate", divisor = "d", bw = bw0)$trans[-1] - S[rp, -1, 1]))
    }
    Gm <- apply(S[seq_len(max(reps, 1)), , , drop = FALSE], c(2, 3), mean)
    res[[length(res) + 1]] <- data.frame(archetype = A, roi_id = roi, source = m$model[k], d = m$d[k], feasible = feas, reps = reps, n_mean = mean(nn),
                                         variant = rep(names(VAR), each = length(R_GRID)), r = R_GRID, g = as.vector(Gm))
  }
  out <- do.call(rbind, res); out$chk_maxabs <- chk; out$bw0 <- bw0; out
}
if (MODE == "test") {
  t0 <- proc.time()[["elapsed"]]; o <- run_unit(11L, "main", nrep = 2L, mods = c("CSR", "PD"))
  print(aggregate(g ~ source + variant, o[o$r %in% c(2, 6, 30), ], function(z) round(z, 3)))
  cat("chk_maxabs vs spatstat::pcf:", unique(o$chk_maxabs), " elapsed", proc.time()[["elapsed"]] - t0, "s\n")
} else {
  units <- expand.grid(i = seq_len(nrow(rois)), block = c("grid", "main"), stringsAsFactors = FALSE)
  units$f <- file.path(SO, sprintf("%02d_%s.rds", units$i, units$block)); units <- units[!file.exists(units$f), ]
  cat(sprintf("[sens] %d unita' da fare\n", nrow(units)))
  invisible(mclapply(seq_len(nrow(units)), function(u) { o <- run_unit(units$i[u], units$block[u]); saveRDS(o, units$f[u]); NULL },
                     mc.cores = NC, mc.preschedule = FALSE))
  cat("[sens] fatto\n")
}
