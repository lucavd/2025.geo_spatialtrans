# tools/S1.2_checkB.R — check B (B1a, B1b, B2) e controprove CP-1, CP-2, CP-3 di S1.2
# Pre-registrazione: results/S1.2/S1.2_preregistration.md (commit 60187fb). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.2_checkB.R [stage]
#   stage = sim (simulazioni, rilanciabile: salta i compiti gia' fatti) | real | b1b | manual | all (default)
# Definizioni operative (fissate prima del run; specificano la pre-registrazione senza cambiarne le soglie):
#  - finestra = maschera di tessuto R2 (results/R2/tissue_masks/<A>_<roi>_valid_ds4.png, 1/4 di risoluzione,
#    pixel = side_um / ncol); stessa finestra per reale e simulato. Regioni per seed_centroids =
#    extract_regions(maschera, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1): i poligoni coincidono
#    con i pixel della maschera.
#  - reale primario: cellpose_rgb in A1, A2, A3, A5; spaceranger in A4, A6; secondario: stardist_he (scale 1, keep).
#  - pcf: spatstat.explore::pcf, correction "translate", divisor "d", r = 0, 0.5, ..., 30 µm; larghezza di banda
#    FISSATA per ROI = quella scelta dalla regola di Stoyan sul pattern reale primario (attr "bw.used" = bw.stoyan), passata a tutte
#    le stime dello stesso ROI (reale secondario e simulati): evita che la densita' diversa cambi il lisciamento.
#  - D = trapezi su r = 0.5..30 (g(0) escluso: non stimabile).
#  - modelli: CSR (d = 0), PD (d = 2/3 eq_r), NUC (d = diametro nucleare B-003...B-038), griglia d = 0..15 passo 0.5;
#    20 repliche per (ROI, modello), densita' = totale manuale; griglia con max_attempts_factor = 50 e arresto alla
#    prima replica con n_failed > 0 (d non fattibile; d maggiori saltati per quel ROI).
#  - B1b: mappe graphclust 8 µm dei ROI (S1.1), extract_regions default (min 100 µm², tol 0.5, stride 1),
#    densita' = n_cells / area totale delle componenti (denominatore di S1.1); seed 42 primario, 10 seed per robustezza.
#  - manuale (CP-3 descrittivo): punti di Luca nelle 24 finestre R2b; pcf per finestra (r 0..15, translate, banda di
#    Stoyan del pattern manuale) e media sulle 4 finestre; stesso calcolo per il segmentatore primario nelle finestre.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel); library(spatstat.geom); library(spatstat.explore)
                                 library(arrow); library(png) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R")
args <- commandArgs(TRUE); STAGE <- if (length(args)) args[1] else "all"
OUT <- "/mnt/micron/geo_spatialtrans/S1.2/checkB"; dir.create(file.path(OUT, "sim"), recursive = TRUE, showWarnings = FALSE)
RES <- "results/S1.2"; dir.create(RES, showWarnings = FALSE)
NCORES <- 40L
R_GRID <- seq(0, 30, by = 0.5)
ARCH <- data.frame(archetype = paste0("A", 1:6),
                   tot = c(12419, 8924, 3100, 28096, 1185, 988), evid = c(10291, 7376, 2311, 26844, 953, 906),
                   d_nuc = c(4.1, 6.3, 6.3, 4.2, 11.3, 4.8),
                   primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   secondary = "stardist_he", stringsAsFactors = FALSE)
ARCH$eq_r <- sqrt(1e6 / (pi * ARCH$tot)); ARCH$d_rule <- 2 / 3 * ARCH$eq_r
write.csv(ARCH, file.path(RES, "S1.2_archetype_params.csv"), row.names = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")
rois <- rois[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px

roi_window <- function(A, roi) {
  side_um <- rois$side_um[rois$archetype == A & rois$roi_id == roi]
  m <- readPNG(sprintf("results/R2/tissue_masks/%s_%s_valid_ds4.png", A, roi))
  if (length(dim(m)) == 3) m <- m[, , 1]
  m <- m > 0.5                                      # riga 1 = alto (y verso il basso)
  stopifnot(nrow(m) == ncol(m))
  px <- side_um / ncol(m)
  W <- owin(mask = m[nrow(m):1, ], xrange = c(0, side_um), yrange = c(0, side_um))   # come tools/R2_spatial.R
  idx <- which(m, arr.ind = TRUE)                   # row = y (giu'), col = x
  clust <- data.frame(x = idx[, 2], y = idx[, 1], cl = 1L)
  reg <- extract_regions(clust, pixel_size_um = px, min_region_area_um2 = 0, cluster_col = "cl",
                         simplify_tol_um = 0, stride = 1L, verbose = FALSE)
  list(W = W, reg = reg, side_um = side_um, px = px, area_um2 = area(W))
}
to_ppp <- function(x, y_down, w) {
  ok <- inside.owin(x, w$side_um - y_down, w$W)
  list(X = ppp(x[ok], w$side_um - y_down[ok], window = w$W, check = FALSE), n_lost = sum(!ok))
}
pcf_fixed <- function(X, bw = NULL, r = R_GRID) {
  g <- if (is.null(bw)) pcf(X, r = r, correction = "translate", divisor = "d")
       else pcf(X, r = r, correction = "translate", divisor = "d", bw = bw)
  list(g = g$trans, bw = attr(g, "bw.used"))
}

# ---- reale: g(r) primario e secondario per ROI ---------------------------------
real_file <- file.path(OUT, "real_pcf.rds")
if (STAGE %in% c("real", "all", "sim") && !file.exists(real_file)) {
  nuc <- read_parquet("results/R2/R2_nuclei_all.parquet", col_select = c("archetype", "roi_id", "method", "scale", "keep", "x_um", "y_um"))
  nuc <- nuc[nuc$scale == 1 & nuc$keep, ]
  jobs <- split(rois[, c("archetype", "roi_id")], seq_len(nrow(rois)))
  real <- mclapply(jobs, function(j) {
    A <- j$archetype; roi <- j$roi_id; w <- roi_window(A, roi); a <- ARCH[ARCH$archetype == A, ]
    out <- list(archetype = A, roi_id = roi, area_mm2 = w$area_um2 / 1e6, n_regions = nrow(w$reg$region_df))
    for (role in c("primary", "secondary")) {
      d <- nuc[nuc$archetype == A & nuc$roi_id == roi & nuc$method == a[[role]], ]
      pp <- to_ppp(d$x_um, d$y_um, w)
      if (role == "primary") { pg <- pcf_fixed(pp$X); out$bw <- pg$bw } else pg <- pcf_fixed(pp$X, bw = out$bw)
      out[[role]] <- list(method = a[[role]], n = pp$X$n, n_lost = pp$n_lost, lambda = pp$X$n / area(w$W) * 1e6, g = pg$g)
    }
    out
  }, mc.cores = min(NCORES, 30L))
  saveRDS(real, real_file)
}

# ---- simulazioni ---------------------------------------------------------------
if (STAGE %in% c("sim", "all")) {
  real <- readRDS(real_file)
  bw_of <- setNames(vapply(real, function(z) z$bw, numeric(1)), vapply(real, function(z) paste(z$archetype, z$roi_id), ""))
  models <- function(a) {
    rbind(data.frame(model = c("CSR", "PD", "NUC"), d = c(0, a$d_rule, a$d_nuc), grid = FALSE),
          data.frame(model = sprintf("G%04.1f", seq(0, 15, 0.5)), d = seq(0, 15, 0.5), grid = TRUE))
  }
  tasks <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) {
    a <- ARCH[ARCH$archetype == rois$archetype[i], ]; m <- models(a)
    data.frame(roi_idx = i, archetype = rois$archetype[i], roi_id = rois$roi_id[i], m, model_idx = seq_len(nrow(m)))
  }))
  # una unita' di lavoro = (ROI, blocco di modelli) per ammortizzare la costruzione della finestra
  units <- split(tasks, paste(tasks$roi_idx, ifelse(tasks$grid, "grid", "main")))
  todo <- units[!file.exists(file.path(OUT, "sim", paste0(gsub(" ", "_", names(units)), ".rds")))]
  cat(sprintf("[sim] %d unita' da fare su %d\n", length(todo), length(units)))
  invisible(mclapply(names(todo), function(u) {
    tk <- todo[[u]]; A <- tk$archetype[1]; roi <- tk$roi_id[1]; a <- ARCH[ARCH$archetype == A, ]
    w <- roi_window(A, roi); bw <- bw_of[[paste(A, roi)]]
    comp <- data.frame(cluster_id = 1L, cell_type = "c", fraction = 1)
    res <- list(); infeasible <- FALSE
    for (k in seq_len(nrow(tk))) {
      if (tk$grid[k] && infeasible) { res[[k]] <- data.frame(tk[k, ], reps = 0L, feasible = FALSE, n_mean = NA, dens_mean = NA, dens_min = NA, dens_max = NA, lost_mean = NA, t(rep(NA_real_, length(R_GRID)))); next }
      G <- matrix(NA_real_, 20, length(R_GRID)); nn <- numeric(20); lost <- numeric(20); feas <- TRUE; reps <- 0L
      for (rp in 1:20) {
        seed <- 1e6 * tk$roi_idx[k] + 1e3 * tk$model_idx[k] + rp
        o <- seed_centroids(w$reg, data.frame(cell_type = "c", density = a$tot, min_dist_um = tk$d[k]), comp,
                            random_seed = seed, max_attempts_factor = if (tk$grid[k]) 50 else 1000, verbose = FALSE)
        if (o$info$n_failed > 0) { feas <- FALSE; if (tk$grid[k]) { infeasible <- TRUE; break } }
        pp <- to_ppp(o$centroids$x, o$centroids$y, w)
        G[rp, ] <- pcf_fixed(pp$X, bw = bw)$g; nn[rp] <- o$info$n_cells; lost[rp] <- pp$n_lost; reps <- reps + 1L
      }
      dens <- nn[seq_len(reps)] / w$area_um2 * 1e6
      res[[k]] <- data.frame(tk[k, ], reps = reps, feasible = feas, n_mean = mean(nn[seq_len(reps)]),
                             dens_mean = mean(dens), dens_min = min(dens), dens_max = max(dens), lost_mean = mean(lost[seq_len(reps)]),
                             t(colMeans(G[seq_len(max(reps, 1)), , drop = FALSE])))
    }
    out <- do.call(rbind, res); names(out)[(ncol(out) - length(R_GRID) + 1):ncol(out)] <- sprintf("g_%04.1f", R_GRID)
    saveRDS(out, file.path(OUT, "sim", paste0(gsub(" ", "_", u), ".rds")))
    NULL
  }, mc.cores = NCORES, mc.preschedule = FALSE))
}

# ---- B1b: mappe graphclust dei ROI ----------------------------------------------
if (STAGE %in% c("b1b", "all")) {
  IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
  b1 <- do.call(rbind, mclapply(seq_len(nrow(rois)), function(i) {
    A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]
    o <- readRDS(file.path(IN, sprintf("roi_%s_%s.rds", A, roi)))
    reg <- extract_regions(o$clust, pixel_size_um = o$pixel_size_um, cluster_col = o$cluster_col, stride = 1L, verbose = FALSE)
    comp <- data.frame(cluster_id = unique(as.character(reg$region_df$cluster_id)), cell_type = "c", fraction = 1)
    dens <- vapply(c(42, 1:9), function(s) seed_centroids(reg, data.frame(cell_type = "c", density = a$tot), comp,
                   random_seed = s, verbose = FALSE)$info$n_cells / reg$info$area_total_um2 * 1e6, numeric(1))
    data.frame(archetype = A, roi_id = roi, frac_excluded = reg$info$frac_area_excluded, area_total_mm2 = reg$info$area_total_um2 / 1e6,
               dens_seed42 = dens[1], dens_min10 = min(dens), dens_max10 = max(dens),
               expected = a$tot * (1 - reg$info$frac_area_excluded), evid = a$evid, tot = a$tot)
  }, mc.cores = 30L))
  write.csv(b1, file.path(RES, "S1.2_B1b.csv"), row.names = FALSE)
}

# ---- manuale R2b (descrittivo CP-3) ----------------------------------------------
# pcf per finestra con ratio = TRUE e banda di Stoyan del pattern manuale della finestra (uguale per le tre fonti),
# poi spatstat pool() sulle 4 finestre (pooling dei rapporti, come da pre-registrazione). In piu', riassunto senza
# lisciamento: K(r)/(pi r^2) a r = 2, 3, 5 µm (Kest, translate) per finestra e pool.
if (STAGE %in% c("manual", "all")) {
  pts <- read.csv("results/R2b/R2b_points.csv"); win <- read.csv("results/R2b/R2b_windows.csv")
  pts <- pts[pts$rater == "Luca", ]
  nuc <- read_parquet("results/R2/R2_nuclei_all.parquet", col_select = c("archetype", "roi_id", "method", "scale", "keep", "x_um", "y_um"))
  nuc <- nuc[nuc$scale == 1 & nuc$keep, ]
  rm_ <- seq(0, 15, 0.5); rk <- c(0, 2, 3, 5)
  pat <- lapply(seq_len(nrow(win)), function(i) {
    wv <- win[i, ]; a <- ARCH[ARCH$archetype == wv$archetype, ]
    W <- owin(c(wv$x0_um, wv$x0_um + wv$side_um), c(wv$y0_um, wv$y0_um + wv$side_um))
    mk <- function(x, y) { ok <- inside.owin(x, y, W); ppp(x[ok], y[ok], window = W, check = FALSE) }
    p <- pts[pts$win_id == wv$win_id, ]
    srcs <- list(manual = mk(p$x_um, p$y_um))
    for (meth in unique(c(a$primary, a$secondary))) { d <- nuc[nuc$archetype == wv$archetype & nuc$roi_id == wv$roi_id & nuc$method == meth, ]; srcs[[meth]] <- mk(d$x_um, d$y_um) }
    bw <- attr(pcf(srcs$manual, r = rm_, correction = "translate", divisor = "d"), "bw.used")
    list(archetype = wv$archetype, win_id = wv$win_id, srcs = srcs, bw = bw,
         g = lapply(srcs, function(X) pcf(X, r = rm_, correction = "translate", divisor = "d", bw = bw, ratio = TRUE)),
         K = lapply(srcs, function(X) Kest(X, r = rk, correction = "translate", ratio = TRUE)))
  })
  man <- list(); kk <- list()
  for (A in ARCH$archetype) {
    pa <- pat[vapply(pat, function(z) z$archetype == A, TRUE)]
    for (src in names(pa[[1]]$srcs)) {
      gp <- do.call(pool, lapply(pa, function(z) z$g[[src]])); Kp <- do.call(pool, lapply(pa, function(z) z$K[[src]]))
      man[[length(man) + 1]] <- data.frame(archetype = A, source = src, n = sum(vapply(pa, function(z) z$srcs[[src]]$n, 0)), r = rm_, g = gp$pooltrans)
      kk[[length(kk) + 1]] <- data.frame(archetype = A, source = src, r = rk[-1], K_over_pir2 = Kp$pooltrans[-1] / (pi * rk[-1]^2))
    }
    for (z in pa) for (src in names(z$srcs)) kk[[length(kk) + 1]] <- data.frame(archetype = A, source = paste0(src, "@", z$win_id), r = rk[-1],
                                                                                 K_over_pir2 = z$K[[src]]$trans[-1] / (pi * rk[-1]^2))
  }
  write.csv(do.call(rbind, man), file.path(RES, "S1.2_manual_pcf.csv"), row.names = FALSE)
  write.csv(do.call(rbind, kk), file.path(RES, "S1.2_manual_K.csv"), row.names = FALSE)
}
cat("[S1.2_checkB] stage", STAGE, "finito\n")
