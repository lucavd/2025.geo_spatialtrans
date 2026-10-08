# R/testing/test_R3.R — check C di R3 (C-R3.1…C-R3.8). Pre-registrazione: results/R3/R3_preregistration.md.
# Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_R3.R synthetic   # C-R3.3a/b, 3.4, 3.1s/3.2s, 3.3c/3.5 (python)
#   ... test_R3.R real [outdir]      # C-R3.1, 3.2, 3.8 sui 30 ROI (dopo tools/R3_run.R real + R3_containment.py all)
#   ... test_R3.R repro <dir1> <dir2> # C-R3.6
# Stampa PASS/FAIL per asserzione; scrive results/R3/R3_test_<mode>[_<mutante>].csv; exit 1 se qualche FAIL.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(arrow) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R")
args <- commandArgs(TRUE); MODE <- if (length(args)) args[1] else "synthetic"
MUT <- r3_mutant(); RES <- "results/R3"; dir.create(RES, showWarnings = FALSE)
rows <- list()
chk <- function(id, case, value, pass, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(check = id, case = case, value = signif(value, 8), result = if (isTRUE(pass)) "PASS" else "FAIL", note = note)
  cat(sprintf("%-4s %-8s %-28s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", id, case, format(signif(value, 6))))
}
sq <- function(x0, y0, s) st_sfc(st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0)))))

if (MODE == "synthetic") {
  # C-R3.3a reticolo quadrato (il quadrato di ritaglio taglia le celle di bordo) ed esagonale
  a <- 5; g <- expand.grid(x = seq(1, 99, by = a), y = seq(1, 99, by = a))
  tv <- r3_tessellate(g$x, g$y, sq(0, 0, 100), sq(0, 0, 100)); i <- tv$interior
  chk("C-R3.3a", "square_area", max(abs(tv$area[i] / a^2 - 1)), sum(i) > 0 && max(abs(tv$area[i] / a^2 - 1)) < 1e-9)
  chk("C-R3.3a", "square_ecc", max(tv$ecc_T[i]), max(tv$ecc_T[i]) < 1e-6)
  chk("C-R3.3a", "square_sides", mean(tv$nsides[i] == 4), all(tv$nsides[i] == 4))
  hx <- expand.grid(i = 0:24, j = 0:28); x <- 1 + a * (hx$i + 0.5 * (hx$j %% 2)); y <- 1 + a * sqrt(3) / 2 * hx$j
  tv <- r3_tessellate(x, y, sq(0, 0, 120), sq(0, 0, 120)); i <- tv$interior
  chk("C-R3.3a", "hex_area", max(abs(tv$area[i] / (sqrt(3) / 2 * a^2) - 1)), sum(i) > 0 && max(abs(tv$area[i] / (sqrt(3) / 2 * a^2) - 1)) < 1e-9)
  chk("C-R3.3a", "hex_ecc", max(tv$ecc_T[i]), max(tv$ecc_T[i]) < 1e-6)
  chk("C-R3.3a", "hex_sides", mean(tv$nsides[i] == 6), all(tv$nsides[i] == 6))
  # C-R3.3b ellisse poligonale (1 024 vertici)
  for (phi_deg in c(0, 30, 60, -30, 75)) {
    phi <- phi_deg * pi / 180; t <- seq(0, 2 * pi, length.out = 1025)[-1]; ea <- 6; eb <- 3
    ex <- 20 + ea * cos(t) * cos(phi) - eb * sin(t) * sin(phi); ey <- 20 + ea * cos(t) * sin(phi) + eb * sin(t) * cos(phi)
    mo <- r3_poly_moments(ex, ey); sh <- r3_shape(mo[["mxx"]], mo[["myy"]], mo[["mxy"]])
    de <- abs(sh$ecc - sqrt(1 - (eb / ea)^2)); dth <- r3_axis_diff(sh$theta, phi) * 180 / pi
    chk("C-R3.3b", sprintf("ellipse_ecc_phi%d", phi_deg), de, de < 1e-3)
    chk("C-R3.3b", sprintf("ellipse_theta_phi%d", phi_deg), dth, dth < 0.5)
  }
  # C-R3.4 CSR 20 000 punti (riferimento: gamma a = b = 7/2, doi:10.1016/j.physa.2007.07.063; Euler)
  set.seed(20261008); n <- 20000; L <- sqrt(n / 5000) * 1000; x <- runif(n, 0, L); y <- runif(n, 0, L)
  tv <- r3_tessellate(x, y, sq(0, 0, L), sq(0, 0, L)); i <- tv$interior; lam <- n / L^2
  chk("C-R3.4", "csr_mean_area_x_lambda", mean(tv$area[i]) * lam, abs(mean(tv$area[i]) * lam - 1) < 0.01)
  vn <- var(tv$area[i] / mean(tv$area[i])); chk("C-R3.4", "csr_var_norm_area", vn, abs(vn - 2 / 7) < 0.015)
  chk("C-R3.4", "csr_mean_nsides", mean(tv$nsides[i]), abs(mean(tv$nsides[i]) - 6) < 0.02)
  # C-R3.1s / C-R3.2s: ritaglio non convesso con buco (anello quadrato) + CSR
  ring <- st_difference(sq(0, 0, 300), sq(100, 100, 100)); set.seed(7)
  p <- st_coordinates(st_sample(ring, 2000)); tv <- r3_tessellate(p[, 1], p[, 2], ring, sq(0, 0, 300))
  chk("C-R3.1", "synthetic_ring_area", abs(sum(tv$area) / as.numeric(st_area(ring)) - 1), abs(sum(tv$area) / as.numeric(st_area(ring)) - 1) < 1e-9)
  chk("C-R3.2", "synthetic_ring_gen_in_own", mean(tv$gen_in_own), all(tv$gen_in_own) && nrow(tv) == 2000)
  # C-R3.3c / C-R3.5 (python)
  py <- if (MUT == "M3") "--mutant M3" else ""
  system(sprintf(".venv/bin/python tools/R3_containment.py selftest %s > /dev/null", py))
  st <- read.csv(file.path(RES, sprintf("R3_selftest_py%s.csv", if (MUT == "M3") "_M3" else "")))
  for (k in seq_len(nrow(st))) chk(st$check[k], st$case[k], st$value[k], isTRUE(as.logical(st$passed[k])))
}

if (MODE == "real") {
  OUT <- if (length(args) > 1) args[2] else "/mnt/micron/geo_spatialtrans/R3"
  suppressPackageStartupMessages({ library(spatstat.geom) })
  rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
  prim <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
  f <- list.files(file.path(OUT, "real"), "_roi\\.rds$", full.names = TRUE)
  rr <- do.call(rbind, lapply(f, readRDS))
  n_prim <- sum(rr$method == prim[rr$archetype])
  chk("C-R3.8", "roi_primary_30", n_prim, n_prim == 30 && nrow(rois) == 30)
  e1 <- abs(rr$sum_area / rr$area_poly_um2 - 1); e2 <- abs(rr$area_poly_um2 / rr$area_px_um2 - 1)
  chk("C-R3.1", "real_sum_vs_poly_max", max(e1), all(e1 < 1e-9), sprintf("%d tassellazioni", nrow(rr)))
  chk("C-R3.1", "real_poly_vs_px_max", max(e2), all(e2 < 1e-6))
  chk("C-R3.2", "real_gen_in_own_min", min(rr$gen_in_own), all(rr$gen_in_own == 1))
  nuc <- as.data.frame(read_parquet("results/R2/R2_nuclei_all.parquet",
                                    col_select = c("archetype", "roi_id", "method", "scale", "keep", "x_um", "y_um")))
  nuc <- nuc[nuc$scale == 1 & nuc$keep, ]
  bad_n <- 0; bad_cells <- 0; dx <- 0; da <- 0; de <- 0; miss <- 0
  for (k in seq_len(nrow(rr))) {
    A <- rr$archetype[k]; roi <- rr$roi_id[k]; me <- rr$method[k]
    side <- rois$side_um[rois$archetype == A & rois$roi_id == roi]
    m <- png::readPNG(sprintf("results/R2/tissue_masks/%s_%s_valid_ds4.png", A, roi)); if (length(dim(m)) == 3) m <- m[, , 1]
    W <- owin(mask = (m > 0.5)[nrow(m):1, ], xrange = c(0, side), yrange = c(0, side))   # percorso indipendente (spatstat)
    d <- nuc[nuc$archetype == A & nuc$roi_id == roi & nuc$method == me, ]
    n_ind <- sum(inside.owin(d$x_um, side - d$y_um, W))
    if (n_ind != rr$n[k]) bad_n <- bad_n + 1
    cells <- read_parquet(file.path(OUT, "real", sprintf("%s_%s_%s_cells.parquet", A, roi, me)))
    if (nrow(cells) != rr$n[k]) bad_cells <- bad_cells + 1
    gen <- read_parquet(file.path(OUT, "real", sprintf("%s_%s_%s_gen.parquet", A, roi, me)))
    fn <- file.path(OUT, "real", sprintf("%s_%s_%s_nuc.parquet", A, roi, me))
    if (!file.exists(fn)) { miss <- miss + 1; next }
    nu <- read_parquet(fn); nu <- nu[match(gen$idx, nu$idx), ]
    dx <- max(dx, abs(nu$cx_px - gen$x), abs(nu$cy_px - gen$y), na.rm = FALSE)
    da <- max(da, abs(nu$area_px_um2 - gen$area_nuc)); de <- max(de, abs(nu$ecc_N - gen$ecc_parquet))
  }
  chk("C-R3.8", "n_used_eq_independent_count", bad_n, bad_n == 0, "inside.owin su maschera spatstat")
  chk("C-R3.2", "n_territories_eq_n", bad_cells, bad_cells == 0)
  chk("C-R3.8", "nuc_files_present", miss, miss == 0)
  chk("C-R3.8", "mask_centroid_vs_parquet_um", dx, is.finite(dx) && dx < 1e-6, "maschere native = maschere di R2")
  chk("C-R3.8", "mask_area_vs_parquet_um2", da, is.finite(da) && da < 1e-9)
  chk("C-R3.8", "mask_ecc_vs_parquet", de, is.finite(de) && de < 1e-6)
  ins <- c("results/R2/R2_nuclei_all.parquet", "results/R2/R2_rois_checked.csv", "results/R2b/R2b_points.csv",
           "results/R2b/R2b_windows.csv", "results/R2b/annotations/R2b_annotations_Luca.json",
           Sys.glob("results/R2/tissue_masks/*_valid_ds4.png"),
           Sys.glob("/mnt/micron/geo_spatialtrans/R2/masks/*_native.npz"))
  writeLines(paste(tools::md5sum(ins), ins), file.path(RES, "R3_inputs_md5.txt"))
  chk("C-R3.8", "inputs_md5_recorded", length(ins), length(ins) > 100)
}

if (MODE == "repro") {
  d1 <- args[2]; d2 <- args[3]
  for (sub in c("real", "null", "calib")) {
    f1 <- sort(list.files(file.path(d1, sub), "(_cells\\.parquet|_nuc\\.parquet|\\.rds|\\.parquet)$"))
    f2 <- sort(list.files(file.path(d2, sub), "(_cells\\.parquet|_nuc\\.parquet|\\.rds|\\.parquet)$"))
    common <- intersect(f1, f2)
    same <- tools::md5sum(file.path(d1, sub, common)) == tools::md5sum(file.path(d2, sub, common))
    chk("C-R3.6", sprintf("%s_md5_identical", sub), sum(!same), length(common) > 0 && all(same), sprintf("%d file confrontati", length(common)))
  }
}

out <- do.call(rbind, rows)
write.csv(out, file.path(RES, sprintf("R3_test_%s%s.csv", MODE, if (MUT != "none") paste0("_", MUT) else "")), row.names = FALSE)
cat(sprintf("== %s%s: %d PASS, %d FAIL\n", MODE, if (MUT != "none") paste0(" [mutante ", MUT, "]") else "", sum(out$result == "PASS"), sum(out$result == "FAIL")))
if (any(out$result == "FAIL")) quit(status = 1)
