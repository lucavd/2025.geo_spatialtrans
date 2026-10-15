.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")

# recalc — punto 4: ricostruzione indipendente di V-b (s = 1) da S1.4 (derive_cell_geometry) per A1, A4, A6 fold 1 (+ altri fold se richiesti)
suppressPackageStartupMessages({ library(sf); library(geos); library(spatstat.geom) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R"); source("tools/S1.4_mutants.R")
E <- dcg_env()
OUTD <- "results/S1.4b/review/recalc"
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px
ksd <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b)))
  max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
args <- commandArgs(TRUE); jobs <- if (length(args)) args else c("A1:1", "A4:1", "A6:1")
meth <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
res <- list()
for (jb in jobs) {
  A <- sub(":.*", "", jb); k <- as.integer(sub(".*:", "", jb)); me <- meth[[A]]; s_idx <- 1L
  L14 <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.4/loo/%s_%s_fold%d.rds", A, me, k))
  Lb  <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.4b/loo/%s_%s_fold%d.rds", A, me, k))
  roi <- L14$roi_test; stopifnot(roi == Lb$roi_test)
  te <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.4/real/%s_%s_%s.rds", A, roi, me)); n <- te$n
  fp <- L14$fits$M2
  ct <- data.frame(cell_type = "c", density = 1e6 / fp[["a_ref"]], nucleus_to_eq_ratio = fp[["rho"]],
                   nucleus_area_slope = fp[["beta"]], nucleus_area_sdlog = fp[["sigma"]])
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = te$cells$x, y = te$cells$y, cell_type = "c")
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = te$cells$a_t))
  seed <- L14$seed + 100000 * s_idx
  o <- E$derive_cell_geometry(cen, tess, ct, random_seed = seed, verbose = FALSE)
  cd <- o$cell_df
  ev <- te$cells$interior & te$cells$matched & is.finite(te$cells$a_n) & te$cells$a_n > 0
  real_an <- te$cells$a_n[ev]; real_nc <- real_an / te$cells$a_t[ev]
  real_sp <- suppressWarnings(cor(real_an, te$cells$a_t[ev], method = "spearman"))
  # controllo C-0 della mia ricostruzione V-c vs S1.4 seeds[s,]
  vc_ksnc <- ksd(cd$nucleus_area[ev] / te$cells$a_t[ev], real_nc)
  # V-b: mio codice
  At <- cd$nucleus_area_target; kp <- sqrt(2 * pi / (32 * sin(2 * pi / 32)))
  R <- kp * sqrt(At / pi)
  X <- te$cells$x; Y <- te$cells$y
  win <- owin(range(X) + c(-1, 1), range(Y) + c(-1, 1))
  P <- ppp(X, Y, window = win, check = FALSE)
  rmax <- 2 * max(R) / 0.99
  cp <- closepairs(P, rmax, what = "indices", twice = TRUE)
  dij <- sqrt((X[cp$i] - X[cp$j])^2 + (Y[cp$i] - Y[cp$j])^2)
  q <- 0.99 * dij / (R[cp$i] + R[cp$j])
  s_pair <- rep(1, n); mq <- tapply(q, cp$i, min); s_pair[as.integer(names(mq))] <- pmin(1, mq)
  w <- r3_roi_window(A, roi, rois)
  bd <- geos_boundary(as_geos_geometry(w$poly))
  dreg <- geos_distance(geos_make_point(X, Y), bd)
  s_reg <- pmin(1, 0.99 * dreg / R)
  s <- pmin(1, s_pair, s_reg)
  Af <- s^2 * At
  # verifica non sovrapposizione
  viol <- sum(s[cp$i] * R[cp$i] + s[cp$j] * R[cp$j] > dij * (1 + 1e-12))
  m1 <- Lb$m[Lb$m$seed_idx == s_idx & Lb$m$variant == "V-b", ]
  mc <- Lb$m[Lb$m$seed_idx == s_idx & Lb$m$variant == "V-c", ]
  b7_mine <- mean(sqrt(Af[ev] / pi) > te$cells$d_gen[ev])
  res[[jb]] <- data.frame(arch = A, fold = k, roi = roi, n = n, n_ev = sum(ev), seed = seed,
    vc_ksnc_mine = vc_ksnc, vc_ksnc_S14 = L14$seeds[s_idx, "ks_nc"], vc_ksnc_S14b = mc$ks_nc,
    vc_trunc_mine = mean(cd$nucleus_clamped), vc_trunc_S14b = mc$trunc_all,
    ks_an_mine = ksd(Af[ev], real_an), ks_an_raw = m1$ks_an,
    ks_nc_mine = ksd(Af[ev] / te$cells$a_t[ev], real_nc), ks_nc_raw = m1$ks_nc,
    d_sp_mine = abs(suppressWarnings(cor(Af[ev], te$cells$a_t[ev], method = "spearman")) - real_sp), d_sp_raw = m1$d_sp,
    trunc_all_mine = mean(s < 1), trunc_all_raw = m1$trunc_all,
    trunc_ev_mine = mean(s[ev] < 1), trunc_ev_raw = m1$trunc_ev,
    loss_all_mine = mean(1 - s^2), loss_all_raw = m1$loss_all,
    b7_mine = b7_mine, b7_raw = m1$b7,
    trunc_pair_mine = mean(s_pair < 1), trunc_pair_raw = Lb$checks$trunc_pair[Lb$checks$seed_idx == s_idx],
    trunc_reg_mine = mean(s_reg < 1), trunc_reg_raw = Lb$checks$trunc_reg[Lb$checks$seed_idx == s_idx],
    n_pairs_mine = length(cp$i) / 2, n_pairs_raw = Lb$checks$n_pairs[Lb$checks$seed_idx == s_idx],
    n_pairs_overlap_before = sum(R[cp$i] + R[cp$j] > 0.99 * dij) / 2,
    viol_after = viol, c3b_mine = NA, c3b_raw = Lb$checks$c3_b_viol[Lb$checks$seed_idx == s_idx])
  # C3 di V-b: nucleo (32-gono, raggio s*R, centro generatore) fuori dal territorio di Voronoi
  # distanza dal bordo del territorio: d_gen; nucleo circoscritto dentro sse s*R <= d_gen (condizione sufficiente), poligono esatto via geos
  th <- seq(0, 2 * pi, length.out = 33)[-33]
  polys <- geos_make_polygon(x = as.vector(outer(th, seq_len(n), function(t, i) X[i] + s[i] * R[i] * cos(t))),
                             y = as.vector(outer(th, seq_len(n), function(t, i) Y[i] + s[i] * R[i] * sin(t))),
                             feature_id = rep(seq_len(n), each = 32))
  res[[jb]]$c3b_mine <- sum(!geos_within(polys, as_geos_geometry(te$terr)))
  cat(jb, "fatto\n")
}
out <- do.call(rbind, res)
fn <- file.path(OUTD, sprintf("recalc_vb_rebuild_%s.csv", paste(gsub(":", "f", jobs), collapse = "_")))
write.csv(out, fn, row.names = FALSE)
print(t(out))
