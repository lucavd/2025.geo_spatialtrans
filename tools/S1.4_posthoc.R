# tools/S1.4_posthoc.R — S1.4: analisi POST HOC (dopo i verdetti; non cambiano i verdetti pre-registrati).
# PH-1 il nucleo reale (cerchio di pari area) entrerebbe nel territorio? PH-2 recupero senza clamp (mutante M1).
# PH-3 B-1/B-2/B-3 di M2 senza clamp (M1) e con clamp_factor = 0.99. Uso: Rscript --vanilla tools/S1.4_posthoc.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
source("tools/S1.4_mutants.R"); source("tools/R3_voronoi_metrics.R")
E0 <- dcg_env(); E1 <- dcg_mutate(dcg_env(), "M1")
OUT <- "/mnt/micron/geo_spatialtrans/S1.4"; RES <- "results/S1.4"; R3 <- "/mnt/micron/geo_spatialtrans/R3"
ARCH <- data.frame(archetype = paste0("A", 1:6), primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b))); max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
spear <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
load_real <- function(A, roi, method) readRDS(file.path(OUT, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
evset <- function(cc) cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
k32 <- sqrt(2 * pi / (32 * sin(2 * pi / 32)))
cls <- function(v, med, mx) if (v <= med) "PASS" else if (v <= mx) "WARN" else "FAIL"
arch_v <- function(o) if (sum(o == "PASS") >= 4) "PASS" else if (sum(o == "FAIL") >= 2) "FAIL" else "WARN"
# PH-1
ph1 <- do.call(rbind, lapply(seq_len(nrow(ARCH)), function(i) { A <- ARCH$archetype[i]; p <- ARCH$primary[i]
  do.call(rbind, lapply(rois$roi_id[rois$archetype == A], function(r) { x <- load_real(A, r, p); cc <- x$cells; e <- evset(cc)
    rr <- sqrt(cc$a_n[e] / pi)
    fo <- if (p == "cellpose_rgb") { nu <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_nuc.parquet", A, r, p)))); mean(nu$frac_out[e] > 0, na.rm = TRUE) } else NA
    data.frame(archetype = A, roi = r, n = sum(e), circ_poi_085 = mean(k32 * rr > 0.85 * cc$d_poi[e]), circ_gen_085 = mean(k32 * rr > 0.85 * cc$d_gen[e]),
               circ_gen_1 = mean(rr > cc$d_gen[e]), circ_poi_1 = mean(rr > cc$d_poi[e]), r3_frac_out_gt0 = fo,
               r_over_dgen_med = median(rr / cc$d_gen[e]), r_over_dpoi_med = median(rr / cc$d_poi[e])) })) }))
write.csv(ph1, file.path(RES, "S1.4_ph1.csv"), row.names = FALSE)
# PH-2 recupero senza clamp
REC <- expand.grid(beta = c(0, 0.5, 1), sigma = c(0.2, 0.4))
ph2 <- do.call(rbind, mclapply(seq_len(nrow(ARCH)), function(ai) { A <- ARCH$archetype[ai]; p <- ARCH$primary[ai]
  rs <- lapply(rois$roi_id[rois$archetype == A], function(r) load_real(A, r, p))
  a_ref <- sum(vapply(rs, function(r) sum(r$cells$a_t), 0)) / sum(vapply(rs, function(r) r$n, 0))
  do.call(rbind, lapply(seq_len(nrow(REC)), function(g) {
    ct <- data.frame(cell_type = "c", density = 1e6 / a_ref, nucleus_to_eq_ratio = 0.3, nucleus_area_slope = REC$beta[g], nucleus_area_sdlog = REC$sigma[g])
    sim <- do.call(rbind, lapply(seq_along(rs), function(i) { r <- rs[[i]]
      o <- E1$derive_cell_geometry(data.frame(cell_id = seq_len(r$n), region_id = 1L, x = r$cells$x, y = r$cells$y, cell_type = "c"),
             list(cell_territories = r$terr, territory_df = data.frame(cell_id = seq_len(r$n), territory_area = r$cells$a_t)), ct,
             random_seed = 20261015 + 50000 + 1000 * ai + 10 * g + i, verbose = FALSE, strict = FALSE)
      data.frame(a_n = o$cell_df$nucleus_area, a_t = r$cells$a_t, interior = r$cells$interior) }))
    s <- sim[sim$interior, ]; ft <- E0$fit_nucleus_size_model(s$a_n, s$a_t, a_ref, "full")
    data.frame(archetype = A, beta = REC$beta[g], sigma = REC$sigma[g], d_beta = abs(ft[["beta"]] - REC$beta[g]),
               d_sig = abs(ft[["sigma"]] - REC$sigma[g]) / REC$sigma[g], d_rho = abs(ft[["rho"]] - 0.3) / 0.3) })) }, mc.cores = 6))
ph2$ok <- ph2$d_beta <= 0.05 & ph2$d_sig <= 0.05 & ph2$d_rho <= 0.02
write.csv(ph2, file.path(RES, "S1.4_ph2.csv"), row.names = FALSE)
# PH-3 B-1..B-3 senza clamp e con clamp_factor 0.99 (parametri e semi dei fold)
jobs <- expand.grid(i = seq_len(nrow(ARCH)), k = 1:5)
ph3 <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(q) { i <- jobs$i[q]; k <- jobs$k[q]; A <- ARCH$archetype[i]; p <- ARCH$primary[i]
  o <- readRDS(file.path(OUT, "loo", sprintf("%s_%s_fold%d.rds", A, p, k))); te <- load_real(A, o$roi_test, p); cc <- te$cells; e <- evset(cc)
  f <- o$fits$M2; ct <- data.frame(cell_type = "c", density = 1e6 / f[["a_ref"]], nucleus_to_eq_ratio = f[["rho"]], nucleus_area_slope = f[["beta"]], nucleus_area_sdlog = f[["sigma"]])
  cen <- data.frame(cell_id = seq_len(te$n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c")
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(te$n), territory_area = cc$a_t))
  rn <- cc$a_n[e] / cc$a_t[e]; ra <- cc$a_n[e]; rs <- spear(ra, cc$a_t[e])
  one <- function(lbl, oo) { a <- oo$cell_df$nucleus_area[e]
    data.frame(archetype = A, fold = k, variant = lbl, ks_nc = ks_d(a / cc$a_t[e], rn), ks_an = ks_d(a, ra), d_sp = abs(spear(a, cc$a_t[e]) - rs),
               clamp = mean(oo$cell_df$nucleus_clamped), not_within = oo$info$n_not_within,
               y_nc_med = median(o$yard$ks_nc), y_nc_max = max(o$yard$ks_nc), y_an_med = median(o$yard$ks_an), y_an_max = max(o$yard$ks_an),
               y_sp_med = median(o$yard$d_sp), y_sp_max = max(o$yard$d_sp)) }
  rbind(one("M2 clamp 0.85 (primaria)", E0$derive_cell_geometry(cen, tess, ct, random_seed = o$seed, verbose = FALSE)),
        one("M2 clamp 0.99", E0$derive_cell_geometry(cen, tess, ct, random_seed = o$seed, verbose = FALSE, clamp_factor = 0.99)),
        one("M2 senza clamp (M1)", E1$derive_cell_geometry(cen, tess, ct, random_seed = o$seed, verbose = FALSE, strict = FALSE))) }, mc.cores = 30))
ph3$v_b1 <- mapply(cls, ph3$ks_nc, ph3$y_nc_med, ph3$y_nc_max); ph3$v_b2 <- mapply(cls, ph3$ks_an, ph3$y_an_med, ph3$y_an_max); ph3$v_b3 <- mapply(cls, ph3$d_sp, ph3$y_sp_med, ph3$y_sp_max)
write.csv(ph3, file.path(RES, "S1.4_ph3.csv"), row.names = FALSE)
cat("PH-1\n"); print(aggregate(cbind(circ_poi_085, circ_gen_085, circ_gen_1, r3_frac_out_gt0, r_over_dgen_med) ~ archetype, ph1, mean, na.action = na.pass), digits = 3)
cat("PH-2 recupero senza clamp:", sum(ph2$ok), "/", nrow(ph2), " max dB", round(max(ph2$d_beta), 4), "\n")
cat("PH-3\n"); s3 <- do.call(rbind, lapply(split(ph3, list(ph3$archetype, ph3$variant)), function(z) data.frame(archetype = z$archetype[1], variant = z$variant[1],
  B1 = arch_v(z$v_b1), B2 = arch_v(z$v_b2), B3 = arch_v(z$v_b3), clamp = round(mean(z$clamp), 3), fuori = sum(z$not_within) / sum(!is.na(z$not_within)))))
print(s3[order(s3$archetype, s3$variant), ], row.names = FALSE)
