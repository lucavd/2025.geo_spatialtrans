# results/S1.4c/review/recalc/r06_c0.R — revisione avversariale S1.4c (recalc): C-0 indipendente.
# K-c = derive_cell_geometry() del pacchetto (R/04b4, invariato) con random_seed = seed_fold + 100000*s, s = 0..10, sui
# fold indicati; metriche (ks_nc, ks_an, d_sp, trunc_all, trunc_ev, loss_all, b7) ricalcolate con codice proprio
# (ks.test di R, cor di Spearman) e confrontate con le righe V-c di S1.4b loo. Uso: Rscript --vanilla r06_c0.R A:k ...
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf) })
E <- new.env(); sys.source("R/04b4_derive_cell_geometry.R", envir = E)
B <- "/mnt/micron/geo_spatialtrans"; ARCH <- c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")
ksD <- function(a, b) unname(suppressWarnings(ks.test(a, b, exact = FALSE)$statistic))
out <- list()
for (arg in commandArgs(TRUE)) {
  A <- strsplit(arg, ":")[[1]][1]; k <- as.integer(strsplit(arg, ":")[[1]][2]); met <- ARCH[[A]]
  o <- readRDS(file.path(B, "S1.4/loo", sprintf("%s_%s_fold%d.rds", A, met, k)))
  ob <- readRDS(file.path(B, "S1.4b/loo", sprintf("%s_%s_fold%d.rds", A, met, k)))
  te <- readRDS(file.path(B, "S1.4/real", sprintf("%s_%s_%s.rds", A, o$roi_test, met))); cc <- te$cells; n <- nrow(cc)
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  f2 <- o$fits$M2
  ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]], nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]])
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c")
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t))
  for (s in 0:10) {
    vc <- E$derive_cell_geometry(cen, tess, ct, random_seed = o$seed + 100000 * s, verbose = FALSE)$cell_df
    a <- vc$nucleus_area; tg <- vc$nucleus_area_target
    mine <- c(ks_nc = ksD(a[ev] / cc$a_t[ev], cc$a_n[ev] / cc$a_t[ev]), ks_an = ksD(a[ev], cc$a_n[ev]),
              d_sp = abs(cor(a[ev], cc$a_t[ev], method = "spearman") - cor(cc$a_n[ev], cc$a_t[ev], method = "spearman")),
              trunc_all = mean(vc$nucleus_clamped), trunc_ev = mean(vc$nucleus_clamped[ev]), loss_all = mean(1 - a / tg),
              b7 = mean(sqrt(a[ev] / pi) > vc$center_dist[ev]))
    ref <- unlist(ob$m[ob$m$seed_idx == s & ob$m$variant == "V-c", names(mine)])
    out[[length(out) + 1]] <- data.frame(archetype = A, fold = k, roi = o$roi_test, seed_idx = s, t(mine), max_abs_diff = max(abs(mine - ref)),
                                         identical_all = identical(unname(mine), unname(ref)))
  }
  cat(A, k, "fatto\n")
}
res <- do.call(rbind, out); write.csv(res, "/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/r06_c0.csv", row.names = FALSE)
cat("C-0 righe", nrow(res), "| max |mio - S1.4b V-c|", format(max(res$max_abs_diff), digits = 3), "| identici bit a bit", sum(res$identical_all), "\n")
print(aggregate(max_abs_diff ~ archetype + fold, res, max))
