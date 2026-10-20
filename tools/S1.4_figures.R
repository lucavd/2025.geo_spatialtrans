# tools/S1.4_figures.R — S1.4: figure del report (un pannello per archetipo, sfondo bianco) + ricalcolo di controllo.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.4_figures.R [outdir]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(ggplot2) })
source("R/04b4_derive_cell_geometry.R")
args <- commandArgs(TRUE); OUT <- if (length(args)) args[1] else "/mnt/micron/geo_spatialtrans/S1.4"
RES <- "results/S1.4"; FIG <- file.path(RES, "fig"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ARCH <- data.frame(archetype = paste0("A", 1:6), primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"),
                   label = c("A1 epitelio", "A2 tumore", "A3 stroma", "A4 linfoide", "A5 cervello", "A6 cuore"), stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id")]
lab <- setNames(ARCH$label, ARCH$archetype)
th <- theme_bw(base_size = 10) + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white", colour = NA),
                                       strip.background = element_rect(fill = "grey95"))
ks_d <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b))); max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
load_real <- function(A, roi, method) readRDS(file.path(OUT, "real", sprintf("%s_%s_%s.rds", A, roi, method)))
evset <- function(cc) cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
cal <- read.csv(file.path(RES, "S1.4_calibration.csv")); folds <- read.csv(file.path(RES, "S1.4_loo_folds.csv"))
set.seed(1)
# F1 dispersione log A_T vs log A_N con le rette dei modelli (5 ROI)
d1 <- do.call(rbind, lapply(seq_len(nrow(ARCH)), function(i) { A <- ARCH$archetype[i]; p <- ARCH$primary[i]
  cc <- do.call(rbind, lapply(rois$roi_id[rois$archetype == A], function(r) { x <- load_real(A, r, p)$cells; x[evset(x), ] }))
  cc <- cc[sample(nrow(cc), min(nrow(cc), 4000)), ]; data.frame(archetype = A, at = cc$a_t, an = cc$a_n) }))
ln <- do.call(rbind, lapply(seq_len(nrow(ARCH)), function(i) { A <- ARCH$archetype[i]; ci <- cal[cal$archetype == A & cal$method == ARCH$primary[i], ]
  xs <- range(log(d1$at[d1$archetype == A])); x <- seq(xs[1], xs[2], length.out = 50)
  cc <- d1[d1$archetype == A, ]; g0 <- mean(log(cc$an) - log(cc$at))
  rbind(data.frame(archetype = A, x = x, y = log(ci$rho^2 * ci$a_ref) + ci$beta * (x - log(ci$a_ref)), model = sprintf("M2 (beta = %.2f)", ci$beta)),
        data.frame(archetype = A, x = x, y = g0 + x, model = "M0 design (beta = 1)")) }))
ln$model <- ifelse(grepl("^M2", ln$model), "M2 stimato", "M0 design (beta = 1)")
p1 <- ggplot(d1, aes(log(at), log(an))) + geom_point(size = 0.3, alpha = 0.25, colour = "grey30") +
  geom_line(data = ln, aes(x, y, colour = model), linewidth = 0.8) + facet_wrap(~archetype, scales = "free", labeller = as_labeller(lab)) +
  scale_colour_manual(values = c("M2 stimato" = "#0072B2", "M0 design (beta = 1)" = "#D55E00")) +
  labs(x = "log area del territorio (µm²)", y = "log area del nucleo reale (µm²)", colour = NULL,
       title = "Taglia del nucleo contro taglia del territorio (celle interne, 5 ROI)") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4_fig1_scatter.png"), p1, width = 9, height = 6.5, dpi = 150)
# F2 ECDF di N/C sul fold 1 (ricalcolo con parametri e seme salvati) + controllo delle KS
chk <- list(); d2 <- list(); dm <- list()
for (i in seq_len(nrow(ARCH))) { A <- ARCH$archetype[i]; p <- ARCH$primary[i]
  o <- readRDS(file.path(OUT, "loo", sprintf("%s_%s_fold1.rds", A, p))); te <- load_real(A, o$roi_test, p); ev <- evset(te$cells)
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(te$n), territory_area = te$cells$a_t))
  cen <- data.frame(cell_id = seq_len(te$n), region_id = 1L, x = te$cells$x, y = te$cells$y, cell_type = "c")
  d2[[A]] <- data.frame(archetype = A, nc = te$cells$a_n[ev] / te$cells$a_t[ev], source = "reale")
  for (m in c("M0", "M2")) { f <- o$fits[[m]]
    ct <- data.frame(cell_type = "c", density = 1e6 / f[["a_ref"]], nucleus_to_eq_ratio = f[["rho"]], nucleus_area_slope = f[["beta"]], nucleus_area_sdlog = f[["sigma"]])
    s <- derive_cell_geometry(nucleus_model = "scale", cen, tess, ct, random_seed = o$seed, verbose = FALSE)
    nc <- s$cell_df$nucleus_area[ev] / te$cells$a_t[ev]
    chk[[length(chk) + 1]] <- data.frame(archetype = A, model = m, ks_recomputed = ks_d(nc, d2[[A]]$nc), ks_stored = o$sims[[m]][["ks_nc"]])
    d2[[paste(A, m)]] <- data.frame(archetype = A, nc = nc, source = if (m == "M0") "M0 design" else "M2 stimato")
    if (m == "M2") dm[[A]] <- list(te = te, s = s) } }
chk <- do.call(rbind, chk); write.csv(chk, file.path(RES, "S1.4_fig_recheck.csv"), row.names = FALSE)
cat("ricalcolo KS fold 1: max |diff| =", max(abs(chk$ks_recomputed - chk$ks_stored)), "\n")
d2 <- do.call(rbind, d2)
p2 <- ggplot(d2, aes(nc, colour = source)) + stat_ecdf(linewidth = 0.7) + facet_wrap(~archetype, scales = "free_x", labeller = as_labeller(lab)) +
  scale_x_log10() + scale_colour_manual(values = c(reale = "black", "M0 design" = "#D55E00", "M2 stimato" = "#0072B2")) +
  labs(x = "N/C = area nucleo / area territorio (scala log)", y = "ECDF", colour = NULL, title = "N/C sul ROI escluso (fold 1): reale contro simulato") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4_fig2_nc_ecdf.png"), p2, width = 9, height = 6.5, dpi = 150)
# F3 KS per fold e modello contro il metro inter-ROI
pr <- folds[folds$method == ARCH$primary[match(folds$archetype, ARCH$archetype)], ]
d3 <- rbind(data.frame(pr[, c("archetype", "fold", "model")], metric = "B-1 KS N/C", ks = pr$ks_nc, med = pr$y_ks_nc_med, mx = pr$y_ks_nc_max),
            data.frame(pr[, c("archetype", "fold", "model")], metric = "B-2 KS A_N", ks = pr$ks_an, med = pr$y_ks_an_med, mx = pr$y_ks_an_max))
ym <- unique(d3[d3$model == "M2", c("archetype", "fold", "metric", "med", "mx")])
p3 <- ggplot(d3, aes(factor(fold), ks, colour = model, shape = model)) +
  geom_linerange(data = ym, aes(x = factor(fold), ymin = 0, ymax = mx), inherit.aes = FALSE, colour = "grey80", linewidth = 3) +
  geom_point(data = ym, aes(x = factor(fold), y = med), inherit.aes = FALSE, shape = 95, size = 6, colour = "grey40") +
  geom_point(position = position_dodge(width = 0.6), size = 1.6) + facet_grid(metric ~ archetype, labeller = labeller(archetype = as_labeller(lab))) +
  scale_colour_manual(values = c(M0 = "#D55E00", M1 = "#E69F00", M2 = "#0072B2", M3 = "#009E73")) +
  labs(x = "fold (ROI escluso)", y = "distanza KS dal ROI escluso", colour = "modello", shape = "modello",
       title = "Distanza simulato-reale per fold; barra grigia = da 0 al massimo inter-ROI, trattino = mediana inter-ROI") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4_fig3_ks.png"), p3, width = 11, height = 5.5, dpi = 150)
# F4 mappe: finestra di 50 µm al centro del ROI del fold 1, territori, nucleo simulato M2, cerchio di area reale al generatore
d4p <- list(); d4n <- list(); d4r <- list()
for (A in ARCH$archetype) { te <- dm[[A]]$te; s <- dm[[A]]$s; cc <- te$cells
  side <- 60; x0 <- median(cc$x) - side / 2; y0 <- median(cc$y) - side / 2
  box <- st_sfc(st_polygon(list(rbind(c(0, 0), c(side, 0), c(side, side), c(0, side), c(0, 0)))))
  sh <- function(g) { g <- g - c(x0, y0); g }
  k <- which(cc$x >= x0 & cc$x <= x0 + side & cc$y >= y0 & cc$y <= y0 + side)
  clipb <- function(g, id) { z <- st_intersection(st_sf(id = id, geometry = sh(g)), st_sf(geometry = box)); z$archetype <- A; z[, c("id", "archetype")] }
  d4p[[A]] <- clipb(te$terr[k], k); d4n[[A]] <- clipb(s$cell_nuclei[k], k)
  rr <- k[cc$matched[k] & is.finite(cc$a_n[k])]
  d4r[[A]] <- clipb(st_buffer(st_sfc(lapply(rr, function(j) st_point(c(cc$x[j], cc$y[j])))), sqrt(cc$a_n[rr] / pi), nQuadSegs = 8), rr) }
d4p <- do.call(rbind, d4p); d4n <- do.call(rbind, d4n); d4r <- do.call(rbind, d4r)
p4 <- suppressWarnings(ggplot() + geom_sf(data = d4p, fill = "white", colour = "grey55", linewidth = 0.2) + geom_sf(data = d4n, fill = "#56B4E9", colour = "#0072B2", linewidth = 0.2, alpha = 0.7) +
  geom_sf(data = d4r, fill = NA, colour = "#D55E00", linewidth = 0.3) + facet_wrap(~archetype, labeller = as_labeller(lab)) + coord_sf(datum = NA) +
  labs(title = "Finestra 60 x 60 µm al centro del ROI del fold 1 (cellule con generatore nella finestra)", subtitle = "territori (grigio), nucleo simulato M2 (blu), cerchio di pari area del nucleo reale al generatore (arancio)") + th +
  theme(axis.text = element_blank(), axis.ticks = element_blank()))
ggsave(file.path(FIG, "S1.4_fig4_maps.png"), p4, width = 10, height = 7, dpi = 150)
# F5 frazione clampata per variante
b5 <- read.csv(file.path(RES, "S1.4_b5.csv")); b5 <- b5[b5$method == ARCH$primary[match(b5$archetype, ARCH$archetype)], ]
m2 <- pr[pr$model == "M2", ]
d5 <- rbind(data.frame(archetype = m2$archetype, variant = "cerchio, POI (primaria)", f = m2$clamp_all),
            data.frame(archetype = b5$archetype, variant = "cerchio, generatore", f = b5$clamp_gen),
            data.frame(archetype = b5$archetype, variant = "ellisse, orient. territorio", f = b5$clamp_ell_terr),
            data.frame(archetype = b5$archetype, variant = "ellisse, orient. casuale", f = b5$clamp_ell_rand))
p5 <- ggplot(d5, aes(variant, f, colour = variant)) + geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.4) + stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "black", linewidth = 0.3) +
  facet_wrap(~archetype, labeller = as_labeller(lab)) + scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "cellule con nucleo clampato", colour = NULL, title = "Frazione clampata per fold (punti) e media (barra); tratteggio = soglia B-4 del 5 %") + th +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
ggsave(file.path(FIG, "S1.4_fig5_clamp.png"), p5, width = 9, height = 6.5, dpi = 150)
cat("figure scritte in", FIG, "\n")
