# tools/S1.4b_figures.R — S1.4b: figure del report (sfondo bianco, un pannello per archetipo). Uso: Rscript --vanilla tools/S1.4b_figures.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(ggplot2); library(sf); library(geos) })
source("R/04b1_extract_regions.R"); source("tools/R3_voronoi_metrics.R"); source("tools/S1.4_mutants.R"); source("tools/S1.4b_laguerre.R")
RES <- "results/S1.4b"; FIG <- file.path(RES, "fig"); dir.create(FIG, showWarnings = FALSE)
OUT <- "/mnt/micron/geo_spatialtrans/S1.4b"; S14 <- "/mnt/micron/geo_spatialtrans/S1.4"
th <- theme_bw(base_size = 11) + theme(plot.background = element_rect(fill = "white", colour = NA), strip.background = element_rect(fill = "grey92"))
col_v <- c("V-c" = "#B2182B", "V-a/V-b" = "#2166AC")
ft <- read.csv(file.path(RES, "S1.4b_folds.csv")); s14 <- lapply(list.files(file.path(S14, "loo"), full.names = TRUE), readRDS)
# yard median per fold
ym <- do.call(rbind, lapply(s14, function(o) data.frame(archetype = o$archetype, fold = o$fold, method = o$method,
  y_nc = median(o$yard$ks_nc), y_an = median(o$yard$ks_an), y_sp = median(o$yard$d_sp))))
ym <- ym[ym$method == c(A1 = "cellpose_rgb", A2 = "cellpose_rgb", A3 = "cellpose_rgb", A4 = "spaceranger", A5 = "cellpose_rgb", A6 = "spaceranger")[ym$archetype], ]
d <- merge(ft, ym, by = c("archetype", "fold"))
long <- rbind(data.frame(d[, c("archetype", "fold", "variant")], check = "B-1 N/C", r = d$ks_nc / d$y_nc),
              data.frame(d[, c("archetype", "fold", "variant")], check = "B-2 A_N", r = d$ks_an / d$y_an),
              data.frame(d[, c("archetype", "fold", "variant")], check = "B-3 Spearman", r = d$d_sp / d$y_sp),
              data.frame(d[, c("archetype", "fold", "variant")], check = "B-7 fuori dal Voronoi", r = d$b7 / pmax(d$y_b7_med, 1e-3)))
p1 <- ggplot(long, aes(check, r, colour = variant)) + geom_hline(yintercept = 1, linetype = 2) +
  geom_point(position = position_dodge(width = 0.5), size = 1.8, alpha = 0.8) + scale_y_log10() + facet_wrap(~archetype, nrow = 2) +
  scale_colour_manual(values = col_v, name = NULL) + labs(x = NULL, y = "distanza simulato–reale / mediana del metro inter-ROI (log)",
  title = "Check B per fold (mediana di 10 semi): sotto la linea = PASS") + th + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(FIG, "S1.4b_fig1_checkB.png"), p1, width = 10, height = 6.5, dpi = 150, bg = "white")
# fig 2: troncamento
ar <- read.csv(file.path(RES, "S1.4b_archetypes.csv"))
p2 <- ggplot(ft, aes(archetype, trunc_all, colour = variant)) + geom_hline(yintercept = 0.05, linetype = 2) +
  geom_point(position = position_dodge(width = 0.5), size = 2) + scale_colour_manual(values = col_v, name = NULL) +
  labs(x = NULL, y = "frazione di nuclei troncati (per fold)", title = "B-4: troncamento dell'area estratta (soglia 5 %)") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4b_fig2_trunc.png"), p2, width = 7, height = 4.5, dpi = 150, bg = "white")
# fig 3: CP-1
cp <- read.csv(file.path(RES, "S1.4b_cp1_raw.csv")); cp <- cp[cp$cond != "check_R3", ]
cp$kind <- ifelse(grepl("^perm", cp$cond), "pesi permutati (20)", ifelse(cp$cond == "vor", "Voronoi", ifelse(cp$cond == "lag", "Laguerre, pesi reali", paste0("Laguerre λ = ", cp$lam))))
cp$kind <- factor(cp$kind, levels = c("Voronoi", "Laguerre λ = 0.5", "Laguerre, pesi reali", "Laguerre λ = 1.5", "Laguerre λ = 2", "pesi permutati (20)"))
p3 <- ggplot(cp, aes(roi, m1, colour = kind)) + geom_point(data = cp[grepl("^perm", cp$cond), ], colour = "grey65", size = 1, position = position_jitter(width = 0.15, height = 0)) +
  geom_point(data = cp[!grepl("^perm", cp$cond), ], size = 2.2, position = position_dodge(width = 0.5)) + facet_wrap(~archetype, scales = "free", nrow = 2) +
  scale_colour_manual(values = c("Voronoi" = "#B2182B", "Laguerre λ = 0.5" = "#92C5DE", "Laguerre, pesi reali" = "#2166AC", "Laguerre λ = 1.5" = "#4393C3", "Laguerre λ = 2" = "#053061", "pesi permutati (20)" = "grey65"), name = NULL) +
  labs(x = "ROI", y = "frazione dei pixel nucleari fuori dalla propria cella (m1)", title = "CP-1 sulle maschere reali: Voronoi, Laguerre con pesi reali, placebo (grigio)") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4b_fig3_cp1.png"), p3, width = 10, height = 6.5, dpi = 150, bg = "white")
# fig 4: mappe (A2, A4; fold 1, seme s = 1): V-c, V-b, V-a in una finestra di 50 µm
E <- dcg_env(); rois <- read.csv("results/R2/R2_rois_checked.csv"); rois$side_um <- rois$side_px * rois$um_per_px
panels <- list()
for (A in c("A2", "A4")) { method <- if (A == "A4") "spaceranger" else "cellpose_rgb"
  o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold1.rds", A, method))); te <- readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, o$roi_test, method)))
  cc <- te$cells; n <- te$n; f2 <- o$fits$M2
  ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]], nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]])
  vc <- E$derive_cell_geometry(nucleus_model = "scale", data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c"),
          list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t)), ct, random_seed = o$seed + 100000, verbose = FALSE)
  w <- r3_roi_window(A, o$roi_test, rois); reg <- geos_set_precision(as_geos_geometry(w$poly), 1e-6)
  vb <- vab_nuclei(cc$x, cc$y, vc$cell_df$nucleus_area_target, geos_boundary(reg)); pd <- power_diagram(cc$x, cc$y, vb$df$Rc^2, reg)
  cx <- median(cc$x[cc$interior]); cy <- median(cc$y[cc$interior]); box <- st_sfc(st_polygon(list(rbind(c(cx, cy), c(cx + 50, cy), c(cx + 50, cy + 50), c(cx, cy + 50), c(cx, cy)))))
  tosf <- function(g) st_as_sfc(geos_write_wkb(g))
  real_c <- tosf(nucleus_polys(cc$x, cc$y, k32 * sqrt(ifelse(is.finite(cc$a_n), cc$a_n, cc$a_own) / pi)))
  sets <- list(`V-c: Voronoi + POI + clamp` = list(t = te$terr, nu = vc$cell_nuclei), `V-b: Voronoi + nuclei al generatore` = list(t = te$terr, nu = tosf(nucleus_polys(cc$x, cc$y, vb$df$Rc))),
               `V-a: Laguerre + stessi nuclei` = list(t = tosf(pd$cells), nu = tosf(nucleus_polys(cc$x, cc$y, vb$df$Rc))), `reale: Voronoi + cerchi di pari area` = list(t = te$terr, nu = real_c))
  for (nm in names(sets)) { s <- sets[[nm]]
    ti <- suppressWarnings(st_intersection(st_sf(geometry = st_make_valid(s$t)), box)); ni <- suppressWarnings(st_intersection(st_sf(geometry = st_make_valid(s$nu)), box))
    sh <- function(g) st_sfc(lapply(g, function(z) z - c(cx, cy)))   # coordinate locali della finestra (0-50 µm)
    panels[[length(panels) + 1]] <- rbind(data.frame(archetype = A, set = nm, layer = "territorio", geometry = sh(st_geometry(ti))), data.frame(archetype = A, set = nm, layer = "nucleo", geometry = sh(st_geometry(ni)))) } }
pm <- st_sf(do.call(rbind, panels)); pm$set <- factor(pm$set, levels = unique(pm$set))
p4 <- ggplot() + geom_sf(data = pm[pm$layer == "territorio", ], fill = "grey97", colour = "grey30", linewidth = 0.25) +
  geom_sf(data = pm[pm$layer == "nucleo", ], fill = "#6A51A3", colour = NA, alpha = 0.75) + facet_grid(archetype ~ set) +
  coord_sf(xlim = c(0, 50), ylim = c(0, 50), expand = FALSE) + labs(title = "Finestre di 50 µm (fold 1, seme s = 1): territori e nuclei per opzione", x = "µm", y = "µm") + th + theme(strip.text = element_text(size = 8))
ggsave(file.path(FIG, "S1.4b_fig4_maps.png"), p4, width = 12, height = 6.8, dpi = 150, bg = "white")
# fig 5: B-9 CV e |dA_T|
b9 <- read.csv(file.path(RES, "S1.4b_b8_b9.csv"))
l9 <- rbind(data.frame(archetype = b9$archetype, k = "Laguerre simulato", v = b9$dat_med), data.frame(archetype = b9$archetype, k = "Laguerre reale (aree reali)", v = b9$dat_real_med))
p5 <- ggplot(l9, aes(archetype, v, colour = k)) + geom_hline(yintercept = 0.10, linetype = 2) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  scale_colour_manual(values = c("Laguerre simulato" = "#2166AC", "Laguerre reale (aree reali)" = "black"), name = NULL) +
  labs(x = NULL, y = "mediana di |A_Laguerre − A_Voronoi| / A_Voronoi (per fold)", title = "B-9: quanto Laguerre cambia i territori (P5: ≤ 0.10)") + th + theme(legend.position = "bottom")
ggsave(file.path(FIG, "S1.4b_fig5_b9.png"), p5, width = 7, height = 4.5, dpi = 150, bg = "white")
cat("figure ok\n")
