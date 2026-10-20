
# review/code/vc_poi.R — revisore "code": (i) ricalcolo indipendente delle metriche V-c e V-a/V-b dai dati (30 fold x 10 semi);
# (ii) V-c con l'operatore di osservazione dichiarato in pre-registrazione: metriche rispetto al Voronoi dei CENTRI dei nuclei
# simulati (POI), non dei generatori; (iii) esporta x, y, a_tgt, d_reg per A2/A4 (s = 1) per l'ottimo di non sovrapposizione.
# Scrive solo in /tmp/code_vc. Uso dalla root del repo.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(parallel); library(arrow) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("R/04b3_tessellate_voronoi.R")
source("tools/R3_voronoi_metrics.R"); source("tools/S1.4_mutants.R")
E <- dcg_env(); L <- new.env(); sys.source("tools/S1.4b_laguerre.R", envir = L)
OD <- "/tmp/code_vc"; dir.create(OD, showWarnings = FALSE)
S14 <- "/mnt/micron/geo_spatialtrans/S1.4"; S14B <- "/mnt/micron/geo_spatialtrans/S1.4b"; GRID <- 1e-6
ARCH <- data.frame(archetype = paste0("A", 1:6), primary = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"))
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]
rois$side_um <- rois$side_px * rois$um_per_px; rois$roi_idx <- seq_len(nrow(rois))
one_region <- function(w) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = w$area_um2), region_polygons = st_sfc(w$poly[[1]]))
ksd <- function(a, b) { a <- sort(a); b <- sort(b); v <- sort(unique(c(a, b))); max(abs(findInterval(v, a) / length(a) - findInterval(v, b) / length(b))) }
spr <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))
one <- function(A, k) {
  method <- ARCH$primary[ARCH$archetype == A]
  f <- file.path(OD, sprintf("%s_fold%d.rds", A, k)); if (file.exists(f)) return("skip")
  o <- readRDS(file.path(S14, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  b <- readRDS(file.path(S14B, "loo", sprintf("%s_%s_fold%d.rds", A, method, k)))
  te <- readRDS(file.path(S14, "real", sprintf("%s_%s_%s.rds", A, o$roi_test, method))); cc <- te$cells; n <- te$n
  ev <- cc$interior & cc$matched & is.finite(cc$a_n) & cc$a_n > 0
  w <- r3_roi_window(A, o$roi_test, rois); reg <- geos_set_precision(as_geos_geometry(w$poly), GRID); bnd <- geos_boundary(reg)
  f2 <- o$fits$M2
  ct <- data.frame(cell_type = "c", density = 1e6 / f2[["a_ref"]], nucleus_to_eq_ratio = f2[["rho"]], nucleus_area_slope = f2[["beta"]], nucleus_area_sdlog = f2[["sigma"]])
  cen <- data.frame(cell_id = seq_len(n), region_id = 1L, x = cc$x, y = cc$y, cell_type = "c")
  tess <- list(cell_territories = te$terr, territory_df = data.frame(cell_id = seq_len(n), territory_area = cc$a_t))
  rnc <- cc$a_n[ev] / cc$a_t[ev]; ran <- cc$a_n[ev]; rsp <- spr(ran, cc$a_t[ev])
  rows <- list(); poi <- NULL; exp1 <- NULL
  for (s in 1:10) {
    vc <- E$derive_cell_geometry(nucleus_model = "scale", cen, tess, ct, random_seed = o$seed + 100000 * s, verbose = FALSE)
    d <- vc$cell_df; an <- d$nucleus_area
    if (is.null(poi)) {
      poi <- cbind(d$nucleus_cx, d$nucleus_cy)
      tvp <- tessellate_voronoi(data.frame(cell_id = seq_len(n), region_id = 1L, x = poi[, 1], y = poi[, 2]), one_region(w), verbose = FALSE)
      at_p <- tvp$territory_df$territory_area; gtp <- as_geos_geometry(tvp$cell_territories)
      d_p <- geos_distance(geos_make_point(poi[, 1], poi[, 2]), geos_boundary(gtp))
      d_p_gen <- geos_distance(geos_make_point(cc$x, cc$y), geos_boundary(gtp))
    } else stopifnot(identical(poi, cbind(d$nucleus_cx, d$nucleus_cy)))
    # V-c come nel banco (operatore = Voronoi dei generatori, distanza dal POI)
    vc_bench <- c(ks_nc = ksd(an[ev] / cc$a_t[ev], rnc), ks_an = ksd(an[ev], ran), d_sp = abs(spr(an[ev], cc$a_t[ev]) - rsp),
                  trunc_all = mean(d$nucleus_clamped), b7 = mean(sqrt(an[ev] / pi) > d$center_dist[ev]))
    # V-c con l'operatore della pre-registrazione: Voronoi dei centri dei nuclei (POI)
    vc_poi <- c(ks_nc = ksd(an[ev] / at_p[ev], rnc), ks_an = ksd(an[ev], ran), d_sp = abs(spr(an[ev], at_p[ev]) - rsp),
                trunc_all = mean(d$nucleus_clamped), b7 = mean(sqrt(an[ev] / pi) > d_p[ev]))
    vb <- L$vab_nuclei(cc$x, cc$y, d$nucleus_area_target, bnd, 0.99)$df
    vab <- c(ks_nc = ksd(vb$a_n[ev] / cc$a_t[ev], rnc), ks_an = ksd(vb$a_n[ev], ran), d_sp = abs(spr(vb$a_n[ev], cc$a_t[ev]) - rsp),
             trunc_all = mean(vb$s < 1), b7 = mean(sqrt(vb$a_n[ev] / pi) > cc$d_gen[ev]), loss_all = mean(1 - vb$s^2))
    sb <- b$m[b$m$seed_idx == s, ]
    rows[[s]] <- rbind(data.frame(seed_idx = s, variant = "V-c_bench", t(vc_bench), loss_all = mean(1 - an / d$nucleus_area_target)),
                       data.frame(seed_idx = s, variant = "V-c_poiVor", t(vc_poi), loss_all = NA),
                       data.frame(seed_idx = s, variant = "V-ab", t(vab)))
    rows[[s]]$stored_equal <- c(isTRUE(all.equal(unname(vc_bench), unname(unlist(sb[sb$variant == "V-c", c("ks_nc", "ks_an", "d_sp", "trunc_all", "b7")])), tolerance = 0)),
                                NA, isTRUE(all.equal(unname(vab[1:5]), unname(unlist(sb[sb$variant == "V-a", c("ks_nc", "ks_an", "d_sp", "trunc_all", "b7")])), tolerance = 0)))
    if (s == 1 && A %in% c("A2", "A4")) {
      d_reg <- geos_distance(geos_make_point(cc$x, cc$y), bnd)
      write_parquet(data.frame(x = cc$x, y = cc$y, a_tgt = d$nucleus_area_target, d_reg = d_reg, ev = ev, s_bench = vb$s, vc_clamped = d$nucleus_clamped,
                               a_t = cc$a_t, a_n_real = cc$a_n, matched = cc$matched), file.path(OD, sprintf("%s_fold%d_s1.parquet", A, k)))
    }
  }
  saveRDS(list(A = A, k = k, roi = o$roi_test, m = do.call(rbind, rows), yard = o$yard, yard_b7 = b$yard_b7, b7_real = b$b7_real,
               poi_shift = summary(sqrt((poi[, 1] - cc$x)^2 + (poi[, 2] - cc$y)^2)), at_p_vs_at = median(abs(at_p[ev] - cc$a_t[ev]) / cc$a_t[ev]),
               gen_in_poi_cell_d = summary(d_p_gen)), f)
  "ok"
}
jobs <- expand.grid(A = ARCH$archetype, k = 1:5, stringsAsFactors = FALSE)
r <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(one(jobs$A[i], jobs$k[i]), error = function(e) paste("ERR", conditionMessage(e))), mc.cores = 30, mc.preschedule = FALSE)
print(table(unlist(r))); print(unlist(r)[grepl("ERR", unlist(r))])
