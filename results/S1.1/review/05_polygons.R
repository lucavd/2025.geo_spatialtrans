# Revisione avversariale S1.1 — 05_polygons.R
# Compito 5: chiama extract_regions() (unico uso del codice in esame) e confronta le aree dei poligoni con
# le aree attese calcolate INDIPENDENTEMENTE in Python (02_components.py: n_px per componente da scipy).
# Abbinamento poligono <-> componente: il centro del primo pixel (ordine di riga) della componente scipy.
.libPaths("/home/user/2025.geo_spatialtrans/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(arrow) })
setwd("/home/user/2025.geo_spatialtrans"); source("R/04b1_extract_regions.R")
IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"; REV <- "/mnt/micron/geo_spatialtrans/S1.1/review"
stride_of <- c(I6_syn6800x6500_c2 = 6)
res <- list()
for (id in commandArgs(trailingOnly = TRUE)) {
  o <- readRDS(file.path(IN, paste0(id, ".rds"))); px <- o$pixel_size_um
  s <- if (id %in% names(stride_of)) stride_of[[id]] else 1
  mine <- as.data.frame(read_parquet(file.path(REV, "py_out", paste0(id, "_components.parquet"))))
  mine$area_exp <- mine$n_px * (s * px)^2
  t0 <- proc.time()[[3]]
  r0 <- extract_regions(o$clust, px, min_region_area_um2 = 0, cluster_col = o$cluster_col, simplify_tol_um = 0, verbose = FALSE)
  rf <- extract_regions(o$clust, px, cluster_col = o$cluster_col, verbose = FALSE)        # uso: min 100, tol 0.5
  el <- proc.time()[[3]] - t0
  a0 <- as.numeric(st_area(r0$region_polygons)); af <- as.numeric(st_area(rf$region_polygons))
  # abbinamento su TUTTE le componenti
  pts <- st_sfc(lapply(seq_len(nrow(mine)), function(k)
           st_point(c((mine$rep_x[k] - 1 + s / 2) * px, (mine$rep_y[k] - 1 + s / 2) * px))))
  hit <- st_intersects(pts, r0$region_polygons); nh <- lengths(hit)
  one <- nh == 1L; pi <- rep(NA_integer_, length(nh)); pi[one] <- unlist(hit[one])
  area_ok <- one & abs(a0[pi] - mine$area_exp) <= 1e-9 * mine$area_exp
  clus_ok <- one & as.character(r0$region_df$cluster_id[pi]) == as.character(mine$cluster)
  bij <- one & !duplicated(pi)
  # rf: componenti >= 100 µm² attese
  keep <- mine$area_exp >= 100
  hitf <- st_intersects(pts[keep], rf$region_polygons); nhf <- lengths(hitf)
  pf <- rep(NA_integer_, length(nhf)); pf[nhf == 1L] <- unlist(hitf[nhf == 1L])
  relf <- abs(af[pf] - mine$area_exp[keep]) / mine$area_exp[keep]
  v0 <- st_is_valid(r0$region_polygons); vf <- st_is_valid(rf$region_polygons)
  set.seed(7); samp <- sample(which(one), min(250L, sum(one)))
  write.csv(data.frame(id = id, comp = samp, cluster = mine$cluster[samp], n_px_scipy = mine$n_px[samp],
                       area_expected = mine$area_exp[samp], polygon = pi[samp], area_polygon_r0 = a0[pi[samp]],
                       valid_r0 = v0[pi[samp]], n_holes = r0$region_df$n_holes[pi[samp]]),
            file.path(REV, paste0("05_sample_", id, ".csv")), row.names = FALSE)
  res[[id]] <- data.frame(id = id, stride_fun = r0$info$stride, n_comp_scipy = nrow(mine), n_poly_r0 = length(r0$region_polygons),
    sum_area_expected = sum(mine$area_exp), sum_area_r0 = sum(a0), relerr_r0 = (sum(a0) - sum(mine$area_exp)) / sum(mine$area_exp),
    sum_area_expected_100 = sum(mine$area_exp[keep]), n_poly_rf = length(rf$region_polygons), n_expected_100 = sum(keep),
    sum_area_rf = sum(af), relerr_rf = (sum(af) - sum(mine$area_exp[keep])) / sum(mine$area_exp[keep]),
    max_rel_region_rf = max(relf, na.rm = TRUE), n_rf_unmatched = sum(nhf != 1L),
    n_pts_one_poly = sum(one), n_pts_multi = sum(nh > 1L), n_pts_none = sum(nh == 0L), n_bijective = sum(bij),
    n_area_ok = sum(area_ok), n_cluster_ok = sum(clus_ok), max_abs_area_diff = max(abs(a0[pi] - mine$area_exp), na.rm = TRUE),
    n_invalid_r0 = sum(!v0), n_invalid_rf = sum(!vf), n_nonpolygon_r0 = sum(st_geometry_type(r0$region_polygons) != "POLYGON"),
    sample_n = length(samp), sample_area_ok = sum(area_ok[samp]), sample_valid = sum(v0[pi[samp]]),
    elapsed_extract_s = el)
  print(res[[id]])
}
write.csv(do.call(rbind, res), file.path(REV, "05_polygons_summary.csv"), row.names = FALSE)

