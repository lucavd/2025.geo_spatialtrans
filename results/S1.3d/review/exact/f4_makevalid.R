# review/exact/f4_makevalid.R — punto (b): st_make_valid (sf/GEOS del progetto) sulle celle F4 invalide; scrive i vertici dell'uscita
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(arrow); library(sf) })
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; O <- file.path(D, "review_exact", "f4")
cat(sf_extSoftVersion()[c("GEOS")], "\n")
R <- read.csv(file.path(O, "f4_invalid.csv")); OUT <- list(); S <- list()
for (q in seq_len(nrow(R))) { key <- R$key[q]; ci <- R$cell[q]
  t <- read_parquet(file.path(D, "geom", paste0(key, "_f4cells.parquet"))); t <- t[t$cell == ci, ]; t <- t[order(t$part, t$ring, t$k), ]
  parts <- lapply(split(t, t$part), function(p) lapply(split(p, p$ring), function(r) { m <- as.matrix(r[, c("x", "y")]); rbind(m, m[1, ]) }))
  g <- if (length(parts) == 1) st_polygon(parts[[1]]) else st_multipolygon(parts)
  v0 <- st_is_valid(g); mv <- st_make_valid(g); v1 <- st_is_valid(mv)
  pol <- st_collection_extract(st_sfc(mv), "POLYGON"); a1 <- sum(st_area(pol)); a0 <- st_area(g)
  S[[q]] <- data.frame(key = key, cell = ci, valid_before = v0, reason_before = st_is_valid(g, reason = TRUE), valid_after = v1, type_after = as.character(st_geometry_type(mv)),
                       n_poly_after = length(st_cast(pol, "POLYGON")), n_rings_after = sum(sapply(st_cast(pol, "POLYGON"), length)), area_sf_before = a0, area_sf_after = a1)
  cc <- st_coordinates(pol); L <- if ("L3" %in% colnames(cc)) cc[, "L3"] else cc[, "L2"] * 0 + 1
  OUT[[q]] <- data.frame(key = key, cell = ci, part = L, ring = cc[, "L1"], k = ave(cc[, "X"], L, cc[, "L1"], FUN = seq_along), x = cc[, "X"], y = cc[, "Y"]) }
S <- do.call(rbind, S); write.csv(S, file.path(O, "f4_makevalid_summary.csv"), row.names = FALSE)
write_parquet(do.call(rbind, OUT), file.path(O, "f4_makevalid_coords.parquet")); print(S)
