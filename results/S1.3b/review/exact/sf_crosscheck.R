
# Controprove in R/sf (GEOS 3.12.1) sulle coordinate esportate in esadecimale (nessuna rigenerazione).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); suppressMessages(library(sf))
E <- "/mnt/micron/geo_spatialtrans/S1.3b/review_exact"
build <- function(d) { d <- d[order(d$part, d$ring, d$k), ]
  parts <- lapply(split(d, d$part), function(dp) lapply(split(dp, dp$ring), function(dr) cbind(dr$X, dr$Y)))
  if (d$type[1] == "MULTIPOLYGON") st_multipolygon(parts) else st_polygon(parts[[1]]) }
rd <- function(f) { d <- read.csv(f, colClasses = c(x = "character", y = "character", xd = "character", yd = "character"))
  d$X <- as.numeric(d$x); d$Y <- as.numeric(d$y)
  stopifnot(identical(d$X, as.numeric(d$xd)), identical(d$Y, as.numeric(d$yd))); d }
d <- rd(file.path(E, "P1_pairs_coords_hex.csv"))
out <- list()
for (cs in unique(d$case)) { gi <- st_sfc(build(d[d$case == cs & d$role == "i", ])); gj <- st_sfc(build(d[d$case == cs & d$role == "j", ]))
  out[[cs]] <- data.frame(case = cs, relate_ij = st_relate(gi, gj)[1, 1], relate_ji = st_relate(gj, gi)[1, 1],
    a_int_ij = sum(as.numeric(st_area(st_intersection(gi, gj)))), a_int_ji = sum(as.numeric(st_area(st_intersection(gj, gi)))),
    a_diff_ij = sum(as.numeric(st_area(st_difference(gi, gj)))), a_diff_ji = sum(as.numeric(st_area(st_difference(gj, gi)))),
    union_deficit = as.numeric(st_area(gi)) + as.numeric(st_area(gj)) - as.numeric(st_area(st_union(gi, gj))),
    a_int_grid1e6 = { p <- st_set_precision(c(gi, gj), 1e6); sum(as.numeric(st_area(st_intersection(p[1], p[2])))) },
    relate_grid1e6 = { p <- st_set_precision(c(gi, gj), 1e6); st_relate(p[1], p[2])[1, 1] }) }
out <- do.call(rbind, out); write.csv(out, file.path(E, "P1_sf_orders.csv"), row.names = FALSE)
print(out, digits = 6, row.names = FALSE)
# A5 r5 RSA 15: quali celle mancano dall'unione?
h <- rd(file.path(E, "A5_r5_RSA_15_coords_hex.csv")); reg <- st_sfc(build(h[h$role == "region", ]))
roles <- sort(unique(h$role[h$role != "region"])); tt <- st_sfc(lapply(roles, function(r) build(h[h$role == r, ])))
a <- as.numeric(st_area(tt)); u <- st_union(tt); cat("n", length(tt), "sum a", sprintf("%.17g", sum(a)), "union", sprintf("%.17g", as.numeric(st_area(u))), "\n")
pos <- st_point_on_surface(tt); inU <- lengths(st_within(pos, u)) > 0
cov <- as.numeric(st_area(st_intersection(tt, u)))   # area di ogni cella coperta dall'unione
miss <- data.frame(cell = which(!inU | abs(cov - a) > 1e-6), stringsAsFactors = FALSE)
miss$area <- a[miss$cell]; miss$pos_in_union <- inU[miss$cell]; miss$covered <- cov[miss$cell]
print(miss); write.csv(miss, file.path(E, "A5_union_missing_sf.csv"), row.names = FALSE)
cat("somma aree mancanti", sum(miss$area - miss$covered), "\n")
# unione cascata vs unione a coppie con la cella mancante
for (k in miss$cell) { nb <- st_intersects(tt[k], tt)[[1]]; nb <- setdiff(nb, k)
  cat("cella", k, "vicini", paste(nb, collapse = " "), "\n")
  for (q in nb) cat("  con", q, "relate", st_relate(tt[k], tt[q])[1, 1], "a_int", sum(as.numeric(st_area(st_intersection(tt[k], tt[q])))), "union_def", a[k] + a[q] - as.numeric(st_area(st_union(tt[k], tt[q]))), "\n") }
