# results/S1.3/review/code/export_inputs.R — revisione del codice S1.3 (RA-code)
# Esporta ingressi e uscite di tessellate_voronoi() per la reimplementazione indipendente in Python.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.3/review/code/export_inputs.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
source("tools/S1.3_roi.R")        # stadio "none": definisce rois, ARCH, one_region, gen_points, r3_roi_window, tessellate_voronoi
D <- "/mnt/micron/geo_spatialtrans/S1.3/review_code_data"; dir.create(D, showWarnings = FALSE, recursive = TRUE)
f17 <- function(v) sprintf("%.17g", v)
wkb <- function(s) vapply(sf::st_as_binary(s, hex = TRUE), as.character, "")
dump_case <- function(lab, cen, reg, tv) {
  write.csv(data.frame(cell_id = cen$cell_id, region_id = cen$region_id, x = f17(cen$x), y = f17(cen$y)),
            file.path(D, paste0(lab, "_cen.csv")), row.names = FALSE, quote = FALSE)
  rw <- t(vapply(seq_len(nrow(reg$region_df)), function(r) .tv_window(reg$region_polygons[r]), numeric(4)))
  write.csv(data.frame(region_id = reg$region_df$region_id, area_sf = f17(as.numeric(st_area(reg$region_polygons))),
                       xmin = f17(rw[, 1]), xmax = f17(rw[, 2]), ymin = f17(rw[, 3]), ymax = f17(rw[, 4]),
                       wkb = wkb(reg$region_polygons)), file.path(D, paste0(lab, "_reg.csv")), row.names = FALSE, quote = FALSE)
  td <- tv$territory_df
  write.csv(data.frame(cell_id = td$cell_id, region_id = td$region_id, territory_area = f17(td$territory_area),
                       tile_area = f17(td$tile_area), clipped = td$clipped, n_pieces_lost = td$n_pieces_lost,
                       n_pieces_gained = td$n_pieces_gained, gtype = as.character(st_geometry_type(tv$cell_territories)),
                       terr_wkb = wkb(tv$cell_territories), tile_wkb = wkb(tv$tiles)),
            file.path(D, paste0(lab, "_tv.csv")), row.names = FALSE, quote = FALSE)
  inf <- tv$info; inf$elapsed_s <- NULL
  write.csv(as.data.frame(lapply(inf, function(z) if (length(z) == 1) z else paste(z, collapse = ";"))),
            file.path(D, paste0(lab, "_info.csv")), row.names = FALSE)
  cat(lab, nrow(cen), "celle;", inf$n_fragments, "frammenti;", inf$n_fragments_reassigned, "riassegnati;", inf$n_snapped, "agganci\n")
}
# (a) sintetici I1-I5 seed 1, esattamente come load_case() + seed_case() di R/testing/test_S1.3.R
src <- readLines("R/testing/test_S1.3.R")
eval(parse(text = src[grep("^IN  <- ", src)]))
eval(parse(text = src[grep("^ct_test <- ", src):(grep("^# ---- verifiche su un caso", src) - 1)]))
syn <- c("I1_syn600_c1", "I2_syn600_c2", "I3_syn600_c3", "I4_syn600_c4", "I5_syn600_c2_labels")
for (id in syn) {
  reg <- load_case(id); cen <- seed_case(reg, 1, adv = FALSE)
  tv <- tessellate_voronoi(cen, reg, keep_tiles = TRUE, verbose = FALSE)
  dump_case(sub("_syn600.*", "", id), cen, reg, tv)
}
# (b) 6 ROI reali: centroidi R3 *_gen.parquet, regione = w$poly (unione della maschera) come one_region()
rr <- data.frame(A = c("A1", "A2", "A3", "A4", "A5", "A6"), roi = c("r1", "r1", "r1", "f1", "r1", "r1"),
                 m = c("cellpose_rgb", "cellpose_rgb", "cellpose_rgb", "spaceranger", "cellpose_rgb", "spaceranger"))
for (k in seq_len(nrow(rr))) {
  w <- r3_roi_window(rr$A[k], rr$roi[k], rois)
  gen <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", rr$A[k], rr$roi[k], rr$m[k]))))
  cen <- data.frame(cell_id = seq_len(nrow(gen)), region_id = 1L, x = gen$x, y = gen$y)
  tv <- tessellate_voronoi(cen, one_region(w), keep_tiles = TRUE, verbose = FALSE)
  lab <- sprintf("%s_%s", rr$A[k], rr$roi[k]); dump_case(lab, cen, one_region(w), tv)
  write.csv(data.frame(side_um = f17(w$side_um), px = f17(w$px), area_um2 = f17(w$area_um2), n_valid_px = w$n_valid_px),
            file.path(D, paste0(lab, "_win.csv")), row.names = FALSE, quote = FALSE)
}
# (c) caso C-2: A3 r3 CSR replica 3 (celle 766 e 1987 di S1.3_c2_overlap_pairs_check.csv)
i <- which(rois$archetype == "A3" & rois$roi_id == "r3"); a <- ARCH[ARCH$archetype == "A3", ]
w <- r3_roi_window("A3", "r3", rois); n_obs <- readRDS(file.path(R3, "real", sprintf("A3_r3_%s_roi.rds", a$primary)))$n
p <- gen_points(i, "CSR", 3, w, n_obs, a)
cen <- data.frame(cell_id = seq_along(p$x), region_id = 1L, x = p$x, y = p$y)
tv <- tessellate_voronoi(cen, one_region(w), keep_tiles = TRUE, verbose = FALSE)
dump_case("A3_r3_CSR03", cen, one_region(w), tv)
t <- tv$cell_territories
for (q in c(766L, 1987L)) cat("cella", q, "tipo", as.character(st_geometry_type(t[q])), "valida", st_is_valid(t[q]), "area", f17(st_area(t[q])),
                              "n_lost", tv$territory_df$n_pieces_lost[q], "n_gained", tv$territory_df$n_pieces_gained[q], "clipped", tv$territory_df$clipped[q], "\n")
g <- st_intersection(t[766], t[1987]); cat("area st_intersection(766,1987) =", f17(sum(as.numeric(st_area(g)))), "\n")
cat("st_relate(766,1987) =", st_relate(t[766], t[1987]), "\n")
cat("sf", as.character(packageVersion("sf")), "GEOS", sf::sf_extSoftVersion()[["GEOS"]], "\n")
