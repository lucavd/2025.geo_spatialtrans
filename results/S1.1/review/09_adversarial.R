# Revisione avversariale S1.1 — 09_adversarial.R: input avversari NUOVI contro extract_regions()
.libPaths("/home/user/2025.geo_spatialtrans/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages(library(sf))
setwd("/home/user/2025.geo_spatialtrans"); source("R/04b1_extract_regions.R")
REV <- "/mnt/micron/geo_spatialtrans/S1.1/review"
R <- list()
rec <- function(test, expected, observed, ok) {
  R[[length(R) + 1L]] <<- data.frame(test = test, expected = expected, observed = observed, outcome = if (isTRUE(ok)) "atteso" else "NON atteso")
  cat(sprintf("[%s] %s | atteso: %s | osservato: %s\n", if (isTRUE(ok)) "ok" else "XX", test, expected, observed))
}
ER <- function(cl, px, ...) extract_regions(cl, px, cluster_col = "k", verbose = FALSE, ...)
fm <- function(m) { i <- which(!is.na(m) & m != "" , arr.ind = TRUE); data.frame(x = i[, 2], y = i[, 1], k = m[i], stringsAsFactors = FALSE) }

# A1 stride: due bin isolati (1,1) e (4,1) a 8 µm -> 2 componenti da 64 µm², entrambe escluse
cl <- data.frame(x = c(1L, 4L), y = c(1L, 1L), k = 1L)
r <- ER(cl, 8, min_region_area_um2 = 0, simplify_tol_um = 0)
rec("A1 stride: bin (1,1),(4,1), px 8", "stride 1, 2 comp, 64+64 µm²",
    sprintf("stride %d, %d comp, aree %s", r$info$stride, r$info$n_components, paste(r$region_df$area_um2, collapse = "+")),
    r$info$stride == 1 && r$info$n_components == 2)
rf <- ER(cl, 8)
rec("A1b stesso input, default (min 100)", "0 regioni (tutto escluso)", sprintf("%d regioni, area %s", rf$info$n_regions, paste(rf$region_df$area_um2, collapse = "+")), rf$info$n_regions == 0)
r1 <- ER(cl, 8, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1)
rec("A1c stesso input, stride = 1 forzato", "2 comp da 64", sprintf("%d comp, aree %s", r1$info$n_components, paste(r1$region_df$area_um2, collapse = "+")), r1$info$n_components == 2)
# A2 una riga con bin alterni (1,3,5) -> 3 componenti
cl <- data.frame(x = c(1L, 3L, 5L), y = 1L, k = 1L)
r <- ER(cl, 1, min_region_area_um2 = 0, simplify_tol_um = 0)
rec("A2 stride: riga x = 1,3,5, px 1", "3 comp da 1 µm²", sprintf("stride %d, %d comp, area tot %g", r$info$stride, r$info$n_components, sum(r$region_df$area_um2)), r$info$n_components == 3)
# A3 due bin lontani (1,1),(5001,5001): stride inferito; e costo con stride forzato
cl <- data.frame(x = c(1L, 5001L), y = c(1L, 5001L), k = 1L)
r <- ER(cl, 8, simplify_tol_um = 0)
rec("A3 bin (1,1),(5001,5001), px 8, default", "0 regioni (2 bin da 64 µm²)",
    sprintf("stride %d, %d regioni, area tot %.3g µm²", r$info$stride, r$info$n_regions, sum(r$region_df$area_um2)), r$info$n_regions == 0)
g0 <- gc(reset = TRUE); t0 <- proc.time()[[3]]
r <- ER(cl, 8, simplify_tol_um = 0, stride = 1)
el <- proc.time()[[3]] - t0; g1 <- gc()
rec("A3b stesso, stride = 1 forzato (griglia 5001x5001 per 2 pixel)", "corretto; costo proporzionale al bounding box",
    sprintf("%d comp, %.1f s, picco R %.0f MB", r$info$n_components, el, sum(g1[, ncol(g1)])), r$info$n_components == 2)
# A4 etichette carattere "1" e "01" (distinte) affiancate
m <- matrix(c("1", "1", "01", "01"), nrow = 2, ncol = 4, byrow = FALSE); m <- matrix(rep(c("1", "1", "01", "01"), each = 2), 2, 4)
r <- ER(fm(m), 10, min_region_area_um2 = 0, simplify_tol_um = 0)
rec("A4 etichette '1' e '01' (carattere) affiancate", "2 cluster, 2 regioni",
    sprintf("%d cluster (%s), %d regioni", length(r$info$cluster_ids), paste(r$info$cluster_ids, collapse = ","), r$info$n_regions), r$info$n_regions == 2)
cl <- fm(m); cl$k <- factor(ifelse(cl$k == "1", "1", "1.0"))
r <- ER(cl, 10, min_region_area_um2 = 0, simplify_tol_um = 0)
rec("A4b fattore con livelli '1' e '1.0'", "2 cluster, 2 regioni",
    sprintf("%d cluster, %d regioni", length(r$info$cluster_ids), r$info$n_regions), r$info$n_regions == 2)
# A5 cluster 0 e negativi
m <- matrix(c(0L, 0L, -1L, -1L, 2L, 2L), nrow = 2); cl <- data.frame(x = rep(1:3, each = 2), y = rep(1:2, 3), k = as.vector(m))
r <- ER(cl, 5, min_region_area_um2 = 0, simplify_tol_um = 0)
rec("A5 cluster 0, -1, 2 (interi)", "3 regioni da 50 µm², id -1,0,2",
    sprintf("%d regioni, id %s, aree %s", r$info$n_regions, paste(r$region_df$cluster_id, collapse = ","), paste(r$region_df$area_um2, collapse = ",")),
    r$info$n_regions == 3 && setequal(r$region_df$cluster_id, c(-1, 0, 2)) && all(r$region_df$area_um2 == 50))
# A6 coordinate negative
cl <- expand.grid(x = -5:4, y = -3:6); cl$k <- 1L
r <- ER(cl, 2, min_region_area_um2 = 0, simplify_tol_um = 0)
bb <- as.numeric(st_bbox(r$region_polygons))
rec("A6 coordinate negative x -5..4, y -3..6, px 2", "bbox -12,-8,8,12; area 400",
    sprintf("bbox %s; area %g", paste(bb, collapse = ","), r$region_df$area_um2), all(bb == c(-12, -8, 8, 12)) && r$region_df$area_um2 == 400)
# A7 stride diverso per asse
cl <- expand.grid(x = seq(1L, 19L, 2L), y = seq(1L, 28L, 3L)); cl$k <- 1L
e <- try(ER(cl, 1), silent = TRUE)
rec("A7 stride x = 2, y = 3", "errore esplicito", if (inherits(e, "try-error")) sub("\n", "", as.character(e)) else "nessun errore", inherits(e, "try-error"))
# A8 due buchi di fondo che si toccano in un vertice, dentro una regione 6x6
m <- matrix(1L, 6, 6); m[3, 3] <- 0L; m[4, 4] <- 0L
cl <- { i <- which(m != 0, arr.ind = TRUE); data.frame(x = i[, 2], y = i[, 1], k = m[i]) }
r0 <- ER(cl, 1, min_region_area_um2 = 0, simplify_tol_um = 0); rs <- ER(cl, 1, min_region_area_um2 = 0)
rec("A8 due buchi diagonali che si toccano", "1 regione, area 34, valida (esatta e semplificata)",
    sprintf("%d reg, area %g / %g, n_holes %d, valid %s/%s, tipo %s", r0$info$n_regions, r0$region_df$area_um2, rs$region_df$area_um2, r0$region_df$n_holes,
            st_is_valid(r0$region_polygons), st_is_valid(rs$region_polygons), as.character(st_geometry_type(r0$region_polygons))),
    r0$info$n_regions == 1 && r0$region_df$area_um2 == 34 && all(st_is_valid(r0$region_polygons)) && all(st_is_valid(rs$region_polygons)))
# A9 semplificazione con pixel piccolo (0.25 µm): disco digitale r = 60 px e linea larga 1 px
d <- expand.grid(x = 1:121, y = 1:121); d <- d[(d$x - 61)^2 + (d$y - 61)^2 <= 60^2, ]; d$k <- 1L
r0 <- ER(d, 0.25, min_region_area_um2 = 0, simplify_tol_um = 0); rs <- ER(d, 0.25)
sd <- as.numeric(st_area(st_sym_difference(r0$region_polygons[[1]], rs$region_polygons[[1]])))
rec("A9 disco r = 60 px a 0.25 µm/px, tol 0.5 µm", "errore area <= 2 % (C1s)",
    sprintf("area esatta %.2f, semplificata %.2f (%.2f %%), diff. simmetrica %.2f µm² (%.1f %%), vertici %d -> %d",
            r0$region_df$area_um2, rs$region_df$area_um2, 100 * (rs$region_df$area_um2 / r0$region_df$area_um2 - 1), sd, 100 * sd / r0$region_df$area_um2,
            nrow(st_coordinates(r0$region_polygons)), nrow(st_coordinates(rs$region_polygons))),
    abs(rs$region_df$area_um2 / r0$region_df$area_um2 - 1) <= 0.02)
ln <- data.frame(x = 1:800, y = 1L, k = 1L)
r0 <- ER(ln, 0.25, min_region_area_um2 = 0, simplify_tol_um = 0); rs <- ER(ln, 0.25, min_region_area_um2 = 0)
rec("A9b linea 800x1 px a 0.25 µm/px (larga 0.25 µm < tol)", "area 50 µm², geometria valida",
    sprintf("area %.3f -> %.3f, valid %s, vuota %s", r0$region_df$area_um2, rs$region_df$area_um2, st_is_valid(rs$region_polygons), st_is_empty(rs$region_polygons)),
    abs(rs$region_df$area_um2 - 50) < 1e-9 && st_is_valid(rs$region_polygons))
# A10 C1s cieco: I2 a 1 µm/px — area per regione invariata ma geometria spostata?
o <- readRDS("/mnt/micron/geo_spatialtrans/S1.1/inputs/I2_syn600_c2.rds")
r0 <- extract_regions(o$clust, 1, min_region_area_um2 = 0, cluster_col = o$cluster_col, simplify_tol_um = 0, verbose = FALSE)
rs <- extract_regions(o$clust, 1, min_region_area_um2 = 0, cluster_col = o$cluster_col, simplify_tol_um = 0.5, verbose = FALSE)
symd <- mapply(function(a, b) as.numeric(st_area(st_sym_difference(a, b))), r0$region_polygons, rs$region_polygons)
da <- abs(as.numeric(st_area(rs$region_polygons)) - as.numeric(st_area(r0$region_polygons)))
nv0 <- nrow(st_coordinates(r0$region_polygons)); nvs <- nrow(st_coordinates(rs$region_polygons))
rec("A10 I2 (1 µm/px) tol 0.5: |Δarea| vs diff. simmetrica per regione", "se Δarea = 0 anche diff. simmetrica = 0",
    sprintf("max |Δarea| %.3g µm², somma diff. simmetrica %.1f µm² (%.4f %% del tessuto), vertici %d -> %d",
            max(da), sum(symd), 100 * sum(symd) / sum(r0$region_df$area_px_um2), nv0, nvs), sum(symd) < 1e-9 || max(da) > 0)
# A11 n_holes / tipo per adv_pinch_hole semplificato (WARN C6)
o <- readRDS("/mnt/micron/geo_spatialtrans/S1.1/inputs/adv_pinch_hole.rds")
r0 <- extract_regions(o$clust, 1, min_region_area_um2 = 0, cluster_col = o$cluster_col, simplify_tol_um = 0, verbose = FALSE)
rs <- extract_regions(o$clust, 1, min_region_area_um2 = 0, cluster_col = o$cluster_col, verbose = FALSE)
rec("A11 adv_pinch_hole esatto vs semplificato", "stessa geometria",
    sprintf("WKT esatto %s | semplificato %s | diff. simm. %.3f", st_as_text(r0$region_polygons[[1]]), st_as_text(rs$region_polygons[[1]]),
            as.numeric(st_area(st_sym_difference(r0$region_polygons[[1]], rs$region_polygons[[1]])))),
    identical(st_as_text(r0$region_polygons[[1]]), st_as_text(rs$region_polygons[[1]])))
TT <- do.call(rbind, R); write.csv(TT, file.path(REV, "09_adversarial_results.csv"), row.names = FALSE)

