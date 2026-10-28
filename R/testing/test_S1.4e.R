# R/testing/test_S1.4e.R - S1.4e: test del ridisegno dei nuclei, delle ellissi orientate e dei centri del banco A (cancello G-1; mutanti G-2).
# Pre-registrazione: docs/steps/S1.4e.md. Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S14E_RSRC=<cartella con R/ e tools/>]
#   Rscript --vanilla R/testing/test_S1.4e.R      I mutanti sono copie modificate del codice (run_S1.4e.sh mutants).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow) })
SRC <- Sys.getenv("S14E_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/R3_voronoi_metrics.R"); source("tools/S1.5c_bench_lib.R"); source("tools/S1.5d_bench_lib.R"); source(file.path(SRC, "tools/S1.4e_bench_lib.R"))
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
wkb <- function(n) lapply(as.list(geos::geos_write_wkb(geos::as_geos_geometry(n))), as.vector)
axis_of <- function(g) {                 # angolo dell'asse maggiore (rad, in (-pi/2, pi/2]) dai momenti del poligono (anello esterno)
  xy <- sf::st_coordinates(g)[, 1:2]; m <- r3_poly_moments(xy[, 1], xy[, 2]); 0.5 * atan2(2 * m[["mxy"]], m[["mxx"]] - m[["myy"]]) }
angdiff <- function(a, b) { d <- (a - b) %% pi; pmin(d, pi - d) }
rect <- function(cx, cy, w, h, al) { p <- rbind(c(-w, -h), c(w, -h), c(w, h), c(-w, h), c(-w, -h)) / 2
  R <- matrix(c(cos(al), sin(al), -sin(al), cos(al)), 2); q <- p %*% t(R); sf::st_polygon(list(cbind(q[, 1] + cx, q[, 2] + cy))) }
ct1 <- function(area, asp, orient, phi = 1) data.frame(cell_type = "c", density = 1e3, nucleus_to_eq_ratio = NA_real_, nucleus_radius = sqrt(area / pi),
  nucleus_area_sdlog = 0, nucleus_aspect = asp, nucleus_orient = orient, nucleus_nc_max = phi, stringsAsFactors = FALSE)
tess_of <- function(polys) { g <- sf::st_sfc(polys); list(cell_territories = g, territory_df = data.frame(cell_id = seq_along(g), territory_area = as.numeric(sf::st_area(g)))) }
cen_of <- function(x, y) data.frame(cell_id = seq_along(x), region_id = 1L, x = x, y = y, cell_type = "c")

# ---- T1: ridisegno con il catalogo del livello = nuclei di place_cells() (3 configurazioni) ----
sq <- expand.grid(x = 1:220, y = 1:220); sq$intensity_cluster <- 1L
pc <- function(preset, ...) suppressWarnings(place_cells(sq, preset, pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
  small_components = "exclude", random_seed = 11, types = "multi", size_arm = "group", keep_tiles = TRUE, strict = FALSE, verbose = FALSE, ...))
cfg <- list(
  "intestino, semina uniforme, residui gaussiani" = function() pc("epithelial", spatial = "uniform", size_resid = "gaussian"),
  "tumore del colon, compartimenti e vuoti RV, residui gaussiani" = function() pc("tumor_microenv", spatial = "compartments", min_dist_mode = "nucleus",
      override_compartments = compart_for_d("A2", "r1"), override_heterogeneity = het_for("A2", "r1", "RV"), quotas = TRUE, size_resid = "gaussian", resid_mode = "decile"),
  "cervello, semina uniforme, residui empirici" = function() pc("brain_cortex", spatial = "uniform", size_resid = "empirical"))
for (nm in names(cfg)) { L <- cfg[[nm]](); o <- redraw_nuclei(L, L$cell_types_resolved)
  ck(paste("T1 ridisegno = place_cells():", nm), identical(wkb(o$cell_nuclei), wkb(L$cell_nuclei)) &&
       isTRUE(all.equal(o$cell_df$nucleus_area, L$cell_df$nucleus_area, tolerance = 0)), sprintf("%d nuclei", nrow(L$cell_df))) }
# ---- T2: "territory" orienta l'asse maggiore lungo il territorio (4 rotazioni, entro 1 grado) ----
al <- c(0, 30, 75, 120) * pi / 180; xs <- seq(0, by = 100, length.out = 4)
o2 <- derive_cell_geometry(cen_of(xs, rep(50, 4)), tess_of(lapply(1:4, function(k) rect(xs[k], 50, 30, 12, al[k]))), ct1(20, 1.6, "territory"), random_seed = 1, verbose = FALSE)
d2 <- vapply(1:4, function(k) angdiff(axis_of(o2$cell_nuclei[[k]]), al[k]), 0) * 180 / pi
ck("T2 asse del nucleo lungo l'asse del territorio (entro 1 grado)", all(d2 < 1), sprintf("scarti %s gradi", paste(round(d2, 3), collapse = " ")))
hr <- sf::st_difference(rect(0, 0, 30, 12, 0.4), sf::st_buffer(sf::st_point(c(10 * cos(0.4), 10 * sin(0.4))), 1.5))   # territorio con un buco lontano dal centro
o2b <- derive_cell_geometry(cen_of(0, 0), tess_of(list(hr)), ct1(20, 1.6, "territory"), random_seed = 1, verbose = FALSE)
d2b <- angdiff(axis_of(sf::st_sfc(sf::st_polygon(list(sf::st_coordinates(o2b$cell_nuclei[[1]])[, 1:2])))[[1]]), 0.4) * 180 / pi
ck("T2b territorio con un buco: asse del nucleo lungo l'asse del territorio (entro 1 grado)", sf::st_geometry_type(hr) == "POLYGON" && length(hr) == 2 && d2b < 1, sprintf("scarto %.3f gradi", d2b))
# ---- T3: "random" da' orientamenti uniformi (2 000 nuclei, Kolmogorov-Smirnov p > 0.01) ----
gx <- rep(0:49, 40) * 25; gy <- rep(0:39, each = 50) * 25
o3 <- derive_cell_geometry(cen_of(gx, gy), tess_of(lapply(seq_along(gx), function(k) rect(gx[k], gy[k], 20, 20, 0))), ct1(20, 1.6, "random"), random_seed = 3, verbose = FALSE)
th3 <- vapply(seq_along(gx), function(k) axis_of(o3$cell_nuclei[[k]]), 0); p3 <- suppressWarnings(stats::ks.test((th3 + pi / 2) / pi, "punif")$p.value)
ck("T3 orientamenti uniformi con random (KS)", p3 > 0.01, sprintf("p = %.3f, n = %d", p3, length(th3)))
# ---- T4: ellissi nel modello clip su territori non convessi e con buco: dentro, area = A*, valide ----
Lsh <- sf::st_polygon(list(rbind(c(0, 0), c(14, 0), c(14, 5), c(5, 5), c(5, 14), c(0, 14), c(0, 0))))
hol <- sf::st_polygon(list(rbind(c(30, 0), c(44, 0), c(44, 12), c(30, 12), c(30, 0)), rbind(c(36, 5), c(38, 5), c(38, 7), c(36, 7), c(36, 5))))
nar <- rect(70, 6, 24, 5, 0.3)
t4 <- tess_of(list(Lsh, hol, nar)); c4 <- cen_of(c(2.5, 33, 70), c(2.5, 6, 6))
for (orient in c("random", "territory")) for (sd in 1:5) {
  o4 <- derive_cell_geometry(c4, t4, ct1(60, 1.7, orient, phi = 0.8), random_seed = sd, verbose = FALSE)
  inv <- nuc_invariants(o4, t4$cell_territories)
  ck(sprintf("T4 ellisse %s seme %d: dentro, nessuna sovrapposizione, area = A*, valida", orient, sd),
     inv$n_out == 0 && inv$n_overlap == 0 && inv$n_area_bad == 0 && inv$n_invalid == 0 && inv$n_empty == 0 && any(o4$cell_df$nucleus_clipped),
     sprintf("max errore relativo %.2e, ritagliati %d/3", inv$max_area_rel, sum(o4$cell_df$nucleus_clipped)))
}
# ---- T7: cataloghi dei bracci ----
ct0 <- data.frame(cell_type = c("a", "b"), nucleus_aspect = 1, nucleus_orient = "random", stringsAsFactors = FALSE)
k2 <- arm_catalog(ct0, "K2", "A1"); er <- arm_catalog(ct0, "ER", "A4"); et <- arm_catalog(ct0, "ET", "A1")
ck("T7 cataloghi dei bracci (K2 invariato; ER, ET allungamento 1/sqrt(1 - e^2) di B-115...B-120)",
   identical(k2, ct0) && all(abs(er$nucleus_aspect - 1.25) < 1e-12) && all(er$nucleus_orient == "random") &&
     all(abs(et$nucleus_aspect - 1 / sqrt(1 - 0.64)) < 1e-12) && all(et$nucleus_orient == "territory"),
   sprintf("ET A1 %.4f, ER A4 %.4f", et$nucleus_aspect[1], er$nucleus_aspect[1]))
# ---- T8: centri del banco A (finestra del tumore del colon e del linfonodo) ----
ann <- jsonlite::fromJSON("results/R2b/annotations/R2b_annotations_Luca.json", simplifyVector = FALSE)$windows
WIN <- read.csv("results/R2b/R2b_windows.csv", stringsAsFactors = FALSE); pts <- read.csv("results/R2b/R2b_points.csv", stringsAsFactors = FALSE)
rois <- read.csv("results/R2/R2_rois_checked.csv")[, c("archetype", "roi_id", "um_per_px", "side_px")]; rois$side_um <- rois$side_px * rois$um_per_px
for (wid in c("A2_w1", "A4_w2", "A6_w4")) {
  wv <- WIN[WIN$win_id == wid, ]; A <- wv$archetype; w <- r3_roi_window(A, wv$roi_id, rois)
  base <- if (A %in% c("A4", "A6")) "/mnt/micron/geo_spatialtrans/R3b/S1.4" else "/mnt/micron/geo_spatialtrans/S1.4"
  te <- readRDS(file.path(base, "real", sprintf("%s_%s_%s.rds", A, wv$roi_id, if (A %in% c("A4", "A6")) "consensus" else "cellpose_rgb")))
  pl <- pts[pts$rater == "Luca" & pts$win_id == wid, ]; cc <- centers_A(wv, te$cells[, c("x", "y")], pl, w, ann)
  s <- wv$side_um; inw <- cc$x >= wv$x0_um & cc$x < wv$x0_um + s & cc$y >= wv$y0_um & cc$y < wv$y0_um + s
  ex <- s14e_excl_um(wv, ann); co_in <- which(cc$src == "CO" & inw)
  te_in <- te$cells$x >= wv$x0_um & te$cells$x < wv$x0_um + s & te$cells$y >= wv$y0_um & te$cells$y < wv$y0_um + s
  n_te_ex <- if (is.null(ex) || !any(te_in)) 0L else sum(lengths(sf::st_intersects(sf::st_as_sf(te$cells[te_in, c("x", "y")], coords = c("x", "y")), ex)) > 0)
  ok_ex <- if (!length(co_in)) TRUE else !is.null(ex) && all(lengths(sf::st_intersects(sf::st_as_sf(cc[co_in, ], coords = c("x", "y")), ex)) > 0)
  ck(sprintf("T8 banco A %s: tutti i punti di Luca nella maschera, nessun centro del segmentatore nella finestra fuori dalle esclusioni", wid),
     sum(cc$src == "manual") == sum(r3_in_mask(pl$x_um, pl$y_um, w$m, w$px)) && sum(cc$src == "manual") > 0 && ok_ex && sum(cc$src == "CO") > 0 &&
       length(co_in) == n_te_ex && (wid != "A6_w4" || n_te_ex > 0),
     sprintf("manuali %d/%d, segmentatore nella finestra (nelle esclusioni) %d, attesi %d", sum(cc$src == "manual"), nrow(pl), length(co_in), n_te_ex))
}
cat(sprintf("RESULT %d/%d PASS\n", sum(unlist(res)), length(res)))
