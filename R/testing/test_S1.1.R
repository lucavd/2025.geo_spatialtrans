# S1.1 — check C, controprove e descrittivo B per extract_regions()
# Pre-registrazione: results/S1.1/S1.1_preregistration.md (approvata 2026-10-07)
# Uso (dalla root del repo, dopo tools/S1.1_inputs.R):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S1.1.R
# Stampa PASS/WARN/FAIL per asserzione; scrive results/S1.1/S1.1_test_results.csv,
# S1.1_per_input.csv, S1.1_checkB_descriptive.csv, S1.1_cp1_moore.csv, S1.1_cp3_null.csv,
# S1.1_cp4_truth.csv; output per input in /mnt/micron/geo_spatialtrans/S1.1/outputs/.
#
# Definizioni operative (fissate prima di eseguire):
#  r0 = extract_regions(min_region_area_um2 = 0, simplify_tol_um = 0)   -> tutte le componenti, esatte
#  rs = extract_regions(min_region_area_um2 = 0, simplify_tol_um = 0.5) -> tutte, semplificate
#  rf = extract_regions(default: min 100 µm², tol 0.5)                  -> output d'uso
#  C1  : |Σ area(r0) − n_px·e²| / (n_px·e²)
#  C1s : idem su rs (globale); per regione |area_s − area_px| / area_px, 95° percentile
#  C4  : griglia ricostruita nel test (indipendente da .er_build_grid) + imager::label(4-vicinato):
#        multinsieme (cluster, n_px) identico a quello di r0
#  C5  : centri di tutte le celle della griglia (fondo incluso) -> cluster del poligono che li
#        contiene (0 se nessuno, NA se > 1) == cluster vero; frazione di celle concordanti
#  C6  : sovrapposizione = Σ area − area(∪); scoperto = area(∪ r0) − area(∪ r0 ∩ ∪ rs);
#        esatto: entrambi / area tessuto < 1e-9; semplificato: (sovr. + scoperto) / area ≤ 1 %
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(imager); library(parallel) })
source("R/04b1_extract_regions.R")
source("tools/S1.1_moore_variant.R")

IN  <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
OUTD <- "/mnt/micron/geo_spatialtrans/S1.1/outputs"
RES <- "results/S1.1"
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
NCORES <- 12L
TOL_SIMPL <- 0.5; MIN_AREA <- 100

inputs <- read.csv(file.path(RES, "S1.1_inputs.csv"), stringsAsFactors = FALSE)
inputs <- inputs[inputs$group != "perf", ]          # C11 in tools/S1.1_perf.R
load_in <- function(id) readRDS(file.path(IN, paste0(id, ".rds")))
er <- function(o, min_area, tol, px = o$pixel_size_um, clust = o$clust)
  extract_regions(clust, px, min_region_area_um2 = min_area, cluster_col = o$cluster_col,
                  simplify_tol_um = tol, verbose = FALSE)

# griglia indipendente: indici di cella e codici di cluster (0 = fondo)
indep_grid <- function(o) {
  x <- round(o$clust$x); y <- round(o$clust$y)
  d <- c(diff(sort(unique(x))), diff(sort(unique(y))))
  s <- if (length(d)) min(d) else 1
  stopifnot(all(d %% s == 0))
  i <- (x - min(x)) / s + 1; j <- (y - min(y)) / s + 1
  cl <- as.character(o$clust[[o$cluster_col]])
  G <- matrix(0L, nrow = max(j), ncol = max(i))       # [riga = y, colonna = x]
  lv <- sort(unique(cl))
  G[cbind(j, i)] <- match(cl, lv)
  list(G = G, s = s, x0 = min(x), y0 = min(y), lv = lv)
}

imager_components <- function(g) {
  G <- rbind(0L, cbind(0L, g$G, 0L), 0L)          # bordo di fondo: as.cimg rifiuta le matrici 1x1
  L <- as.matrix(imager::label(as.cimg(G), high_connectivity = FALSE))
  fg <- G != 0L
  tab <- aggregate(n_px ~ lab + code, data = data.frame(lab = L[fg], code = G[fg], n_px = 1L), FUN = sum)
  data.frame(cluster = g$lv[tab$code], n_px = tab$n_px)
}

comp_key <- function(cluster, n_px) {
  k <- data.frame(cluster = as.character(cluster), n_px = as.integer(n_px))
  k[order(k$cluster, k$n_px), , drop = FALSE]
}

raster_roundtrip <- function(g, polys, cl_of_poly, px) {
  ny <- nrow(g$G); nx <- ncol(g$G)
  ii <- rep(seq_len(nx), each = ny); jj <- rep(seq_len(ny), times = nx)   # ordine column-major di G
  cx <- (g$x0 - 1 + (ii - 1) * g$s + g$s / 2) * px
  cy <- (g$y0 - 1 + (jj - 1) * g$s + g$s / 2) * px
  pts <- sf::st_as_sf(data.frame(cx = cx, cy = cy), coords = c("cx", "cy"))
  hit <- sf::st_intersects(pts, polys)
  nh <- lengths(hit)
  pred <- rep(0L, length(nh))
  one <- nh == 1L
  pred[one] <- match(as.character(cl_of_poly[unlist(hit[one])]), g$lv)
  pred[nh > 1L] <- NA_integer_
  truth <- as.vector(g$G)
  c(match = mean(!is.na(pred) & pred == truth), n_multi = sum(nh > 1L), n_cells = length(nh))
}

area_union <- function(p) if (length(p)) as.numeric(sf::st_area(sf::st_union(p))) else 0

run_one <- function(id) {
  o <- load_in(id); e <- NULL
  t0 <- proc.time()[["elapsed"]]
  r0 <- er(o, 0, 0); rs <- er(o, 0, TOL_SIMPL); rf <- er(o, MIN_AREA, TOL_SIMPL)
  e2 <- r0$info$effective_pixel_um^2
  tissue <- r0$info$n_pixels * e2
  g <- indep_grid(o)
  big <- o$group %in% c("real_full", "null_full")
  # C4
  ours <- comp_key(r0$region_df$cluster_id, r0$region_df$n_px)
  ref  <- if (o$group == "null_full") NULL else imager_components(g)
  c4 <- if (is.null(ref)) NA else isTRUE(all.equal(ours, comp_key(ref$cluster, ref$n_px), check.attributes = FALSE))
  # C5 (non sui nulli: milioni di micro-componenti, verificati da C1/C4 sulle mappe reali)
  rt0 <- rts <- c(match = NA, n_multi = NA, n_cells = NA)
  if (o$group != "null_full") {
    rt0 <- raster_roundtrip(g, r0$region_polygons, r0$region_df$cluster_id, o$pixel_size_um)
    rts <- raster_roundtrip(g, rs$region_polygons, rs$region_df$cluster_id, o$pixel_size_um)
  }
  # C6
  ov0 <- unc0 <- ovs <- uncs <- NA
  if (o$group != "null_full") {
    u0 <- sf::st_union(r0$region_polygons); us <- sf::st_union(rs$region_polygons)
    a_u0 <- as.numeric(sf::st_area(u0)); a_us <- as.numeric(sf::st_area(us))
    ov0  <- (sum(r0$region_df$area_um2) - a_u0) / tissue
    unc0 <- (tissue - a_u0) / tissue
    ovs  <- (sum(rs$region_df$area_um2) - a_us) / tissue
    uncs <- (a_u0 - as.numeric(sf::st_area(sf::st_intersection(u0, us)))) / tissue
  }
  relreg <- abs(rs$region_df$area_um2 - rs$region_df$area_px_um2) / rs$region_df$area_px_um2
  # C8 (riproducibilita' e invarianza all'ordine delle righe)
  c8 <- NA
  if (o$group %in% c("adv", "syn") || id == "real_full_A1") {
    rf2 <- er(o, MIN_AREA, TOL_SIMPL)
    set.seed(1); sh <- o$clust[sample(nrow(o$clust)), ]
    rf3 <- er(o, MIN_AREA, TOL_SIMPL, clust = sh)
    strip <- function(r) { r$info$elapsed_s <- NULL; r }
    c8 <- identical(strip(rf), strip(rf2)) && identical(strip(rf), strip(rf3))
  }
  # C7 scala (px = 2) su avversari e I1
  c7s <- NA
  if (o$group == "adv" || id == "I1_syn600_c1") {
    r2 <- er(o, 0, 0, px = 2 * o$pixel_size_um)
    c7s <- isTRUE(all.equal(r2$region_df$area_um2, 4 * r0$region_df$area_um2, tolerance = 1e-12)) &&
           isTRUE(all.equal(r2$region_df$perimeter_um, 2 * r0$region_df$perimeter_um, tolerance = 1e-12))
  }
  # C7 filtro
  c7f <- all(rf$region_df$area_px_um2 >= MIN_AREA) && all(rf$excluded_df$area_px_um2 < MIN_AREA) &&
         (rf$info$n_regions + rf$info$n_excluded == r0$info$n_components) &&
         isTRUE(all.equal(comp_key(rf$region_df$cluster_id, rf$region_df$n_px),
                          comp_key(r0$region_df$cluster_id[r0$region_df$area_px_um2 >= MIN_AREA],
                                   r0$region_df$n_px[r0$region_df$area_px_um2 >= MIN_AREA]),
                          check.attributes = FALSE))
  # C3
  k_big <- unique(r0$region_df$cluster_id[r0$region_df$area_px_um2 >= MIN_AREA])
  c3 <- rf$info$n_regions >= length(k_big) && all(k_big %in% rf$region_df$cluster_id)
  # CP-1 (Moore sui centri) su I1-I4 e blob
  moore <- NULL
  if (grepl("^I[1-4]_", id) || id == "adv_blobs_9_10") {
    G <- rbind(0L, cbind(0L, g$G, 0L), 0L)
    L <- as.matrix(imager::label(as.cimg(G), high_connectivity = FALSE)); fg <- G != 0L
    labs <- unique(L[fg])
    moore <- do.call(rbind, lapply(labs, function(l) {
      w <- which(L == l, arr.ind = TRUE)
      B <- matrix(FALSE, diff(range(w[, 1])) + 1, diff(range(w[, 2])) + 1)
      B[cbind(w[, 1] - min(w[, 1]) + 1, w[, 2] - min(w[, 2]) + 1)] <- TRUE
      data.frame(id = id, n_px = nrow(w), area_moore_px = moore_centre_area(B))
    }))
  }
  out <- list(
    per_input = data.frame(
      id = id, group = o$group, n_px = r0$info$n_pixels, stride = r0$info$stride,
      stride_indep = g$s, n_components = r0$info$n_components,
      n_components_imager = if (is.null(ref)) NA else nrow(ref), c4_identical = c4,
      c1_relerr = abs(sum(r0$region_df$area_um2) - tissue) / tissue,
      c1s_relerr = abs(sum(rs$region_df$area_um2) - tissue) / tissue,
      c1s_p95_region = unname(quantile(relreg, 0.95)), c1s_max_region = max(relreg),
      c2_valid_exact = all(sf::st_is_valid(r0$region_polygons)),
      c2_valid_simpl = all(sf::st_is_valid(rs$region_polygons)) && all(sf::st_is_valid(rf$region_polygons)),
      c2_polygon_only = all(sf::st_geometry_type(r0$region_polygons) == "POLYGON") &&
                        all(sf::st_geometry_type(rs$region_polygons) == "POLYGON"),
      c2_nonempty = !any(sf::st_is_empty(r0$region_polygons)) && !any(sf::st_is_empty(rs$region_polygons)),
      c3 = c3, c5_match_exact = rt0[["match"]], c5_match_simpl = rts[["match"]],
      c5_multi_simpl = rts[["n_multi"]], c6_overlap_exact = ov0, c6_uncovered_exact = unc0,
      c6_overlap_simpl = ovs, c6_uncovered_simpl = uncs, c7_filter = c7f, c7_scale = c7s, c8 = c8,
      n_regions_100 = rf$info$n_regions, n_excluded_100 = rf$info$n_excluded,
      frac_area_excluded_100 = rf$info$frac_area_excluded, n_clusters = length(r0$info$cluster_ids),
      elapsed_total_s = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE),
    regions = data.frame(id = rep(id, nrow(rf$region_df)), group = rep(o$group, nrow(rf$region_df)),
                         cluster_id = rf$region_df$cluster_id,
                         area_um2 = rf$region_df$area_um2, area_px_um2 = rf$region_df$area_px_um2),
    comps = data.frame(id = id, n_px = r0$region_df$n_px, area_px_um2 = r0$region_df$area_px_um2),
    moore = moore, meta = o$meta)
  saveRDS(list(rf = rf, r0_df = r0$region_df, meta = o$meta, group = o$group,
               pixel_size_um = o$pixel_size_um), file.path(OUTD, paste0(id, ".rds")))
  out
}

cat(sprintf("[test_S1.1] %d input, %d core\n", nrow(inputs), NCORES))
ord <- inputs$id[order(-inputs$n_px)]                 # grandi per primi
res <- mclapply(ord, function(id) tryCatch(run_one(id), error = function(e) list(error = conditionMessage(e), id = id)),
                mc.cores = NCORES, mc.preschedule = FALSE)
names(res) <- ord
errs <- Filter(function(z) !is.null(z$error), res)
if (length(errs)) { for (z in errs) cat("ERRORE", z$id, z$error, "\n"); stop("input falliti") }
per <- do.call(rbind, lapply(res, `[[`, "per_input")); per <- per[match(inputs$id, per$id), ]
write.csv(per, file.path(RES, "S1.1_per_input.csv"), row.names = FALSE)
regions <- do.call(rbind, lapply(res, `[[`, "regions"))
comps <- do.call(rbind, lapply(res, `[[`, "comps"))
saveRDS(list(regions = regions, comps = comps), file.path(OUTD, "S1.1_regions_all.rds"))

# ---- tabella delle asserzioni ---------------------------------------------------
T <- list()
add <- function(check, input, metric, value, threshold, status)
  T[[length(T) + 1L]] <<- data.frame(check = check, input = input, metric = metric,
                                     value = format(value, digits = 6), threshold = threshold, status = status)
pf <- function(ok) ifelse(isTRUE(ok), "PASS", "FAIL")
pw <- function(ok) ifelse(isTRUE(ok), "PASS", "WARN")
for (k in seq_len(nrow(per))) {
  p <- per[k, ]; id <- p$id
  add("C1", id, "errore relativo area (esatta)", p$c1_relerr, "< 1e-9", pf(p$c1_relerr < 1e-9))
  add("C1s", id, "errore relativo area globale (semplificata)", p$c1s_relerr, "<= 0.02", pf(p$c1s_relerr <= 0.02))
  add("C1s", id, "95° percentile errore per regione (semplificata)", p$c1s_p95_region, "<= 0.02 (WARN)", pw(p$c1s_p95_region <= 0.02))
  add("C2", id, "valide + non vuote + solo POLYGON", paste(p$c2_valid_exact, p$c2_valid_simpl, p$c2_nonempty, p$c2_polygon_only),
      "tutte TRUE", pf(p$c2_valid_exact && p$c2_valid_simpl && p$c2_nonempty && p$c2_polygon_only))
  add("C3", id, "n_regioni >= k e ogni cluster con componente >= soglia presente", p$c3, "TRUE", pf(p$c3))
  if (!is.na(p$c4_identical))
    add("C4", id, sprintf("multinsieme (cluster, n_px) vs imager::label (%d vs %d)", p$n_components, p$n_components_imager),
        p$c4_identical, "identico", pf(p$c4_identical))
  add("C4b", id, "stride rilevato = stride ricostruito nel test", paste(p$stride, p$stride_indep), "uguali", pf(p$stride == p$stride_indep))
  if (!is.na(p$c5_match_exact)) {
    add("C5", id, "concordanza raster (esatta)", p$c5_match_exact, "= 1", pf(p$c5_match_exact == 1))
    add("C5", id, "concordanza raster (semplificata)", p$c5_match_simpl, ">= 0.98", pf(p$c5_match_simpl >= 0.98))
    add("C6", id, "sovrapposizione + scoperto (esatta)", p$c6_overlap_exact + p$c6_uncovered_exact, "< 1e-9",
        pf(abs(p$c6_overlap_exact) < 1e-9 && abs(p$c6_uncovered_exact) < 1e-9))
    add("C6", id, "sovrapposizione + scoperto (semplificata)", p$c6_overlap_simpl + p$c6_uncovered_simpl, "<= 0.01 (WARN)",
        pw(p$c6_overlap_simpl + p$c6_uncovered_simpl <= 0.01))
  }
  add("C7", id, "filtro a 100 µm²", p$c7_filter, "esatto", pf(p$c7_filter))
  if (!is.na(p$c7_scale)) add("C7", id, "px = 2: aree x4, perimetri x2", p$c7_scale, "esatto", pf(p$c7_scale))
  if (!is.na(p$c8)) add("C8", id, "run ripetuto e righe permutate -> identico", p$c8, "identico", pf(p$c8))
}
# attese sugli avversari
gp <- function(id) per[per$id == id, ]
r0_of <- function(id) readRDS(file.path(OUTD, paste0(id, ".rds")))
adv_exp <- list(adv_checker = 2500, adv_corner = 2, adv_pinch_hole = 1, adv_line = 3, adv_border = 7)
for (id in names(adv_exp)) add("ADV", id, "n componenti", gp(id)$n_components, adv_exp[[id]], pf(gp(id)$n_components == adv_exp[[id]]))
d <- r0_of("adv_donut")$r0_df
add("ADV", "adv_donut", "area / buchi", paste(d$area_um2, d$n_holes), "84 / 1", pf(nrow(d) == 1 && d$area_um2 == 84 && d$n_holes == 1))
d <- r0_of("adv_donut_island")$r0_df
add("ADV", "adv_donut_island", "regioni (ciambella + isola)", nrow(d), "2 (3 con il buco come fondo)", pf(nrow(d) == 2 && all(sort(d$area_um2) == c(4, 84))))
d <- r0_of("adv_pinch_hole")$r0_df
add("ADV", "adv_pinch_hole", "area, buco che tocca il guscio in un vertice", d$area_um2, "16, valido", pf(nrow(d) == 1 && d$area_um2 == 16 && gp("adv_pinch_hole")$c2_valid_exact))
z <- r0_of("adv_single_cluster")
add("ADV/CP-2", "adv_single_cluster", "regioni / area", paste(z$rf$info$n_regions, z$rf$region_df$area_um2), "1 / 8000",
    pf(z$rf$info$n_regions == 1 && z$rf$region_df$area_um2 == 8000))
z <- r0_of("adv_blobs_9_10")
add("ADV/CP-2", "adv_blobs_9_10", "9x9 escluso (81) in excluded_df, 10x10 tenuto (100)",
    paste(z$rf$excluded_df$area_px_um2, z$rf$region_df$area_px_um2), "81 escluso / 100 tenuto",
    pf(identical(z$rf$excluded_df$area_px_um2, 81) && identical(z$rf$region_df$area_px_um2, 100)))
# C9
o <- load_in("adv_pixel35")
b1 <- as.numeric(sf::st_bbox(er(o, 0, 0)$region_polygons)); b2 <- as.numeric(sf::st_bbox(er(o, 0, 0, px = 2)$region_polygons))
add("C9", "adv_pixel35", "bbox px=1 / px=2", paste(paste(b1, collapse = ","), "/", paste(b2, collapse = ",")),
    "2,4,3,5 / 4,8,6,10", pf(all(b1 == c(2, 4, 3, 5)) && all(b2 == c(4, 8, 6, 10))))
o <- load_in("adv_L"); pl <- er(o, 0, 0)$region_polygons
bl <- as.numeric(sf::st_bbox(pl))
inn <- lengths(sf::st_intersects(sf::st_sfc(sf::st_point(c(4.5, 0.5)), sf::st_point(c(0.5, 2.5)), sf::st_point(c(0.5, 4.5))), pl))
add("C9", "adv_L", "bbox; (4.5,0.5) e (0.5,2.5) dentro, (0.5,4.5) fuori", paste(paste(bl, collapse = ","), paste(inn, collapse = "")),
    "0,0,5,3; 110", pf(all(bl == c(0, 0, 5, 3)) && identical(as.integer(inn), c(1L, 1L, 0L))))
# C10
o <- load_in("adv_corner"); cl <- o$clust
errs10 <- c(
  missing_px = inherits(try(extract_regions(cl, cluster_col = "intensity_cluster", verbose = FALSE), silent = TRUE), "try-error"),
  px_le0     = inherits(try(extract_regions(cl, 0, cluster_col = "intensity_cluster", verbose = FALSE), silent = TRUE), "try-error") &&
               inherits(try(extract_regions(cl, -1, cluster_col = "intensity_cluster", verbose = FALSE), silent = TRUE), "try-error"),
  no_col     = inherits(try(extract_regions(cl[, c("x", "y")], 1, verbose = FALSE), silent = TRUE), "try-error"),
  non_int    = inherits(try(extract_regions(transform(cl, x = x + 0.5), 1, cluster_col = "intensity_cluster", verbose = FALSE), silent = TRUE), "try-error"))
add("C10", "adv_corner", "errori espliciti (px mancante, px<=0, colonna mancante, coord. non intere)",
    paste(names(errs10)[errs10], collapse = ","), "4/4", pf(all(errs10)))
# I6 stride atteso
add("C4b", "I6_syn6800x6500_c2", "stride atteso 6", gp("I6_syn6800x6500_c2")$stride, "6", pf(gp("I6_syn6800x6500_c2")$stride == 6))

# CP-1 Moore
moore <- do.call(rbind, lapply(res, `[[`, "moore"))
write.csv(moore, file.path(RES, "S1.1_cp1_moore.csv"), row.names = FALSE)
cp1 <- aggregate(cbind(n_px, area_moore_px) ~ id, data = moore, FUN = sum)
cp1$relerr <- (cp1$area_moore_px - cp1$n_px) / cp1$n_px
for (k in seq_len(nrow(cp1))) {
  if (cp1$id[k] == "adv_blobs_9_10") {
    a10 <- moore$area_moore_px[moore$id == "adv_blobs_9_10" & moore$n_px == 100]
    add("CP-1", "adv_blobs_9_10", "Moore sui centri: area del 10x10", a10, "81 (fallisce C1)", pf(a10 == 81))
  } else {
    add("CP-1", cp1$id[k], "Moore sui centri: errore relativo area globale", cp1$relerr[k], "|err| > 0.02 (fallisce C1)", pf(abs(cp1$relerr[k]) > 0.02))
  }
}
# CP-3 nullo
cp3 <- do.call(rbind, lapply(c("A1", "A2A3", "A4", "A5", "A6"), function(a) {
  r <- gp(paste0("real_full_", a)); n <- gp(paste0("null_full_", a))
  data.frame(dataset = a, n_regions_real = r$n_regions_100, n_regions_null = n$n_regions_100,
             ratio_regions = n$n_regions_100 / r$n_regions_100,
             frac_excl_real = r$frac_area_excluded_100, frac_excl_null = n$frac_area_excluded_100,
             n_components_real = r$n_components, n_components_null = n$n_components)
}))
write.csv(cp3, file.path(RES, "S1.1_cp3_null.csv"), row.names = FALSE)
ok3 <- cp3$n_regions_real < cp3$n_regions_null & cp3$frac_excl_real < cp3$frac_excl_null
add("CP-3", "real_full vs null_full", "n_regioni reale < nullo e frazione scartata reale < nullo",
    sprintf("%d/5", sum(ok3)), "5/5", pf(all(ok3)))
# CP-4 verita' nota (I5)
o <- load_in("I5_syn600_c2_labels"); z <- r0_of("I5_syn600_c2_labels")$r0_df
rec <- aggregate(n_px ~ cluster_id, data = z, FUN = sum)
tr <- o$truth; tr$patch <- as.numeric(tr$patch)
mm <- merge(tr, rec, by.x = "patch", by.y = "cluster_id", all = TRUE)
write.csv(transform(mm, n_components = as.integer(table(z$cluster_id)[as.character(mm$patch)])),
          file.path(RES, "S1.1_cp4_truth.csv"), row.names = FALSE)
add("CP-4", "I5_syn600_c2_labels", sprintf("aree delle patch recuperate (%d patch, %d componenti)", nrow(tr), nrow(z)),
    sum(mm$Freq == mm$n_px, na.rm = TRUE), sprintf("%d/%d esatte", nrow(tr), nrow(tr)),
    pf(nrow(mm) == nrow(tr) && all(mm$Freq == mm$n_px)))

# ---- descrittivo B -------------------------------------------------------------
wmed <- function(a) { o <- order(a); a <- a[o]; a[which(cumsum(a) >= sum(a) / 2)[1]] }
desc <- do.call(rbind, lapply(split(regions, regions$id), function(r) {
  p <- per[per$id == r$id[1], ]; cp <- comps[comps$id == r$id[1], ]
  data.frame(id = r$id[1], group = r$group[1], n_regions = nrow(r),
             area_median_um2 = median(r$area_px_um2), area_q25 = unname(quantile(r$area_px_um2, .25)),
             area_q75 = unname(quantile(r$area_px_um2, .75)), area_p90 = unname(quantile(r$area_px_um2, .9)),
             area_weighted_median_um2 = wmed(r$area_px_um2),
             frac_area_excluded = p$frac_area_excluded_100,
             frac_components_single_px = mean(cp$n_px == 1), n_clusters = p$n_clusters)
}))
meta <- lapply(res, `[[`, "meta")
desc$archetype <- vapply(desc$id, function(i) { m <- meta[[i]]; if (!is.null(m$archetype)) m$archetype else NA_character_ }, "")
write.csv(desc, file.path(RES, "S1.1_checkB_descriptive.csv"), row.names = FALSE)

TT <- do.call(rbind, T)
write.csv(TT, file.path(RES, "S1.1_test_results.csv"), row.names = FALSE)
for (k in seq_len(nrow(TT))) cat(sprintf("%-4s %-9s %-24s %s = %s\n", TT$status[k], TT$check[k], TT$input[k], TT$metric[k], TT$value[k]))
cat(sprintf("\n[test_S1.1] PASS %d · WARN %d · FAIL %d (asserzioni: %d)\n",
            sum(TT$status == "PASS"), sum(TT$status == "WARN"), sum(TT$status == "FAIL"), nrow(TT)))
