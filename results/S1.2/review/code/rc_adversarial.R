# rc_adversarial.R — revisione avversariale S1.2, track "code" (scritto da zero dal revisore)
# Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/code/rc_adversarial.R <modo>
#   modo: main | mem | locale | perf_bigd | perf_ref
# Usa seed_centroids() SOLO per generare gli output da verificare; tutte le verifiche sono in rc_verify.py.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel) })
MODE <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "main"
OUT <- "results/S1.2/review/code/out"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SC <- "R/04b2_seed_centroids.R"
newE <- function() { E <- new.env(parent = globalenv()); sys.source(SC, envir = E); E }
E <- newE()
NC <- 20L

sq <- function(x0, y0, w, h = w) st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
mkreg <- function(polys, cluster_id, region_id = seq_along(polys), area = NULL, crs = NULL) {
  sfc <- if (is.null(crs)) st_sfc(polys) else st_sfc(polys, crs = crs)
  if (is.null(area)) area <- as.numeric(st_area(sfc))
  list(region_df = data.frame(region_id = region_id, cluster_id = cluster_id, area_um2 = area, stringsAsFactors = FALSE),
       region_polygons = sfc)
}
f17 <- function(v) sprintf("%.17g", v)
export_polys <- function(reg, d) {
  rows <- list()
  for (k in seq_along(reg$region_polygons)) {
    g <- reg$region_polygons[[k]]
    parts <- if (inherits(g, "MULTIPOLYGON")) unclass(g) else list(unclass(g))
    for (p in seq_along(parts)) for (r in seq_along(parts[[p]])) {
      m <- parts[[p]][[r]]
      if (!length(m)) next
      rows[[length(rows) + 1]] <- data.frame(feat = k, region_id = as.character(reg$region_df$region_id[k]), part = p, ring = r,
                                             vertex = seq_len(nrow(m)), x = f17(m[, 1]), y = f17(m[, 2]))
    }
  }
  if (length(rows)) write.csv(do.call(rbind, rows), file.path(d, "in_polygons.csv"), row.names = FALSE)
}
write_out <- function(o, d) {
  cdf <- o$centroids
  write.csv(data.frame(cell_id = cdf$cell_id, region_id = as.character(cdf$region_id), x = f17(cdf$x), y = f17(cdf$y),
                       cell_type = as.character(cdf$cell_type), cluster_id = as.character(cdf$cluster_id)),
            file.path(d, "centroids.csv"), row.names = FALSE)
  rd <- o$region_df; rd$region_id <- as.character(rd$region_id); rd$cluster_id <- as.character(rd$cluster_id)
  write.csv(rd, file.path(d, "region_df.csv"), row.names = FALSE)
  write.csv(o$type_df, file.path(d, "type_df.csv"), row.names = FALSE)
}
LOG <- list()
run_case <- function(name, reg, ct, comp, seed = 1, E_ = E, write = TRUE, ...) {
  d <- file.path(OUT, name); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  o <- tryCatch(withCallingHandlers(E_$seed_centroids(reg, ct, comp, random_seed = seed, verbose = FALSE, ...),
                  warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
                error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  isE <- inherits(o, "error")
  ctw <- as.data.frame(lapply(ct, function(v) if (is.factor(v)) as.character(v) else v), stringsAsFactors = FALSE)
  write.csv(ctw, file.path(d, "in_cell_types.csv"), row.names = FALSE)
  cpw <- as.data.frame(lapply(comp, function(v) as.character(v)), stringsAsFactors = FALSE)
  write.csv(cpw, file.path(d, "in_comp.csv"), row.names = FALSE)
  rdw <- reg$region_df; rdw$region_id <- as.character(rdw$region_id); rdw$cluster_id <- as.character(rdw$cluster_id)
  write.csv(rdw, file.path(d, "in_region_df.csv"), row.names = FALSE)
  export_polys(reg, d)
  if (!isE && write) write_out(o, d)
  inf <- if (!isE) o$info else list()
  g <- function(k) if (is.null(inf[[k]])) NA else inf[[k]]
  meta <- data.frame(case = name, seed = seed, elapsed_s = el, error = if (isE) conditionMessage(o) else "",
                     warnings = paste(warns, collapse = " | "), n_cells = g("n_cells"), n_target = g("n_target"),
                     n_failed = g("n_failed"), n_attempts = g("n_attempts"), grid_cell_um = g("grid_cell_um"),
                     hardcore = g("hardcore"),
                     class_region_id_in = class(reg$region_df$region_id)[1],
                     class_region_id_out = if (!isE) class(o$centroids$region_id)[1] else NA)
  write.csv(meta, file.path(d, "meta.csv"), row.names = FALSE)
  LOG[[length(LOG) + 1]] <<- meta
  cat(sprintf("[%s] %.2fs err='%s' warn=%d n=%s fail=%s\n", name, el, meta$error, length(warns), meta$n_cells, meta$n_failed))
  invisible(o)
}
mutant_env <- function(fun, from, to) {
  E2 <- newE(); src <- paste(deparse(E2[[fun]], width.cutoff = 500L), collapse = "\n")
  stopifnot(grepl(from, src, fixed = TRUE))
  f <- eval(parse(text = sub(from, to, src, fixed = TRUE))); environment(f) <- E2; assign(fun, f, envir = E2)
  E2
}

if (MODE == "main") {
  ## A01 multipoligono con buco e isola, L concava, scheggia, area 0, coordinate negative, region_id non consecutivi,
  ##     cluster_id factor con livelli in ordine non alfabetico, tipi condivisi fra cluster, tipo d utente, tipo inutilizzato
  off <- c(-5000, -3000)
  donut <- st_polygon(list(rbind(c(0, 0), c(200, 0), c(200, 200), c(0, 200), c(0, 0)),
                           rbind(c(50, 50), c(50, 150), c(150, 150), c(150, 50), c(50, 50))))
  r1 <- st_multipolygon(list(unclass(donut), unclass(sq(80, 80, 40))))
  r2 <- st_polygon(list(rbind(c(1000, 0), c(1300, 0), c(1300, 40), c(1040, 40), c(1040, 300), c(1000, 300), c(1000, 0))))
  r3 <- sq(200, 0, 200)
  r4 <- sq(2000, 2000, 5)
  r5 <- st_polygon(list(rbind(c(500, 500), c(900, 505), c(500, 510), c(500, 500))))
  polys <- lapply(list(r1, r2, r3, r4, r5), function(p) p + off)
  reg1 <- mkreg(polys, cluster_id = factor(c("a", "b", "a", "z", "b"), levels = c("z", "b", "a")), region_id = c(11L, 3L, 7L, 100L, 5L))
  reg1$region_df$area_um2[4] <- 0
  ct1 <- data.frame(cell_type = factor(c("Thi", "Tmid", "Tlo", "Tusr", "Tunused")), density = c(20000, 8000, 1000, 500, 0),
                    min_dist_um = c(NA, NA, NA, 6, NA))
  comp1 <- data.frame(cluster_id = c("a", "a", "a", "a", "b", "b", "b", "z"),
                      cell_type = c("Thi", "Tmid", "Tlo", "Tusr", "Tmid", "Tlo", "Tusr", "Thi"),
                      fraction = c(.6, .25, .1, .05, .5, .3, .2, 1))
  for (s in 1:5) run_case(sprintf("A01_multi_s%d", s), reg1, ct1, comp1, seed = s)

  ## A02 invarianza alla griglia: stesse regioni 1-2, regione 3 vicina (griglia fine) o a 200 mm (griglia 50 um)
  ct_a1 <- data.frame(cell_type = "a1", density = 12419); cp_a1 <- data.frame(cluster_id = 1, cell_type = "a1", fraction = 1)
  run_case("A02_near", mkreg(list(sq(0, 0, 300), sq(400, 0, 300), sq(1000, 0, 30)), 1L, area = c(9e4, 9e4, 900)), ct_a1, cp_a1, seed = 2)
  run_case("A02_far",  mkreg(list(sq(0, 0, 300), sq(400, 0, 300), sq(200000, 0, 30)), 1L, area = c(9e4, 9e4, 900)), ct_a1, cp_a1, seed = 2)

  ## A04 geometrie degeneri
  run_case("A04a_empty_geom", mkreg(list(sq(0, 0, 100), st_polygon()), 1L, area = c(1e4, 100)), ct_a1, cp_a1)
  run_case("A04b_zero_area_poly", mkreg(list(sq(0, 0, 100), st_polygon(list(rbind(c(0, 200), c(100, 200), c(50, 200), c(0, 200))))), 1L,
                                        area = c(1e4, 100)), ct_a1, cp_a1)
  run_case("A04c_only_empty", mkreg(list(st_polygon()), 1L, area = 100), ct_a1, cp_a1)
  run_case("A04d_area0_regions", mkreg(list(sq(0, 0, 100), sq(200, 0, 100)), 1L, area = c(1e4, 0)), ct_a1, cp_a1)

  ## A05 tipo di region_id
  base5 <- list(sq(0, 0, 100), sq(100, 0, 100)); ct5 <- data.frame(cell_type = "a", density = 5000); cp5 <- data.frame(cluster_id = 1, cell_type = "a", fraction = 1)
  run_case("A05a_rid_factor", mkreg(base5, 1L, region_id = factor(c("R10", "R2"))), ct5, cp5)
  run_case("A05b_rid_char", mkreg(base5, 1L, region_id = c("R10", "R2")), ct5, cp5)
  run_case("A05c_rid_double", mkreg(base5, 1L, region_id = c(1.5, 2.5)), ct5, cp5)
  run_case("A05d_rid_big", mkreg(base5, 1L, region_id = c(1e10, 2e10)), ct5, cp5)
  run_case("A05e_rid_factor_revlev", mkreg(list(sq(0, 0, 100), sq(100, 0, 300)), 1L, region_id = factor(c("2", "1"), levels = c("2", "1"))), ct5, cp5)
  run_case("A05f_rid_int_ref", mkreg(list(sq(0, 0, 100), sq(100, 0, 300)), 1L, region_id = c(2L, 1L)), ct5, cp5)

  ## A06 colonne numeriche come factor
  reg6 <- mkreg(list(sq(0, 0, 1000), sq(1000, 0, 1000)), c("c1", "c2"))
  run_case("A06a_density_factor", reg6, data.frame(cell_type = c("x", "y"), density = factor(c("12419", "988"))),
           data.frame(cluster_id = c("c1", "c2"), cell_type = c("x", "y"), fraction = 1))
  run_case("A06b_fraction_factor", mkreg(list(sq(0, 0, 1000)), "c1"), data.frame(cell_type = c("x", "y"), density = c(12419, 988)),
           data.frame(cluster_id = "c1", cell_type = c("x", "y"), fraction = factor(c("0.3", "0.7"))))
  run_case("A06c_density_char", reg6, data.frame(cell_type = c("x", "y"), density = c("12419", "988"), stringsAsFactors = FALSE),
           data.frame(cluster_id = c("c1", "c2"), cell_type = c("x", "y"), fraction = 1))
  run_case("A06d_mindist_factor", mkreg(list(sq(0, 0, 300)), "c1"), data.frame(cell_type = "x", density = 3000, min_dist_um = factor("8")),
           data.frame(cluster_id = "c1", cell_type = "x", fraction = 1))

  ## A07 cluster_id
  reg7 <- mkreg(list(sq(0, 0, 200), sq(200, 0, 200)), c(1e5, 2))
  run_case("A07a_cluster_1e5", reg7, ct5, data.frame(cluster_id = c("100000", "2"), cell_type = "a", fraction = 1))
  run_case("A07b_cluster_factor_unused", mkreg(list(sq(0, 0, 200), sq(200, 0, 200)), factor(c("k1", "k2"), levels = c("k3", "k2", "k1"))), ct5,
           data.frame(cluster_id = c("k1", "k2", "k9"), cell_type = "a", fraction = 1))
  run_case("A07c_cluster_int_vs_char", mkreg(list(sq(0, 0, 200)), 1L), ct5, data.frame(cluster_id = "01", cell_type = "a", fraction = 1))

  ## A09 RNG
  reg9 <- mkreg(list(sq(0, 0, 200)), 1L); rng <- list()
  RNGkind("L'Ecuyer-CMRG"); set.seed(5); s0 <- .Random.seed; k0 <- RNGkind()
  o_lec <- E$seed_centroids(reg9, ct5, cp5, random_seed = 9, verbose = FALSE)
  rng$lecuyer_state_restored <- identical(.Random.seed, s0); rng$lecuyer_kind_restored <- identical(RNGkind(), k0)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(5)
  o_mt <- E$seed_centroids(reg9, ct5, cp5, random_seed = 9, verbose = FALSE)
  rng$output_indep_of_caller_kind <- identical(o_lec$centroids, o_mt$centroids)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding")); set.seed(5); s1 <- .Random.seed; k1 <- RNGkind()
  wr <- character(0)
  o_ro <- withCallingHandlers(E$seed_centroids(reg9, ct5, cp5, random_seed = 9, verbose = FALSE),
                              warning = function(w) { wr <<- c(wr, conditionMessage(w)); invokeRestart("muffleWarning") })
  rng$rounding_state_restored <- identical(.Random.seed, s1); rng$rounding_kind_restored <- identical(RNGkind(), k1)
  rng$rounding_output_same <- identical(o_ro$centroids, o_mt$centroids); rng$rounding_n_warnings <- length(wr)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(77); s2 <- .Random.seed
  e1 <- tryCatch(E$seed_centroids(mkreg(list(sq(0, 0, 200)), 1L, crs = 32632), ct5, cp5, random_seed = 9, verbose = FALSE), error = function(e) e)
  rng$crs_error <- if (inherits(e1, "error")) conditionMessage(e1) else "nessun errore"
  rng$error_path_state_restored <- identical(.Random.seed, s2)
  oa <- E$seed_centroids(reg9, ct5, cp5, random_seed = 1, verbose = FALSE); ob <- E$seed_centroids(reg9, ct5, cp5, random_seed = 1.7, verbose = FALSE)
  rng$seed_1_equals_1.7 <- identical(oa$centroids, ob$centroids)
  rng$seed_2pow31_error <- inherits(tryCatch(E$seed_centroids(reg9, ct5, cp5, random_seed = 2^31, verbose = FALSE), error = function(e) e), "error")
  rng$state_after_2pow31 <- identical(.Random.seed, s2)
  write.csv(data.frame(item = names(rng), value = vapply(rng, as.character, "")), file.path(OUT, "A09_rng.csv"), row.names = FALSE)
  print(unlist(rng))

  ## A10 tipo raro grande in miscela densa (ordine per densita' decrescente) e controfattuale con ordine invertito
  reg10 <- mkreg(list(sq(0, 0, 1000)), 1L)
  ct10 <- data.frame(cell_type = c("lym", "big"), density = c(20000, 150))
  cp10 <- data.frame(cluster_id = 1, cell_type = c("lym", "big"), fraction = c(.99, .01))
  for (s in 1:3) run_case(sprintf("A10_bigrare_s%d", s), reg10, ct10, cp10, seed = s)
  EX3 <- mutant_env(".sc_type_table", "order(-used$density, used$cell_type)", "order(used$density, used$cell_type)")
  for (s in 1:3) run_case(sprintf("A10_bigfirst_X3_s%d", s), reg10, ct10, cp10, seed = s, E_ = EX3)

  ## A11 effetti di bordo: quadrato isolato e due rettangoli adiacenti (ordine di riga invertito nel caso 'swap')
  regsq <- mkreg(list(sq(0, 0, 1000)), 1L)
  reg_lr <- mkreg(list(sq(0, 0, 500, 1000), sq(500, 0, 500, 1000)), 1L)
  reg_rl <- mkreg(list(sq(500, 0, 500, 1000), sq(0, 0, 500, 1000)), 1L)
  for (nm in c("A11_square", "A11_adjLR", "A11_adjRL")) {
    rg <- switch(nm, A11_square = regsq, A11_adjLR = reg_lr, A11_adjRL = reg_rl)
    L <- mclapply(1:20, function(s) { o <- E$seed_centroids(rg, ct_a1, cp_a1, random_seed = s, verbose = FALSE)
      data.frame(seed = s, region_row = match(o$centroids$region_id, rg$region_df$region_id), x = f17(o$centroids$x), y = f17(o$centroids$y)) }, mc.cores = NC)
    dir.create(file.path(OUT, nm), showWarnings = FALSE); write.csv(do.call(rbind, L), file.path(OUT, nm, "centroids_20seeds.csv"), row.names = FALSE)
    cat(nm, "done\n")
  }

  ## A12 C-sat come pre-registrato (1 x 1000 um) e ruotato di 45 gradi
  ct12 <- data.frame(cell_type = "s", density = 2e5, min_dist_um = 10); cp12 <- data.frame(cluster_id = 1, cell_type = "s", fraction = 1)
  run_case("A12a_strip_1x1000", mkreg(list(sq(0, 0, 1000, 1)), 1L), ct12, cp12)
  u <- c(1, 1) / sqrt(2); v <- c(-1, 1) / sqrt(2)
  rot <- st_polygon(list(rbind(c(0, 0), 1000 * u, 1000 * u + v, v, c(0, 0))))
  run_case("A12b_strip_rot45", mkreg(list(rot), 1L), ct12, cp12)

  ## A13 area_um2 incoerente con il poligono (nessun controllo) + coerenza area_um2/st_area sugli input reali di S1.1
  run_case("A13a_area_x10", mkreg(list(sq(0, 0, 100)), 1L, area = 1e5), data.frame(cell_type = "a", density = 1000, min_dist_um = 0), cp5)
  source("R/04b1_extract_regions.R", local = (ER <- new.env()))
  IN <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
  ar <- list()
  for (id in c("I1_syn600_c1", "I2_syn600_c2", "I3_syn600_c3", "I4_syn600_c4", "roi_A1_r1", "roi_A4_f1", "adv_donut_island", "adv_pinch_hole")) {
    ob <- readRDS(file.path(IN, paste0(id, ".rds")))
    st <- if (ob$group == "real_roi") 1L else NULL
    rg <- ER$extract_regions(ob$clust, pixel_size_um = ob$pixel_size_um, cluster_col = ob$cluster_col, stride = st, verbose = FALSE)
    a_geom <- as.numeric(st_area(rg$region_polygons))
    ar[[id]] <- data.frame(input = id, n_regions = nrow(rg$region_df), sum_area_um2 = sum(rg$region_df$area_um2), sum_st_area = sum(a_geom),
                           max_abs_rel_diff = max(abs(rg$region_df$area_um2 - a_geom) / pmax(a_geom, 1e-12)),
                           geom_types = paste(unique(as.character(st_geometry_type(rg$region_polygons))), collapse = "/"),
                           n_invalid = sum(!st_is_valid(rg$region_polygons)))
    if (id == "I1_syn600_c1") I1 <- rg
  }
  write.csv(do.call(rbind, ar), file.path(OUT, "A13b_area_consistency.csv"), row.names = FALSE)

  ## A15 coordinate enormi (1e8 um)
  run_case("A15_offset1e8", mkreg(list(sq(1e8, 1e8, 300), sq(1e8 + 300, 1e8, 300)), c("a", "b")),
           data.frame(cell_type = c("p", "q"), density = c(12419, 3100)),
           data.frame(cluster_id = c("a", "a", "b"), cell_type = c("p", "q", "q"), fraction = c(.7, .3, 1)), seed = 4)

  ## A16 Madow: diretto sull'helper e attraverso la funzione completa (lambda = 7 esatto, 4 tipi, d = 0)
  set.seed(2026); f4 <- c(.37, .29, .21, .13)
  M7 <- t(vapply(1:200000, function(i) E$.sc_allocate_types(7L, f4), integer(4)))
  M1q <- t(vapply(1:200000, function(i) E$.sc_allocate_types(1L, rep(.25, 4)), integer(4)))
  M10 <- t(vapply(1:20000, function(i) E$.sc_allocate_types(10L, c(.1, .2, .3, .4)), integer(4)))
  write.csv(data.frame(case = rep(c("n7", "n1q", "n10"), c(200000, 200000, 20000)), rbind(M7, M1q, M10)), file.path(OUT, "A16a_madow_direct.csv"), row.names = FALSE)
  side <- sqrt(700); nreg <- 2000; nc <- ceiling(sqrt(nreg)); gg <- expand.grid(i = 0:(nc - 1), j = 0:(nc - 1))[1:nreg, ]
  reg16 <- mkreg(lapply(1:nreg, function(k) sq(gg$i[k] * 40, gg$j[k] * 40, side)), 1L, area = rep(700, nreg))
  ct16 <- data.frame(cell_type = c("w", "x", "y", "z"), density = 10000, min_dist_um = 0)
  cp16 <- data.frame(cluster_id = 1, cell_type = c("w", "x", "y", "z"), fraction = f4)
  L16 <- mclapply(1:10, function(s) { o <- E$seed_centroids(reg16, ct16, cp16, random_seed = s, verbose = FALSE)
    tb <- table(factor(o$centroids$region_id, levels = 1:nreg), factor(o$centroids$cell_type, levels = c("w", "x", "y", "z")))
    data.frame(seed = s, region = 1:nreg, unclass(tb)[, 1:4], n_target = o$region_df$n_target) }, mc.cores = NC)
  write.csv(do.call(rbind, L16), file.path(OUT, "A16b_madow_full.csv"), row.names = FALSE)

  ## A17 C1b ricalcolato con 1 000 seed (2 000 quadrati da 200 um2, passo 30 um, 1 185/mm2)
  s17 <- sqrt(200); gg <- expand.grid(i = 0:(nc - 1), j = 0:(nc - 1))[1:2000, ]
  reg17 <- mkreg(lapply(1:2000, function(k) sq(gg$i[k] * 30, gg$j[k] * 30, s17)), 1L, area = rep(200, 2000))
  L17 <- mclapply(1:1000, function(s) { o <- E$seed_centroids(reg17, data.frame(cell_type = "a5", density = 1185),
                  data.frame(cluster_id = 1, cell_type = "a5", fraction = 1), random_seed = s, verbose = FALSE)
    c(seed = s, n_cells = o$info$n_cells, n_target = o$info$n_target, n_failed = o$info$n_failed) }, mc.cores = NC)
  write.csv(do.call(rbind, L17), file.path(OUT, "A17_c1b_1000seeds.csv"), row.names = FALSE)

  ## X: mutanti aggiuntivi del revisore su I1 (stesso catalogo/composizione del test) + M1/M2/M4 del progetto
  ct_test <- data.frame(cell_type = c("T1", "T2", "T3", "T4"), density = c(8000, 3000, 20000, 1000), min_dist_um = c(NA, NA, NA, 15))
  ids <- sort(unique(as.character(I1$region_df$cluster_id)))
  base <- list(data.frame(cell_type = c("T1", "T2", "T3"), fraction = c(.92, .05, .03)), data.frame(cell_type = c("T2", "T3"), fraction = c(.95, .05)),
               data.frame(cell_type = c("T4", "T2"), fraction = c(.85, .15)), data.frame(cell_type = "T3", fraction = 1))
  compI1 <- do.call(rbind, lapply(seq_along(ids), function(k) data.frame(cluster_id = ids[k], if (k <= 4) base[[k]] else data.frame(cell_type = "T2", fraction = 1))))
  run_case("X0_orig_I1", I1, ct_test, compI1, seed = 42)
  run_case("X1_dmin_min_I1", I1, ct_test, compI1, seed = 42, E_ = mutant_env(".sc_place_all", "thr <- (dt + D[ids]) * 0.5", "thr <- pmin(dt, D[ids])"))
  run_case("X2_dmin_max_I1", I1, ct_test, compI1, seed = 42, E_ = mutant_env(".sc_place_all", "thr <- (dt + D[ids]) * 0.5", "thr <- pmax(dt, D[ids])"))
  run_case("X3_order_asc_I1", I1, ct_test, compI1, seed = 42, E_ = EX3)
  run_case("X4_bernoulli_half_I1", I1, ct_test, compI1, seed = 42,
           E_ = mutant_env(".sc_round_regions", "(lambda - fl)", "(lambda - fl)/2"))
  source("tools/S1.2_mutants.R", local = (MU <- new.env()))
  run_case("M1_proj_I1", I1, ct_test, compI1, seed = 42, E_ = MU$sc_mutate(MU$sc_env(), "M1"))
  run_case("M2_proj_I1", I1, ct_test, compI1, seed = 42, E_ = MU$sc_mutate(MU$sc_env(), "M2"))
  run_case("M4_proj_unif", mkreg(list(sq(0, 0, 1000)), 1L, area = 1e6), data.frame(cell_type = "a1", density = 12419),
           data.frame(cluster_id = 1, cell_type = "a1", fraction = 1), seed = 1, E_ = MU$sc_mutate(MU$sc_env(), "M4"))
  run_case("M4_ref_unif", mkreg(list(sq(0, 0, 1000)), 1L, area = 1e6), data.frame(cell_type = "a1", density = 12419),
           data.frame(cluster_id = 1, cell_type = "a1", fraction = 1), seed = 1)
  write.csv(do.call(rbind, LOG), file.path(OUT, "LOG_main.csv"), row.names = FALSE)
}

if (MODE == "mem") {   # regioni lontanissime in diagonale: griglia 250 um, 16 M celle (eseguire sotto ulimit -v)
  run_case("A03_far_diagonal", mkreg(list(sq(0, 0, 300), sq(1e6, 1e6, 300)), 1L), data.frame(cell_type = "a4", density = 28096),
           data.frame(cluster_id = 1, cell_type = "a4", fraction = 1), seed = 1)
  write.csv(do.call(rbind, LOG), file.path(OUT, "LOG_mem.csv"), row.names = FALSE)
}
if (MODE == "locale") {  # parita' di densita': ordine dei tipi dipende dalla collazione?
  lc <- Sys.getlocale("LC_COLLATE")
  run_case(paste0("A08_tie_", gsub("[^A-Za-z0-9]", "", lc)), mkreg(list(sq(0, 0, 400)), 1L),
           data.frame(cell_type = c("a", "B"), density = 5000, min_dist_um = c(2, 6)),
           data.frame(cluster_id = 1, cell_type = c("a", "B"), fraction = .5), seed = 3)
}
if (MODE %in% c("perf_bigd", "perf_ref")) {   # C-perf con un tipo a d grande nella meta' destra (h = 40 um) vs riferimento
  s <- 4000
  if (MODE == "perf_bigd") {
    reg <- mkreg(list(sq(0, 0, 2000, s), sq(2000, 0, 2000, s)), c(1L, 2L))
    ct <- data.frame(cell_type = c("a4", "big"), density = c(28096, 200), min_dist_um = c(NA, 40))
    cp <- data.frame(cluster_id = c(1, 2), cell_type = c("a4", "big"), fraction = 1)
  } else {
    reg <- mkreg(list(sq(0, 0, 2000, s), sq(2000, 0, 2000, s)), c(1L, 1L))
    ct <- data.frame(cell_type = "a4", density = 28096); cp <- data.frame(cluster_id = 1, cell_type = "a4", fraction = 1)
    reg$region_df$area_um2[2] <- 0
  }
  run_case(paste0("A18_", MODE), reg, ct, cp, seed = 42)
}

if (MODE == "a11big") {   # A11 con 200 seed; si salvano solo i punti entro 30 um dal confine condiviso (x = 500) o dal bordo esterno
  ct_a1 <- data.frame(cell_type = "a1", density = 12419); cp_a1 <- data.frame(cluster_id = 1, cell_type = "a1", fraction = 1)
  cfg <- list(A11b_square = mkreg(list(sq(0, 0, 1000)), 1L),
              A11b_adjLR = mkreg(list(sq(0, 0, 500, 1000), sq(500, 0, 500, 1000)), 1L),
              A11b_adjRL = mkreg(list(sq(500, 0, 500, 1000), sq(0, 0, 500, 1000)), 1L))
  for (nm in names(cfg)) {
    rg <- cfg[[nm]]
    L <- mclapply(1:200, function(s) { o <- E$seed_centroids(rg, ct_a1, cp_a1, random_seed = s, verbose = FALSE)
      x <- o$centroids$x; y <- o$centroids$y
      keep <- if (nm == "A11b_square") pmin(x, 1000 - x, y, 1000 - y) < 30 else abs(x - 500) < 30
      data.frame(seed = s, region_row = match(o$centroids$region_id, rg$region_df$region_id)[keep], x = f17(x[keep]), y = f17(y[keep])) }, mc.cores = NC)
    dir.create(file.path(OUT, nm), showWarnings = FALSE); write.csv(do.call(rbind, L), file.path(OUT, nm, "centroids_200seeds_near_edges.csv"), row.names = FALSE)
  }
}
