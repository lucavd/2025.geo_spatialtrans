# S1.3 — check C di tessellate_voronoi() con mutanti (BL-048)
# Pre-registrazione: results/S1.3/S1.3_preregistration.md (approvata da Luca 2026-10-08, commit e2f32fc)
# Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S13_BACKEND=geos|deldir (variante C-10)] [S13_MUTANT=none|M1..M6] \
#     Rscript --vanilla R/testing/test_S1.3.R
# Stampa PASS/FAIL per asserzione; scrive results/S1.3/S1.3_test_<backend>_<mutante>.csv.
# C-5, C-6 (ROI reali, nulli R3) sono in tools/S1.3_roi.R; C-perf in tools/S1.3_perf.R.
#
# Definizioni operative (fissate prima di eseguire; specificano la pre-registrazione senza cambiarne le soglie):
#  Ingressi: sintetici S1.1 I1–I5 (extract_regions di default) con seed 1, 2, 3; avversari S1.1 adv_* (11) con
#    min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1, seed 1. Catalogo e composizione di TEST di
#    R/testing/test_S1.2.R (T1 8000, T2 3000, T3 20000, T4 1000/mm² con min_dist_um 15), non biologici.
#  C-1 : regioni con >= 1 cellula: max |area − Σ territori| / area <= 1e-9
#  C-2 : per regione: (Σ aree − area(unione)) / area <= 1e-9 (sovrapposizioni) e area(unione Δ regione) / area <= 1e-9
#  C-3 : tutti POLYGON, validi, non vuoti; cell_id in ordine crescente; generatore dentro il proprio territorio;
#        info$n_multipart = 0, info$n_repaired = 0
#  C-4 : casi analitici (sotto), entrambe le soglie 1e-9 relative
#  C-7 : cs = 0.1, 1/3, 2/3, 1 → n_iter 1, 1, 2, 3; per cellula area(smussato − originale) <= 1e-9 · area originale;
#        area > 0; Σ lacune non decrescente in n_iter; area(regione) − area(unione smussati) = Σ lacune (rel. 1e-9)
#  C-8 : secondo run identico (serialize); righe di centroids permutate → stessi territori per cell_id (identical)
#  C-9 : errori con messaggio atteso (grepl) per 7 ingressi non validi
#  Aggiunte dopo la revisione avversariale (RA-code-03/04/12, dichiarate nel report): C-4b regola dei frammenti su
#    geometrie costruite (due destinatari, contatto puntiforme, catena in due ordini); C-8c rietichettatura dei cell_id;
#    C-3 anche con 0 agganci falliti; mutante M7 (confine piu' corto) → C-4.
#  Correzioni dello sviluppo (prima dell'esecuzione ufficiale, dichiarate nel report): C-4 esagonale con interne = punti
#    del reticolo con 6 vicini (non "non ritagliate"); avversari: 1 centroide (st_point_on_surface) nelle regioni senza cellule.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel) })
source("R/04b1_extract_regions.R"); source("R/04b2_seed_centroids.R"); source("tools/S1.3_mutants.R"); source("tools/S1.3_variants.R")
MUT <- Sys.getenv("S13_MUTANT", "none"); BACK <- Sys.getenv("S13_BACKEND", "geos")   # motore: variante del confronto C-10
E <- tv_mutate(tv_variant(tv_env(), BACK), MUT)
stopifnot(E$.tv_engine() == BACK)
TV <- function(cen, reg, cs = 0, keep_tiles = FALSE) E$tessellate_voronoi(cen, reg, corner_smoothing = cs,
                                                                         keep_tiles = keep_tiles, verbose = FALSE)
IN  <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
RES <- "results/S1.3"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)
NCORES <- 16L; TOL <- 1e-9
rows <- list()
rec <- function(check, input, metric, value, threshold, status) {
  rows[[length(rows) + 1]] <<- data.frame(check = check, input = input, metric = metric, value = as.character(value),
                                          threshold = threshold, status = status)
  cat(sprintf("%-5s %-6s %-26s %-30s %s (soglia %s)\n", status, check, input, metric, value, threshold))
}
pf <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"
fmt <- function(v) formatC(v, format = "e", digits = 2)

ct_test <- data.frame(cell_type = c("T1", "T2", "T3", "T4"), density = c(8000, 3000, 20000, 1000), min_dist_um = c(NA, NA, NA, 15))
comp_for <- function(cluster_ids) {
  base <- list(`1` = data.frame(cell_type = c("T1", "T2", "T3"), fraction = c(.92, .05, .03)),
               `2` = data.frame(cell_type = c("T2", "T3"), fraction = c(.95, .05)),
               `3` = data.frame(cell_type = c("T4", "T2"), fraction = c(.85, .15)),
               `4` = data.frame(cell_type = "T3", fraction = 1))
  ids <- sort(unique(as.character(cluster_ids)))
  do.call(rbind, lapply(seq_along(ids), function(k) {
    b <- if (k <= 4) base[[k]] else data.frame(cell_type = "T2", fraction = 1)
    data.frame(cluster_id = ids[k], b)
  }))
}
load_case <- function(id) {
  inp <- readRDS(file.path(IN, paste0(id, ".rds")))
  reg <- if (grepl("^adv_", id)) extract_regions(inp$clust, inp$pixel_size_um, min_region_area_um2 = 0, cluster_col = inp$cluster_col,
                                                 simplify_tol_um = 0, stride = 1L, verbose = FALSE) else
    extract_regions(inp$clust, inp$pixel_size_um, cluster_col = inp$cluster_col, verbose = FALSE)
  reg
}
seed_case <- function(reg, seed, adv = FALSE) {
  sc <- suppressWarnings(seed_centroids(reg, ct_test, comp_for(reg$region_df$cluster_id), random_seed = seed, verbose = FALSE))
  cen <- sc$centroids[, c("cell_id", "region_id", "x", "y")]
  if (adv) {                                       # avversari: 1 centroide (punto interno) nelle regioni rimaste vuote
    empty <- which(!reg$region_df$region_id %in% cen$region_id)
    if (length(empty)) {
      pos <- st_coordinates(st_point_on_surface(reg$region_polygons[empty]))
      cen <- rbind(cen, data.frame(cell_id = nrow(cen) + seq_along(empty), region_id = reg$region_df$region_id[empty],
                                   x = pos[, 1], y = pos[, 2]))
    }
  }
  cen
}

# ---- verifiche su un caso ---------------------------------------------------------------------
check_case <- function(id, reg, cen, do_c8 = FALSE) {
  out <- list(); add <- function(...) out[[length(out) + 1]] <<- list(...)
  tv <- tryCatch(TV(cen, reg), error = function(e) e)
  if (inherits(tv, "error")) { add("C-1", id, "errore", conditionMessage(tv), "-", "FAIL"); return(out) }
  rc <- tv$region_check[tv$region_check$n_cells > 0, ]
  c1 <- if (nrow(rc)) max(abs(rc$area_um2 - rc$sum_territory) / rc$area_um2) else 0
  add("C-1", id, "max rel |area-Σterr|", fmt(c1), "<= 1e-9", pf(c1 <= TOL))
  terr <- tv$cell_territories; td <- tv$territory_df
  ov <- 0; sd <- 0
  for (r in which(reg$region_df$region_id %in% rc$region_id)) {
    k <- which(td$region_id == reg$region_df$region_id[r]); A <- as.numeric(st_area(reg$region_polygons[r]))
    u <- st_union(terr[k]); ov <- max(ov, (sum(td$territory_area[k]) - as.numeric(st_area(u))) / A)
    sd <- max(sd, sum(as.numeric(st_area(st_sym_difference(u, reg$region_polygons[r])))) / A)
  }
  add("C-2", id, "max rel sovrapposizioni", fmt(ov), "<= 1e-9", pf(ov <= TOL))
  add("C-2", id, "max rel unione Δ regione", fmt(sd), "<= 1e-9", pf(sd <= TOL))
  gt <- as.character(st_geometry_type(terr))
  pts <- st_sfc(lapply(seq_len(nrow(cen)), function(i) st_point(c(cen$x[i], cen$y[i]))))
  cs <- cen[order(cen$cell_id), ]; pts <- pts[order(cen$cell_id)]
  own <- st_intersects(pts, terr)
  gin <- mean(vapply(seq_along(own), function(i) i %in% own[[i]], TRUE))
  ok3 <- all(gt == "POLYGON") && all(st_is_valid(terr)) && !any(st_is_empty(terr)) && identical(td$cell_id, cs$cell_id) &&
    tv$info$n_multipart == 0 && tv$info$n_repaired == 0 && tv$info$n_snap_failed == 0
  add("C-3", id, "POLYGON valido, ordine, 0 rip.", sprintf("%s; multi %d; rip %d", paste(unique(gt), collapse = "/"), tv$info$n_multipart, tv$info$n_repaired), "100 %", pf(ok3))
  add("C-3", id, "generatore nel territorio", sprintf("%.6f", gin), "1", pf(gin == 1))
  add("ROB", id, "agganci / multiparte / riparazioni", sprintf("%d/%d/%d", tv$info$n_snapped, tv$info$n_multipart, tv$info$n_repaired), "info", "INFO")
  add("ROB", id, "agganci falliti / isolati / rip. smussatura", sprintf("%d/%d/-", tv$info$n_snap_failed, tv$info$n_isolated_cells), "info", "INFO")
  # C-7 smussatura
  prev_gap <- -Inf; ok_n <- TRUE; worst_out <- 0; min_a <- Inf; worst_gap <- 0; mono <- TRUE
  for (csv in c(0.1, 1/3, 2/3, 1)) {
    s <- tryCatch(TV(cen, reg, cs = csv), error = function(e) e)
    if (inherits(s, "error")) { add("C-7", id, sprintf("errore cs=%.2f", csv), conditionMessage(s), "-", "FAIL"); next }
    ok_n <- ok_n && s$info$n_iter == c(1L, 1L, 2L, 3L)[match(csv, c(0.1, 1/3, 2/3, 1))]
    a0 <- td$territory_area
    outa <- vapply(seq_along(terr), function(i) {
      d <- st_difference(s$cell_territories[[i]], terr[[i]]); if (st_is_empty(d)) 0 else as.numeric(st_area(st_sfc(d)))
    }, 0) / a0
    worst_out <- max(worst_out, outa); min_a <- min(min_a, s$territory_df$territory_area)
    g <- sum(s$region_check$gap_area, na.rm = TRUE)
    if (g < prev_gap - TOL * sum(rc$area_um2)) mono <- FALSE
    prev_gap <- g
    for (r in which(reg$region_df$region_id %in% rc$region_id)) {
      k <- which(s$territory_df$region_id == reg$region_df$region_id[r]); A <- as.numeric(st_area(reg$region_polygons[r]))
      gu <- A - as.numeric(st_area(st_union(s$cell_territories[k])))
      worst_gap <- max(worst_gap, abs(gu - s$region_check$gap_area[r]) / A)
    }
  }
  add("C-7", id, "n_iter secondo la formula", ok_n, "1,1,2,3", pf(ok_n))
  add("C-7", id, "max rel area fuori dall'orig.", fmt(worst_out), "<= 1e-9", pf(worst_out <= TOL))
  add("C-7", id, "min area smussata", sprintf("%.3g", min_a), "> 0", pf(min_a > 0))
  add("C-7", id, "lacune monotone in n_iter", mono, "TRUE", pf(mono))
  add("C-7", id, "lacuna(unione) = Σ lacune", fmt(worst_gap), "<= 1e-9", pf(worst_gap <= TOL))
  if (do_c8) {
    tv2 <- TV(cen, reg)
    same <- identical(serialize(tv[c("cell_territories", "territory_df", "region_check")], NULL),
                      serialize(tv2[c("cell_territories", "territory_df", "region_check")], NULL))
    add("C-8", id, "secondo run identico", same, "TRUE", pf(same))
    set.seed(99); p <- sample.int(nrow(cen)); tv3 <- TV(cen[p, ], reg)
    perm <- identical(unclass(tv3$cell_territories), unclass(tv$cell_territories)) && identical(tv3$territory_df, tv$territory_df)
    add("C-8", id, "invarianza alla permutazione", perm, "TRUE", pf(perm))
    # C-8c (dopo la revisione, RA-code-12a): cell_id rietichettati rispetto ai punti → il territorio segue il punto
    cen_r <- cen; cen_r$cell_id <- rev(cen$cell_id); tv4 <- TV(cen_r, reg)
    m1 <- match(cen$cell_id, tv$territory_df$cell_id); m4 <- match(cen_r$cell_id, tv4$territory_df$cell_id)
    da <- max(abs(tv$territory_df$territory_area[m1] - tv4$territory_df$territory_area[m4]) / tv$territory_df$territory_area[m1])
    own4 <- st_intersects(st_sfc(lapply(seq_len(nrow(cen_r)), function(i) st_point(c(cen_r$x[i], cen_r$y[i])))), tv4$cell_territories)
    gin4 <- mean(vapply(seq_along(own4), function(i) m4[i] %in% own4[[i]], TRUE))
    add("C-8", id, "rietichettatura: territorio segue il punto", sprintf("max rel area %s; gen %.4f", fmt(da), gin4), "<= 1e-9; 1", pf(da <= TOL && gin4 == 1))
  }
  out
}

syn <- c("I1_syn600_c1", "I2_syn600_c2", "I3_syn600_c3", "I4_syn600_c4", "I5_syn600_c2_labels")
adv <- sub("\\.rds$", "", list.files(IN, pattern = "^adv_.*\\.rds$"))
jobs <- c(lapply(syn, function(id) lapply(1:3, function(s) list(id = id, seed = s))), lapply(adv, function(id) list(list(id = id, seed = 1))))
jobs <- unlist(jobs, recursive = FALSE)
res <- mclapply(jobs, function(j) {
  reg <- load_case(j$id); cen <- seed_case(reg, j$seed, adv = grepl("^adv_", j$id))
  lab <- sprintf("%s_s%d", sub("_syn600", "", j$id), j$seed)
  if (nrow(cen) == 0) return(list(list("C-1", lab, "nessuna cellula", 0, "-", "SKIP")))
  tryCatch(check_case(lab, reg, cen, do_c8 = j$id %in% syn && j$seed == 1),
           error = function(e) list(list("C-1", lab, "errore", conditionMessage(e), "-", "FAIL")))
}, mc.cores = NCORES, mc.preschedule = FALSE)
for (r in unlist(res, recursive = FALSE)) if (is.list(r)) do.call(rec, r) else rec("C-1", "?", "errore worker", as.character(r), "-", "FAIL")

# ---- C-4 casi analitici --------------------------------------------------------------------------
sq <- function(x0, y0, x1, y1) st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
mkreg <- function(g) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = as.numeric(st_area(st_sfc(g)))),
                          region_polygons = st_sfc(g))
mkcen <- function(x, y) data.frame(cell_id = seq_along(x), region_id = 1L, x = x, y = y)
n_sides <- function(p) {
  m <- unclass(p)[[1]]; m <- m[-nrow(m), , drop = FALSE]
  keep <- c(TRUE, rowSums(abs(diff(m))) > 1e-9); m <- m[keep, , drop = FALSE]
  if (nrow(m) > 1 && sum(abs(m[1, ] - m[nrow(m), ])) <= 1e-9) m <- m[-nrow(m), , drop = FALSE]
  k <- nrow(m); if (k < 3) return(k)
  pr <- m[c(k, 1:(k - 1)), ]; nx <- m[c(2:k, 1), ]
  cr <- (m[, 1] - pr[, 1]) * (nx[, 2] - m[, 2]) - (m[, 2] - pr[, 2]) * (nx[, 1] - m[, 1])
  sum(abs(cr) > 1e-9 * max(abs(m)))
}
c4 <- function(lab, expr) tryCatch(expr, error = function(e) rec("C-4", lab, "errore", conditionMessage(e), "-", "FAIL"))
c4("quadrato", {
  g <- expand.grid(i = 0:9, j = 0:9); tv <- TV(mkcen(5 + 10 * g$i, 5 + 10 * g$j), mkreg(sq(0, 0, 100, 100)))
  e <- max(abs(tv$territory_df$territory_area - 100) / 100); s <- vapply(tv$cell_territories, n_sides, 1L)
  rec("C-4", "quadrato 10x10", "max rel |area − a²|", fmt(e), "<= 1e-9", pf(e <= TOL))
  rec("C-4", "quadrato 10x10", "lati = 4", sprintf("%d/%d", sum(s == 4), length(s)), "100 %", pf(all(s == 4)))
})
c4("esagonale", {
  a <- 10; g <- expand.grid(i = 0:11, j = 0:11)
  x <- g$i * a + (g$j %% 2) * a / 2 + 5; y <- g$j * a * sqrt(3) / 2 + 5
  tv <- TV(mkcen(x, y), mkreg(sq(0, 0, 130, 105)))
  k <- which(g$i %in% 1:10 & g$j %in% 1:10)       # interne del reticolo: tutti e 6 i vicini presenti (def. corretta in sviluppo)
  e <- max(abs(tv$territory_df$territory_area[k] - sqrt(3) / 2 * a^2) / (sqrt(3) / 2 * a^2))
  s <- vapply(tv$cell_territories[k], n_sides, 1L)
  rec("C-4", "esagonale (interne)", sprintf("max rel |area − √3/2·a²| (n=%d)", length(k)), fmt(e), "<= 1e-9", pf(length(k) > 50 && e <= TOL))
  rec("C-4", "esagonale (interne)", "lati = 6", sprintf("%d/%d", sum(s == 6), length(s)), "100 %", pf(all(s == 6)))
})
c4("bisettrice", {
  tv <- TV(mkcen(c(20, 60), c(50, 50)), mkreg(sq(0, 0, 100, 100)))
  e <- max(abs(tv$territory_df$territory_area - c(4000, 6000)) / c(4000, 6000))
  rec("C-4", "2 punti, bisettrice x=40", "max rel area (4000, 6000)", fmt(e), "<= 1e-9", pf(e <= TOL))
})
c4("buco", {
  g <- st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)), rbind(c(40, 40), c(40, 60), c(60, 60), c(60, 40), c(40, 40))))
  tv <- TV(mkcen(c(20, 80, 20, 80), c(20, 20, 80, 80)), mkreg(g))
  e <- max(abs(tv$territory_df$territory_area - 2400) / 2400)
  rec("C-4", "quadrato con buco, 4 punti", "max rel area (2400)", fmt(e), "<= 1e-9", pf(e <= TOL))
})
c4("uno", {
  g <- st_polygon(list(rbind(c(0, 0), c(50, 0), c(50, 10), c(10, 10), c(10, 40), c(0, 40), c(0, 0))))
  tv <- TV(mkcen(5, 5), mkreg(g))
  d <- sum(as.numeric(st_area(st_sym_difference(tv$cell_territories, st_sfc(g))))) / as.numeric(st_area(st_sfc(g)))
  rec("C-4", "regione a L, 1 centroide", "rel area(territorio Δ regione)", fmt(d), "<= 1e-9", pf(d <= TOL))
})
c4("U", {
  U <- st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(80, 100), c(80, 20), c(20, 20), c(20, 100), c(0, 100), c(0, 0))))
  P <- c(50, 10); Q <- c(18, 60)
  half <- function(a, b, L = 1e4) {          # semipiano dei punti piu' vicini ad a che a b
    m <- (a + b) / 2; nn <- (a - b) / sqrt(sum((a - b)^2)); t <- c(-nn[2], nn[1])
    st_polygon(list(rbind(m + L * t, m + L * t + L * nn, m - L * t + L * nn, m - L * t, m + L * t)))
  }
  hq <- st_cast(st_intersection(st_sfc(half(Q, P)), st_sfc(U)), "POLYGON")
  qmain <- hq[lengths(st_intersects(hq, st_sfc(st_point(Q)))) > 0]
  exp_q <- as.numeric(st_area(qmain)); exp_p <- as.numeric(st_area(st_sfc(U))) - exp_q
  tv <- TV(mkcen(c(P[1], Q[1]), c(P[2], Q[2])), mkreg(U))
  e <- max(abs(tv$territory_df$territory_area - c(exp_p, exp_q)) / c(exp_p, exp_q))
  rec("C-4", "U, orfano nel braccio destro", sprintf("max rel area attesa (%.1f, %.1f)", exp_p, exp_q), fmt(e), "<= 1e-9", pf(e <= TOL))
  okf <- tv$info$n_fragments == 1 && tv$info$n_fragments_reassigned == 1 && all(as.character(st_geometry_type(tv$cell_territories)) == "POLYGON")
  rec("C-4", "U, orfano nel braccio destro", "1 frammento riassegnato, POLYGON", sprintf("%d/%d", tv$info$n_fragments_reassigned, tv$info$n_fragments), "1/1", pf(okf))
})

# ---- C-4b regola dei frammenti su geometrie costruite (dopo la revisione, RA-code-03/04/12) -------------
FR <- function(geoms, x, y) E$.tv_fragments(geoms, x, y, rep(1L, length(geoms)))
sqp <- function(x0, y0, x1, y1) sq(x0, y0, x1, y1)
c4("frammenti", {
  # (a) due destinatari: confine 2 (cella 2) e 1 (cella 3); cella 4 tocca l'orfano in un solo punto
  O <- sqp(10, 0, 12, 1)
  g <- list(st_multipolygon(list(unclass(sqp(0, 0, 1, 1)), unclass(O))), sqp(10, 1, 12, 2), sqp(12, 0, 13, 1), sqp(9, -1, 10, 0))
  f <- FR(g, c(0.5, 11, 12.5, 9.5), c(0.5, 1.5, 0.5, -0.5)); a <- vapply(f$geoms, function(z) as.numeric(st_area(st_sfc(z))), 0)
  rec("C-4", "frammenti: due destinatari", "orfano al confine piu' lungo (cella 2)", sprintf("aree %s", paste(a, collapse = "/")), "1/4/1/1", pf(isTRUE(all.equal(a, c(1, 4, 1, 1)))))
  # (b) solo contatto in un punto: l'orfano resta al donatore (multiparte, deviazione 3)
  g <- list(st_multipolygon(list(unclass(sqp(0, 0, 1, 1)), unclass(O))), sqp(9, -1, 10, 0), sqp(30, 30, 31, 31))
  f <- FR(g, c(0.5, 9.5, 30.5), c(0.5, -0.5, 30.5))
  rec("C-4", "frammenti: contatto puntiforme", "nessuna riassegnazione, 1 isolato", sprintf("riass. %d, isolati %d", f$n_reassigned, f$n_isolated_cells), "0, 1", pf(f$n_reassigned == 0 && f$n_isolated_cells == 1))
  # (c) catena: O1 confina con la cella B (1) e con O2 (1); O2 confina con O1 (1) e con la cella D (0.5) → O2 a D in entrambi gli ordini
  O1 <- sqp(10, 0, 11, 1); O2 <- sqp(10, 1, 11, 2)
  gA <- list(st_multipolygon(list(unclass(sqp(0, 0, 1, 1)), unclass(O1))), sqp(11, 0, 12, 1),
             st_multipolygon(list(unclass(sqp(20, 0, 21, 1)), unclass(O2))), sqp(9, 1, 10, 1.5))
  xs <- c(0.5, 11.5, 20.5, 9.5); ys <- c(0.5, 0.5, 0.5, 1.25)
  fa <- FR(gA, xs, ys); aa <- vapply(fa$geoms, function(z) as.numeric(st_area(st_sfc(z))), 0)
  p <- 4:1; fb <- FR(gA[p], xs[p], ys[p]); ab <- vapply(fb$geoms, function(z) as.numeric(st_area(st_sfc(z))), 0)[order(p)]
  rec("C-4", "frammenti: catena, due ordini", "aree identiche; O2 alla cella D", sprintf("%s | %s", paste(aa, collapse = "/"), paste(ab, collapse = "/")), "1/2/1/1.5",
      pf(isTRUE(all.equal(aa, ab)) && isTRUE(all.equal(aa, c(1, 2, 1, 1.5)))))
})

# ---- C-9 errori ----------------------------------------------------------------------------------
R1 <- mkreg(sq(0, 0, 100, 100)); C1 <- mkcen(c(10, 50, 90), c(10, 50, 90))
err_cases <- list(
  list("punti duplicati", function() TV(mkcen(c(10, 10, 50), c(10, 10, 50)), R1), "duplicati"),
  list("centroide fuori regione", function() TV(mkcen(c(10, 150, 50), c(10, 50, 50)), R1), "fuori dalla propria regione"),
  list("region_id assente", function() TV(transform(C1, region_id = c(1L, 1L, 7L)), R1), "assente"),
  list("colonna mancante", function() TV(C1[, c("cell_id", "x", "y")], R1), "senza colonne"),
  list("cs = 1.5", function() TV(C1, R1, cs = 1.5), "corner_smoothing"),
  list("cs = -0.1", function() TV(C1, R1, cs = -0.1), "corner_smoothing"),
  list("region_id factor", function() TV(transform(C1, region_id = factor(region_id)), R1), "factor"))
for (ec in err_cases) {
  m <- tryCatch({ ec[[2]](); "nessun errore" }, error = function(e) conditionMessage(e))
  rec("C-9", ec[[1]], "messaggio", substr(m, 1, 60), ec[[3]], pf(grepl(ec[[3]], m, fixed = TRUE)))
}

tab <- do.call(rbind, rows)
tab$backend <- BACK; tab$mutant <- MUT
f <- file.path(RES, sprintf("S1.3_test_%s_%s.csv", BACK, MUT)); write.csv(tab, f, row.names = FALSE)
cat(sprintf("\nTOTALE %s/%s: %d PASS, %d FAIL, %d SKIP su %d -> %s\n", BACK, MUT, sum(tab$status == "PASS"), sum(tab$status == "FAIL"),
            sum(tab$status == "SKIP"), nrow(tab), f))
