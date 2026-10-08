# S1.2 — check C di seed_centroids() con mutanti (BL-048)
# Pre-registrazione: results/S1.2/S1.2_preregistration.md (approvata da Luca 2026-10-08, commit 60187fb)
# Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S1.2.R
# Stampa PASS/WARN/FAIL per asserzione; scrive results/S1.2/S1.2_test_results.csv e
# results/S1.2/S1.2_mutants.csv. C-perf e' in tools/S1.2_perf.R (processo separato, /usr/bin/time).
#
# Definizioni operative (fissate prima di eseguire, specificano la pre-registrazione senza cambiarne le soglie):
#  Catalogo di TEST (non biologico): T1 8000/mm², T2 3000, T3 20000, T4 1000 con min_dist_um = 15 (utente);
#    composizione per cluster dei sintetici I1–I4: c1 = T1 .92 / T2 .05 / T3 .03; c2 = T2 .95 / T3 .05;
#    c3 = T4 .85 / T2 .15; c4 = T3 1.0 (cluster oltre il 4°: T2 1.0).
#  C1a  : I1–I4 (extract_regions default), n_cells == n_target in ogni regione, n_failed = 0
#  C1b  : 2 000 quadrati da 200 µm² (passo 30 µm), un tipo a 1 185/mm² (d regola 10.93 µm), 200 seed;
#         |media(totale) − 474| < 3·sqrt(Σ p(1−p))/sqrt(200)
#  C4a  : I1–I4 e input C4b: |n_i − n·f_i| < 1 per regione e tipo
#  C4c  : I1–I4: per cluster |n_i/Σn − f_i| / f_i ≤ 10 % per i tipi con atteso ≥ 50
#  C4b  : 2 000 quadrati da 400 µm², miscela T1 .90 / T2 .08 / T3 .02, 200 seed; per tipo |media − atteso| < 3 SE (SE empirico)
#  C-in : ray casting even-odd scritto qui (indipendente da GEOS) su I1–I4, adv_donut, adv_donut_island, adv_L, roi_A1_r1
#  C-d  : spatstat closepairs (indipendente dalla griglia) su I1–I4 e roi_A4_f1: d >= (d_i+d_j)/2 − 1e-9 per tutte le coppie,
#         incluse quelle fra regioni diverse
#  C-unif: quadrato 1 000 × 1 000 µm, densità A1 12 419, d regola; quadrati 10×10, chisq.test p > 0.001 in ≥ 19/20 seed
#  C-mix: regione 1 mm², f = .5/.5, ρ = 10 000/1 000 → target 1/(.5/1e4 + .5/1e3) = 1 818.18; n_target ∈ {floor, ceil}
#  C-seed: identical(centroids, region_df) a parita' di seed; diverso a seed diverso; .Random.seed invariato; assente resta assente
#  C-sat: striscia 1 × 100 µm, 200 000/mm², min_dist_um 10 → n target 20, al piu' 11 posizionabili: warning, n_failed >= 9, < 60 s
#  Mutanti: M1 → C-in, M2 → C-d, M3 → C1b, M3b → C4b, M4 → C-unif, M5 → C-seed, M6 → C-mix; ognuno DEVE far fallire il proprio check
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(parallel); library(spatstat.geom) })
source("R/04b1_extract_regions.R")
source("tools/S1.2_mutants.R")
E0 <- sc_env()
IN  <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
RES <- "results/S1.2"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)
NCORES <- 24L
rows <- list()
rec <- function(check, input, metric, value, threshold, status) {
  rows[[length(rows) + 1]] <<- data.frame(check = check, input = input, metric = metric,
                                          value = as.character(value), threshold = threshold, status = status)
  cat(sprintf("%-5s %-8s %-22s %-28s %s (soglia %s)\n", status, check, input, metric, value, threshold))
}
pf <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"

# ---- input -------------------------------------------------------------------
ct_test <- data.frame(cell_type = c("T1", "T2", "T3", "T4"), density = c(8000, 3000, 20000, 1000),
                      min_dist_um = c(NA, NA, NA, 15))
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
load_regions <- function(id, ...) {
  o <- readRDS(file.path(IN, paste0(id, ".rds")))
  st <- if (o$group == "real_roi") 1L else NULL
  extract_regions(o$clust, pixel_size_um = o$pixel_size_um, cluster_col = o$cluster_col, stride = st, verbose = FALSE, ...)
}
sq <- function(x0, y0, s) st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
squares <- function(n, area, step = 30) {
  s <- sqrt(area); nc <- ceiling(sqrt(n))
  g <- expand.grid(i = 0:(nc - 1), j = 0:(nc - 1))[seq_len(n), ]
  list(region_df = data.frame(region_id = seq_len(n), cluster_id = 1L, area_um2 = area),
       region_polygons = st_sfc(lapply(seq_len(n), function(k) sq(g$i[k] * step, g$j[k] * step, s))))
}
syn_ids <- c("I1_syn600_c1", "I2_syn600_c2", "I3_syn600_c3", "I4_syn600_c4")
syn <- setNames(lapply(syn_ids, load_regions), syn_ids)

# ray casting even-odd su tutti gli anelli di un POLYGON/MULTIPOLYGON (indipendente da GEOS)
pip_evenodd <- function(x, y, g) {
  rings <- if (inherits(g, "MULTIPOLYGON")) unlist(unclass(g), recursive = FALSE) else unclass(g)
  inside <- logical(length(x))
  for (R in rings) {
    n <- nrow(R)
    for (e in seq_len(n - 1)) {
      x1 <- R[e, 1]; y1 <- R[e, 2]; x2 <- R[e + 1, 1]; y2 <- R[e + 1, 2]
      cr <- ((y1 > y) != (y2 > y)) & (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1)
      inside <- xor(inside, cr)
    }
  }
  inside
}
check_inside <- function(o, reg) {
  ok <- logical(nrow(o$centroids))
  sp <- split(seq_len(nrow(o$centroids)), o$centroids$region_id)
  for (rid in names(sp)) {
    ii <- sp[[rid]]; g <- reg$region_polygons[[match(as.integer(rid), reg$region_df$region_id)]]
    ok[ii] <- pip_evenodd(o$centroids$x[ii], o$centroids$y[ii], g)
  }
  ok
}
check_dmin <- function(o) {
  cdf <- o$centroids
  if (nrow(cdf) < 2) return(list(viol = 0L, viol_cross = 0L, n_pairs = 0L, min_ratio = NA))
  dt <- o$type_df$min_dist_um[match(as.character(cdf$cell_type), o$type_df$cell_type)]
  W <- owin(range(cdf$x) + c(-1, 1), range(cdf$y) + c(-1, 1))
  X <- ppp(cdf$x, cdf$y, window = W, check = FALSE)
  cp <- closepairs(X, rmax = 1.5 * max(dt), what = "indices", twice = FALSE)   # 1.5x: conta anche le coppie vicine (non solo le violazioni)
  dd <- sqrt((cdf$x[cp$i] - cdf$x[cp$j])^2 + (cdf$y[cp$i] - cdf$y[cp$j])^2)
  thr <- (dt[cp$i] + dt[cp$j]) / 2
  v <- dd < thr - 1e-9
  cross <- cdf$region_id[cp$i] != cdf$region_id[cp$j]
  list(viol = sum(v), viol_cross = sum(v & cross), n_pairs = length(dd), n_cross = sum(cross),
       min_ratio = if (length(dd)) min(dd / thr) else NA)
}

# ---- C1a / C4a / C4c / C-in / C-d su I1–I4 -----------------------------------
for (id in syn_ids) {
  reg <- syn[[id]]; comp <- comp_for(reg$region_df$cluster_id)
  o <- E0$seed_centroids(reg, ct_test, comp, random_seed = 42, verbose = FALSE)
  rd <- o$region_df
  rec("C1a", id, "regioni con n != target", sum(rd$n_cells != rd$n_target), "0", pf(all(rd$n_cells == rd$n_target)))
  rec("C1a", id, "n_failed", sum(rd$n_failed), "0", pf(sum(rd$n_failed) == 0))
  # C4a
  tab <- table(factor(o$centroids$region_id, levels = rd$region_id), factor(o$centroids$cell_type, levels = o$type_df$cell_type))
  dev <- 0
  for (i in seq_len(nrow(rd))) {
    cc <- comp[comp$cluster_id == as.character(rd$cluster_id[i]), ]
    ni <- tab[i, cc$cell_type]; dev <- max(dev, abs(ni - rd$n_cells[i] * cc$fraction))
  }
  rec("C4a", id, "max |n_i - n f_i|", sprintf("%.3f", dev), "< 1", pf(dev < 1))
  # C4c
  worst <- 0; nchk <- 0
  for (cl in unique(as.character(rd$cluster_id))) {
    cc <- comp[comp$cluster_id == cl, ]; rr <- rd$region_id[as.character(rd$cluster_id) == cl]
    sel <- o$centroids$region_id %in% rr; ntot <- sum(sel)
    for (k in seq_len(nrow(cc))) {
      expct <- ntot * cc$fraction[k]; if (expct < 50) next
      nchk <- nchk + 1
      worst <- max(worst, abs(sum(o$centroids$cell_type[sel] == cc$cell_type[k]) / ntot - cc$fraction[k]) / cc$fraction[k])
    }
  }
  rec("C4c", id, sprintf("max scarto rel. (%d tipi)", nchk), sprintf("%.4f", worst), "<= 0.10", pf(worst <= 0.10))
  ins <- check_inside(o, reg)
  rec("C-in", id, "centroidi fuori regione", sum(!ins), "0", pf(all(ins)))
  dm <- check_dmin(o)
  rec("C-d", id, sprintf("violazioni (%d coppie)", dm$n_pairs), dm$viol, "0", pf(dm$viol == 0))
  rec("C-d", id, sprintf("violaz. fra regioni (%d coppie)", dm$n_cross), dm$viol_cross, "0", pf(dm$viol_cross == 0))
}

# ---- C-in su input avversari e reali; C-d su un ROI reale denso ------------
ct_a <- data.frame(cell_type = "a", density = 12419)
for (id in c("adv_donut", "adv_donut_island", "adv_L", "roi_A1_r1")) {
  reg <- load_regions(id, min_region_area_um2 = 0)
  comp <- data.frame(cluster_id = unique(as.character(reg$region_df$cluster_id)), cell_type = "a", fraction = 1)
  dens <- if (grepl("^adv", id)) 2e5 else 12419             # adversari piccoli: densita' alta per avere punti
  o <- E0$seed_centroids(reg, data.frame(cell_type = "a", density = dens, min_dist_um = if (grepl("^adv", id)) 0.5 else NA),
                         comp, random_seed = 42, verbose = FALSE)
  ins <- check_inside(o, reg)
  rec("C-in", id, sprintf("fuori regione (n=%d)", nrow(o$centroids)), sum(!ins), "0", pf(all(ins) && nrow(o$centroids) > 0))
}
reg_a4 <- load_regions("roi_A4_f1")
comp_a4 <- data.frame(cluster_id = unique(as.character(reg_a4$region_df$cluster_id)), cell_type = "a4", fraction = 1)
o_a4 <- E0$seed_centroids(reg_a4, data.frame(cell_type = "a4", density = 28096), comp_a4, random_seed = 42, verbose = FALSE)
dm <- check_dmin(o_a4)
rec("C-d", "roi_A4_f1", sprintf("violazioni (%d coppie)", dm$n_pairs), dm$viol, "0", pf(dm$viol == 0))
rec("C-d", "roi_A4_f1", sprintf("violaz. fra regioni (%d)", dm$n_cross), dm$viol_cross, "0", pf(dm$viol_cross == 0))

# ---- C1b: non distorsione su regioni piccole --------------------------------
run_c1b <- function(E, seeds) {
  reg <- squares(2000, 200)
  unlist(mclapply(seeds, function(s) E$seed_centroids(reg, data.frame(cell_type = "a5", density = 1185),
         data.frame(cluster_id = 1, cell_type = "a5", fraction = 1), random_seed = s, verbose = FALSE)$info$n_cells,
         mc.cores = NCORES))
}
lam <- 1185 * 200 / 1e6
se_th <- sqrt(2000 * lam * (1 - lam)) / sqrt(200)
tot <- run_c1b(E0, 1:200)
rec("C1b", "2000x200um2", sprintf("media totale (SD %.1f)", sd(tot)), sprintf("%.2f", mean(tot)),
    sprintf("474 ± %.2f", 3 * se_th), pf(abs(mean(tot) - 474) < 3 * se_th))

# ---- C4b: non distorsione per tipo su regioni piccole -----------------------
run_c4b <- function(E, seeds) {
  reg <- squares(2000, 400)
  ctm <- data.frame(cell_type = c("T1", "T2", "T3"), density = c(8000, 3000, 20000))
  cm <- data.frame(cluster_id = 1, cell_type = c("T1", "T2", "T3"), fraction = c(.90, .08, .02))
  do.call(rbind, mclapply(seeds, function(s) {
    o <- E$seed_centroids(reg, ctm, cm, random_seed = s, verbose = FALSE)
    as.numeric(table(factor(o$centroids$cell_type, levels = c("T1", "T2", "T3"))))
  }, mc.cores = NCORES))
}
rho_mix <- 1 / sum(c(.90, .08, .02) / c(8000, 3000, 20000))
lam4 <- 2000 * rho_mix * 400 / 1e6
M <- run_c4b(E0, 1:200)
for (k in 1:3) {
  ex <- lam4 * c(.90, .08, .02)[k]; se <- sd(M[, k]) / sqrt(nrow(M))
  rec("C4b", "2000x400um2", sprintf("media T%d (atteso %.1f)", k, ex), sprintf("%.2f", mean(M[, k])),
      sprintf("± %.2f (3 SE)", 3 * se), pf(abs(mean(M[, k]) - ex) < 3 * se))
}
o4b <- E0$seed_centroids(squares(2000, 400), data.frame(cell_type = c("T1", "T2", "T3"), density = c(8000, 3000, 20000)),
                         data.frame(cluster_id = 1, cell_type = c("T1", "T2", "T3"), fraction = c(.90, .08, .02)), random_seed = 1, verbose = FALSE)
tb <- table(factor(o4b$centroids$region_id, levels = 1:2000), factor(o4b$centroids$cell_type, levels = c("T1", "T2", "T3")))
dev4 <- max(abs(unclass(tb) - outer(o4b$region_df$n_cells, c(.90, .08, .02))))
rec("C4a", "2000x400um2", "max |n_i - n f_i|", sprintf("%.3f", dev4), "< 1", pf(dev4 < 1))

# ---- C-unif ------------------------------------------------------------------
reg_sq <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = 1e6), region_polygons = st_sfc(sq(0, 0, 1000)))
run_unif <- function(E, seeds) unlist(mclapply(seeds, function(s) {
  o <- E$seed_centroids(reg_sq, data.frame(cell_type = "a1", density = 12419), data.frame(cluster_id = 1, cell_type = "a1", fraction = 1),
                        random_seed = s, verbose = FALSE)
  q <- table(factor(pmin(floor(o$centroids$x / 100), 9) + 10 * pmin(floor(o$centroids$y / 100), 9), levels = 0:99))
  suppressWarnings(chisq.test(as.numeric(q))$p.value)
}, mc.cores = NCORES))
pu <- run_unif(E0, 1:20)
rec("C-unif", "1mm2 A1", "seed con p > 0.001", sprintf("%d/20 (p min %.3g)", sum(pu > 0.001), min(pu)), ">= 19/20", pf(sum(pu > 0.001) >= 19))

# ---- C-mix -------------------------------------------------------------------
run_mix <- function(E) {
  E$seed_centroids(reg_sq, data.frame(cell_type = c("hi", "lo"), density = c(10000, 1000)),
                   data.frame(cluster_id = 1, cell_type = c("hi", "lo"), fraction = c(.5, .5)), random_seed = 3, verbose = FALSE)
}
om <- run_mix(E0); tgt <- 1 / (.5 / 1e4 + .5 / 1e3)
okm <- abs(om$region_df$target_density_weighted - tgt) < 1e-9 && om$region_df$n_target %in% c(floor(tgt), ceiling(tgt))
rec("C-mix", "1mm2 50/50", "densita' target; n_target", sprintf("%.2f; %d", om$region_df$target_density_weighted, om$region_df$n_target),
    sprintf("%.2f; {%d,%d}", tgt, floor(tgt), ceiling(tgt)), pf(okm))

# ---- C-seed ------------------------------------------------------------------
run_seed <- function(E) {
  reg <- syn[[1]]; comp <- comp_for(reg$region_df$cluster_id)
  a <- E$seed_centroids(reg, ct_test, comp, random_seed = 11, verbose = FALSE)
  b <- E$seed_centroids(reg, ct_test, comp, random_seed = 11, verbose = FALSE)
  d <- E$seed_centroids(reg, ct_test, comp, random_seed = 12, verbose = FALSE)
  set.seed(999); s0 <- .Random.seed
  invisible(E$seed_centroids(reg, ct_test, comp, random_seed = 11, verbose = FALSE))
  same_state <- identical(s0, .Random.seed)
  rm(".Random.seed", envir = globalenv())
  invisible(E$seed_centroids(reg, ct_test, comp, random_seed = 11, verbose = FALSE))
  absent <- !exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  c(same = identical(a$centroids, b$centroids) && identical(a$region_df, b$region_df),
    differ = !identical(a$centroids, d$centroids), state = same_state, absent = absent)
}
rs <- run_seed(E0)
for (nm in names(rs)) rec("C-seed", "I1", nm, rs[[nm]], "TRUE", pf(rs[[nm]]))

# ---- C-sat -------------------------------------------------------------------
strip <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = 100),
              region_polygons = st_sfc(st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 1), c(0, 1), c(0, 0))))))
w <- NULL
t_sat <- system.time(os <- withCallingHandlers(
  E0$seed_centroids(strip, data.frame(cell_type = "s", density = 2e5, min_dist_um = 10),
                    data.frame(cluster_id = 1, cell_type = "s", fraction = 1), random_seed = 1, verbose = FALSE),
  warning = function(cw) { w <<- conditionMessage(cw); invokeRestart("muffleWarning") }))[["elapsed"]]
rec("C-sat", "striscia 1x100", "n_target; n_cells; n_failed", sprintf("%d; %d; %d", os$info$n_target, os$info$n_cells, os$info$n_failed),
    "20; <= 11; >= 9", pf(os$info$n_target == 20 && os$info$n_cells <= 11 && os$info$n_failed >= 9))
rec("C-sat", "striscia 1x100", "warning emesso; tempo (s)", sprintf("%s; %.1f", !is.null(w), t_sat), "TRUE; < 60", pf(!is.null(w) && t_sat < 60))

# ---- C-schema e input non validi ---------------------------------------------
o1 <- E0$seed_centroids(syn[[1]], ct_test, comp_for(syn[[1]]$region_df$cluster_id), random_seed = 1, verbose = FALSE)
rec("C-schema", "I1", "colonne", paste(names(o1$centroids), collapse = ","), "cell_id,region_id,x,y,cell_type,cluster_id",
    pf(identical(names(o1$centroids), c("cell_id", "region_id", "x", "y", "cell_type", "cluster_id"))))
rec("C-schema", "I1", "cell_id == 1..N", identical(o1$centroids$cell_id, seq_len(nrow(o1$centroids))), "TRUE",
    pf(identical(o1$centroids$cell_id, seq_len(nrow(o1$centroids)))))
rec("C-schema", "I1", "cluster_id coerente", all(o1$centroids$cluster_id == o1$region_df$cluster_id[match(o1$centroids$region_id, o1$region_df$region_id)]),
    "TRUE", pf(all(o1$centroids$cluster_id == o1$region_df$cluster_id[match(o1$centroids$region_id, o1$region_df$region_id)])))
err <- function(expr) inherits(tryCatch({ suppressWarnings(expr); NULL }, error = function(e) e), "error")
r1 <- reg_sq
rec("C-schema", "1mm2", "tipo ignoto -> errore", err(E0$seed_centroids(r1, ct_a, data.frame(cluster_id = 1, cell_type = "zz", fraction = 1))), "TRUE",
    pf(err(E0$seed_centroids(r1, ct_a, data.frame(cluster_id = 1, cell_type = "zz", fraction = 1)))))
z <- err(E0$seed_centroids(r1, data.frame(cell_type = c("a", "ecm"), density = c(1000, 0)), data.frame(cluster_id = 1, cell_type = c("a", "ecm"), fraction = c(.5, .5))))
rec("C-schema", "1mm2", "densita' 0 in composizione -> errore", z, "TRUE", pf(z))
z <- err(E0$seed_centroids(r1, data.frame(cell_type = "a", density = -5), data.frame(cluster_id = 1, cell_type = "a", fraction = 1)))
rec("C-schema", "1mm2", "densita' < 0 -> errore", z, "TRUE", pf(z))
z <- err(E0$seed_centroids(r1, ct_a, data.frame(cluster_id = 2, cell_type = "a", fraction = 1)))
rec("C-schema", "1mm2", "cluster senza composizione -> errore", z, "TRUE", pf(z))
wn <- NULL
on <- withCallingHandlers(E0$seed_centroids(r1, data.frame(cell_type = c("a", "b"), density = c(1000, 2000)),
                                            data.frame(cluster_id = 1, cell_type = c("a", "b"), fraction = c(.45, .45)), random_seed = 1, verbose = FALSE),
                          warning = function(cw) { wn <<- conditionMessage(cw); invokeRestart("muffleWarning") })
okn <- !is.null(wn) && abs(on$region_df$target_density_weighted - 1 / (.5 / 1000 + .5 / 2000)) < 1e-9
rec("C-schema", "1mm2", "frazioni .45/.45 -> warning + normalizzate", okn, "TRUE", pf(okn))
empty <- list(region_df = data.frame(region_id = integer(0), cluster_id = integer(0), area_um2 = numeric(0)), region_polygons = st_sfc())
oe <- tryCatch(E0$seed_centroids(empty, ct_a, data.frame(cluster_id = 1, cell_type = "a", fraction = 1), verbose = FALSE), error = function(e) e)
rec("C-schema", "vuoto", "0 regioni -> 0 righe senza errore", !inherits(oe, "error") && nrow(oe$centroids) == 0, "TRUE",
    pf(!inherits(oe, "error") && nrow(oe$centroids) == 0))

# ---- mutanti: ognuno deve far fallire il proprio check -----------------------
mrows <- list()
mrec <- function(id, check, value, fails) {
  st <- if (fails) "PASS" else "FAIL"     # PASS = il mutante e' stato rilevato
  mrows[[length(mrows) + 1]] <<- data.frame(mutant = id, check = check, value = as.character(value), detected = fails, status = st)
  rec("MUT", id, paste("rileva", check), value, "check fallisce", st)
}
# M1 -> C-in (I2)
E <- sc_mutate(E0, "M1"); reg <- syn[[2]]
o <- E$seed_centroids(reg, ct_test, comp_for(reg$region_df$cluster_id), random_seed = 42, verbose = FALSE)
nout <- sum(!check_inside(o, reg)); mrec("M1", "C-in", sprintf("%d fuori", nout), nout > 0)
# M2 -> C-d (I1)
E <- sc_mutate(E0, "M2"); reg <- syn[[1]]
o <- E$seed_centroids(reg, ct_test, comp_for(reg$region_df$cluster_id), random_seed = 42, verbose = FALSE)
v <- check_dmin(o)$viol; mrec("M2", "C-d", sprintf("%d violazioni", v), v > 0)
# M3 -> C1b
tm <- run_c1b(sc_mutate(E0, "M3"), 1:20); mrec("M3", "C1b", sprintf("media %.1f", mean(tm)), abs(mean(tm) - 474) >= 3 * se_th)
# M3b -> C4b (T3 raro)
Mm <- run_c4b(sc_mutate(E0, "M3b"), 1:20)
mrec("M3b", "C4b", sprintf("T3 media %.1f vs %.1f", mean(Mm[, 3]), lam4 * .02), abs(mean(Mm[, 3]) - lam4 * .02) >= 3 * sd(M[, 3]) / sqrt(nrow(M)))
# M4 -> C-unif
pm <- run_unif(sc_mutate(E0, "M4"), 1:20); mrec("M4", "C-unif", sprintf("%d/20 p > 0.001", sum(pm > 0.001)), sum(pm > 0.001) < 19)
# M5 -> C-seed
rm5 <- run_seed(sc_mutate(E0, "M5")); mrec("M5", "C-seed", paste(names(rm5), rm5, sep = "=", collapse = " "), !all(rm5))
# M6 -> C-mix
o6 <- run_mix(sc_mutate(E0, "M6"))
mrec("M6", "C-mix", sprintf("%.1f", o6$region_df$target_density_weighted), abs(o6$region_df$target_density_weighted - tgt) > 1e-9)

res <- do.call(rbind, rows)
write.csv(res, file.path(RES, "S1.2_test_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, mrows), file.path(RES, "S1.2_mutants.csv"), row.names = FALSE)
cat(sprintf("\nTOTALE: %d asserzioni, %d PASS, %d WARN, %d FAIL\n", nrow(res), sum(res$status == "PASS"), sum(res$status == "WARN"), sum(res$status == "FAIL")))
