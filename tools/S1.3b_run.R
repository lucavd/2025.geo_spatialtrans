# tools/S1.3b_run.R — S1.3b (BL-069): run ufficiale. Pre-registrazione results/S1.3b/S1.3b_preregistration.md (d96fea1).
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.3b_run.R <stage> [outdir]
#   stage = md5 | null | grid | cp1 | cp2 | mutants ; outdir default /mnt/micron/geo_spatialtrans/S1.3b ; stadi rilanciabili.
# Definizioni operative (dichiarate qui, prima del run):
#  - nulli: rigenerati con gen_points()/s13_tessellate() di tools/S1.3_roi.R (stessi seed); regione = w$poly (una regione).
#  - C-3: 1e5 punti uniformi nella regione V0 (st_sample, seed 20261010 + indice del job); stessi punti per tutte le varianti.
#  - C-5/B-1 rispetto alle aree V0 del pacchetto (territory_area); interne ricalcolate con le aree della variante (regola R3).
#  - C-6: griglia di quadrati da 2 µm allineata a multipli di 2 µm; pezzi = intersezione cella x bin con l'operazione della
#    variante; tolleranza = g * (somma dei perimetri dei pezzi [+ perimetro di bin ∩ regione]) con g = 1e-6 anche per V0
#    (stessa scala di V1/V2); bin interni (predicato st_within) = area del bin; bin di bordo = bin ∩ regione della variante.
#  - CP-1b placebo: tanti spostamenti (stesse ampiezze) e inserimenti collineari quanti in CP-1a, su vertici/lati a più di
#    1e-6 µm dall'altro poligono.
#  - CP-2: perturbazione U(-δ, δ) per coordinata, indipendente per poligono; regione non perturbata.
REAL_ARGS <- commandArgs(TRUE)
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
source("tools/S1.3_roi.R")          # rois, ARCH, gen_points, s13_tessellate, c123, one_region, R3, TOL, BASE_SEED, NREP
rm(commandArgs)
source("tools/S1.3b_overlay.R")
STAGE <- REAL_ARGS[1]; OUTB <- if (length(REAL_ARGS) > 1) REAL_ARGS[2] else "/mnt/micron/geo_spatialtrans/S1.3b"
OUT13 <- "/mnt/micron/geo_spatialtrans/S1.3"; SEL_SEED <- 20261010
NC <- as.integer(Sys.getenv("S13B_CORES", "40")); MODELS <- c("CSR", "RSA", "RSArule"); VARS <- c("V0", "V0g", "V1", "V2", "V3")
for (d in c("null", "grid", "cp1", "cp2", "mutants")) dir.create(file.path(OUTB, d), recursive = TRUE, showWarnings = FALSE)

regen <- function(i, model, rp) {
  A <- rois$archetype[i]; roi <- rois$roi_id[i]; a <- ARCH[ARCH$archetype == A, ]
  w <- r3_roi_window(A, roi, rois)
  n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  p <- gen_points(i, model, rp, w, n_obs, a)
  s <- s13_tessellate(p$x, p$y, w, backend = "geos")
  list(A = A, roi = roi, w = w, p = p, s = s, terr = s$tv$cell_territories, reg = st_sfc(w$poly[[1]]))
}
reg_features <- function(reg) { g <- reg[[1]]; parts <- if (inherits(g, "MULTIPOLYGON")) unclass(g) else list(unclass(g))
  c(n_parts = length(parts), n_holes = sum(lengths(parts) - 1L), n_vertices = sum(vapply(parts, function(p) sum(vapply(p, nrow, 1L)), 1))) }
job_index <- function(i, model, rp) (i - 1L) * 60L + (match(model, MODELS) - 1L) * 20L + rp
fname <- function(A, roi, model, rp) sprintf("%s_%s_%s_%02d.rds", A, roi, model, rp)
par_run <- function(jobs, FUN, cores = NC) {
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) tryCatch(FUN(j), error = function(e) paste("ERRORE:", conditionMessage(e))),
                  mc.cores = cores, mc.preschedule = FALSE)
  err <- vapply(res, function(r) is.character(r) && length(r) == 1 && startsWith(r, "ERRORE"), TRUE)
  for (k in which(err)) cat(sprintf("[%s] job %d: %s\n", STAGE, k, res[[k]]))
  cat(sprintf("[%s] errori %d su %d\n", STAGE, sum(err), length(res))); if (any(err)) quit(status = 1)
}

# ---- md5 degli ingressi ---------------------------------------------------------------------------------------------
if (STAGE == "md5") { st <- system("md5sum -c --quiet results/R3/R3_inputs_md5.txt"); cat("md5 exit", st, "\n"); quit(status = st) }

# ---- stadio null: C-0..C-5, B-1, B-2 (riassunti per variante) ------------------------------------------------------
null_one <- function(j) {
  i <- j$i; model <- j$model; rp <- j$rp; A <- rois$archetype[i]; roi <- rois$roi_id[i]
  f <- file.path(OUTB, "null", fname(A, roi, model, rp)); if (file.exists(f)) return(NULL)
  t0 <- proc.time()[[3]]; G <- regen(i, model, rp); t_tess <- proc.time()[[3]] - t0
  s <- G$s; tv <- s$tv; w <- G$w; old <- readRDS(file.path(OUT13, "null", fname(A, roi, model, rp)))
  c0 <- c123(tv, w, G$p$x, G$p$y)
  P <- s3b_perim(G$terr); Preg <- s3b_perim(G$reg)
  set.seed(SEL_SEED + job_index(i, model, rp)); pts <- st_sample(G$reg, 1e5)
  in_frame <- lengths(st_within(tv$tiles, w$frame)) > 0; li <- r3_local_intensity(G$p$x, G$p$y, w)
  res <- list()
  for (v in VARS) {
    t1 <- proc.time()[[3]]; V <- s3b_variant(G$terr, G$reg, v); t_build <- proc.time()[[3]] - t1
    m <- s3b_measure(V, P, Preg, pts = if (v != "V0g") pts else NULL, full = v != "V0g"); m$time["build"] <- t_build
    d <- s$df; d$area <- m$area; d$interior <- !d$clipped | (in_frame & d$area / d$area_tile > 0.999)
    m$summary <- r3_summary(d, li$lambda_loc); m$interior_flips <- sum(d$interior != s$df$interior)
    dA <- abs(m$area - s$df$area); gk <- G_CHECK[[v]]
    m$c5 <- c(n_over = sum(dA > gk * P), max_ratio = if (gk > 0) max(dA / (gk * P)) else NA, max_abs = max(dA))
    big <- s$df$area >= 1; rel <- dA[big] / s$df$area[big]
    m$b1 <- c(n_cells = sum(big), max_rel = max(rel), n_over = sum(rel > 1e-6))
    if (rp <= 2) m$cells <- d[, c("x", "y", "area", "interior", "nsides", "ecc_T", "theta_T")]
    m$info <- V$info; m$area <- NULL; res[[v]] <- m
  }
  saveRDS(list(archetype = A, roi_id = roi, model = model, rep = rp, n = length(G$terr), A = w$area_um2, Preg = Preg, sumP = sum(P),
               reg = reg_features(G$reg), c0 = c0, old_c123 = old$c123, old_summary = old$summary, t_tess = t_tess, res = res), f)
  NULL
}
if (STAGE == "null") {
  jobs <- expand.grid(i = seq_len(nrow(rois)), model = MODELS, rp = seq_len(NREP), stringsAsFactors = FALSE)
  if (nzchar(Sys.getenv("S13B_ONLY"))) { k <- strsplit(Sys.getenv("S13B_ONLY"), ",")[[1]]; jobs <- jobs[paste(rois$archetype[jobs$i], rois$roi_id[jobs$i], jobs$model, jobs$rp) %in% k, ] }
  jobs <- jobs[order(-jobs$i), ]
  par_run(jobs, null_one)
}

# ---- selezioni (fissate dal seed) -------------------------------------------------------------------------------------
ns13 <- read.csv("results/S1.3/S1.3_null_summary.csv"); ns13$fail <- ns13$c2_overlap > TOL | ns13$c2_symdiff > TOL
sel_grid <- function() { set.seed(SEL_SEED); ok <- ns13[!ns13$fail, ]
  pick <- do.call(rbind, lapply(split(ok, ok$archetype), function(d) d[sample.int(nrow(d), 2), ]))
  rest <- ok[!paste(ok$archetype, ok$roi_id, ok$model, ok$rep) %in% paste(pick$archetype, pick$roi_id, pick$model, pick$rep), ]
  pick <- rbind(pick, rest[sample.int(nrow(rest), 3), ])
  rbind(cbind(ns13[ns13$fail, c("archetype", "roi_id", "model", "rep")], set = "fail"), cbind(pick[, c("archetype", "roi_id", "model", "rep")], set = "pass")) }
sel_cp2 <- function() { set.seed(SEL_SEED + 1); ok <- ns13[!ns13$fail, ]
  do.call(rbind, lapply(split(ok, ok$archetype), function(d) d[sample.int(nrow(d), 10), c("archetype", "roi_id", "model", "rep")])) }
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)

# ---- stadio grid: C-6 ---------------------------------------------------------------------------------------------
grid_pieces <- function(V, gridV) {
  if (V$kind == "sf") {
    pc <- st_intersection(V$terr, gridV); ix <- attr(pc, "idx")
    list(i = ix[, 1], j = ix[, 2], area = as.numeric(st_area(pc)), perim = s3b_perim(pc))
  } else {
    im <- geos_intersects_matrix(V$terr, gridV); ii <- rep(seq_along(im), lengths(im)); jj <- unlist(im)
    pc <- geos_intersection_prec(V$terr[ii], gridV[jj], V$g)
    list(i = ii, j = jj, area = geos_area(pc), perim = geos_length(pc))
  }
}
grid_one <- function(j) {
  f <- file.path(OUTB, "grid", fname(j$archetype, j$roi_id, j$model, j$rep)); if (file.exists(f)) return(NULL)
  G <- regen(roi_i(j$archetype, j$roi_id), j$model, j$rep)
  bb <- st_bbox(G$reg); x0 <- floor(bb[["xmin"]] / 2) * 2; y0 <- floor(bb[["ymin"]] / 2) * 2
  nx <- ceiling((bb[["xmax"]] - x0) / 2); ny <- ceiling((bb[["ymax"]] - y0) / 2)
  grid <- st_make_grid(cellsize = 2, offset = c(x0, y0), n = c(nx, ny), what = "polygons",
                       x = st_as_sfc(st_bbox(c(xmin = x0, ymin = y0, xmax = x0 + 2 * nx, ymax = y0 + 2 * ny))))
  out <- list(set = j$set, archetype = j$archetype, roi_id = j$roi_id, model = j$model, rep = j$rep, n = length(G$terr), n_bins = length(grid))
  for (v in c("V0", "V1", "V2", "V3")) {
    t0 <- proc.time()[[3]]; V <- s3b_variant(G$terr, G$reg, v); g <- if (v == "V3") 1e-7 else 1e-6
    gridV <- if (V$kind == "sf") grid else as_geos_geometry(grid)
    pc <- grid_pieces(V, gridV); aV <- s3b_area(V, V$terr)
    sc <- tapply(pc$area, factor(as.integer(pc$i), levels = seq_along(aV)), sum); sc[is.na(sc)] <- 0
    pp <- tapply(pc$perim, factor(as.integer(pc$i), levels = seq_along(aV)), sum); pp[is.na(pp)] <- 0
    dc <- abs(as.numeric(sc) - aV); tc <- g * as.numeric(pp)
    rs <- s3b_sf(V, V$reg)                    # predicati con la regione come primo argomento (preparata da sf)
    within <- logical(length(grid)); within[unlist(st_contains(rs, grid))] <- TRUE
    touch <- logical(length(grid)); touch[unlist(st_intersects(rs, grid))] <- TRUE
    truth <- numeric(length(grid)); tperim <- numeric(length(grid))
    truth[within] <- as.numeric(st_area(grid[within])); tperim[within] <- 8
    bnd <- which(touch & !within)
    if (length(bnd)) {
      if (V$kind == "sf") { bi <- st_intersection(grid[bnd], rs); ixb <- attr(bi, "idx")[, 1]
        truth[bnd[ixb]] <- as.numeric(st_area(bi)); tperim[bnd[ixb]] <- s3b_perim(bi)
      } else { bi <- geos_intersection_prec(gridV[bnd], V$reg, V$g); truth[bnd] <- geos_area(bi); tperim[bnd] <- geos_length(bi) }
    }
    sb <- tapply(pc$area, factor(as.integer(pc$j), levels = seq_along(grid)), sum); sb[is.na(sb)] <- 0
    pb <- tapply(pc$perim, factor(as.integer(pc$j), levels = seq_along(grid)), sum); pb[is.na(pb)] <- 0
    db <- abs(as.numeric(sb) - truth); tb <- g * (as.numeric(pb) + tperim)
    out[[v]] <- c(n_pieces = length(pc$area), cell_viol = sum(dc > tc), cell_max = max(dc), cell_sum = sum(dc),
                  bin_viol = sum(db > tb), bin_max = max(db), bin_sum = sum(db), n_bnd = length(bnd), n_within = sum(within),
                  elapsed_s = proc.time()[[3]] - t0)
  }
  saveRDS(out, f); NULL
}
if (STAGE == "grid") { jobs <- sel_grid(); write.csv(jobs, file.path(OUTB, "grid", "selection.csv"), row.names = FALSE)
  if (nzchar(Sys.getenv("S13B_ONLY"))) { k <- strsplit(Sys.getenv("S13B_ONLY"), ",")[[1]]; jobs <- jobs[paste(jobs$archetype, jobs$roi_id, jobs$model, jobs$rep) %in% k, ] }
  par_run(jobs, grid_one, cores = min(NC, 30L)) }

# ---- stadio cp1: D-0, CP-1a, CP-1b, CP-1c ------------------------------------------------------------------------------
nonshared_vertices <- function(gk, other, tol) {   # vertici di gk a più di tol dal bordo di other
  v <- s3b_geoms_df(list(gk))
  d <- as.numeric(st_distance(st_sfc(lapply(seq_len(nrow(v)), function(k) st_point(c(v$x[k], v$y[k])))), st_boundary(st_sfc(other))))
  list(v = v, far = d > tol)
}
if (STAGE == "cp1") {
  f <- file.path(OUTB, "cp1", "cp1.rds")
  G <- regen(roi_i("A3", "r3"), "CSR", 3L); t1 <- G$terr[[766]]; t2 <- G$terr[[1987]]
  a_v0 <- sum(as.numeric(st_area(st_intersection(st_sfc(t1), st_sfc(t2)))))
  rel <- st_relate(st_sfc(t1), st_sfc(t2))[1, 1]
  d0 <- s3b_node(list(t1, t2), 1e-6)
  v12 <- s3b_geoms_df(list(t1, t2)); fu <- s3b_fuse(v12$x, v12$y, 1e-6); cr <- v12$g[fu$I] != v12$g[fu$J] & fu$d > 0
  d_reg <- if (any(cr)) as.numeric(st_distance(st_sfc(lapply(fu$I[cr], function(k) st_point(c(v12$x[k], v12$y[k])))), st_boundary(G$reg))) else numeric(0)
  a_1a <- sum(as.numeric(st_area(st_intersection(st_sfc(d0$geoms[[1]]), st_sfc(d0$geoms[[2]])))))
  k_move <- d0$info[["n_moved"]]; amp <- d0$d
  ulp0 <- function(z) 2^(floor(log2(abs(z))) - 52)
  mv0 <- which(v12$g[fu$I] != v12$g[fu$J] & fu$d > 0)      # coppia di vertici diversi: spostamento in ulp della coordinata
  k_ulp <- if (length(mv0)) max(round(abs(v12$x[fu$I[mv0]] - v12$x[fu$J[mv0]]) / ulp0(v12$x[fu$I[mv0]])), round(abs(v12$y[fu$I[mv0]] - v12$y[fu$J[mv0]]) / ulp0(v12$y[fu$I[mv0]]))) else 0
  k_ins1 <- sum(d0$tj_g == 1); k_ins2 <- sum(d0$tj_g == 2)
  ns2 <- nonshared_vertices(t2, t1, 1e-6); ns1 <- nonshared_vertices(t1, t2, 1e-6)
  ins <- function(df, far, k) { n <- nrow(df); nx <- c(2:n, 1); ok <- which(far & far[nx]); if (!k || !length(ok)) return(df[0, ])
    e <- ok[sample.int(length(ok), min(k, length(ok)))]; t <- runif(length(e), 0.1, 0.9)
    data.frame(g = df$g[e], p = df$p[e], r = df$r[e], i = df$i[e] + t, x = df$x[e] + t * (df$x[nx[e]] - df$x[e]), y = df$y[e] + t * (df$y[nx[e]] - df$y[e])) }
  placebo <- do.call(rbind, lapply(1:20, function(sd) {
    set.seed(SEL_SEED + 100 + sd)
    vv <- ns2$v; cand <- which(ns2$far); km <- min(k_move, length(cand)); mv <- if (km) cand[sample.int(length(cand), km)] else integer(0)
    # CORREZIONE dopo la revisione (RC-01): la versione pre-run sommava 8.9e-16 µm a coordinate ~1e3 (ulp 1.1e-13) -> nessuna coordinata
    # cambiava (placebo nullo). Ora: stesso numero di ulp dello spostamento di CP-1a, nella coordinata e col segno casuali; si conta n_changed.
    ulp <- function(z) 2^(floor(log2(abs(z))) - 52)
    if (km) { cc <- sample(c("x", "y"), km, TRUE); sg <- sample(c(-1, 1), km, TRUE)
      for (q in seq_len(km)) { vv[[cc[q]]][mv[q]] <- vv[[cc[q]]][mv[q]] + sg[q] * k_ulp * ulp(vv[[cc[q]]][mv[q]]) } }
    n_changed <- sum(vv$x != ns2$v$x | vv$y != ns2$v$y)
    n2 <- rbind(vv, ins(ns2$v, ns2$far, k_ins2)); n1 <- rbind(ns1$v, ins(ns1$v, ns1$far, k_ins1))
    g2p <- s3b_rebuild(n2, class(t2)[2])$geoms[[1]]; g1p <- s3b_rebuild(n1, class(t1)[2])$geoms[[1]]
    data.frame(seed = sd, moved = km, k_ulp = k_ulp, n_changed = n_changed, ins1 = k_ins1, ins2 = k_ins2, area = sum(as.numeric(st_area(st_intersection(st_sfc(g1p), st_sfc(g2p))))),
               valid = all(st_is_valid(st_sfc(g1p, g2p))))
  }))
  case1 <- list(a_v0 = a_v0, relate = rel, d0_info = d0$info, d0_d = d0$d, d0_dtj = d0$d_tj, d0_tj_g = d0$tj_g, diff_dist_to_region = d_reg,
                a_cp1a = a_1a, cp1a_valid = all(st_is_valid(st_sfc(d0$geoms))), placebo = placebo, n_far2 = sum(ns2$far), n_far1 = sum(ns1$far),
                clipped = G$s$df$clipped[c(766, 1987)], types = c(class(t1)[2], class(t2)[2]))
  H <- regen(roi_i("A5", "r5"), "RSA", 15L); A <- H$w$area_um2; aH <- as.numeric(st_area(H$terr))
  u0 <- st_union(H$terr); un0 <- sum(aH) - as.numeric(st_area(u0)); sd0 <- as.numeric(st_area(st_sym_difference(u0, H$reg)))
  d0g <- s3b_node(c(list(H$reg[[1]]), s3b_sfg_list(H$terr)), 1e-6)
  vH <- s3b_geoms_df(c(list(H$reg[[1]]), s3b_sfg_list(H$terr))); fuH <- s3b_fuse(vH$x, vH$y, 1e-6)
  crH <- vH$g[fuH$I] != vH$g[fuH$J] & fuH$d > 0
  on_reg <- if (any(crH)) as.numeric(st_distance(st_sfc(lapply(fuH$I[crH], function(k) st_point(c(vH$x[k], vH$y[k])))), st_boundary(H$reg))) <= 1e-6 else logical(0)
  V3 <- s3b_variant(H$terr, H$reg, "V3"); a3 <- as.numeric(st_area(V3$terr)); u3 <- st_union(V3$terr)
  case2 <- list(A = A, union_deficit_v0 = un0, symd_v0 = sd0, d0_info = d0g$info, d0_d = d0g$d, d0_dtj = d0g$d_tj, diff_on_region_boundary = on_reg,
                cp1c_rel_ov = abs(sum(a3) - as.numeric(st_area(u3))) / A, cp1c_rel_symd = as.numeric(st_area(st_sym_difference(u3, V3$reg))) / A,
                v3_info = V3$info)
  saveRDS(list(case1 = case1, case2 = case2), f)
  print(case1[c("a_v0", "relate", "a_cp1a", "clipped", "types")]); print(case1$d0_info); print(summary(placebo$area))
  print(case2[c("union_deficit_v0", "symd_v0", "cp1c_rel_ov", "cp1c_rel_symd")]); print(case2$d0_info); print(table(on_reg))
}

# ---- stadio cp2: stress ---------------------------------------------------------------------------------------------
cp2_one <- function(j) {
  f <- file.path(OUTB, "cp2", fname(j$archetype, j$roi_id, j$model, j$rep)); if (file.exists(f)) return(NULL)
  i <- roi_i(j$archetype, j$roi_id); G <- regen(i, j$model, j$rep); P <- s3b_perim(G$terr); Preg <- s3b_perim(G$reg)
  out <- list(archetype = j$archetype, roi_id = j$roi_id, model = j$model, rep = j$rep, A = G$w$area_um2, Preg = Preg, sumP = sum(P), rows = list())
  for (dl in c(1e-12, 1e-10, 1e-8)) {
    tp <- s3b_perturb(G$terr, dl, SEL_SEED + 7 * job_index(i, j$model, j$rep) + round(-log10(dl)))
    n_inv <- sum(!st_is_valid(tp))
    for (v in c("V0", "V1", "V2", "V3")) {   # eccezione GEOS registrata (post hoc, dichiarato: conta come FAIL della variante)
      m <- tryCatch({ V <- s3b_variant(tp, G$reg, v); s3b_measure(V, P, Preg, pts = NULL, full = FALSE) }, error = function(e) conditionMessage(e))
      out$rows[[length(out$rows) + 1]] <- if (is.character(m))
        data.frame(delta = dl, variant = v, n_invalid_input = n_inv, ov = NA, symd = NA, tau = NA, n_pairs = NA, n_pairs_over = NA, pair_max = NA, error = substr(m, 1, 200))
      else data.frame(delta = dl, variant = v, n_invalid_input = n_inv, ov = m$ov, symd = m$symd, tau = m$tau,
                      n_pairs = m$n_pairs, n_pairs_over = m$n_pairs_over, pair_max = if (m$n_pairs) max(m$pairs$area) else 0, error = NA_character_)
    }
  }
  out$rows <- do.call(rbind, out$rows); saveRDS(out, f); NULL
}
if (STAGE == "cp2") { jobs <- sel_cp2(); write.csv(jobs, file.path(OUTB, "cp2", "selection.csv"), row.names = FALSE)
  if (nzchar(Sys.getenv("S13B_ONLY"))) { k <- strsplit(Sys.getenv("S13B_ONLY"), ",")[[1]]; jobs <- jobs[paste(jobs$archetype, jobs$roi_id, jobs$model, jobs$rep) %in% k, ] }
  par_run(jobs, cp2_one) }

# ---- stadio mutants: M1, M2, M3 -----------------------------------------------------------------------------------
if (STAGE == "mutants") {
  rows <- list()
  G <- regen(roi_i("A1", "r1"), "CSR", 1L); P <- s3b_perim(G$terr)
  V <- s3b_variant(G$terr, G$reg, "V1", g = 1e-3); a <- s3b_area(V, V$terr); a0 <- as.numeric(st_area(G$terr))
  rows[[1]] <- data.frame(mutant = "M1", case = "A1 r1 CSR 1", target = "C-5", n_fail = sum(abs(a - a0) > 1e-6 * P))
  for (cs in list(list("A3", "r3", "CSR", 3L), list("A5", "r5", "RSA", 15L))) {
    H <- regen(roi_i(cs[[1]], cs[[2]]), cs[[3]], cs[[4]]); P <- s3b_perim(H$terr); Preg <- s3b_perim(H$reg); lab <- paste(cs, collapse = " ")
    V <- s3b_variant(H$terr, H$reg, "V3", tol = 0); m <- s3b_measure(V, P, Preg, full = FALSE, gchk = G_CHECK[["V3"]])
    rows[[length(rows) + 1]] <- data.frame(mutant = "M2", case = lab, target = "C-1/C-2", n_fail = (abs(m$ov) > m$tau | m$symd > m$tau) + m$n_pairs_over)
    V <- s3b_variant(H$terr, H$reg, "V1", round_terr = FALSE); m <- s3b_measure(V, P, Preg, full = FALSE, gchk = G_CHECK[["V1"]])
    rows[[length(rows) + 1]] <- data.frame(mutant = "M3", case = lab, target = "C-1/C-2", n_fail = (abs(m$ov) > m$tau | m$symd > m$tau) + m$n_pairs_over)
  }
  out <- do.call(rbind, rows); saveRDS(out, file.path(OUTB, "mutants", "mutants.rds")); print(out)
}
