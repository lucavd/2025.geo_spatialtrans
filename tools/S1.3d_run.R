# tools/S1.3d_run.R — S1.3d (BL-084): stadi del banco di prova degli overlay.
# Pre-registrazione results/S1.3d/S1.3d_preregistration.md (ee6208b); definizioni operative results/S1.3d/S1.3d_dev_notes.md.
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.3d_run.R <stage> <set> [outdir]
#   stage = md5 | export | cands ; set = N | K | O2 | H | CP1 | U | S | <lista di chiavi separate da virgola>
#   outdir default /mnt/micron/geo_spatialtrans/S1.3d ; stadi rilanciabili (saltano i file già scritti).
# Ordine: export -> arbitro (tools/S1.3d_arbiter_batch.py) -> cands (candidati + confronto con la verità) -> coverage (python).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
REAL_ARGS <- commandArgs(TRUE)
source("tools/S1.3d_common.R"); source("tools/S1.3d_candidates.R"); source("tools/S1.3d_f4.R"); source("tools/S1.3d_compare.R")
suppressPackageStartupMessages(library(jsonlite))
STAGE <- REAL_ARGS[1]; SET <- REAL_ARGS[2]; OUTD <- if (length(REAL_ARGS) > 2) REAL_ARGS[3] else "/mnt/micron/geo_spatialtrans/S1.3d"
NC <- as.integer(Sys.getenv("S13D_CORES", "40"))
GEOM <- file.path(OUTD, "geom"); TRUTH <- file.path(OUTD, "truth"); RES <- file.path(OUTD, "res")
for (d in c(GEOM, TRUTH, RES)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
ns13 <- read.csv("results/S1.3/S1.3_null_summary.csv"); ns13$fail <- ns13$c2_overlap > TOL | ns13$c2_symdiff > TOL
key_of <- function(A, roi, model, rp) sprintf("%s_%s_%s_%02d", A, roi, model, rp)

# ---- insiemi (fissati dal seed) ------------------------------------------------------------------------------------
set_N <- function() { j <- expand.grid(i = seq_len(nrow(rois)), model = MODELS, rp = seq_len(NREP), stringsAsFactors = FALSE)
  j$A <- rois$archetype[j$i]; j$roi <- rois$roi_id[j$i]; j$key <- key_of(j$A, j$roi, j$model, j$rp); j }
set_K <- function() {                       # nulli che contengono i 26 casi noti (15 unioni + 10 coppie P1 + 71/175)
  p1 <- read.csv("results/S1.3b/S1.3b_posthoc_P1_pairs_v0.csv")
  k <- unique(c(with(ns13[ns13$fail, ], key_of(archetype, roi_id, model, rep)), with(p1, key_of(archetype, roi_id, model, rep)), key_of("A5", "r5", "RSA", 15)))
  N <- set_N(); N[N$key %in% k, ] }
set_O2 <- function() { N <- set_N(); K <- set_K(); set.seed(SEL_SEED + 3)
  rest <- N[!N$key %in% K$key, ]
  add <- do.call(rbind, lapply(split(rest, rest$A), function(d) { nk <- sum(K$A == d$A[1]); d[sample.int(nrow(d), max(0, 20 - nk)), ] }))
  rbind(K, add) }
set_H <- function() { j <- rbind(expand.grid(i = seq_len(nrow(rois)), model = MODELS, rp = 61:63, stringsAsFactors = FALSE),
                                data.frame(i = seq_len(nrow(rois)), model = "CSR", rp = 64L))
  j$A <- rois$archetype[j$i]; j$roi <- rois$roi_id[j$i]; j$key <- key_of(j$A, j$roi, j$model, j$rp); j }
get_set <- function(S) switch(S, N = set_N(), K = set_K(), O2 = set_O2(), H = set_H(),
  { N <- rbind(set_N(), set_H()); k <- strsplit(S, ",")[[1]]; N[N$key %in% k, ] })

par_run <- function(jobs, FUN, cores = NC) {
  res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) tryCatch(FUN(j), error = function(e) paste("ERRORE:", conditionMessage(e))),
                  mc.cores = cores, mc.preschedule = FALSE)
  err <- vapply(res, function(r) is.character(r) && length(r) == 1 && startsWith(r, "ERRORE"), TRUE)
  for (k in which(err)) cat(sprintf("[%s] job %s: %s\n", STAGE, jobs$key[k], res[[k]]))
  cat(sprintf("[%s %s] errori %d su %d\n", STAGE, SET, sum(err), length(res))); if (any(err)) quit(status = 1)
}

# ---- export: geometria bit-esatta + metadati -------------------------------------------------------------------------
write_case <- function(key, terr, reg, pts, extra = list()) {
  n <- length(terr)
  write_parquet(geom_to_df(terr), file.path(GEOM, paste0(key, "_cells.parquet")))
  write_parquet(geom_to_df(reg, ids = 0L), file.path(GEOM, paste0(key, "_region.parquet")))
  pr <- st_intersects(terr, terr)
  pairs <- do.call(rbind, lapply(seq_along(pr), function(q) { j <- pr[[q]]; j <- j[j > q]; if (length(j)) cbind(q, j) }))
  if (is.null(pairs)) pairs <- matrix(integer(0), 0, 2)
  write_parquet(data.frame(i = as.integer(pairs[, 1]), j = as.integer(pairs[, 2])), file.path(GEOM, paste0(key, "_pairs.parquet")))
  if (!is.null(pts)) write_parquet(data.frame(x = pts$x, y = pts$y), file.path(GEOM, paste0(key, "_pts.parquet")))
  bb <- st_bbox(c(st_geometry(reg), terr))
  meta <- c(list(key = key, n = n, A = as.numeric(st_area(reg)), Preg = s3b_perim(reg), P = s3b_perim(terr), gs = grid_spec(st_bbox(terr)),
                 max_abs_over_g = max(abs(as.numeric(bb))) / G_GRID, has_pts = !is.null(pts)), extra)
  write_json(meta, file.path(GEOM, paste0(key, "_meta.json")), auto_unbox = TRUE, digits = NA)
  # controllo bit-esatto del giro parquet
  back <- df_to_sfc(as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_cells.parquet")))), n)
  identical(st_as_binary(back), st_as_binary(terr))
}
export_null <- function(j) {
  if (file.exists(file.path(GEOM, paste0(j$key, "_meta.json")))) return(NULL)
  G <- regen(j$i, j$model, j$rp); c0 <- c123(G$s$tv, G$w, G$p$x, G$p$y)
  # tabella per cella di S1.3 (per B-1: verdetti B ricalcolati con le aree dei candidati)
  dfc <- G$s$df; dfc$in_frame <- lengths(st_within(G$s$tv$tiles, G$w$frame)) > 0; dfc$lambda_loc <- r3_local_intensity(G$p$x, G$p$y, G$w)$lambda_loc
  write_parquet(dfc, file.path(GEOM, paste0(j$key, "_df.parquet")))
  ok <- write_case(j$key, G$terr, G$reg, G$p, extra = list(c0 = as.list(c0[c("c1", "c2_overlap", "c2_symdiff")]), clipped = G$s$df$clipped,
                                                            n_lost = G$s$df$n_lost, n_gained = G$s$df$n_gained,
                                                            archetype = j$A, roi = j$roi, model = j$model, rep = j$rp))
  if (!ok) stop("giro parquet non bit-esatto")
  NULL
}

# ---- cands: candidati + confronto ------------------------------------------------------------------------------------
read_case <- function(key) {
  meta <- fromJSON(file.path(GEOM, paste0(key, "_meta.json")), simplifyVector = TRUE)
  terr <- df_to_sfc(as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_cells.parquet")))), meta$n)
  reg <- df_to_sfc(transform(as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_region.parquet")))), cell = 1L), 1L)
  pairs <- as.matrix(as.data.frame(read_parquet(file.path(GEOM, paste0(key, "_pairs.parquet")))))
  pf <- file.path(GEOM, paste0(key, "_pts.parquet")); pts <- if (file.exists(pf)) as.data.frame(read_parquet(pf)) else NULL
  list(meta = meta, terr = terr, reg = reg, pairs = pairs, pts = pts)
}
cands_one <- function(key, do_o2) {
  f <- file.path(RES, paste0(key, if (do_o2) "_o12" else "_o1", ".rds")); if (file.exists(f)) return(NULL)
  D <- read_case(key); m <- D$meta; P <- m$P; out <- list(key = key, n = m$n, o1 = list(), o2 = list(), viol = list())
  tr <- fromJSON(file.path(TRUTH, paste0(key, "_o1.json")))$o1
  tp <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o1_pairs.parquet"))))
  tc <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o1_cells.parquet"))))
  stopifnot(identical(as.integer(tp$i), as.integer(D$pairs[, 1])), identical(as.integer(tp$j), as.integer(D$pairs[, 2])))
  f4 <- NULL
  if (!is.null(D$pts)) {
    f4 <- tryCatch(f4_build(D$pts$x, D$pts$y, D$reg), error = function(e) conditionMessage(e))
    if (is.list(f4)) {
      write_parquet(geom_to_df(f4$terr), file.path(GEOM, paste0(key, "_f4cells.parquet")))
      a0 <- as.numeric(st_area(D$terr)); a4 <- as.numeric(st_area(f4$terr))
      out$c5 <- list(n_over = sum(abs(a4 - a0) > G_GRID * P), max_ratio = max(abs(a4 - a0) / (G_GRID * P)), info = f4$info[setdiff(names(f4$info), c("lost", "gained"))],
                     same_lost = identical(as.integer(f4$info$lost), as.integer(m$n_lost)), same_gained = identical(as.integer(f4$info$gained), as.integer(m$n_gained)))
    } else out$c5 <- list(error = f4)
  }
  for (cd in CANDS_O1) {
    if (cd == "F4" && !is.list(f4)) { out$o1[[cd]] <- data.frame(cand = cd, applicable = !is.null(D$pts), fail = !is.null(D$pts)); next }
    o <- tryCatch(o1_run(cand_prepare(cd, D$terr, D$reg, f4 = if (is.list(f4)) f4$terr), D$pairs), error = function(e) conditionMessage(e))
    if (is.character(o)) { out$o1[[cd]] <- data.frame(cand = cd, applicable = TRUE, fail = TRUE, n_exceptions = 1, error = substr(o, 1, 200)); next }
    cm <- cmp_o1(o, tr, tp, tc, P, m$Preg); cm$summary$applicable <- TRUE
    dff <- file.path(GEOM, paste0(key, "_df.parquet"))
    if (file.exists(dff)) { d <- as.data.frame(read_parquet(dff)); li <- d$lambda_loc; d$area <- o$area
      d$interior <- !d$clipped | (d$in_frame & d$area / d$area_tile > 0.999)
      out$b1[[cd]] <- list(summary = r3_summary(d[, setdiff(names(d), c("in_frame", "lambda_loc"))], li),
                           cells = if (isTRUE(m$rep <= 2)) d[, c("x", "y", "area", "interior", "nsides", "ecc_T", "theta_T")] else NULL) }
    if (cd == "F2") { cm$summary$n_invalid_rounded <- length(o$info$invalid_idx); cm$summary$still_invalid <- o$info$still_invalid }
    out$o1[[cd]] <- cm$summary; if (!is.null(cm$viol)) out$viol[[paste0("o1_", cd)]] <- cm$viol
  }
  if (do_o2) {
    t2 <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o2_s2.parquet")))); t8 <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o2_s8.parquet"))))
    k8 <- paste(t2$i, t2$ix %/% 4, t2$iy %/% 4); tau3 <- tapply(G_GRID * t2$perim, k8, sum)[paste(t8$i, t8$ix, t8$iy)]
    gs <- m$gs; grids <- list(`2` = grid_sfc(gs, 2), `8` = grid_sfc(gs, 8))
    for (cd in CANDS_O2) {
      if (cd == "F4" && !is.list(f4)) next
      C <- tryCatch(cand_prepare(cd, D$terr, D$reg, f4 = if (is.list(f4)) f4$terr), error = function(e) NULL)
      for (s in c(2, 8)) {
        o <- if (is.null(C)) list(pieces = NULL, time = c(o2 = NA)) else o2_run(C, gs, s, grid = grids[[as.character(s)]])
        cm <- cmp_o2(o$pieces, if (s == 2) t2 else t8, s, cd); cm$summary$t_o2 <- o$time[["o2"]]; cm$summary$op <- "O2"
        out$o2[[paste(cd, s)]] <- cm$summary; if (!is.null(cm$viol)) out$viol[[paste0("o2_", cd, "_", s)]] <- cm$viol
        if (s == 2) { cm3 <- cmp_o2(o3_aggregate(o$pieces), t8, 8, cd, tau_truth = as.numeric(tau3)); cm3$summary$op <- "O3"; cm3$summary$t_o2 <- NA
          out$o2[[paste(cd, "O3")]] <- cm3$summary; if (!is.null(cm3$viol)) out$viol[[paste0("o3_", cd)]] <- cm3$viol }
      }
    }
  }
  out$o1 <- rbind_fill(out$o1)
  if (length(out$o2)) out$o2 <- do.call(rbind, out$o2)
  saveRDS(out, f); NULL
}
rbind_fill <- function(L) { L <- Filter(Negate(is.null), L); nm <- unique(unlist(lapply(L, names)))
  do.call(rbind, lapply(L, function(d) { for (k in setdiff(nm, names(d))) d[[k]] <- NA; d[, nm, drop = FALSE] })) }

if (STAGE == "md5") { st <- system("md5sum -c --quiet results/R3/R3_inputs_md5.txt"); cat("md5 exit", st, "\n"); quit(status = st) }
if (STAGE == "list") { J <- get_set(SET); cat(J$key, sep = "\n"); quit(status = 0) }
if (STAGE == "export") { J <- get_set(SET); J <- J[order(-J$i), ]; par_run(J, export_null) }
if (STAGE == "cands") { J <- get_set(SET); o2keys <- set_O2()$key
  do2 <- Sys.getenv("S13D_O2", "auto"); J$do_o2 <- if (do2 == "auto") J$key %in% o2keys else do2 == "1"
  par_run(J, function(j) cands_one(j$key, j$do_o2)) }

# ---- export degli insiemi costruiti (U, CP1, S) ---------------------------------------------------------------------
export_built <- function(cs) {
  if (file.exists(file.path(GEOM, paste0(cs$key, "_meta.json")))) return(NULL)
  ok <- write_case(cs$key, cs$terr, cs$reg, cs$pts, extra = c(list(note = cs$note, changed = cs$changed, n_invalid_input = sum(!st_is_valid(cs$terr)),
                                                                 types = paste(unique(as.character(st_geometry_type(cs$terr))), collapse = ",")), cs$extra))
  if (!ok) stop("giro parquet non bit-esatto: ", cs$key)
  NULL
}
u5_tessellate <- function(u) {
  cen <- data.frame(cell_id = seq_len(nrow(u$pts)), region_id = 1L, x = u$pts$x, y = u$pts$y)
  rg <- list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = as.numeric(st_area(st_sfc(u$reg)))), region_polygons = st_sfc(u$reg))
  tv <- tryCatch(tessellate_voronoi(cen, rg, verbose = FALSE), error = function(e) conditionMessage(e))
  if (is.character(tv)) { write_json(list(key = u$key, error = tv), file.path(GEOM, paste0(u$key, "_FAILED.json")), auto_unbox = TRUE); return(NULL) }
  mk(u$key, tv$cell_territories, u$reg, pts = u$pts, note = u$note, extra = list(n_lost = tv$territory_df$n_pieces_lost, n_gained = tv$territory_df$n_pieces_gained))
}
built_keys_file <- function(S) file.path(OUTD, paste0("keys_", S, ".txt"))
if (STAGE == "export_built") {
  source("tools/S1.3d_cases.R")
  if (SET == "U") {
    L <- c(u1_cases(), u3_cases(), u4_cases(), Filter(Negate(is.null), lapply(u5_point_sets(), u5_tessellate)), u2_cases())
    res <- mclapply(L, function(cs) tryCatch(export_built(cs), error = function(e) paste("ERRORE", cs$key, conditionMessage(e))), mc.cores = min(NC, 20))
    err <- Filter(is.character, res); for (e in err) cat(e, "\n")
    writeLines(vapply(L, `[[`, "", "key"), built_keys_file("U")); cat("U: casi", length(L), "errori", length(err), "\n"); if (length(err)) quit(status = 1)
  }
  if (SET == "CP1") {
    sel <- cp1_selection(ns13)
    res <- mclapply(seq_len(nrow(sel)), function(q) tryCatch({ L <- cp1_cases_one(sel$archetype[q], sel$roi_id[q], sel$model[q], sel$rep[q]); for (cs in L) export_built(cs); vapply(L, `[[`, "", "key") },
                    error = function(e) paste("ERRORE", q, conditionMessage(e))), mc.cores = NC, mc.preschedule = FALSE)
    err <- Filter(function(r) length(r) == 1 && startsWith(r, "ERRORE"), res); for (e in err) cat(e, "\n")
    writeLines(unlist(res[!vapply(res, function(r) length(r) == 1 && startsWith(r, "ERRORE"), TRUE)]), built_keys_file("CP1")); cat("CP1: errori", length(err), "\n"); if (length(err)) quit(status = 1)
  }
  if (SET == "S") { cs <- s_case(); export_built(cs); writeLines(cs$key, built_keys_file("S")); cat("S esportato\n") }
}
if (STAGE == "cands_built") {          # candidati sugli insiemi costruiti; O2 per U3, U4, S (vedi dev_notes)
  keys <- readLines(built_keys_file(SET)); keys <- keys[file.exists(file.path(GEOM, paste0(keys, "_meta.json")))]
  J <- data.frame(key = keys, do_o2 = grepl("^(U1|U3|U4|U5|CP1)", keys), stringsAsFactors = FALSE)
  par_run(J, function(j) cands_one(j$key, j$do_o2), cores = if (SET == "S") 1L else NC)
}

# ---- mutanti C-3 (errore pianificato ~2 tau) su un nullo di N con verità esatta -------------------------------------
if (STAGE == "mutants") {
  key <- SET; D <- read_case(key); m <- D$meta; P <- m$P
  tr <- fromJSON(file.path(TRUTH, paste0(key, "_o1.json")))$o1
  tp <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o1_pairs.parquet")))); tc <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o1_cells.parquet"))))
  t2 <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o2_s2.parquet")))); t8 <- as.data.frame(read_parquet(file.path(TRUTH, paste0(key, "_o2_s8.parquet"))))
  rows <- list(); tauU <- G_GRID * (m$Preg + sum(P))
  # M-O1: la cella più grande ristretta attorno al baricentro fino a perdere 2 tau_U
  k <- which.max(tc$a); g <- D$terr[[k]]; ctr <- st_coordinates(st_centroid(st_sfc(g)))[1, 1:2]; f <- sqrt(1 - 2 * tauU / tc$a[k])
  Tm <- D$terr; Tm[[k]] <- (g - ctr) * f + ctr
  o <- o1_run(cand_prepare("F1", Tm, D$reg), D$pairs); s <- cmp_o1(o, tr, tp, tc, P, m$Preg)$summary
  rows[[1]] <- data.frame(mutant = "M-O1", planted = 2 * tauU, detected = s$union_viol, detail = sprintf("union_err/tauU = %.3f", s$union_err / tauU))
  # M-O2a: tolto il pezzo più grande dell'uscita F1 a 2 µm
  C <- cand_prepare("F1", D$terr, D$reg); o2 <- o2_run(C, m$gs, 2); pc <- o2$pieces; q <- which.max(pc$area)
  s <- cmp_o2(pc[-q, ], t2, 2, "F1")$summary; rows[[2]] <- data.frame(mutant = "M-O2a", planted = pc$area[q], detected = s$n_missing >= 1, detail = sprintf("n_missing %d", s$n_missing))
  # M-O2b: griglia spostata di 1e-3 µm (stessi indici)
  gsx <- m$gs; gsx$X0 <- gsx$X0 + 1e-3; s <- cmp_o2(o2_run(C, gsx, 2)$pieces, t2, 2, "F1")$summary
  rows[[3]] <- data.frame(mutant = "M-O2b", planted = 1e-3, detected = s$n_viol > 0, detail = sprintf("n_viol %d su %d", s$n_viol, s$n_truth))
  # M-O3: annidamento sbagliato ((ix + 1) %/% 4)
  pw <- pc; pw$ix <- pw$ix + 1L; a3 <- o3_aggregate(pw); a3$ix <- a3$ix
  k8 <- paste(t2$i, t2$ix %/% 4, t2$iy %/% 4); tau3 <- tapply(G_GRID * t2$perim, k8, sum)[paste(t8$i, t8$ix, t8$iy)]
  s <- cmp_o2(a3, t8, 8, "F1", tau_truth = as.numeric(tau3))$summary; rows[[4]] <- data.frame(mutant = "M-O3", planted = NA, detected = s$n_viol > 0, detail = sprintf("n_viol %d", s$n_viol))
  out <- do.call(rbind, rows); out$key <- key; print(out); saveRDS(out, file.path(RES, "mutants.rds"))
}
