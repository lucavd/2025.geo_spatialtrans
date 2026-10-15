# results/S1.3d/review/recalc/recalc_independent.R — revisione «ricalcolo» S1.3d: ricalcolo indipendente dei confronti con la verità.
# Codice proprio (non sorgente tools/S1.3d_*.R): lettore parquet, V0 (sf, GEOS di sf) e F1 (geos *_prec, g = 1e-6) riscritti dalla definizione
# della pre-registrazione; O2 con enumerazione dei bin dal riquadro della cella (non dal predicato), chiavi intere, tolleranze dalla pre-registrazione:
#   tau_U = g (P_reg + sum P); tau_ij = g (P_i + P_j); tau_i = g P_i; tau_kb = g * perimetro esatto del pezzo; spurio: g * perimetro del pezzo candidato;
#   O3: tau = somma dei tau dei sotto-pezzi esatti da 2 µm.
# Uso: Rscript --vanilla recalc_independent.R <modo> ; modo = nulls (O1 V0+F1 su N ∪ H; O2/O3 F1 su O2 ∪ H) | S_F1 | S_V0
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow); library(jsonlite); library(parallel) })
D <- "/mnt/micron/geo_spatialtrans/S1.3d"; HEAVY <- file.path(D, "review_recalc"); OUTI <- file.path(HEAVY, "indep"); dir.create(OUTI, showWarnings = FALSE)
g <- 1e-6; MODE <- commandArgs(TRUE)[1]; if (is.na(MODE)) MODE <- "none"; NC <- as.integer(Sys.getenv("RC_CORES", "26"))
read_poly <- function(f, n) { t <- as.data.frame(read_parquet(f)); t <- t[order(t$cell, t$part, t$ring, t$k), ]
  by <- split(seq_len(nrow(t)), t$cell); out <- vector("list", n)
  for (q in seq_len(n)) { ix <- by[[as.character(q)]]
    if (is.null(ix)) { out[[q]] <- st_polygon(); next }
    d <- t[ix, ]; parts <- lapply(split(d, d$part), function(dp) lapply(split(dp, dp$ring), function(dr) { m <- cbind(dr$x, dr$y); rbind(m, m[1, ]) }))
    out[[q]] <- if (length(parts) == 1) st_polygon(parts[[1]]) else st_multipolygon(unname(parts)) }
  st_sfc(out) }
read_reg <- function(k) { t <- as.data.frame(read_parquet(file.path(D, "geom", paste0(k, "_region.parquet")))); t$cell <- 1L; read_poly_df(t) }
read_poly_df <- function(t) { t <- t[order(t$part, t$ring, t$k), ]; parts <- lapply(split(t, t$part), function(dp) lapply(split(dp, dp$ring), function(dr) { m <- cbind(dr$x, dr$y); rbind(m, m[1, ]) }))
  st_sfc(if (length(parts) == 1) st_polygon(parts[[1]]) else st_multipolygon(unname(parts))) }
pq <- function(...) as.data.frame(read_parquet(file.path(...)))
safe <- function(expr) tryCatch(expr, error = function(e) structure(NA_real_, err = conditionMessage(e)))

o1_cmp <- function(cand, a, ua, sd, pij, pji, nexc, m, tr, tp, tc) {
  P <- m$P; tauU <- g * (m$Preg + sum(P)); tij <- g * (P[tp$i] + P[tp$j]); ti <- g * P
  eU <- abs(ua - tr$union); eS <- abs(sd - tr$symd); ac <- a[tc$cell]
  vij <- is.na(pij) | abs(pij - tp$a) > tij; vji <- is.na(pji) | abs(pji - tp$a) > tij; vc <- is.na(ac) | abs(ac - tc$a) > ti[tc$cell]
  uv <- is.na(eU) || eU > tauU; sv <- is.na(eS) || eS > tauU; viol <- uv || sv || any(vij) || any(vji) || any(vc)
  data.frame(cand = cand, union_err = eU, symd_err = eS, tauU = tauU, uv = uv, sv = sv, n_ij = sum(vij), n_ji = sum(vji), n_any = sum(vij | vji), n_cell = sum(vc),
             pmax = max(c(abs(pij - tp$a), abs(pji - tp$a)) / tij, na.rm = TRUE), n_exc = nexc, silent = viol && nexc == 0, fail = viol || nexc > 0,
             viol_pairs = paste(head(paste(tp$i[vij | vji], tp$j[vij | vji], sep = "-"), 50), collapse = " "), n_tc = nrow(tc), tc_cells_ok = setequal(tc$cell, seq_along(P)))
}
run_o1 <- function(k, cands = c("V0", "F1"), check_pairs = TRUE) {
  m <- fromJSON(file.path(D, "geom", paste0(k, "_meta.json"))); t <- read_poly(file.path(D, "geom", paste0(k, "_cells.parquet")), m$n); reg <- read_reg(k)
  pr <- pq(D, "geom", paste0(k, "_pairs.parquet")); tr <- fromJSON(file.path(D, "truth", paste0(k, "_o1.json")))$o1
  tp <- pq(D, "truth", paste0(k, "_o1_pairs.parquet")); tc <- pq(D, "truth", paste0(k, "_o1_cells.parquet"))
  # coppie della verità = coppie della geometria (come insieme), poi allineate
  kp <- paste(pr$i, pr$j); kt <- paste(tp$i, tp$j); stopifnot(setequal(kp, kt), !anyDuplicated(kt)); tp <- tp[match(kp, kt), ]
  pc_chk <- NA
  if (check_pairs) { it <- st_intersects(t, t); q <- do.call(rbind, lapply(seq_along(it), function(a) { b <- it[[a]]; b <- b[b > a]; if (length(b)) cbind(a, b) }))
    pc_chk <- length(setdiff(paste(q[, 1], q[, 2]), kp)) + length(setdiff(kp, paste(q[, 1], q[, 2]))) }
  out <- list()
  if ("V0" %in% cands) { ne <- 0; a <- as.numeric(st_area(t)); U <- safe(st_union(t)); ua <- if (inherits(U, "sfc")) sum(as.numeric(st_area(U))) else { ne <- ne + 1; NA }
    sd <- if (inherits(U, "sfc")) safe(sum(as.numeric(st_area(st_sym_difference(U, reg))))) else NA; if (!is.null(attr(sd, "err"))) ne <- ne + 1
    pa <- function(I, J) vapply(seq_along(I), function(q) tryCatch(sum(as.numeric(st_area(st_intersection(t[I[q]], t[J[q]])))), error = function(e) NA_real_), 0)
    pij <- pa(pr$i, pr$j); pji <- pa(pr$j, pr$i); ne <- ne + sum(is.na(pij)) + sum(is.na(pji))
    out$V0 <- o1_cmp("V0", a, ua, sd, pij, pji, ne, m, tr, tp, tc) }
  if ("F1" %in% cands) { ne <- 0; tg <- geos_set_precision(as_geos_geometry(t), g); rg <- geos_set_precision(as_geos_geometry(reg), g); a <- geos_area(tg)
    U <- tryCatch(geos_unary_union_prec(geos_make_collection(tg), g), error = function(e) NULL); ua <- if (is.null(U)) { ne <- ne + 1; NA } else geos_area(U)
    sd <- if (is.null(U)) NA else tryCatch(geos_area(geos_sym_difference_prec(U, rg, g)), error = function(e) { ne <<- ne + 1; NA })
    pij <- tryCatch(geos_area(geos_intersection_prec(tg[pr$i], tg[pr$j], g)), error = function(e) { ne <<- ne + 1; rep(NA_real_, nrow(pr)) })
    pji <- tryCatch(geos_area(geos_intersection_prec(tg[pr$j], tg[pr$i], g)), error = function(e) { ne <<- ne + 1; rep(NA_real_, nrow(pr)) })
    out$F1 <- o1_cmp("F1", a, ua, sd, pij, pji, ne, m, tr, tp, tc); out$F1$n_invalid_cells <- sum(!geos_is_valid(tg)) }
  r <- do.call(rbind, lapply(out, function(d) { if (is.null(d$n_invalid_cells)) d$n_invalid_cells <- NA; d })); r$key <- k; r$pairs_mismatch <- pc_chk; r
}
# ---- O2/O3 F1 (o V0) con bin dal riquadro --------------------------------------------------------------------------------
pieces <- function(cand, t, m, s) {
  gs <- m$gs; nx <- gs$nx8 * 8 / s; ny <- gs$ny8 * 8 / s
  if (cand == "F1") { tg <- geos_set_precision(as_geos_geometry(t), g); ex <- geos_extent(tg) } else { tg <- t; bb <- lapply(t, st_bbox); ex <- data.frame(xmin = sapply(bb, `[[`, 1), ymin = sapply(bb, `[[`, 2), xmax = sapply(bb, `[[`, 3), ymax = sapply(bb, `[[`, 4)) }
  ok <- is.finite(ex$xmin)
  x0 <- pmax(0, floor((ex$xmin - gs$X0) / s)); x1 <- pmin(nx - 1, floor((ex$xmax - gs$X0) / s)); y0 <- pmax(0, floor((ex$ymin - gs$Y0) / s)); y1 <- pmin(ny - 1, floor((ex$ymax - gs$Y0) / s))
  nb <- ifelse(ok, (x1 - x0 + 1) * (y1 - y0 + 1), 0); ii <- rep(seq_along(nb), nb)
  off <- sequence(nb) - 1L; w <- (x1 - x0 + 1)[ii]; bx <- x0[ii] + off %% w; by <- y0[ii] + off %/% w
  if (cand == "F1") { R <- geos_create_rectangle(gs$X0 + s * bx, gs$Y0 + s * by, gs$X0 + s * (bx + 1), gs$Y0 + s * (by + 1))
    pc <- geos_intersection_prec(tg[ii], R, g); ar <- geos_area(pc); pe <- geos_length(pc)
    keep <- ar > 0; return(data.frame(i = as.integer(ii[keep]), ix = as.integer(bx[keep]), iy = as.integer(by[keep]), area = ar[keep], perim = pe[keep]))
  }
  # V0 come sf: st_intersection(celle, bin necessari) (stessa operazione GEOS di sf); indici dei bin propri (bx, by) dalla tabella dei bin unici
  bk <- bx + by * nx; ub <- sort(unique(bk)); ubx <- ub %% nx; uby <- ub %/% nx
  R <- st_sfc(lapply(seq_along(ub), function(q) { xa <- gs$X0 + s * ubx[q]; ya <- gs$Y0 + s * uby[q]; st_polygon(list(rbind(c(xa, ya), c(xa + s, ya), c(xa + s, ya + s), c(xa, ya + s), c(xa, ya)))) }))
  pc <- st_intersection(t, R); id <- attr(pc, "idx"); ar <- as.numeric(st_area(pc)); pe <- geos_length(as_geos_geometry(pc))
  keep <- ar > 0; data.frame(i = as.integer(id[keep, 1]), ix = as.integer(ubx[id[keep, 2]]), iy = as.integer(uby[id[keep, 2]]), area = ar[keep], perim = pe[keep])
}
o2_cmp <- function(pc, tt, tau_t = NULL, s) {
  kt <- sprintf("%d_%d_%d", as.integer(tt$i), as.integer(tt$ix), as.integer(tt$iy)); kc <- sprintf("%d_%d_%d", pc$i, pc$ix, pc$iy)
  stopifnot(!anyDuplicated(kt), !anyDuplicated(kc))
  ka <- union(kt, kc); at <- tt$area[match(ka, kt)]; ac <- pc$area[match(ka, kc)]
  tau <- if (is.null(tau_t)) g * tt$perim[match(ka, kt)] else tau_t[match(ka, kt)]
  spur <- is.na(at); pcp <- pc$perim[match(ka, kc)]; tau[spur] <- g * ifelse(is.na(pcp[spur]), 4 * s, pcp[spur])
  at[is.na(at)] <- 0; mis <- is.na(ac); ac[is.na(ac)] <- 0; e <- abs(ac - at); v <- e > tau
  data.frame(n_truth = nrow(tt), n_cand = nrow(pc), n_viol = sum(v), n_missing = sum(v & mis), n_spurious = sum(v & spur), max_ratio = max(e / tau), n_missing_any = sum(mis), n_spurious_any = sum(spur))
}
run_o2 <- function(k, cand = "F1") {
  m <- fromJSON(file.path(D, "geom", paste0(k, "_meta.json"))); t <- read_poly(file.path(D, "geom", paste0(k, "_cells.parquet")), m$n)
  t2 <- pq(D, "truth", paste0(k, "_o2_s2.parquet")); t8 <- pq(D, "truth", paste0(k, "_o2_s8.parquet"))
  ne <- 0; p2 <- tryCatch(pieces(cand, t, m, 2), error = function(e) { ne <<- ne + 1; NULL }); p8 <- tryCatch(pieces(cand, t, m, 8), error = function(e) { ne <<- ne + 1; NULL })
  r2 <- if (is.null(p2)) NULL else o2_cmp(p2, t2, NULL, 2); r8 <- if (is.null(p8)) NULL else o2_cmp(p8, t8, NULL, 8)
  # O3: somma dei sotto-pezzi (chiavi intere), tau = somma dei tau esatti dei sotto-pezzi
  k8t <- sprintf("%d_%d_%d", as.integer(t8$i), as.integer(t8$ix), as.integer(t8$iy))
  g2 <- sprintf("%d_%d_%d", as.integer(t2$i), as.integer(t2$ix) %/% 4L, as.integer(t2$iy) %/% 4L)
  tau3 <- tapply(g * t2$perim, g2, sum); s3 <- tapply(t2$area, g2, sum)
  ident <- c(missing_groups = sum(!k8t %in% names(s3)), extra_groups = sum(!names(s3) %in% k8t), max_rel = max(abs(s3[k8t] - t8$area) / pmax(t8$area, 1e-300), na.rm = TRUE))
  r3 <- NULL
  if (!is.null(p2)) { a3 <- tapply(p2$area, sprintf("%d_%d_%d", p2$i, p2$ix %/% 4L, p2$iy %/% 4L), sum); sp <- do.call(rbind, strsplit(names(a3), "_"))
    p3 <- data.frame(i = as.integer(sp[, 1]), ix = as.integer(sp[, 2]), iy = as.integer(sp[, 3]), area = as.numeric(a3), perim = NA_real_)
    r3 <- o2_cmp(p3, t8, as.numeric(tau3[k8t]), 8) }
  rr <- rbind(if (!is.null(r2)) cbind(op = "O2s2", r2), if (!is.null(r8)) cbind(op = "O2s8", r8), if (!is.null(r3)) cbind(op = "O3", r3))
  rr$key <- k; rr$cand <- cand; rr$n_exc <- ne; rr$o3_ident_missing <- ident[["missing_groups"]]; rr$o3_ident_extra <- ident[["extra_groups"]]; rr$o3_ident_maxrel <- ident[["max_rel"]]; rr
}
dojob <- function(k, what) { f <- file.path(OUTI, sprintf("%s_%s.rds", k, what)); if (file.exists(f)) return(NULL)
  r <- switch(what, o1 = run_o1(k), o2F1 = run_o2(k, "F1"), o2V0 = run_o2(k, "V0"), o1F1 = run_o1(k, "F1", check_pairs = FALSE)); saveRDS(r, f); NULL }
if (MODE == "nulls") {
  N <- readLines(file.path(D, "keys_N.txt")); H <- readLines(file.path(D, "keys_H.txt")); O2 <- readLines(file.path(D, "keys_O2.txt"))
  J <- rbind(data.frame(k = c(H, O2), w = "o2F1"), data.frame(k = c(H, N), w = "o1"))
  set.seed(1); J <- J[c(which(J$w == "o2F1"), sample(which(J$w == "o1"))), ]
  res <- mclapply(seq_len(nrow(J)), function(q) tryCatch(dojob(J$k[q], J$w[q]), error = function(e) paste("ERRORE", J$k[q], J$w[q], conditionMessage(e))), mc.cores = NC, mc.preschedule = FALSE)
  e <- Filter(is.character, res); for (x in e) cat(x, "\n"); cat("nulls: errori", length(e), "su", nrow(J), "\n")
}
if (MODE == "S_F1") { t0 <- proc.time()[[3]]; dojob("S_cperf", "o2F1"); cat("S F1 O2/O3", proc.time()[[3]] - t0, "s\n"); t0 <- proc.time()[[3]]; dojob("S_cperf", "o1F1"); cat("S F1 O1", proc.time()[[3]] - t0, "s\n") }
if (MODE == "S_V0") { t0 <- proc.time()[[3]]; dojob("S_cperf", "o2V0"); cat("S V0 O2/O3", proc.time()[[3]] - t0, "s\n") }
cat("FINE", MODE, "\n")
