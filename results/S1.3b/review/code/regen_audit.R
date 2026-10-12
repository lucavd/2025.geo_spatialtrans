
# review/code/regen_audit.R — rigenera i 10 nulli con coppia V0 sbagliata: firma del vertice, CP-1a/placebo per coppia, C-6 V0/V2 indipendente.
REAL_ARGS <- commandArgs(TRUE)
RV <- "/mnt/micron/geo_spatialtrans/S1.3b/review_code"; dir.create(RV, showWarnings = FALSE)
args <- c("none", file.path(RV, "s13out")); commandArgs <- function(trailingOnly = TRUE) args
suppressMessages(source("tools/S1.3_roi.R")); rm(commandArgs); source("tools/S1.3b_overlay.R")
src <- readLines("tools/S1.3b_run.R"); eval(parse(text = src[grep("^regen <- function", src):(grep("^reg_features <- function", src) - 1)]))
suppressPackageStartupMessages(library(parallel))
OUTR <- "results/S1.3b/review/code"
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)
P1 <- read.csv("results/S1.3b/S1.3b_posthoc_P1_pairs_v0.csv")
ulp <- function(x) 2^(floor(log2(pmax(abs(x), 2^-1022))) - 52)
verts <- function(g) { v <- s3b_geoms_df(list(g)); v }
inter_area <- function(a, b) sum(as.numeric(st_area(st_intersection(st_sfc(a), st_sfc(b)))))
setv <- function(g, df) s3b_rebuild(df, class(g)[2])$geoms[[1]]
DO_GRID <- !("nogrid" %in% REAL_ARGS)
one <- function(q) {
  cs <- P1[q, ]; G <- regen(roi_i(cs$archetype, cs$roi_id), cs$model, cs$rep); ti <- G$terr[[cs$i]]; tj <- G$terr[[cs$j]]
  vi <- verts(ti); vj <- verts(tj)
  # cross pairs within 1e-6
  D <- outer(vi$x, vj$x, "-")^2 + outer(vi$y, vj$y, "-")^2; D <- sqrt(D); w <- which(D <= 1e-6, arr.ind = TRUE)
  ident <- sum(D[w] == 0); dif <- w[D[w] > 0, , drop = FALSE]
  rb <- st_boundary(G$reg); regv <- st_coordinates(G$reg)[, 1:2]
  info <- lapply(seq_len(nrow(dif)), function(k) { a <- dif[k, 1]; b <- dif[k, 2]
    xa <- vi$x[a]; ya <- vi$y[a]; xb <- vj$x[b]; yb <- vj$y[b]
    data.frame(xi = sprintf("%.17g", xa), yi = sprintf("%.17g", ya), xj = sprintf("%.17g", xb), yj = sprintf("%.17g", yb),
      ulp_x = abs(xa - xb) / ulp(xa), ulp_y = abs(ya - yb) / ulp(ya), d = D[a, b],
      dist_reg_i = as.numeric(st_distance(st_sfc(st_point(c(xa, ya))), rb)), dist_reg_j = as.numeric(st_distance(st_sfc(st_point(c(xb, yb))), rb)),
      i_is_regvertex = any(regv[, 1] == xa & regv[, 2] == ya), j_is_regvertex = any(regv[, 1] == xb & regv[, 2] == yb)) })
  info <- do.call(rbind, info)
  a0 <- inter_area(ti, tj)
  # CP-1a per coppia: vertice di j portato su quello di i (solo i vertici diversi)
  vj2 <- vj; for (k in seq_len(nrow(dif))) { vj2$x[dif[k, 2]] <- vi$x[dif[k, 1]]; vj2$y[dif[k, 2]] <- vi$y[dif[k, 1]] }
  a_fix <- inter_area(ti, setv(tj, vj2))
  # e nel verso opposto (i portato su j)
  vi2 <- vi; for (k in seq_len(nrow(dif))) { vi2$x[dif[k, 1]] <- vj$x[dif[k, 2]]; vi2$y[dif[k, 1]] <- vj$y[dif[k, 2]] }
  a_fix_rev <- inter_area(setv(ti, vi2), tj)
  # placebo ufficiale replicato (solo per 766/1987 e' identico al codice ufficiale; qui applicato a ogni coppia con la stessa logica)
  far_j <- as.numeric(st_distance(st_sfc(lapply(seq_len(nrow(vj)), function(k) st_point(c(vj$x[k], vj$y[k])))), st_boundary(st_sfc(ti)))) > 1e-6
  amp <- if (nrow(info)) info$d else 0; km <- nrow(dif); cand <- which(far_j)
  pl_off <- do.call(rbind, lapply(1:20, function(sd) { set.seed(20261010 + 100 + sd); vv <- vj
    mv <- cand[sample.int(length(cand), km)]; th <- runif(km, 0, 2 * pi); am <- rep_len(amp, km)
    nx <- vv$x[mv] + am * cos(th); ny <- vv$y[mv] + am * sin(th); changed <- sum(nx != vv$x[mv] | ny != vv$y[mv])
    vv$x[mv] <- nx; vv$y[mv] <- ny
    data.frame(seed = sd, changed = changed, coord_mag = max(abs(c(vj$x[mv], vj$y[mv]))), area = inter_area(ti, setv(tj, vv))) }))
  # placebo corretto: stesso numero di ulp (della coordinata stessa) sullo stesso asse, su vertici non condivisi
  ux <- if (nrow(info)) round(info$ulp_x[1]) else 0; uy <- if (nrow(info)) round(info$ulp_y[1]) else 0
  sx <- if (nrow(info)) sign(vj$x[dif[1, 2]] - vi$x[dif[1, 1]]) else 0; sy <- if (nrow(info)) sign(vj$y[dif[1, 2]] - vi$y[dif[1, 1]]) else 0
  pl_ulp <- do.call(rbind, lapply(1:20, function(sd) { set.seed(777 + sd); vv <- vj; mv <- cand[sample.int(length(cand), km)]
    s <- sample(c(-1, 1), 2, TRUE)
    nx <- vv$x[mv] + s[1] * ux * ulp(vv$x[mv]); ny <- vv$y[mv] + s[2] * uy * ulp(vv$y[mv]); changed <- sum(nx != vv$x[mv] | ny != vv$y[mv])
    vv$x[mv] <- nx; vv$y[mv] <- ny; g2 <- setv(tj, vv)
    data.frame(seed = sd, changed = changed, area = inter_area(ti, g2), valid = st_is_valid(st_sfc(g2))) }))
  # controllo inverso: introdurre un vertice condiviso non identico (1 ulp) su un vertice condiviso IDENTICO della stessa coppia
  idw <- w[D[w] == 0, , drop = FALSE]
  inv <- if (nrow(idw)) do.call(rbind, lapply(seq_len(nrow(idw)), function(k) { vv <- vj; b <- idw[k, 2]
      vv$y[b] <- vv$y[b] + ulp(vv$y[b]); data.frame(k = k, dist_reg = as.numeric(st_distance(st_sfc(st_point(c(vj$x[b], vj$y[b]))), rb)), area = inter_area(ti, setv(tj, vv))) })) else NULL
  out <- list(case = cs, n_ident = ident, n_diff = nrow(dif), info = info, a0 = a0, a_fix = a_fix, a_fix_rev = a_fix_rev,
              n_far_j = sum(far_j), pl_off = pl_off, pl_ulp = pl_ulp, inv = inv, clipped = G$s$df$clipped[c(cs$i, cs$j)])
  # C-6 indipendente su questo nullo, V0 e V2
  if (DO_GRID) {
    bb <- st_bbox(G$reg); x0 <- floor(bb[["xmin"]] / 2) * 2; y0 <- floor(bb[["ymin"]] / 2) * 2
    nx <- ceiling((bb[["xmax"]] - x0) / 2); ny <- ceiling((bb[["ymax"]] - y0) / 2)
    xs <- x0 + 2 * (0:(nx - 1)); ys <- y0 + 2 * (0:(ny - 1)); gx <- rep(xs, times = ny); gy <- rep(ys, each = nx)
    grid <- st_sfc(lapply(seq_along(gx), function(k) st_polygon(list(rbind(c(gx[k], gy[k]), c(gx[k] + 2, gy[k]), c(gx[k] + 2, gy[k] + 2), c(gx[k], gy[k] + 2), c(gx[k], gy[k]))))))
    terr <- G$terr; aT <- as.numeric(st_area(terr)); P <- s3b_perim(terr)
    # V0
    pc <- st_intersection(terr, grid); ix <- attr(pc, "idx"); ap <- as.numeric(st_area(pc))
    sc <- numeric(length(terr)); s <- tapply(ap, ix[, 1], sum); sc[as.integer(names(s))] <- s
    dcell <- abs(sc - aT)
    pij <- c(cs$i, cs$j)
    reg_g <- G$reg; bi <- st_intersection(grid, reg_g); bix <- attr(bi, "idx")[, 1]; truth <- numeric(length(grid)); truth[bix] <- as.numeric(st_area(bi))
    sb <- numeric(length(grid)); s <- tapply(ap, ix[, 2], sum); sb[as.integer(names(s))] <- s; dbin <- abs(sb - truth)
    # V2
    tg <- geos_set_precision(as_geos_geometry(terr), 1e-6); gg <- as_geos_geometry(grid); rg <- geos_set_precision(as_geos_geometry(reg_g), 1e-6)
    im <- geos_intersects_matrix(tg, gg); ii <- rep(seq_along(im), lengths(im)); jj <- unlist(im)
    a2 <- geos_area(geos_intersection_prec(tg[ii], gg[jj], 1e-6)); aV2 <- geos_area(tg)
    sc2 <- numeric(length(terr)); s <- tapply(a2, ii, sum); sc2[as.integer(names(s))] <- s; dcell2 <- abs(sc2 - aV2)
    tr2 <- geos_area(geos_intersection_prec(gg, rg, 1e-6)); sb2 <- numeric(length(grid)); s <- tapply(a2, jj, sum); sb2[as.integer(names(s))] <- s; dbin2 <- abs(sb2 - tr2)
    out$grid <- c(n_bins = length(grid), n_pieces_v0 = length(ap), cell_viol_v0 = sum(dcell > 1e-6 * P), cell_max_v0 = max(dcell), pair_cells_max_v0 = max(dcell[pij]),
                  bin_max_v0 = max(dbin), bin_over_1e9 = sum(dbin > 1e-9), n_pieces_v2 = length(a2), cell_viol_v2 = sum(dcell2 > 1e-6 * P), cell_max_v2 = max(dcell2),
                  bin_max_v2 = max(dbin2), total_v0 = sum(ap) - sum(truth), total_v2 = sum(a2) - sum(tr2))
  }
  out
}
res <- mclapply(seq_len(nrow(P1)), function(q) tryCatch(one(q), error = function(e) paste("ERR", conditionMessage(e))), mc.cores = 10, mc.preschedule = FALSE)
saveRDS(res, file.path(RV, "regen_audit.rds"))
for (q in seq_along(res)) { r <- res[[q]]; if (is.character(r)) { cat(q, r, "\n"); next }
  cs <- r$case; cat(sprintf("\n== %s %s %s %d  %d/%d  clipped %s  ident %d diff %d  a0 %.4f  a_fix %.3g  a_fix_rev %.3g  n_far_j %d\n", cs$archetype, cs$roi_id, cs$model, cs$rep, cs$i, cs$j,
    paste(r$clipped, collapse = ","), r$n_ident, r$n_diff, r$a0, r$a_fix, r$a_fix_rev, r$n_far_j))
  print(r$info, row.names = FALSE)
  cat(sprintf("placebo ufficiale-like: coord changed in %d/20 seeds (sum changed %d); area intact (|d|<=0.01) %d/20; coord mag max %.1f\n",
      sum(r$pl_off$changed > 0), sum(r$pl_off$changed), sum(abs(r$pl_off$area - r$a0) <= 0.01), max(r$pl_off$coord_mag)))
  cat(sprintf("placebo ulp: coord changed %d/20; area intact %d/20; area==0 %d/20; valid %d/20\n", sum(r$pl_ulp$changed > 0), sum(abs(r$pl_ulp$area - r$a0) <= 0.01), sum(r$pl_ulp$area < 1e-6), sum(r$pl_ulp$valid)))
  if (!is.null(r$inv)) { cat("inverse (1 ulp on identical shared vertex):\n"); print(r$inv, row.names = FALSE) }
  if (!is.null(r$grid)) print(signif(r$grid, 4))
}
