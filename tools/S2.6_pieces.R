# tools/S2.6_pieces.R - S2.6: simulazione in un pezzo solo e a pezzi con le funzioni di oggi (nessuna modifica a R/).
# Definizioni: docs/steps/S2.6.md. Il pezzo e' un livello cellulare ritagliato dal livello dell'area intera: cellule con il centro nel
# riquadro del pezzo (area utile allargata del bordo, tagliata al riquadro dell'area intera), territori e nuclei tagliati al riquadro,
# regioni intersecate con il riquadro, piu' quattro "ancore" (quadratini di 1e-4 um agli angoli del riquadro) perche' il riquadro dei
# territori e delle regioni, da cui simulate_visium_hd() costruisce la griglia e il riquadro della riflessione, sia esattamente il
# riquadro del pezzo (griglia allineata a quella dell'area intera). Si tengono i quadratini di 2 um con il centro nell'area utile.
# Uso (root del repo):
#   Rscript --vanilla tools/S2.6_pieces.R one    <layer.rds> <seed> <outdir> <tag> <stripes 0|1> [frame_open 0|1] [workers]
#   Rscript --vanilla tools/S2.6_pieces.R pieces <layer.rds> <tile_um> <halo_um> <seed> <outdir> <tag> <stripes 0|1> [workers]
#   Rscript --vanilla tools/S2.6_pieces.R cut    <layer.rds> <tile_um> <halo_um> <p1,p2,...> <outdir>   livelli di alcuni pezzi
# Uscite per esecuzione in <outdir>/<tag>/: u2.bin (UMI per quadratino di 2 um sulla griglia dell'area intera, int32, ordine di bin_id),
# c8.rds (matrice geni x blocchi di 8 um), own8.bin / dup8.bin (trascritti catturati per blocco da cellule del pezzo / da cellule di un
# altro pezzo, int32), run.json (tempi, picco di memoria, semi, pezzi).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
peak_mb <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }

# riquadro dell'area intera e griglia di riferimento (come build_hd_grid())
s26_global_grid <- function(L, margin_um = 280) {
  bb <- sf::st_bbox(c(sf::st_set_crs(L$cell_territories, NA), sf::st_set_crs(L$region_polygons, NA)))
  list(bbox = bb, grid = hd_grid_geometry(bb, margin_um = margin_um, bin_um = 2))
}
.s26_box_sfc <- function(b) sf::st_sfc(sf::st_polygon(list(rbind(c(b[1], b[2]), c(b[3], b[2]), c(b[3], b[4]), c(b[1], b[4]), c(b[1], b[2])))))
# solo la parte poligonale di un'intersezione (un territorio tangente al lato del riquadro puo' dare poligono + segmento)
.s26_polygonal <- function(x) {
  ty <- as.character(sf::st_geometry_type(x))
  for (i in which(ty == "GEOMETRYCOLLECTION")) x[i] <- sf::st_union(sf::st_collection_extract(x[i], "POLYGON", warn = FALSE))
  x
}
# livello del pezzo: cellule con il centro in box (semiaperto a destra e in alto, salvo il lato dell'area intera), geometrie tagliate
s26_piece_layer <- function(L, box, gbb, anchor = 1e-4) {
  cd <- L$cell_df
  inx <- (cd$x >= box[1] & cd$x < box[3]) | (box[3] >= gbb[["xmax"]] & cd$x == box[3])
  iny <- (cd$y >= box[2] & cd$y < box[4]) | (box[4] >= gbb[["ymax"]] & cd$y == box[4])
  keep <- which(inx & iny)
  B <- .s26_box_sfc(box)
  Tt <- sf::st_set_crs(L$cell_territories[keep], NA); Nn <- sf::st_set_crs(L$cell_nuclei[keep], NA)
  gT <- geos::as_geos_geometry(Tt); gN <- geos::as_geos_geometry(Nn); gB <- geos::as_geos_geometry(B)
  cross <- which(!geos::geos_within(gT, gB))
  if (length(cross)) {
    Tt[cross] <- .s26_polygonal(sf::st_as_sfc(geos::geos_write_wkb(geos::geos_intersection(gT[cross], gB))))
    ncr <- cross[!geos::geos_within(gN[cross], gB)]
    if (length(ncr)) Nn[ncr] <- .s26_polygonal(sf::st_as_sfc(geos::geos_write_wkb(geos::geos_intersection(gN[ncr], gB))))
  }
  if (all(sf::st_geometry_type(Tt) == "POLYGON")) Tt <- sf::st_cast(Tt, "POLYGON")
  if (all(sf::st_geometry_type(Nn) == "POLYGON")) Nn <- sf::st_cast(Nn, "POLYGON")
  R <- sf::st_set_crs(L$region_polygons, NA)
  gR <- geos::geos_make_valid(geos::geos_intersection(geos::as_geos_geometry(R), gB))
  Rp <- sf::st_as_sfc(geos::geos_write_wkb(gR[!geos::geos_is_empty(gR)]))
  Rp <- suppressWarnings(sf::st_collection_extract(Rp, "POLYGON"))
  Rp <- sf::st_cast(Rp, "POLYGON", warn = FALSE)
  a <- anchor
  anc <- sf::st_sfc(lapply(list(c(box[1], box[2]), c(box[3] - a, box[2]), c(box[1], box[4] - a), c(box[3] - a, box[4] - a)),
                           function(p) sf::st_polygon(list(rbind(p, p + c(a, 0), p + c(a, a), p + c(0, a), p)))))
  P <- L
  P$cell_df <- cd[keep, , drop = FALSE]; rownames(P$cell_df) <- NULL
  P$cell_territories <- Tt; P$cell_nuclei <- Nn
  P$region_polygons <- c(Rp, anc)
  P$metadata$n_cells <- length(keep)
  attr(P, "n_clipped") <- length(cross)
  P
}
# indici globali dei quadratini di 2 um di una griglia di pezzo con il centro in U
s26_keep_bins <- function(gp, gg, U) {
  dj <- (gp$x0 - gg$x0) / gp$bin_um; di <- (gp$y0 - gg$y0) / gp$bin_um
  if (abs(dj - round(dj)) > 1e-6 || abs(di - round(di)) > 1e-6) stop("griglia del pezzo non allineata", call. = FALSE)
  dj <- as.integer(round(dj)); di <- as.integer(round(di))
  jj <- seq_len(gp$n_col) - 1L; ii <- seq_len(gp$n_row) - 1L
  xc <- gp$x0 + (jj + 0.5) * gp$bin_um; yc <- gp$y0 + (ii + 0.5) * gp$bin_um
  kj <- jj[xc >= U[1] & xc < U[3]]; ki <- ii[yc >= U[2] & yc < U[4]]
  li <- rep(ki, each = length(kj)); lj <- rep(kj, times = length(ki))
  gi <- li + di; gj <- lj + dj
  ok <- gi >= 0 & gi < gg$n_row & gj >= 0 & gj < gg$n_col
  list(local = (li * gp$n_col + lj + 1L)[ok], global = (gi * gg$n_col + gj + 1L)[ok])
}
b8_of <- function(id, g) { i <- (id - 1L) %/% g$n_col; j <- (id - 1L) %% g$n_col; as.integer((i %/% g$block) * (g$n_col %/% g$block) + j %/% g$block + 1L) }
s26_layers <- function(stripes) list(stripes = stripes, background = TRUE, field = TRUE, second_tail = TRUE)
s26_seeds <- function(seed, n) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); sample.int(1e8, n) }

save_run <- function(dir, u2, c8, own8, dup8, info) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeBin(as.integer(u2), file.path(dir, "u2.bin"), size = 4L); saveRDS(c8, file.path(dir, "c8.rds"))
  writeBin(as.integer(own8), file.path(dir, "own8.bin"), size = 4L); writeBin(as.integer(dup8), file.path(dir, "dup8.bin"), size = 4L)
  writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, digits = NA), file.path(dir, "run.json.tmp")); file.rename(file.path(dir, "run.json.tmp"), file.path(dir, "run.json"))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE); MODE <- a[1]
  if (MODE == "one") {
    L <- readRDS(a[2]); seed <- as.integer(a[3]); OUT <- a[4]; tag <- a[5]; st <- a[6] == "1"
    fo <- if (length(a) >= 7) a[7] == "1" else TRUE; W <- if (length(a) >= 8) as.integer(a[8]) else 1L
    dir <- file.path(OUT, tag); if (file.exists(file.path(dir, "run.json"))) quit(status = 0)
    G <- s26_global_grid(L); gg <- G$grid
    gc(); t0 <- proc.time()[["elapsed"]]
    out <- simulate_visium_hd(L, seed = seed, frame_open = fo, layers = s26_layers(st), workers = W)
    te <- proc.time()[["elapsed"]] - t0
    stopifnot(isTRUE(all.equal(out$metadata$grid[c("x0", "y0", "n_row", "n_col")], gg[c("x0", "y0", "n_row", "n_col")], tolerance = 0)))
    u2 <- as.integer(round(Matrix::colSums(out$counts_2um)))
    og <- out$truth$origin_8um; n8 <- (gg$n_row %/% 4L) * (gg$n_col %/% 4L)
    own8 <- tabulate(rep(og$bin8_id, og$n), n8)
    info <- list(mode = "one", tag = tag, layer = basename(a[2]), seed = seed, stripes = st, frame_open = fo, workers = W,
                 archetype = out$metadata$archetype, n_cells = nrow(L$cell_df), grid = gg[c("x0", "y0", "n_row", "n_col")],
                 n_generated = out$metadata$n_generated, n_captured = out$metadata$n_captured, n_background = out$metadata$n_background,
                 M = out$metadata$overgeneration, chunk_size = out$metadata$chunk_size,
                 elapsed = as.list(out$metadata$elapsed_s), elapsed_s = te, peak_mb = peak_mb(), load = readLines("/proc/loadavg"))
    save_run(dir, u2, out$counts_8um, own8, integer(n8), info)
    cat(sprintf("[one %s] %.3g generati, %.0f s, %.0f MB\n", tag, info$n_generated, te, info$peak_mb))
  } else if (MODE == "pieces") {
    L <- readRDS(a[2]); Tu <- as.numeric(a[3]); H <- as.numeric(a[4]); seed <- as.integer(a[5]); OUT <- a[6]; tag <- a[7]; st <- a[8] == "1"
    W <- if (length(a) >= 9) as.integer(a[9]) else 1L
    dir <- file.path(OUT, tag); if (file.exists(file.path(dir, "run.json"))) quit(status = 0)
    pdir <- file.path(dir, "pieces"); dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    G <- s26_global_grid(L); gg <- G$grid; gbb <- G$bbox
    xs <- seq(gbb[["xmin"]], gbb[["xmax"]], by = Tu); if (tail(xs, 1) < gbb[["xmax"]] - 1e-9) xs <- c(xs, gbb[["xmax"]])
    ys <- seq(gbb[["ymin"]], gbb[["ymax"]], by = Tu); if (tail(ys, 1) < gbb[["ymax"]] - 1e-9) ys <- c(ys, gbb[["ymax"]])
    P <- expand.grid(ix = seq_len(length(xs) - 1L), iy = seq_len(length(ys) - 1L))
    sd_p <- s26_seeds(seed, nrow(P))
    one_piece <- function(p) {
      fp <- file.path(pdir, sprintf("p%03d.rds", p)); if (file.exists(fp)) return(fp)
      U <- c(xs[P$ix[p]], ys[P$iy[p]], xs[P$ix[p] + 1L], ys[P$iy[p] + 1L])
      box <- c(max(gbb[["xmin"]], U[1] - H), max(gbb[["ymin"]], U[2] - H), min(gbb[["xmax"]], U[3] + H), min(gbb[["ymax"]], U[4] + H))
      t0 <- proc.time()[["elapsed"]]
      Lp <- s26_piece_layer(L, box, gbb)
      res <- list(p = p, U = U, box = box, n_cells = nrow(Lp$cell_df), n_clipped = attr(Lp, "n_clipped"), seed = sd_p[p])
      if (nrow(Lp$cell_df) == 0L) { res$empty <- TRUE; saveRDS(res, fp); return(fp) }
      rss0 <- { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmRSS", s)])) / 1024 }
      out <- simulate_visium_hd(Lp, seed = sd_p[p], frame_open = TRUE, layers = s26_layers(st), workers = 1L)
      cl <- out$truth$cells
      res$cons_ok <- all(cl$n_generated == cl$n_captured + cl$n_thinned + cl$n_lost); res$rss_start_mb <- rss0
      Uk <- U                                               # sui lati dell'area intera si tiene anche il margine di vetro
      if (U[1] <= gbb[["xmin"]]) Uk[1] <- gg$x0; if (U[2] <= gbb[["ymin"]]) Uk[2] <- gg$y0
      if (U[3] >= gbb[["xmax"]]) Uk[3] <- gg$x0 + gg$n_col * gg$bin_um; if (U[4] >= gbb[["ymax"]]) Uk[4] <- gg$y0 + gg$n_row * gg$bin_um
      kb <- s26_keep_bins(out$metadata$grid, gg, Uk)
      X <- out$counts_2um[, kb$local, drop = FALSE]
      og <- out$truth$origin; m <- match(og$bin_id, kb$local); og <- og[!is.na(m), ]; og$gbin <- kb$global[m[!is.na(m)]]
      cdp <- Lp$cell_df; own_ids <- cdp$cell_id[cdp$x >= U[1] & (cdp$x < U[3] | (U[3] >= gbb[["xmax"]] & cdp$x == U[3])) &
                                                 cdp$y >= U[2] & (cdp$y < U[4] | (U[4] >= gbb[["ymax"]] & cdp$y == U[4]))]
      og$own <- og$cell_id %in% own_ids
      res$gbin <- kb$global; res$X <- X; res$og <- og[, c("gbin", "n", "own")]
      res$n_generated <- out$metadata$n_generated; res$M <- out$metadata$overgeneration
      res$elapsed_s <- proc.time()[["elapsed"]] - t0; res$peak_mb <- peak_mb(); res$sim_elapsed <- as.list(out$metadata$elapsed_s)
      saveRDS(res, fp); fp
    }
    t0 <- proc.time()[["elapsed"]]
    r <- parallel::mclapply(seq_len(nrow(P)), function(p) tryCatch(one_piece(p), error = function(e) paste("ERRORE", p, conditionMessage(e))),
                            mc.cores = W, mc.preschedule = FALSE)
    err <- unlist(r)[grepl("^ERRORE", unlist(r))]; if (length(err)) { print(err); quit(status = 1) }
    tw <- proc.time()[["elapsed"]] - t0
    nb <- gg$n_row * gg$n_col; n8 <- (gg$n_row %/% 4L) * (gg$n_col %/% 4L)
    u2 <- integer(nb); own8 <- integer(n8); dup8 <- integer(n8); cover <- integer(nb); Xs <- list(); ids <- list(); pinfo <- list()
    for (p in seq_len(nrow(P))) {
      z <- readRDS(file.path(pdir, sprintf("p%03d.rds", p)))
      pinfo[[p]] <- z[intersect(names(z), c("p", "U", "box", "n_cells", "n_clipped", "seed", "n_generated", "M", "elapsed_s", "peak_mb", "rss_start_mb", "cons_ok", "empty"))]
      if (isTRUE(z$empty)) next
      cover[z$gbin] <- cover[z$gbin] + 1L
      u2[z$gbin] <- as.integer(round(Matrix::colSums(z$X)))
      b8 <- b8_of(z$og$gbin, gg)
      own8 <- own8 + tabulate(rep(b8[z$og$own], z$og$n[z$og$own]), n8); dup8 <- dup8 + tabulate(rep(b8[!z$og$own], z$og$n[!z$og$own]), n8)
      Xs[[length(Xs) + 1L]] <- z$X; ids[[length(ids) + 1L]] <- z$gbin
    }
    X2 <- do.call(cbind, Xs); gid <- unlist(ids)
    # somma a 8 um sulla griglia dell'area intera (come sum_to_8um, con le colonne nell'ordine di gid)
    b8 <- b8_of(gid, gg)
    Sm <- Matrix::sparseMatrix(i = seq_along(gid), j = b8, x = 1, dims = c(length(gid), n8))
    c8 <- X2 %*% Sm
    info <- list(mode = "pieces", tag = tag, layer = basename(a[2]), seed = seed, tile_um = Tu, halo_um = H, stripes = st, workers = W,
                 n_pieces = nrow(P), grid = gg[c("x0", "y0", "n_row", "n_col")], cover_max = max(cover),
                 n_bins_covered = sum(cover > 0), pieces = pinfo, elapsed_wall_s = tw,
                 elapsed_process_s = sum(vapply(pinfo, function(z) if (is.null(z$elapsed_s)) 0 else z$elapsed_s, 0)),
                 peak_mb = peak_mb(), load = readLines("/proc/loadavg"))
    save_run(dir, u2, c8, own8, dup8, info)
    cat(sprintf("[pieces %s] %d pezzi, %.0f s di orologio, %.0f s di processo\n", tag, nrow(P), tw, info$elapsed_process_s))
  } else if (MODE == "cut") {
    # livelli di alcuni pezzi, per simularli in processi separati (`one` sul livello del pezzo = simulazione del pezzo)
    L <- readRDS(a[2]); Tu <- as.numeric(a[3]); H <- as.numeric(a[4]); only <- as.integer(strsplit(a[5], ",")[[1]]); OUT <- a[6]
    dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
    G <- s26_global_grid(L); gbb <- G$bbox
    xs <- seq(gbb[["xmin"]], gbb[["xmax"]], by = Tu); if (tail(xs, 1) < gbb[["xmax"]] - 1e-9) xs <- c(xs, gbb[["xmax"]])
    ys <- seq(gbb[["ymin"]], gbb[["ymax"]], by = Tu); if (tail(ys, 1) < gbb[["ymax"]] - 1e-9) ys <- c(ys, gbb[["ymax"]])
    P <- expand.grid(ix = seq_len(length(xs) - 1L), iy = seq_len(length(ys) - 1L)); info <- list()
    for (p in only) {
      U <- c(xs[P$ix[p]], ys[P$iy[p]], xs[P$ix[p] + 1L], ys[P$iy[p] + 1L])
      box <- c(max(gbb[["xmin"]], U[1] - H), max(gbb[["ymin"]], U[2] - H), min(gbb[["xmax"]], U[3] + H), min(gbb[["ymax"]], U[4] + H))
      Lp <- s26_piece_layer(L, box, gbb); saveRDS(Lp, file.path(OUT, sprintf("piece_p%03d.rds", p)))
      info[[length(info) + 1]] <- list(p = p, U = U, box = box, useful_mm2 = prod(U[3:4] - U[1:2]) / 1e6, box_mm2 = prod(box[3:4] - box[1:2]) / 1e6,
                                       n_cells = nrow(Lp$cell_df), n_clipped = attr(Lp, "n_clipped"))
    }
    writeLines(jsonlite::toJSON(list(layer = basename(a[2]), tile_um = Tu, halo_um = H, pieces = info), auto_unbox = TRUE, digits = NA), file.path(OUT, "cut.json"))
    cat(sprintf("[cut] %d pezzi\n", length(only)))
  } else stop("modo sconosciuto")
}
