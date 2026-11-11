# tools/S2.8_run.R - S2.8: esecuzioni del calcolo (docs/steps/S2.8.md): tempo e memoria dei due motori, uguaglianza con S2.7 e S2.5,
# controprova CP-2, descrittiva di BL-237. Uso (root del repo):
#   Rscript --vanilla tools/S2.8_run.R pieces <tag> <seed> <piece_um|Inf> <outdir> <k> [default|window] [ceiling_q]
#   Rscript --vanilla tools/S2.8_run.R ref    <tag> <seed> <k> <outdir>
# Livelli: /mnt/micron/geo_spatialtrans/S1.6/layers/<tag>_U6.rds (come P-1 di S2.7); frame_open = TRUE; un processo (workers = 1).
# "ref" = motore di S2.5 che calcola da se' griglia e verita' geometrica (nessun `hd` dato), blocchi di 4e6 come il pannello di S2.5.
# "window" = riquadro della finestra di R2 dell'area (0, 0, lato, lato), con lato = massimo di xmax e ymax di regioni e territori.
# Uscite in <outdir>/<tag>_k<k>: .bin, .g2.bin, .u8.bin, .g8.bin (come tools/S2.7_run.R), .tis.bin (in_tissue per quadratino, int8),
# .json (tempi, memoria, metadati), .md5 (impronte delle uscite come S2.7, piu' layers_s27 = impronta degli strati riportati ai nomi di S2.7).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
peak_mb <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }
md5_of <- function(x) { f <- tempfile(); saveRDS(x, f, compress = FALSE); on.exit(unlink(f)); unname(tools::md5sum(f)) }
write_bins <- function(fo, X2, X8) {
  writeBin(as.integer(round(Matrix::colSums(X2))), paste0(fo, ".bin"), size = 4L); writeBin(as.integer(diff(X2@p)), paste0(fo, ".g2.bin"), size = 4L)
  X8 <- as(X8, "CsparseMatrix")
  writeBin(as.integer(round(Matrix::colSums(X8))), paste0(fo, ".u8.bin"), size = 4L); writeBin(as.integer(diff(X8@p)), paste0(fo, ".g8.bin"), size = 4L)
}
wjson <- function(x, f) { writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, null = "null"), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f) }
a <- commandArgs(TRUE); MODE <- a[1]; tag <- a[2]; seed <- as.numeric(a[3]); A <- strsplit(tag, "_")[[1]][1]
layf <- sprintf("/mnt/micron/geo_spatialtrans/S1.6/layers/%s_U6.rds", tag)
grid_json <- function(g) list(x0 = g$x0, y0 = g$y0, n_row = g$n_row, n_col = g$n_col, bin_um = g$bin_um, block = g$block)
if (MODE == "pieces") {
  pu <- as.numeric(a[4]); OUT <- a[5]; k <- as.integer(a[6]); fr <- if (length(a) >= 7) a[7] else "default"; cq <- if (length(a) >= 8) as.numeric(a[8]) else 1e-5
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE); base <- file.path(OUT, sprintf("%s_k%d", tag, k))
  if (file.exists(paste0(base, ".json"))) quit(status = 0L)
  frame <- NULL
  if (fr == "window") { L <- readRDS(layf); bb <- sf::st_bbox(c(sf::st_set_crs(L$cell_territories, NA), sf::st_set_crs(L$region_polygons, NA)))
    side <- max(bb[["xmax"]], bb[["ymax"]]); frame <- c(xmin = 0, ymin = 0, xmax = side, ymax = side); rm(L); invisible(gc(FALSE)) }
  t0 <- proc.time()[["elapsed"]]
  out <- simulate_visium_hd_pieces(layf, seed = seed, piece_um = pu, frame_open = TRUE, frame = frame, ceiling_q = cq, out_dir = paste0(base, "_pieces"),
                                   workers = 1L, assemble = "full", verbose = TRUE)
  te <- proc.time()[["elapsed"]] - t0
  m <- out$metadata; ps <- m$piece_stats
  j <- list(tag = tag, archetype = m$archetype, k = k, mode = "pieces", seed = seed, piece_um = m$piece_um, n_pieces = m$n_pieces, frame_kind = fr,
            frame = as.list(m$frame), ceiling_q = cq, grid = grid_json(m$grid), n_cells = m$n_cells, overgeneration = m$overgeneration,
            n_generated = m$n_generated, n_captured = m$n_captured, n_sent = m$n_sent_between_pieces, checks = m$checks,
            elapsed = as.list(m$elapsed_s), elapsed_total_s = te, main_peak_mb = peak_mb(), piece_peak_max_mb = max(ps$peak_mb, na.rm = TRUE),
            stage_time_sum_s = as.list(tapply(ps$elapsed_s, ps$stage, sum, na.rm = TRUE)), n_reused = sum(ps$reused),
            package_version = m$package_version, code_hash = m$code_hash, profiles_hash = m$profiles_hash,
            load_end = paste(readLines("/proc/loadavg"), collapse = ""))
  j$peak_mb <- max(j$main_peak_mb, j$piece_peak_max_mb)
  write_bins(base, out$counts_2um, out$counts_8um); writeBin(as.integer(out$bins$in_tissue), paste0(base, ".tis.bin"), size = 1L)
  tr <- out$truth; ly27 <- tr$layers
  if (!is.null(ly27$stripes)) { ly27$stripes$M <- ly27$stripes$M_raw; ly27$stripes$M_raw <- NULL }
  h <- list(counts_2um = md5_of(out$counts_2um), counts_8um = md5_of(out$counts_8um), bins = md5_of(out$bins), bins_8um = md5_of(out$bins_8um),
            overlap = md5_of(tr$overlap), overlap_8um = md5_of(tr$overlap_8um), origin = md5_of(tr$origin), origin_8um = md5_of(tr$origin_8um),
            background = md5_of(tr$background), background_8um = md5_of(tr$background_8um), arrived_2um = md5_of(tr$arrived_2um),
            cells = md5_of(tr$cells), genes = md5_of(tr$genes), layers = md5_of(tr$layers), layers_s27 = md5_of(ly27))
  wjson(h, paste0(base, ".md5")); wjson(j, paste0(base, ".json")); unlink(paste0(base, "_pieces"), recursive = TRUE)
  cat(sprintf("[%s] %s pezzi %s, %.0f s, picco %.0f MB\n", basename(base), A, format(pu), te, j$peak_mb))
} else if (MODE == "ref") {
  k <- as.integer(a[4]); OUT <- a[5]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE); base <- file.path(OUT, sprintf("%s_k%d", tag, k))
  if (file.exists(paste0(base, ".json"))) quit(status = 0L)
  t0 <- proc.time()[["elapsed"]]
  L <- readRDS(layf)
  out <- simulate_visium_hd(L, seed = seed, archetype = A, frame_open = TRUE, workers = 1L, chunk_transcripts = 4e6)
  te <- proc.time()[["elapsed"]] - t0
  m <- out$metadata
  write_bins(base, out$counts_2um, out$counts_8um); writeBin(as.integer(out$bins$in_tissue), paste0(base, ".tis.bin"), size = 1L)
  wjson(list(tag = tag, archetype = A, k = k, mode = "ref", seed = seed, grid = grid_json(m$grid), overgeneration = m$overgeneration,
             n_captured = m$n_captured, elapsed = as.list(m$elapsed_s), elapsed_total_s = te, peak_mb = peak_mb(),
             package_version = m$package_version, code_hash = m$code_hash, load_end = paste(readLines("/proc/loadavg"), collapse = "")), paste0(base, ".json"))
  cat(sprintf("[%s ref] %.0f s, picco %.0f MB\n", basename(base), te, peak_mb()))
} else stop("modo sconosciuto")
